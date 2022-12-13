/**
 * @file B_STS.h
 * @brief Header file for B_STS.c
 */
#ifndef B_STS_H
#define B_STS_H
#include "../ascot5.h"
#include "../linint/linint.h"
#include "../spline/interp.h"

/**
 * @brief stellarator magnetic field parameters on the host
 */
typedef struct {
    int psigrid_n_r;       /**< Number of R grid points in psi data           */
    int psigrid_n_z;       /**< Number of z grid points in psi data           */
    int psigrid_n_phi;     /**< Number of phi grid points in psi data         */
    real psigrid_r_min;    /**< Minimum R grid point in psi data [m]          */
    real psigrid_r_max;    /**< Maximum R grid point in psi data [m]          */
    real psigrid_z_min;    /**< Minimum z grid point in psi data [m]          */
    real psigrid_z_max;    /**< Maximum z grid point in psi data [m]          */
    real psigrid_phi_min;  /**< Minimum phi grid point in psi data [rad]      */
    real psigrid_phi_max;  /**< Maximum phi grid point in psi data [rad]      */

    int Bgrid_n_r;       /**< Number of R grid points in B data               */
    int Bgrid_n_z;       /**< Number of z grid points in B data               */
    int Bgrid_n_phi;     /**< Number of phi grid points in B data             */
    real Bgrid_r_min;    /**< Minimum R coordinate in the grid in B data [m]  */
    real Bgrid_r_max;    /**< Maximum R coordinate in the grid in B data [m]  */
    real Bgrid_z_min;    /**< Minimum z coordinate in the grid in B data [m]  */
    real Bgrid_z_max;    /**< Maximum z coordinate in the grid in B data [m]  */
    real Bgrid_phi_min;  /**< Minimum phi grid point in B data [rad]          */
    real Bgrid_phi_max;  /**< Maximum phi grid point in B data [rad]          */

    real psi0;           /**< Poloidal flux value at magnetic axis [V*s*m^-1] */
    real psi1;           /**< Poloidal flux value at separatrix [V*s*m^-1]    */
    int offload_array_length; /**< Number of elements in offload_array        */

    int n_axis;          /**< Number of phi grid points in axis data          */
    real axis_min;       /**< Minimum phi grid point in axis data [rad]       */
    real axis_max;       /**< Maximum phi grid point in axis data [rad]       */
    real axis_grid;      /**< phi grid interval in axis data [rad]            */
} B_STS_offload_data;

/**
 * @brief stellarator magnetic field parameters on the target
 */
typedef struct {
    real psi0;           /**< Poloidal flux value at magnetic axis [v*s*m^-1] */
    real psi1;           /**< Poloidal flux value at separatrix [V*s*m^-1]    */
    linint1D_data axis_r;/**< 1D axis r-value interpolation data struct       */
    linint1D_data axis_z;/**< 1D axis z-value interpolation data struct       */
    interp3D_data psi;   /**< 3D psi interpolation data struct                */
    interp3D_data B_r;   /**< 3D B_r interpolation data struct                */
    interp3D_data B_phi ;/**< 3D B_phi interpolation data struct              */
    interp3D_data B_z;   /**< 3D B_z interpolation data struct                */
} B_STS_data;

int B_STS_init_offload(B_STS_offload_data* offload_data, real** offload_array);
void B_STS_free_offload(B_STS_offload_data* offload_data, real** offload_array);

DECLARE_TARGET
void B_STS_init(B_STS_data* Bdata, B_STS_offload_data* offload_data,
               real* offload_array);
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
a5err B_STS_eval_psi(real* psi, real r, real phi, real z, B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    interperr += interp3Dcomp_eval_f(&psi[0], &Bdata->psi, r, phi, z);

    /* Test for psi interpolation error */
    if(interperr) {
        return error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );
    }

    return err;
}//;
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
a5err B_STS_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                          B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real psi_dpsi_temp[10];

    interperr += interp3Dcomp_eval_df(psi_dpsi_temp, &Bdata->psi, r, phi, z);

    psi_dpsi[0] = psi_dpsi_temp[0];
    psi_dpsi[1] = psi_dpsi_temp[1];
    psi_dpsi[2] = psi_dpsi_temp[2];
    psi_dpsi[3] = psi_dpsi_temp[3];

    if(interperr) {
        return error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );
    }

    return err;
}
//;
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
a5err B_STS_eval_rho(real* rho, real psi, B_STS_data* Bdata) {
    a5err err = 0;

    /* Check that the values seem valid */
    real delta = (Bdata->psi1 - Bdata->psi0);
    if( (psi - Bdata->psi0) / delta < 0 ) {
         return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_STS );
    }

    rho[0] = sqrt( (psi - Bdata->psi0) / delta );

    return err;
}
//;
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
a5err B_STS_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                          B_STS_data* Bdata) {
    a5err err = 0;
    real psi_dpsi[4];

    B_STS_eval_psi_dpsi(psi_dpsi, r, phi, z, Bdata);

    /* Check that the values seem valid */
    real delta = Bdata->psi1 - Bdata->psi0;
    if( (psi_dpsi[0] - Bdata->psi0) / delta < 0 ) {
         return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_STS );
    }

    /* Normalize psi to get rho */
    rho_drho[0] = sqrt( (psi_dpsi[0] - Bdata->psi0) / delta );

    rho_drho[1] = psi_dpsi[1] / (2*delta*rho_drho[0]);
    rho_drho[2] = psi_dpsi[2] / (2*delta*rho_drho[0]);
    rho_drho[3] = psi_dpsi[3] / (2*delta*rho_drho[0]);

    return err;
}
//;
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
a5err B_STS_eval_B(real B[3], real r, real phi, real z, B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    interperr += interp3Dcomp_eval_f(&B[0], &Bdata->B_r, r, phi, z);
    interperr += interp3Dcomp_eval_f(&B[1], &Bdata->B_phi, r, phi, z);
    interperr += interp3Dcomp_eval_f(&B[2], &Bdata->B_z, r, phi, z);

    /* Test for B field interpolation error */
    if(interperr) {
        return error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );
    }

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) == 0);
    if(!err && check) {
        return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_STS );
    }

    return err;

}
//;
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
a5err B_STS_eval_B_dB(real B_dB[12], real r, real phi, real z,
                      B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real B_dB_temp[10];

    interperr += interp3Dcomp_eval_df(B_dB_temp, &Bdata->B_r, r, phi, z);

    B_dB[0] = B_dB_temp[0];
    B_dB[1] = B_dB_temp[1];
    B_dB[2] = B_dB_temp[2];
    B_dB[3] = B_dB_temp[3];

    interperr += interp3Dcomp_eval_df(B_dB_temp, &Bdata->B_phi, r, phi, z);

    B_dB[4] = B_dB_temp[0];
    B_dB[5] = B_dB_temp[1];
    B_dB[6] = B_dB_temp[2];
    B_dB[7] = B_dB_temp[3];

    interperr += interp3Dcomp_eval_df(B_dB_temp, &Bdata->B_z, r, phi, z);

    B_dB[8] = B_dB_temp[0];
    B_dB[9] = B_dB_temp[1];
    B_dB[10] = B_dB_temp[2];
    B_dB[11] = B_dB_temp[3];

    /* Test for B field interpolation error */
    if(interperr) {
        return error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );
    }

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B_dB[0]*B_dB[0] + B_dB[4]*B_dB[4] + B_dB[8]*B_dB[8]) == 0);
    if(!err && check) {
        return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_STS );
    }

    return 0;
}//;
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
a5err B_STS_get_axis_r(real* axis_r, B_STS_data* Bdata, real phi) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    interperr += linint1D_eval_f(axis_r, &Bdata->axis_r, phi);
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );
    }
    return err;
}//;
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
a5err B_STS_get_axis_z(real* axis_z, B_STS_data* Bdata, real phi) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    interperr += linint1D_eval_f(axis_z, &Bdata->axis_z, phi);
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );
    }
    return err;
}//;
DECLARE_TARGET_END
#endif
