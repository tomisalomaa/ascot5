/**
 * @file B_2DS.h
 * @brief Header file for B_2DS.c
 *
 * Contains declaration of B_2DS_offload_data and B_2DS_data structs.
 */
#ifndef B_2DS_H
#define B_2DS_H
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp.h"

/**
 * @brief 2D magnetic field parameters that will be offloaded to target
 */
typedef struct {
    int n_r;                  /**< Number of r grid points                    */
    int n_z;                  /**< Number of z grid points                    */
    real r_min;               /**< Minimum R coordinate in the grid [m]       */
    real r_max;               /**< Maximum R coordinate in the grid [m]       */
    real z_min;               /**< Minimum z coordinate in the grid [m]       */
    real z_max;               /**< Maximum z coordinate in the grid [m]       */
    real psi0;                /**< Poloidal flux at magnetic axis [V*s*m^-1]  */
    real psi1;                /**< Poloidal flux at separatrix [V*s*m^-1]     */
    real axis_r;              /**< R coordinate of magnetic axis [m]          */
    real axis_z;              /**< z coordinate of magnetic axis [m]          */
    int offload_array_length; /**< Number of elements in offload_array        */
} B_2DS_offload_data;

/**
 * @brief 2D magnetic field parameters on the target
 */
typedef struct {
    real psi0;           /**< Poloidal flux value at magnetic axis [V*s*m^-1] */
    real psi1;           /**< Poloidal flux value at separatrix [V*s*m^-1]    */
    real axis_r;         /**< R coordinate of magnetic axis [m]               */
    real axis_z;         /**< z coordinate of magnetic axis [m]               */
    interp2D_data psi;   /**< psi interpolation 2D spline struct              */
    interp2D_data B_r;   /**< B_R interpolation 2D spline struct              */
    interp2D_data B_phi; /**< B_phi interpolation 2D spline struct            */
    interp2D_data B_z;   /**< B_z interpolation 2D spline struct              */
} B_2DS_data;

int B_2DS_init_offload(B_2DS_offload_data* offload_data, real** offload_array);
void B_2DS_free_offload(B_2DS_offload_data* offload_data, real** offload_array);

DECLARE_TARGET
void B_2DS_init(B_2DS_data* Bdata, B_2DS_offload_data* offload_data,
                real* offload_array);
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
a5err B_2DS_eval_psi(real* psi, real r, real phi, real z, B_2DS_data* Bdata) {

    int interperr = 0;
    interperr += interp2Dcomp_eval_f(&psi[0], &Bdata->psi, r, z);

    a5err err = 0;
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
    }
    return err;
}
//;
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
a5err B_2DS_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                          B_2DS_data* Bdata) {

    int interperr = 0;
    real psi_dpsi_temp[6];

    interperr += interp2Dcomp_eval_df(psi_dpsi_temp, &Bdata->psi, r, z);
    psi_dpsi[0] = psi_dpsi_temp[0];
    psi_dpsi[1] = psi_dpsi_temp[1];
    psi_dpsi[2] = 0;
    psi_dpsi[3] = psi_dpsi_temp[2];

     a5err err = 0;
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
    }

    return err;
}
//;
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
a5err B_2DS_eval_rho(real* rho, real psi, B_2DS_data* Bdata) {

    /* Check that the values seem valid */
    real delta = (Bdata->psi1 - Bdata->psi0);
    if( (psi - Bdata->psi0) / delta < 0 ) {
         return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_2DS );
    }

    rho[0] = sqrt( (psi - Bdata->psi0) / delta );
    return 0;
}
//;
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
a5err B_2DS_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                          B_2DS_data* Bdata) {
    int interperr = 0;
    real psi_dpsi[6];

    interperr += interp2Dcomp_eval_df(psi_dpsi, &Bdata->psi, r, z);

    a5err err = 0;
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
    }

    /* Check that the values seem valid */
    real delta = Bdata->psi1 - Bdata->psi0;
    if( !err && (psi_dpsi[0] - Bdata->psi0) / delta < 0 ) {
         err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_2DS );
    }

    /* Normalize psi to get rho */
    rho_drho[0] = sqrt(fabs((psi_dpsi[0] - Bdata->psi0) / delta));

    rho_drho[1] = psi_dpsi[1] / (2*delta*rho_drho[0]);
    rho_drho[2] = 0;
    rho_drho[3] = psi_dpsi[2] / (2*delta*rho_drho[0]);

    return err;
}//;
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
a5err B_2DS_eval_B(real B[3], real r, real phi, real z, B_2DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0;

    interperr += interp2Dcomp_eval_f(&B[0], &Bdata->B_r, r, z);
    interperr += interp2Dcomp_eval_f(&B[1], &Bdata->B_phi, r, z);
    interperr += interp2Dcomp_eval_f(&B[2], &Bdata->B_z, r, z);

    /* Test for B field interpolation error */
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
    }

    if(!err) {
        real psi_dpsi[6];
        interperr += interp2Dcomp_eval_df(psi_dpsi, &Bdata->psi, r, z);
        B[0] = B[0] - psi_dpsi[2]/r;
        B[2] = B[2] + psi_dpsi[1]/r;

        /* Test for psi interpolation error */
        if(interperr) {
            err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
        }
    }

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) == 0);
    if(!err && check) {
        err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_2DS );
    }

    return err;
}
//;
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
a5err B_2DS_eval_B_dB(real B_dB[12], real r, real phi, real z, B_2DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0;
    real B_dB_temp[6];

    interperr += interp2Dcomp_eval_df(B_dB_temp, &Bdata->B_r, r, z);

    B_dB[0] = B_dB_temp[0];
    B_dB[1] = B_dB_temp[1];
    B_dB[2] = 0;
    B_dB[3] = B_dB_temp[2];

    interperr += interp2Dcomp_eval_df(B_dB_temp, &Bdata->B_phi, r, z);

    B_dB[4] = B_dB_temp[0];
    B_dB[5] = B_dB_temp[1];
    B_dB[6] = 0;
    B_dB[7] = B_dB_temp[2];

    interperr += interp2Dcomp_eval_df(B_dB_temp, &Bdata->B_z, r, z);

    B_dB[8] = B_dB_temp[0];
    B_dB[9] = B_dB_temp[1];
    B_dB[10] = 0;
    B_dB[11] = B_dB_temp[2];

    /* Test for B field interpolation error */
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
    }


    real psi_dpsi[6];

    if(!err) {
        interperr += interp2Dcomp_eval_df(psi_dpsi, &Bdata->psi, r, z);

        B_dB[0] = B_dB[0] - psi_dpsi[2]/r;
        B_dB[1] = B_dB[1] + psi_dpsi[2]/(r*r)-psi_dpsi[5]/r;
        B_dB[3] = B_dB[3] - psi_dpsi[4]/r;
        B_dB[8] = B_dB[8] + psi_dpsi[1]/r;
        B_dB[9] = B_dB[9] - psi_dpsi[1]/(r*r) + psi_dpsi[3]/r;
        B_dB[11] = B_dB[11] + psi_dpsi[5]/r;

        /* Test for psi interpolation error */
        if(interperr) {
            err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
        }
    }

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B_dB[0]*B_dB[0] + B_dB[4]*B_dB[4] + B_dB[8]*B_dB[8]) == 0);
    if(!err && check) {
        err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_2DS );
    }

    return err;
}
//;
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
real B_2DS_get_axis_r(B_2DS_data* Bdata) {
    return Bdata->axis_r;
}
//;
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
inline
__attribute__((always_inline))
real B_2DS_get_axis_z(B_2DS_data* Bdata) {
    return Bdata->axis_z;
}//;
DECLARE_TARGET_END
#endif
