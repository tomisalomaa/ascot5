/**
 * @file B_field.h
 * @brief Header file for B_field.c
 *
 * Contains a list declaring all B_field_types, and declaration of
 * B_field_offload_data and B_field_data structs.
 */
#ifndef B_FIELD_H
#define B_FIELD_H

#include "ascot5.h"
#include "error.h"
#include "Bfield/B_GS.h"
#include "Bfield/B_2DS.h"
#include "Bfield/B_3DS.h"
#include "Bfield/B_STS.h"
#include "Bfield/B_TC.h"

/**
 * @brief Magnetic field types
 *
 * Magnetic field types are used in the magnetic field interface (B_field.c) to
 * direct function calls to correct magnetic field instances. Each magnetic
 * field instance must have a corresponding type.
 */
typedef enum B_field_type {
    B_field_type_GS,  /**< Analytic magnetic field                          */
    B_field_type_2DS, /**< Spline-interpolated axisymmetric  magnetic field */
    B_field_type_3DS, /**< Spline-interpolated 3D magnetic field            */
    B_field_type_STS, /**< Spline-interpolated stellarator magnetic field   */
    B_field_type_TC   /**< Trivial Cartesian magnetic field                 */
} B_field_type;

/**
 * @brief Magnetic field offload data
 *
 * This struct holds data necessary for offloading. The struct is initialized in
 * B_field_init_offload().
 *
 * The intended usage is that only single offload data is used at the time, and
 * the type of the data is declared with the "type" field.
 */
typedef struct {
    B_field_type type;        /**< Magnetic field type wrapped by this struct */
    B_GS_offload_data BGS;    /**< GS field or NULL if not active             */
    B_2DS_offload_data B2DS;  /**< 2DS field or NULL if not active            */
    B_3DS_offload_data B3DS;  /**< 3DS field or NULL if not active            */
    B_STS_offload_data BSTS;  /**< STS field or NULL if not active            */
    B_TC_offload_data BTC;    /**< TC field or NULL if not active             */
    int offload_array_length; /**< Allocated offload array length             */
} B_field_offload_data;

/**
 * @brief Magnetic field simulation data
 *
 * This struct holds data necessary for simulation. The struct is initialized
 * from the B_field_offload_data in B_field_init().
 *
 * The intended usage is that only single B_field_data is used at the time, and
 * the type of the data is declared with the "type" field.
 */
typedef struct {
    B_field_type type; /**< Magnetic field type wrapped by this struct */
    B_GS_data BGS;     /**< GS field or NULL if not active             */
    B_2DS_data B2DS;   /**< 2DS field or NULL if not active            */
    B_3DS_data B3DS;   /**< 3DS field or NULL if not active            */
    B_STS_data BSTS;   /**< STS field or NULL if not active            */
    B_TC_data BTC;     /**< TC field or NULL if not active             */
} B_field_data;

int B_field_init_offload(B_field_offload_data* offload_data,
                         real** offload_array);
void B_field_free_offload(B_field_offload_data* offload_data,
                          real** offload_array);

DECLARE_TARGET
int B_field_init(
    B_field_data* Bdata, B_field_offload_data* offload_data,
    real* offload_array);
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
DECLARE_TARGET
inline
__attribute__((always_inline))
a5err B_field_eval_psi(
    real* psi, real r, real phi, real z, real t, B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
            err = B_GS_eval_psi(psi, r, phi, z, &(Bdata->BGS));
            break;

        case B_field_type_2DS:
            err = B_2DS_eval_psi(psi, r, phi, z, &(Bdata->B2DS));
            break;

        case B_field_type_3DS:
            err = B_3DS_eval_psi(psi, r, phi, z, &(Bdata->B3DS));
            break;

        case B_field_type_STS:
            err = B_STS_eval_psi(psi, r, phi, z, &(Bdata->BSTS));
            break;

        case B_field_type_TC:
            err = B_TC_eval_psi(psi, r, phi, z, &(Bdata->BTC));
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD );
            break;
    }

    if(err) {
        /* In case of error, return some reasonable value to avoid further
           complications */
        psi[0] = 1;
    }

    return err;
}//;
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
DECLARE_TARGET
a5err B_field_eval_psi_dpsi(
    real psi_dpsi[4], real r, real phi, real z, real t, B_field_data* Bdata);
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
DECLARE_TARGET
inline
__attribute__((always_inline))
a5err B_field_eval_rho(real* rho, real psi, B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
            err = B_GS_eval_rho(rho, psi, &(Bdata->BGS));
            break;

        case B_field_type_2DS:
            err = B_2DS_eval_rho(rho, psi, &(Bdata->B2DS));
            break;

        case B_field_type_3DS:
            err = B_3DS_eval_rho(rho, psi, &(Bdata->B3DS));
            break;

        case B_field_type_STS:
            err = B_STS_eval_rho(rho, psi, &(Bdata->BSTS));
            break;

        case B_field_type_TC:
            err = B_TC_eval_rho(rho, psi, &(Bdata->BTC));
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD );
            break;
    }

    if(err) {
        /* In case of error, return some reasonable value to avoid further
           complications */
        rho[0] = 1;
    }

    return err;
}//;
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
DECLARE_TARGET
a5err B_field_eval_rho_drho(
    real rho_drho[4], real r, real phi, real z, B_field_data* Bdata);
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
DECLARE_TARGET
inline
__attribute__((always_inline))
a5err B_field_eval_B(real B[3], real r, real phi, real z, real t,
                     B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
            err = B_GS_eval_B(B, r, phi, z, &(Bdata->BGS));
            break;

        case B_field_type_2DS:
            err = B_2DS_eval_B(B, r, phi, z, &(Bdata->B2DS));
            break;

        case B_field_type_3DS:
            err = B_3DS_eval_B(B, r, phi, z, &(Bdata->B3DS));
            break;

        case B_field_type_STS:
            err = B_STS_eval_B(B, r, phi, z, &(Bdata->BSTS));
            break;

        case B_field_type_TC:
            err = B_TC_eval_B(B, r, phi, z, &(Bdata->BTC));
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD );
            break;
    }

    if(err) {
        /* In case of error, return some reasonable values to avoid further
           complications */
        B[0] = 1;
        for(int k=1; k<3; k++) {B[k] = 0;}
    }

    return err;
}//;
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
DECLARE_TARGET
inline
__attribute__((always_inline))
a5err B_field_eval_B_dB(
    real B_dB[15], real r, real phi, real z, real t, B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
            err = B_GS_eval_B_dB(B_dB, r, phi, z, &(Bdata->BGS));
            break;

        case B_field_type_2DS:
            err = B_2DS_eval_B_dB(B_dB, r, phi, z, &(Bdata->B2DS));
            break;

        case B_field_type_3DS:
            err = B_3DS_eval_B_dB(B_dB, r, phi, z, &(Bdata->B3DS));
            break;

        case B_field_type_STS:
            err = B_STS_eval_B_dB(B_dB, r, phi, z, &(Bdata->BSTS));
            break;

        case B_field_type_TC:
            err = B_TC_eval_B_dB(B_dB, r, phi, z, &(Bdata->BTC));
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD );
            break;
    }

    if(err) {
        /* In case of error, return some reasonable values to avoid further
           complications */
        B_dB[0] = 1;
        for(int k=1; k<12; k++) {B_dB[k] = 0;}
    }

    return err;
}//;
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
DECLARE_TARGET
inline
__attribute__((always_inline))
real B_field_get_axis_r(B_field_data* Bdata, real phi) {
    a5err err = 0;
    real axis_r = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
            axis_r = B_GS_get_axis_r(&(Bdata->BGS));
            break;

        case B_field_type_2DS:
            axis_r = B_2DS_get_axis_r(&(Bdata->B2DS));
            break;

        case B_field_type_3DS:
            axis_r = B_3DS_get_axis_r(&(Bdata->B3DS));
            break;

        case B_field_type_STS:
            err = B_STS_get_axis_r(&axis_r, &(Bdata->BSTS), phi);
            if(err) {
                /* In case of error, return some reasonable value to avoid
                   further complications */
                axis_r = 5;
            }
            break;

        case B_field_type_TC:
            axis_r = B_TC_get_axis_r(&(Bdata->BTC));
            break;
    }

    return axis_r;
}//;
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
DECLARE_TARGET
inline
__attribute__((always_inline))
real B_field_get_axis_z(B_field_data* Bdata, real phi) {
    a5err err = 0;
    real axis_z = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
            axis_z = B_GS_get_axis_z(&(Bdata->BGS));
            break;

        case B_field_type_2DS:
            axis_z = B_2DS_get_axis_z(&(Bdata->B2DS));
            break;

        case B_field_type_3DS:
            axis_z = B_3DS_get_axis_z(&(Bdata->B3DS));
            break;

        case B_field_type_STS:
            err = B_STS_get_axis_z(&axis_z, &(Bdata->BSTS), phi);
            if(err) {
                /* In case of error, return some reasonable value to avoid
                   further complications */
                axis_z = 0;
            }
            break;

        case B_field_type_TC:
            axis_z = B_TC_get_axis_z(&(Bdata->BTC));
            break;
    }

    return axis_z;
}
//;
DECLARE_TARGET_END

#endif
