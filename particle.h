/**
 * @file particle.h
 * @brief Header file for particle.c
 *
 * The relationship between the seven different marker structs is:
 *
 *    particle >--+            +--> particle_simd_fo
 *                |            |
 * particle_gc >--particle_state--> particle_simd_gc
 *                |            |
 * particle_ml >--+            +--> particle_simd_ml
 *
 * which is explained in particle.c. This file contains the definitions
 * of these structs as well as input_particle and particle_queue structs. Former
 * is a wrapper for particle, particle_gc, particle_ml, and particle_state
 * structs while latter is a queue from where markers are chosen when simulation
 * begins and updated when simulation ends.
 */
#ifndef PARTICLE_H
#define PARTICLE_H

#include "ascot5.h"
#include "B_field.h"
#include "E_field.h"
#include "error.h"
//#include "consts.h"

/**
 * @brief General representation of a marker
 *
 * This struct is a self-consistent representation of a marker which can be
 * constructed from any input and from which any simulation marker struct can
 * be constructed. This structure is not intended to be used during the
 * simulation, but before and after for storing marker data, and making it easy
 * to switch between different marker types and construct markers from inputs.
 *
 * Whenever marker properties are changed, it is the responsibility of the
 * function performing the change to make sure all parameters remain consistent.
 * For example, magnetic field must be updated when marker position changes.
 */
typedef struct {
    real r;           /**< Guiding center R coordinate [m]                 */
    real phi;         /**< Guiding center phi coordinate [rad]             */
    real z;           /**< Guiding center z coordinate [m]                 */
    real vpar;        /**< Parallel velocity [m/s]                         */
    real mu;          /**< Magnetic moment [J/T]                           */
    real zeta;        /**< Gyroangle [rad]                                 */
    real rprt;        /**< Particle R coordinate [m]                       */
    real phiprt;      /**< Particle phi coordinate [phi]                   */
    real zprt;        /**< Particle z coordinate [m]                       */
    real rdot;        /**< dr/dt [m/s]                                     */
    real phidot;      /**< dphi/dt [rad/s]                                 */
    real zdot;        /**< dz/dt [m/s]                                     */
    real mass;        /**< Mass [kg]                                       */
    real charge;      /**< Charge [C]                                      */
    int  anum;        /**< Atomic mass number of marker species            */
    int  znum;        /**< Charge number of marker species                 */
    real weight;      /**< Marker weight                                   */
    real time;        /**< Marker simulation time [s]                      */
    real cputime;     /**< Marker wall-clock time [s]                      */
    real rho;         /**< Marker rho coordinate                           */
    real theta;       /**< Marker poloidal coordinate [rad]                */
    integer id;       /**< Arbitrary but unique ID for the marker          */
    integer endcond;  /**< Marker end condition                            */
    integer walltile; /**< ID of walltile if marker has hit the wall       */
    real B_r;         /**< Magnetic field R component at (r, phi, z) [T]   */
    real B_phi;       /**< Magnetic field phi component at (r, phi, z) [T] */
    real B_z;         /**< Magnetic field z component at (r, phi, z) [T]   */
    real B_r_dr;      /**< dB_R/dR at (r, phi, z) [T/m]                    */
    real B_phi_dr;    /**< dB_phi/dR at (r, phi, z) [T/m]                  */
    real B_z_dr;      /**< dB_z/dR at (r, phi, z) [T/m]                    */
    real B_r_dphi;    /**< dB_R/dphi at (r, phi, z) [T/m]                  */
    real B_phi_dphi;  /**< dB_phi/dphi at (r, phi, z) [T/m]                */
    real B_z_dphi;    /**< dB_z/dphi at (r, phi, z) [T/m]                  */
    real B_r_dz;      /**< dB_R/dz at (r, phi, z) [T/m]                    */
    real B_phi_dz;    /**< dB_phi/dz at (r, phi, z) [T/m]                  */
    real B_z_dz;      /**< dB_z/dz at (r, phi, z) [T/m]                    */

    a5err err;        /**< error flag */
} particle_state;

/**
 * @brief Particle input
 *
 * When particle marker data is read, this struct is created and filled. This
 * struct is then converted to a particle_state struct.
 */
typedef struct {
    real r;      /**< R coordinate [m]                    */
    real phi;    /**< phi coordinate [rad]                */
    real z;      /**< z coordinate [m]                    */
    real v_r;    /**< Velocity R-component [m/s]          */
    real v_phi;  /**< Velocity phi-component [m/s]        */
    real v_z;    /**< Velocity z-component [m/s]          */
    real mass;   /**< Mass [kg]                           */
    real charge; /**< Charge [C]                          */
    int  anum;   /**< Atomic mass number [1]              */
    int  znum;   /**< Charge number [1]                   */
    real weight; /**< Particle marker weight              */
    real time;   /**< Particle marker simulation time [s] */
    integer id;  /**< Unique ID for the particle marker   */
} particle;

/**
 * @brief Guiding center input
 *
 * When guiding center marker data is read, this struct is created and filled.
 * This struct is then converted to a particle_state struct.
 */
typedef struct {
    real r;      /**< R coordinate [m]                          */
    real phi;    /**< phi coordinate [rad]                      */
    real z;      /**< z coordinate [m]                          */
    real energy; /**< Kinetic energy [J]                        */
    real pitch;  /**< Pitch                                     */
    real zeta;   /**< Gyroangle [rad]                           */
    real mass;   /**< Mass [kg]                                 */
    real charge; /**< Charge [C]                                */
    int  anum;   /**< Atomic mass number [1]                    */
    int  znum;   /**< Charge number [1]                         */
    real weight; /**< Guiding center marker weight              */
    real time;   /**< Guiding center marker simulation time [s] */
    integer id;  /**< Unique ID for the guiding center marker   */
} particle_gc;

/**
 * @brief Field line input
 *
 * When field line marker data is read, this struct is created and filled. This
 * struct is then converted to a particle_state struct.
 */
typedef struct {
    real r;      /**< R coordinate [m]                      */
    real phi;    /**< phi coordinate [rad]                  */
    real z;      /**< z coordinate [m]                      */
    real pitch;  /**< Direction                             */
    real weight; /**< Field line marker weight              */
    real time;   /**< Field line marker simulation time [s] */
    integer id;  /**< Unique ID for the field line marker   */
} particle_ml;

/**
 * @brief Marker queue
 *
 * Each time a marker has finished simulation, a new marker is chosen from this
 * queue and the old marker's data is updated. Markers are never removed from
 * the queue but an index is kept to mark where the next not yet simulated
 * marker is found. Markers are represented by particle_state struct when they
 * are stored in the queue.
 *
 * Note: The queue can and is accessed by several threads, so make sure each
 * access is thread-safe.
 */
typedef struct {
    int n;                 /**< Total number of markers in this queue        */
    particle_state** p;    /**< Pointer to an array storing pointers to all
                                markers within this queue.                   */
    volatile int next;     /**< Index where next unsimulated marker is found */
    volatile int finished; /**< Number of markers who have finished
                                simulation                                   */
} particle_queue;

/**
 * @brief Marker types enum
 *
 * Used to indicate what marker type is stored in input_particle wrapper.
 */
typedef enum input_particle_type {
    input_particle_type_p,  /**< Type corresponding to particle struct       */
    input_particle_type_gc, /**< Type corresponding to particle_gc struct    */
    input_particle_type_ml, /**< Type corresponding to particle_ml struct    */
    input_particle_type_s   /**< Type corresponding to particle_state struct */
} input_particle_type;

/**
 * @brief Wrapper for marker structs
 *
 * This struct wraps particle_state struct and all input structs. Reason is
 * because input data can have several marker types and with this wrapper only
 * a single array is required to store them. The same array can be used when
 * the input markers are turned into marker states.
 *
 * Only a single type is stored here indicated by the "type" field. Storing a
 * a new marker struct removes the old marker struct.
 */
typedef struct {
    input_particle_type type; /**< Type of data currently stored */
    union {
        particle p;           /**< Particle input                */
        particle_gc p_gc;     /**< Guiding center input          */
        particle_ml p_ml;     /**< Field line tracer input       */
        particle_state p_s;   /**< Marker state                  */
    };
} input_particle;

/**
 * @brief Struct representing NSIMD particle markers
 *
 * This struct is used in simulation when the simulation loop in
 * simulate_fo_fixed.c is used.
 *
 * It contains physical and simulation parameters necessary for the simulation
 * If a function makes changes to any of these parameters, it is that function's
 * responsibility to make sure all fields remain consistent, i.e., if position
 * changes then the magnetic field should be updated. Each field is a memory
 * aligned array with length NSIMD, so this struct represents NSIMD markers
 * (they can be dummy or markers whose simulation has been terminated) and so it
 * can be used within SIMD loops.
 *
 * The fields are aligned to 64 bit with __memalign__ (see ascot5.h).
 */
typedef struct {
    /* Physical coordinates and parameters */
    real r[NSIMD] __memalign__;       /**< Particle R coordinate [m]          */
    real phi[NSIMD] __memalign__;     /**< Particle phi coordinate [phi]      */
    real z[NSIMD] __memalign__;       /**< Particle z coordinate [m]          */
    real rdot[NSIMD] __memalign__;    /**< dr/dt [m/s]                        */
    real phidot[NSIMD] __memalign__;  /**< dphi/dt [rad/s]                    */
    real zdot[NSIMD] __memalign__;    /**< dz/dt [m/s]                        */
    real mass[NSIMD] __memalign__;    /**< Mass [kg]                          */
    real charge[NSIMD] __memalign__;  /**< Charge [C]                         */
    real time[NSIMD] __memalign__;    /**< Marker simulation time [s]         */

    /* Magnetic field data */
    real B_r[NSIMD] __memalign__;        /**< Magnetic field R component at
                                              marker position [T]             */
    real B_phi[NSIMD] __memalign__;      /**< Magnetic field phi component at
                                              marker position [T]             */
    real B_z[NSIMD] __memalign__;        /**< Magnetic field z component at
                                              marker position [T]             */

    real B_r_dr[NSIMD] __memalign__;     /**< dB_R/dR at marker position [T/m]     */
    real B_phi_dr[NSIMD] __memalign__;   /**< dB_phi/dR at marker position [T/m]   */
    real B_z_dr[NSIMD] __memalign__;     /**< dB_z/dR at marker position [T/m]     */
    real B_r_dphi[NSIMD] __memalign__;   /**< dB_R/dphi at marker position [T/m]   */
    real B_phi_dphi[NSIMD] __memalign__; /**< dB_phi/dphi at marker position [T/m] */
    real B_z_dphi[NSIMD] __memalign__;   /**< dB_z/dphi at marker position [T/m]   */
    real B_r_dz[NSIMD] __memalign__;     /**< dB_R/dz at marker position [T/m]     */
    real B_phi_dz[NSIMD] __memalign__;   /**< dB_phi/dz at marker position [T/m]   */
    real B_z_dz[NSIMD] __memalign__;     /**< dB_z/dz at marker position [T/m]     */

    /* Quantities used in diagnostics */
    real weight[NSIMD] __memalign__;  /**< Marker weight                      */
    real cputime[NSIMD] __memalign__; /**< Marker wall-clock time [s]         */
    real rho[NSIMD] __memalign__;     /**< Marker rho coordinate              */
    real theta[NSIMD] __memalign__;   /**< Marker poloidal coordinate [rad]   */

    integer id[NSIMD] __memalign__;       /**< Unique ID for the marker       */
    integer endcond[NSIMD] __memalign__;  /**< Marker end condition           */
    integer walltile[NSIMD] __memalign__; /**< ID of walltile if marker has
                                               hit the wall                   */

    /* Meta data */
    integer running[NSIMD] __memalign__; /**< Indicates whether this marker is
                                              currently simulated (1) or not  */
    a5err err[NSIMD] __memalign__;       /**< Error flag, zero if no error    */
    integer index[NSIMD] __memalign__;   /**< Marker index at marker queue    */
} particle_simd_fo;

/**
 * @brief Struct representing NSIMD guiding center markers
 *
 * This struct is used in simulation when the simulation loop in
 * simulate_gc_fixed.c or simulate_gc_adaptive.c is used.
 *
 * It contains physical and simulation parameters necessary for the simulation
 * If a function makes changes to any of these parameters, it is that function's
 * responsibility to make sure all fields remain consistent, i.e., if position
 * changes then the magnetic field should be updated. Each field is a memory
 * aligned array with length NSIMD, so this struct represents NSIMD markers
 * (they can be dummy or markers whose simulation has been terminated) and so it
 * can be used within SIMD loops.
 *
 * The fields are aligned to 64 bit with __memalign__ (see ascot5.h).
 */
typedef struct {
    /* Physical coordinates and parameters */
    real r[NSIMD] __memalign__;      /**< Guiding center R coordinate [m]     */
    real phi[NSIMD] __memalign__;    /**< Guiding center phi coordinate [phi] */
    real z[NSIMD] __memalign__;      /**< Guiding center z coordinate [m]     */
    real vpar[NSIMD] __memalign__;   /**< Parallel velocity [m/s]             */
    real mu[NSIMD] __memalign__;     /**< Magnetic moment [J/T]               */
    real zeta[NSIMD] __memalign__;   /**< Gyroangle [rad]                     */
    real mass[NSIMD] __memalign__;   /**< Mass [kg]                           */
    real charge[NSIMD] __memalign__; /**< Charge [C]                          */
    real time[NSIMD] __memalign__;   /**< Marker simulation time [s]          */

    /* Magnetic field data */
    real B_r[NSIMD] __memalign__;        /**< Magnetic field R component at
                                              marker position [T]             */
    real B_phi[NSIMD] __memalign__;      /**< Magnetic field phi component at
                                              marker position [T]             */
    real B_z[NSIMD] __memalign__;        /**< Magnetic field z component at
                                              marker position [T]             */

    real B_r_dr[NSIMD] __memalign__;     /**< dB_R/dR at marker position [T/m]     */
    real B_phi_dr[NSIMD] __memalign__;   /**< dB_phi/dR at marker position [T/m]   */
    real B_z_dr[NSIMD] __memalign__;     /**< dB_z/dR at marker position [T/m]     */
    real B_r_dphi[NSIMD] __memalign__;   /**< dB_R/dphi at marker position [T/m]   */
    real B_phi_dphi[NSIMD] __memalign__; /**< dB_phi/dphi at marker position [T/m] */
    real B_z_dphi[NSIMD] __memalign__;   /**< dB_z/dphi at marker position [T/m]   */
    real B_r_dz[NSIMD] __memalign__;     /**< dB_R/dz at marker position [T/m]     */
    real B_phi_dz[NSIMD] __memalign__;   /**< dB_phi/dz at marker position [T/m]   */
    real B_z_dz[NSIMD] __memalign__;     /**< dB_z/dz at marker position [T/m]     */

    /* Quantities used in diagnostics */
    real weight[NSIMD] __memalign__;  /**< Marker weight                      */
    real cputime[NSIMD] __memalign__; /**< Marker wall-clock time [s]         */
    real rho[NSIMD] __memalign__;     /**< Marker rho coordinate              */
    real theta[NSIMD] __memalign__;   /**< Marker poloidal coordinate [rad]   */

    integer id[NSIMD] __memalign__;       /**< Unique ID for the marker       */
    integer endcond[NSIMD] __memalign__;  /**< Marker end condition           */
    integer walltile[NSIMD] __memalign__; /**< ID of walltile if marker has
                                               hit the wall                   */

    /* Meta data */
    integer running[NSIMD] __memalign__; /**< Indicates whether this marker is
                                              currently simulated (1) or not  */
    a5err err[NSIMD] __memalign__;       /**< Error flag, zero if no error    */
    integer index[NSIMD] __memalign__;   /**< Marker index at marker queue    */
} particle_simd_gc;

/**
 * @brief Struct representing NSIMD field line markers
 *
 * This struct is used in simulation when the simulation loop in
 * simulate_ml_adaptive.c is used.
 *
 * It contains physical and simulation parameters necessary for the simulation
 * If a function makes changes to any of these parameters, it is that function's
 * responsibility to make sure all fields remain consistent, i.e., if position
 * changes then the magnetic field should be updated. Each field is a memory
 * aligned array with length NSIMD, so this struct represents NSIMD markers
 * (they can be dummy or markers whose simulation has been terminated) and so it
 * can be used within SIMD loops.
 *
 * The fields are aligned to 64 bit with __memalign__ (see ascot5.h).
 */
typedef struct {
    /* Physical coordinates and parameters */
    real r[NSIMD] __memalign__;     /**< Field line R coordinate [m]          */
    real phi[NSIMD] __memalign__;   /**< Field line phi coordinate [phi]      */
    real z[NSIMD] __memalign__;     /**< Field line z coordinate [m]          */
    real pitch[NSIMD] __memalign__; /**< Field line direction: along (1) or
                                         against (-1) magnetic field vector   */
    real time[NSIMD] __memalign__;  /**< Field line simulation "time" i.e.
                                         (distance / speed of light) [m]      */

    /* Magnetic field data */
    real B_r[NSIMD] __memalign__;        /**< Magnetic field R component at
                                              marker position [T]             */
    real B_phi[NSIMD] __memalign__;      /**< Magnetic field phi component at
                                              marker position [T]             */
    real B_z[NSIMD] __memalign__;        /**< Magnetic field z component at
                                              marker position [T]             */

    real B_r_dr[NSIMD] __memalign__;     /**< dB_R/dR at marker position [T/m]     */
    real B_phi_dr[NSIMD] __memalign__;   /**< dB_phi/dR at marker position [T/m]   */
    real B_z_dr[NSIMD] __memalign__;     /**< dB_z/dR at marker position [T/m]     */
    real B_r_dphi[NSIMD] __memalign__;   /**< dB_R/dphi at marker position [T/m]   */
    real B_phi_dphi[NSIMD] __memalign__; /**< dB_phi/dphi at marker position [T/m] */
    real B_z_dphi[NSIMD] __memalign__;   /**< dB_z/dphi at marker position [T/m]   */
    real B_r_dz[NSIMD] __memalign__;     /**< dB_R/dz at marker position [T/m]     */
    real B_phi_dz[NSIMD] __memalign__;   /**< dB_phi/dz at marker position [T/m]   */
    real B_z_dz[NSIMD] __memalign__;     /**< dB_z/dz at marker position [T/m]     */

    /* Quantities used in diagnostics */
    real weight[NSIMD] __memalign__;  /**< Marker weight                      */
    real cputime[NSIMD] __memalign__; /**< Marker wall-clock time [s]         */
    real rho[NSIMD] __memalign__;     /**< Marker rho coordinate              */
    real theta[NSIMD] __memalign__;   /**< Marker poloidal coordinate [rad]   */

    integer id[NSIMD] __memalign__;       /**< Unique ID for the marker       */
    integer endcond[NSIMD] __memalign__;  /**< Marker end condition           */
    integer walltile[NSIMD] __memalign__; /**< ID of walltile if marker has
                                               hit the wall                   */

    /* Meta data */
    integer running[NSIMD] __memalign__; /**< Indicates whether this marker is
                                              currently simulated (1) or not  */
    a5err err[NSIMD] __memalign__;       /**< Error flag, zero if no error    */
    integer index[NSIMD] __memalign__;   /**< Marker index at marker queue    */
} particle_simd_ml;




DECLARE_TARGET
void particle_to_fo_dummy(particle_simd_fo* p_fo, int j);
DECLARE_TARGET_END
DECLARE_TARGET
void particle_to_gc_dummy(particle_simd_gc* p_gc, int j);
DECLARE_TARGET_END
DECLARE_TARGET
void particle_to_ml_dummy(particle_simd_ml* p_ml, int j);
DECLARE_TARGET_END

DECLARE_TARGET
int particle_cycle_fo(particle_queue* q, particle_simd_fo* p,
                      B_field_data* Bdata, int* cycle);
DECLARE_TARGET_END
DECLARE_TARGET
int particle_cycle_gc(particle_queue* q, particle_simd_gc* p,
                      B_field_data* Bdata, int* cycle);
DECLARE_TARGET_END
DECLARE_TARGET
int particle_cycle_ml(particle_queue* q, particle_simd_ml* p,
                      B_field_data* Bdata, int* cycle);

DECLARE_TARGET_END
DECLARE_TARGET
void particle_input_to_state(input_particle* p, particle_state* ps,
                             B_field_data* Bdata);
DECLARE_TARGET_END

#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
DECLARE_TARGET
inline
__attribute__((always_inline))
a5err particle_state_to_fo(particle_state* p, int i, particle_simd_fo* p_fo,
                           int j, B_field_data* Bdata)
{
    a5err err = p->err;

    if(!err) {
        p_fo->r[j]          = p->rprt;
        p_fo->phi[j]        = p->phiprt;
        p_fo->z[j]          = p->zprt;
        p_fo->rdot[j]       = p->rdot;
        p_fo->phidot[j]     = p->phidot;
        p_fo->zdot[j]       = p->zdot;

        p_fo->mass[j]       = p->mass;
        p_fo->charge[j]     = p->charge;
        p_fo->weight[j]     = p->weight;
        p_fo->time[j]       = p->time;
        p_fo->theta[j]      = p->theta;
        p_fo->id[j]         = p->id;
        p_fo->endcond[j]    = p->endcond;
        p_fo->walltile[j]   = p->walltile;
    }

    /* Magnetic field stored in state is for the gc position */
    real B_dB[15], psi[1], rho[1];
    if(!err) {
        err = B_field_eval_B_dB(B_dB, p->rprt, p->phiprt, p->zprt, p->time,
                                Bdata);
    }
    if(!err) {
        err = B_field_eval_psi(psi, p->rprt, p->phiprt, p->zprt, p->time,
                               Bdata);
    }
    if(!err) {
        err = B_field_eval_rho(rho, psi[0], Bdata);
    }

    if(!err) {
        p_fo->rho[j]        = rho[0];

        p_fo->B_r[j]        = B_dB[0];
        p_fo->B_r_dr[j]     = B_dB[1];
        p_fo->B_r_dphi[j]   = B_dB[2];
        p_fo->B_r_dz[j]     = B_dB[3];

        p_fo->B_phi[j]      = B_dB[4];
        p_fo->B_phi_dr[j]   = B_dB[5];
        p_fo->B_phi_dphi[j] = B_dB[6];
        p_fo->B_phi_dz[j]   = B_dB[7];

        p_fo->B_z[j]        = B_dB[8];
        p_fo->B_z_dr[j]     = B_dB[9];
        p_fo->B_z_dphi[j]   = B_dB[10];
        p_fo->B_z_dz[j]     = B_dB[11];

        p_fo->running[j] = 1;
        if(p->endcond) {
            p_fo->running[j] = 0;
        }
        p_fo->cputime[j] = p->cputime;
        p_fo->index[j]   = i;

        p_fo->err[j] = 0;
    }

    return err;
}
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
DECLARE_TARGET
void particle_fo_to_state(particle_simd_fo* p_fo, int j, particle_state* p,
                          B_field_data* Bdata);
DECLARE_TARGET_END
/* inline */
/* void particle_fo_to_state_inline(particle_simd_fo* p_fo, int j, particle_state* p, */
/*                           B_field_data* Bdata) */
/* { */
/*     a5err err = 0; */

/*     p->rprt       = p_fo->r[j]; */
/*     p->phiprt     = p_fo->phi[j]; */
/*     p->zprt       = p_fo->z[j]; */
/*     p->rdot       = p_fo->rdot[j]; */
/*     p->phidot     = p_fo->phidot[j]; */
/*     p->zdot       = p_fo->zdot[j]; */

/*     p->mass       = p_fo->mass[j]; */
/*     p->charge     = p_fo->charge[j]; */
/*     p->weight     = p_fo->weight[j]; */
/*     p->time       = p_fo->time[j]; */
/*     p->theta      = p_fo->theta[j]; */
/*     p->id         = p_fo->id[j]; */
/*     p->endcond    = p_fo->endcond[j]; */
/*     p->walltile   = p_fo->walltile[j]; */
/*     p->cputime    = p_fo->cputime[j]; */

/*     /\* Particle to guiding center *\/ */
/*     real B_dB[15], psi[1], rho[1]; */
/*     rho[0]        = p_fo->rho[j]; */
/*     B_dB[0]       = p_fo->B_r[j]; */
/*     B_dB[1]       = p_fo->B_r_dr[j]; */
/*     B_dB[2]       = p_fo->B_r_dphi[j]; */
/*     B_dB[3]       = p_fo->B_r_dz[j]; */
/*     B_dB[4]       = p_fo->B_phi[j]; */
/*     B_dB[5]       = p_fo->B_phi_dr[j]; */
/*     B_dB[6]       = p_fo->B_phi_dphi[j]; */
/*     B_dB[7]       = p_fo->B_phi_dz[j]; */
/*     B_dB[8]       = p_fo->B_z[j]; */
/*     B_dB[9]       = p_fo->B_z_dr[j]; */
/*     B_dB[10]      = p_fo->B_z_dphi[j]; */
/*     B_dB[11]      = p_fo->B_z_dz[j]; */

/*     /\* Guiding center transformation *\/ */
/*     real vpar; */
/*     if(!err) { */
/*         real vR   = p->rdot; */
/*         real vphi = p->phidot * p->rprt; */
/*         real vz   = p->zdot; */

/*         gctransform_particle2guidingcenter( */
/*             p->mass, p->charge, B_dB, */
/*             p->rprt, p->phiprt, p->zprt, vR , vphi, vz, */
/*             &p->r, &p->phi, &p->z, &vpar, &p->mu, &p->zeta); */
/*     } */
/*     if(!err && p->r <= 0)  {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);} */
/*     if(!err && p->mu < 0)  {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);} */

/*     if(!err) { */
/*         err = B_field_eval_B_dB(B_dB, p->r, p->phi, p->z, p->time, Bdata); */
/*     } */
/*     if(!err) { */
/*         err = B_field_eval_psi(psi, p->r, p->phi, p->z, p->time, Bdata); */
/*     } */
/*     if(!err) { */
/*         err = B_field_eval_rho(rho, psi[0], Bdata); */
/*     } */

/*     if(!err) { */
/*         p->vpar = vpar; */
/*     } */
/*     if(!err && p->vpar >= CONST_C) { */
/*         err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE); */
/*     } */

/*     /\* Normally magnetic field data at gc position is stored here */
/*      * but, if gc transformation fails, field at particle position is */
/*      * stored instead. *\/ */
/*     p->rho        = rho[0]; */

/*     p->B_r        = B_dB[0]; */
/*     p->B_r_dr     = B_dB[1]; */
/*     p->B_r_dphi   = B_dB[2]; */
/*     p->B_r_dz     = B_dB[3]; */

/*     p->B_phi      = B_dB[4]; */
/*     p->B_phi_dr   = B_dB[5]; */
/*     p->B_phi_dphi = B_dB[6]; */
/*     p->B_phi_dz   = B_dB[7]; */

/*     p->B_z        = B_dB[8]; */
/*     p->B_z_dr     = B_dB[9]; */
/*     p->B_z_dphi   = B_dB[10]; */
/*     p->B_z_dz     = B_dB[11]; */

/*     /\* If marker already has error flag, make sure it is not overwritten here *\/ */
/*     a5err simerr  = p_fo->err[j]; */
/*     if(simerr) { */
/*         p->err = simerr; */
/*     } */
/*     else { */
/*         p->err = err; */
/*     } */
/* }//; */
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
DECLARE_TARGET
a5err particle_state_to_gc(particle_state* p, int i, particle_simd_gc* p_gc,
                           int j, B_field_data* Bdata);
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
DECLARE_TARGET
void particle_gc_to_state(particle_simd_gc* p_gc, int j, particle_state* p,
                          B_field_data* Bdata);
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
DECLARE_TARGET
a5err particle_state_to_ml(particle_state* p, int i, particle_simd_ml* p_ml,
                           int j, B_field_data* Bdata);
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd uniform(Bdata)
#endif
DECLARE_TARGET
void particle_ml_to_state(particle_simd_ml* p_ml, int j, particle_state* p,
                          B_field_data* Bdata);
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd uniform(p_fo,Bdata)
#endif
DECLARE_TARGET
int particle_fo_to_gc(particle_simd_fo* p_fo, int j, particle_simd_gc* p_gc,
                      B_field_data* Bdata);
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd
#endif
DECLARE_TARGET
inline
__attribute__((always_inline))
void particle_copy_fo(particle_simd_fo* p1, int i, particle_simd_fo* p2, int j)
{
        p2->r[j]          = p1->r[i];
        p2->phi[j]        = p1->phi[i];
        p2->z[j]          = p1->z[i];
        p2->rdot[j]       = p1->rdot[i];
        p2->phidot[j]     = p1->phidot[i];
        p2->zdot[j]       = p1->zdot[i];

        p2->time[j]       = p1->time[i];
        p2->cputime[j]    = p1->cputime[i];
        p2->rho[j]        = p1->rho[i];
        p2->weight[j]     = p1->weight[i];
        p2->cputime[j]    = p1->cputime[i];
        p2->rho[j]        = p1->rho[i];
        p2->theta[j]      = p1->theta[i];

        p2->mass[j]       = p1->mass[i];
        p2->charge[j]     = p1->charge[i];

        p2->id[j]         = p1->id[i];
        p2->running[j]    = p1->running[i];
        p2->endcond[j]    = p1->endcond[i];
        p2->walltile[j]   = p1->walltile[i];

        p2->B_r[j]        = p1->B_r[i];
        p2->B_phi[j]      = p1->B_phi[i];
        p2->B_z[j]        = p1->B_z[i];

        p2->B_r_dr[j]     = p1->B_r_dr[i];
        p2->B_r_dphi[j]   = p1->B_r_dphi[i];
        p2->B_r_dz[j]     = p1->B_r_dz[i];

        p2->B_phi_dr[j]   = p1->B_phi_dr[i];
        p2->B_phi_dphi[j] = p1->B_phi_dphi[i];
        p2->B_phi_dz[j]   = p1->B_phi_dz[i];

        p2->B_z_dr[j]     = p1->B_z_dr[i];
        p2->B_z_dphi[j]   = p1->B_z_dphi[i];
        p2->B_z_dz[j]     = p1->B_z_dz[i];
}
DECLARE_TARGET_END

#ifdef SIMD
#pragma omp declare simd
#endif
DECLARE_TARGET
void particle_copy_gc(particle_simd_gc* p1, int i, particle_simd_gc* p2, int j);
DECLARE_TARGET_END
#ifdef SIMD
#pragma omp declare simd
#endif
DECLARE_TARGET
void particle_copy_ml(particle_simd_ml* p1, int i, particle_simd_ml* p2, int j);
DECLARE_TARGET_END

#endif
