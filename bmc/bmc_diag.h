#ifndef BMCDIAG_H
#define BMCDIAG_H

#include "../diag.h"
#include "../diag/dist_6D.h"
#include "../diag/dist_5D.h"
#include "../simulate.h"
#include "../endcond.h"
#include "../consts.h"

#ifdef MPI
#include "mpi.h"
#endif

#include "bmc_wall.h"
#include <string.h>
#include "../print.h"
#include "../math.h"
#include "../physlib.h"

typedef struct {
    int index[NSIMD * 32] __memalign__;
    int target_hit[NSIMD * 32] __memalign__;
    real weight[NSIMD * 32] __memalign__;
} particle_deposit_weights;

void particle_state_to_particle_gc(
    particle_state* ps,
    particle_gc* p_gc
);

void diag_move_distribution(sim_offload_data* sim, diag_data* diag_dest, diag_data* diag_src, int* updated, int* nloss, int* n_err);

int bmc_dist5D_gc_indexes_from_coordinates(int* indexes, real* weights, int* target_hit, real r, real phi, real z, real ppar, real pperp, dist_5D_data* dist, wall_data* wallData);

real fmc_compute_signal_from_states(
    int n_particles,
    particle_state* ps,
    int* n_updated,
    int* n_loss,
    int* n_err
);

int bmc_update_distr5D(
        dist_5D_data* dist1,
        dist_5D_data* dist0,
        int* p0_index,
        particle_simd_gc* p1,
        particle_simd_gc* p0,
        int n_simd_particles,
        wall_data* wallData 
    );

int fmc_update_distr5D_from_states(
        dist_5D_data* dist1,
        int* p0_indexes,
        particle_state* p1,
        int n_particles,
        wall_2d_data* w2d,
        int *n_loss,
        int *n_err
    );


int bmc_dist5D_state_indexes(int* indexes, real* weights, int* target_hit, particle_state* ps, dist_5D_data* dist, wall_data* wdata);
int bmc_dist5D_gc_indexes(int* indexes, real* weights, int* target_hit, particle_simd_gc* p, int i, dist_5D_data* dist, wall_data* wallData);
int bmc_dist6D_fo_index(particle_state* ps, dist_6D_data* dist);

void compute_5d_indexes_from_hist_index(int i, int* i_x, dist_5D_offload_data* dist5D);
void compute_5d_coordinates_from_hist_index(int i, int* i_x, real* r, real* phi, real* z, real* ppara, real* pperp, dist_5D_offload_data* dist5D);
void compute_element_5d_coordinates(int* i_x_new, real* r, real* phi, real* z, real* ppara, real* pperp, dist_5D_offload_data* dist);
void bmc_5D_to_particle_state(
        B_field_data* Bdata,
        real r, real phi, real z,
        real ppara, real pperp,
        real t,
        int id,
        particle_state* ps,
        real m,
        real q,
        int rk4_subcycles
    );

void write_probability_distribution(
    sim_offload_data* sim_offload,
    diag_data* distr,
    real* distr_array,
    int mpi_rank,
    int write_for_importance_sampling
);

void bmc_update_distr5D_from_weights(
        particle_deposit_weights *p1_weights,
        dist_5D_data* dist1,
        dist_5D_data* dist0,
        particle_simd_gc* p1,
        int n_simd_particles,
        int* p0_indexes
    );

void bmc_compute_prob_weights(particle_deposit_weights *p1_weightsIndexes,
                              int n_simd_particles, particle_simd_gc *p1, particle_simd_gc *p0,
                              dist_5D_data *dist1, dist_5D_data *dist0, wall_data *wallData,
                              int* p0_indexes, int debugExitVelocitySpace);

#endif
