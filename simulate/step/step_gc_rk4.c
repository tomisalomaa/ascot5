/**
 * @file step_gc_rk4.c
 * @brief Calculate a guiding center step for a struct of particles with RK4
 **/
#include <math.h>
#include "../../ascot5.h"
#include "../../consts.h"
#include "step_gc_rk4.h"
#include "../../phys_orbit.h"
#include "../../B_field.h"
#include "../../E_field.h"
#include "../../math.h"
#include "../../particle.h"

/**
 * @brief Integrate a guiding center step for a struct of particles with RK4
 *
 * This function calculates a guiding center step for a struct of NSIMD 
 * particles with RK4 simultaneously using SIMD instructions. All arrays in the 
 * function are of NSIMD length so vectorization can be performed directly 
 * without gather and scatter operations.
 *
 * @param p simd_gc struct that will be updated
 * @param h pointer to array containing time steps
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 */
void step_gc_rk4(particle_simd_gc* p, real* h, B_field_data* Bdata, E_field_data* Edata) {

    int i;
    /* Following loop will be executed simultaneously for all i */
    #pragma omp simd 
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            real k1[5];
            real k2[5];
            real k3[5];
            real k4[5];
            real tempy[5];
            real yprev[5];
            real y[5];

            real mass;
            real charge;

            real B[3];
            real B_dB[12];
	    real E[3];

	    real R0   = p->r[i];
	    real z0   = p->z[i];

            /* Coordinates are copied from the struct into an array to make 
             * passing parameters easier */
            yprev[0] = p->r[i];
            yprev[1] = p->phi[i];
            yprev[2] = p->z[i];
            yprev[3] = p->vpar[i];
            yprev[4] = p->mu[i];
            mass = p->mass[i];
            charge = p->charge[i];

	    /* Magnetic field at initial position already known */
	    B_dB[0] = p->B_r[i];
	    B_dB[1] = p->B_r_dr[i];
	    B_dB[2] = p->B_r_dphi[i];
	    B_dB[3] = p->B_r_dz[i];

	    B_dB[4] = p->B_phi[i];
	    B_dB[5] = p->B_phi_dr[i];
	    B_dB[6] = p->B_phi_dphi[i];
	    B_dB[7] = p->B_phi_dz[i];

	    B_dB[8] = p->B_z[i];
	    B_dB[9] = p->B_z_dr[i];
	    B_dB[10] = p->B_z_dphi[i];
	    B_dB[11] = p->B_z_dz[i];

	    E_field_eval_E(E, yprev[0], yprev[1], yprev[2], Edata, Bdata);
            phys_eomgc(k1, yprev, mass, charge, B_dB, E);
            int j;
            /* particle coordinates for the subsequent ydot evaluations are
             * stored in tempy */
            for(j = 0; j < 5; j++) {
                tempy[j] = yprev[j] + h[i]/2.0*k1[j];
            }

            B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
            E_field_eval_E(E, tempy[0], tempy[1], tempy[2], Edata, Bdata);
            phys_eomgc(k2, tempy, mass, charge, B_dB, E);
            for(j = 0; j < 5; j++) {
                tempy[j] = yprev[j] + h[i]/2.0*k2[j];
            }

            B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);;
	    E_field_eval_E(E, tempy[0], tempy[1], tempy[2], Edata, Bdata);
            phys_eomgc(k3, tempy, mass, charge, B_dB, E);
            for(j = 0; j < 5; j++) {
                tempy[j] = yprev[j] + h[i]*k3[j];
            }

            B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
	    E_field_eval_E(E, tempy[0], tempy[1], tempy[2], Edata, Bdata);
            phys_eomgc(k4, tempy, mass, charge, B_dB, E);
            for(j = 0; j < 5; j++) {
                y[j] = yprev[j]
                    + h[i]/6.0 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
            } 

            p->r[i] = y[0];
            p->phi[i] = y[1];
            p->z[i] = y[2];
            p->vpar[i] = y[3];

	    /* Evaluate magnetic field (and gradient) and rho at new position */
	    B_field_eval_B_dB(B_dB, p->r[i], p->phi[i], p->z[i], Bdata);
	    p->B_r[i]        = B_dB[0];
	    p->B_r_dr[i]     = B_dB[1];
	    p->B_r_dphi[i]   = B_dB[2];
	    p->B_r_dz[i]     = B_dB[3];

	    p->B_phi[i]      = B_dB[4];
	    p->B_phi_dr[i]   = B_dB[5];
	    p->B_phi_dphi[i] = B_dB[6];
	    p->B_phi_dz[i]   = B_dB[7];

	    p->B_z[i]        = B_dB[8];
	    p->B_z_dr[i]     = B_dB[9];
	    p->B_z_dphi[i]   = B_dB[10];
	    p->B_z_dz[i]     = B_dB[11];

	    
	    real psi[1];
	    real rho[1];
	    B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i], Bdata);
	    B_field_eval_rho(rho, psi[0], Bdata);
	    p->rho[i] = rho[0];

	    
	    /* Evaluate pol angle so that it is cumulative */
	    real axis_r = B_field_get_axis_r(Bdata);
	    real axis_z = B_field_get_axis_z(Bdata);
	    p->pol[i] += atan2( (R0-axis_r) * (p->z[i]-axis_z) - (z0-axis_z) * (p->r[i]-axis_r), 
				(R0-axis_r) * (p->r[i]-axis_r) + (z0-axis_z) * (p->z[i]-axis_z) );
	    
        }
    }
}
