/**
 * @file simulate_fo_fixed.h
 * @brief Header file for simulate_fo_fixed.c
 */
#ifndef SIMULATE_FO_FIXED_H
#define SIMULATE_FO_FIXED_H

#include "../ascot5.h"
#include "../simulate.h"
#include "../particle.h"

//DECLARE_TARGET
#pragma acc routine gang
void simulate_fo_fixed(particle_queue* pq, sim_data* sim);
//DECLARE_TARGET_END

#endif
