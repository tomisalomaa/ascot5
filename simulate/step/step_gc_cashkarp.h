/**
 * @file step_gc_rk4.h
 * @brief Header file for step_gc_rk4.c
 */
#ifndef STEP_GC_CASHKARP_H
#define STEP_GC_CASHKARP_H

#include "../../B_field.h"
#include "../../E_field.h"
#include "../../particle.h"

DECLARE_TARGET
void step_gc_cashkarp(particle_simd_gc* p, real* h, real* hnext,
                      real tol, B_field_data* Bdata, E_field_data* Edata);
DECLARE_TARGET_END

#endif
