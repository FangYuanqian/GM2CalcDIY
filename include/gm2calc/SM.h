/* ====================================================================
 * This file is part of GM2Calc.
 *
 * GM2Calc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * GM2Calc is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GM2Calc.  If not, see
 * <http://www.gnu.org/licenses/>.
 * ==================================================================== */

#ifndef GM2_SM_H
#define GM2_SM_H

/**
 * @file SM.h
 * @brief contains declarations of C interface functions for the SM
 *
 * This file contains the declarations for the C interface functions
 * used to set and retrieve the SM parameters and masses.
 *
 * @todo adding getters and setters for SM parameters
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct gm2calc_SM gm2calc_SM;

/** allocate new SM */
gm2calc_SM* gm2calc_generalthdm_sm_new();

/** delete SM */
void gm2calc_generalthdm_sm_free(gm2calc_SM*);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
