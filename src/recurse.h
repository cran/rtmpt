// copied from GSL

/* linalg/recurse.h
 *
 * Copyright (C) 2016, 2017, 2018, 2019 Patrick Alken
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 3, or (at your option) any
 * later version.
 *
 * This source is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * This module contains code to invert triangular matrices
 */

// modified for our purposes

/* define how a problem is split recursively */
#define GSL_LINALG_SPLIT(n)         ((n >= 16) ? ((n + 8) / 16) * 8 : n / 2)
#define GSL_LINALG_SPLIT_COMPLEX(n) ((n >= 8) ? ((n + 4) / 8) * 4 : n / 2)

/* matrix size for crossover to Level 2 algorithms */
#define CROSSOVER              24
#define CROSSOVER_LU           CROSSOVER
#define CROSSOVER_CHOLESKY     CROSSOVER
#define CROSSOVER_INVTRI       CROSSOVER
#define CROSSOVER_TRIMULT      CROSSOVER

int gsl_linalg_tri_lower_invert_dings(gsl_matrix* T);
