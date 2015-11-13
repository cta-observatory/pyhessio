/* ============================================================================

   Copyright (C) 1995, 2001, 2006, 2007, 2010  Konrad Bernloehr

   This file is part of the eventio/hessio library.

   The eventio/hessio library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this library. If not, see <http://www.gnu.org/licenses/>.

============================================================================ */

/** @file moments.c
 *  @short Calculate mean, rms, skewness, and kurtosis of data.
 *
 *  @author  Konrad Bernloehr
 *  @date    1995 to 2010
 *  $Date: 2011/02/28 09:56:42 $
 *  $Revision: 1.3 $
 */


#include "histogram.h"

/* ------------------------- alloc_moments ---------------------- */
/**
 *    Allocate a structure for sums of powers of data.
 *    Returns NULL if no structure could be allocated.
 *
 *    @param low Lower limit of range for truncation
 *    @param high Upper limit of range for truncation
 *
 *    @return Pointer to allocated structure or NULL.
 *
 */

MOMENTS *alloc_moments (HISTVALUE_REAL low, HISTVALUE_REAL high)
{
   MOMENTS *tmom;

   if ( low > high )
      return ((MOMENTS *) NULL);
   if ( (tmom = (MOMENTS *) malloc(sizeof(MOMENTS))) == (MOMENTS *) NULL )
      return ((MOMENTS *) NULL);
   tmom->lower_limit = low;
   tmom->upper_limit = high;
   clear_moments(tmom);
   return(tmom);
}

/* ---------------------- clear_moments --------------------------- */
/**
 *    Initialize an existing moments structure (except for its
 *    range limits).
 *
 *    @param   mom  Pointer to moments structure
 *
 */

void clear_moments (MOMENTS *mom)
{
   if ( mom != (MOMENTS *) NULL )
   {
      mom->sum = mom->sum2 = mom->sum3 = mom->sum4 =
         mom->tsum = mom->tsum2 = mom->tsum3 = mom->tsum4 = 0.;
      mom->entries = mom->tentries = 0;
   }
   mom->level = 0;
}

/* ------------------------ free_moments ------------------------- */
/**
 *    Deallocates memory previously allocated to a moments structure.
 *
 *    @param mom Pointer to previously allocated structure
 *
 */

void free_moments (MOMENTS *mom)
{
   if ( mom == (MOMENTS *) NULL )
      return;
   free(mom);
}

/* ---------------------------- fill_moments -------------------------- */
/**
 *    Add up those things needed to compute
 *         mean,
 *         standard deviation,
 *         skewness, and
 *         kurtosis
 *      (both for all data and separately for data in a
 *      range defined in alloc_moments().
 *
 *    @param  mom   Pointer to previously allocated MOMENTS structure.
 *    @param  value One measurement value
 *
 */

void fill_moments (MOMENTS *mom, HISTVALUE_REAL value)
{
   double x2, x3, x4;

   if ( mom == (MOMENTS *) NULL )
      return; /* Histogram has not been allocated */

   x2 = value*value;
   x3 = x2*value;
   x4 = x2*x2;

   mom->sum  += (HISTSUM_REAL) value;
   mom->sum2 += x2;
   mom->sum3 += x3;
   mom->sum4 += x4;
   mom->entries++;

   if ( value >= mom->lower_limit && value < mom->upper_limit )
   {
      mom->tsum += (HISTSUM_REAL) value;
      mom->tsum2 += x2;
      mom->tsum3 += x3;
      mom->tsum4 += x4;
      mom->tentries++;
   }
   mom->level = 4;
}

/* ---------------------- fill_mean_and_sigma ------------------------ */
/**
 *    Add up those things needed to compute
 *        -- mean,
 *        -- standard deviation,
 *      (both for all data and separately for data in a
 *      range defined in alloc_moments().
 *
 *    @param mom   Pointer to previously allocated MOMENTS structure.
 *    @param value One measurement value
 *
 */

void fill_mean_and_sigma (MOMENTS *mom, HISTVALUE_REAL value)
{
   HISTSUM_REAL x2;

   if ( mom == (MOMENTS *) NULL )
      return; /* Histogram has not been allocated */

   x2 = value*value;

   mom->sum  += (HISTSUM_REAL) value;
   mom->sum2 += x2;
   mom->entries++;

   if ( value >= mom->lower_limit && value < mom->upper_limit )
   {
      mom->tsum += (HISTSUM_REAL) value;
      mom->tsum2 += x2;
      mom->tentries++;
   }
   mom->level = 2;
}

/* -------------------------- fill_mean ----------------------------- */
/**
 *    Add up those things needed to compute
 *        -- mean,
 *      (both for all data and separately for data in a
 *      range defined in alloc_moments().
 *
 *    @param mom   Pointer to previously allocated MOMENTS structure.
 *    @param value One measurement value
 *
 */

void fill_mean (MOMENTS *mom, HISTVALUE_REAL value)
{
   if ( mom == (MOMENTS *) NULL )
      return; /* Histogram has not been allocated */

   mom->sum  += (HISTSUM_REAL) value;
   mom->entries++;

   if ( value >= mom->lower_limit && value < mom->upper_limit )
   {
      mom->tsum += (HISTSUM_REAL) value;
      mom->tentries++;
   }
}

/* ------------------------- fill_real_moments ----------------------- */
/**
 *    Add up those things needed to compute
 *        -- mean,
 *        -- standard deviation,
 *        -- skewness, and
 *        -- kurtosis
 *      (both for all data and separately for data in a
 *      range defined in alloc_moments().
 *
 *    @param  mom    Pointer to previously allocated MOMENTS structure.
 *    @param  value  One measurement value
 *    @param  weight Weighting factor of this value
 *
 */

void fill_real_moments (MOMENTS *mom, HISTVALUE_REAL value, double weight)
{
   HISTSUM_REAL x1, x2, x3, x4;

   if ( mom == (MOMENTS *) NULL )
      return; /* Histogram has not been allocated */

   x1 = value*weight;
   x2 = x1*value;
   x3 = x2*value;
   x4 = x3*value;

   mom->sum  += x1;
   mom->sum2 += x2;
   mom->sum3 += x3;
   mom->sum4 += x4;
   mom->entries++;

   if ( value >= mom->lower_limit && value < mom->upper_limit )
   {
      mom->tsum  += x1;
      mom->tsum2 += x2;
      mom->tsum3 += x3;
      mom->tsum4 += x4;
      mom->tentries++;
   }
   mom->level = 4;
}

/* -------------------- fill_real_mean_and_sigma ---------------------- */
/**
 *    Add up those things needed to compute
 *        -- mean,
 *        -- standard deviation,
 *      (both for all data and separately for data in a
 *      range defined in alloc_moments().
 *
 *    @param  mom    Pointer to previously allocated MOMENTS structure.
 *    @param  value  One measurement value
 *    @param  weight Weighting factor of this value
 *
 */

void fill_real_mean_and_sigma (MOMENTS *mom, HISTVALUE_REAL value, double weight)
{
   HISTSUM_REAL x1, x2;

   if ( mom == (MOMENTS *) NULL )
      return; /* Histogram has not been allocated */

   x1 = (HISTSUM_REAL) (value*weight);
   x2 = (HISTSUM_REAL) (value*value*weight);

   mom->sum  += x1;
   mom->sum2 += x2;
   mom->entries++;

   if ( value >= mom->lower_limit && value < mom->upper_limit )
   {
      mom->tsum += x1;
      mom->tsum2 += x2;
      mom->tentries++;
   }
   mom->level = 2;
}

/* ----------------------- fill_real_mean --------------------------- */
/**
 *    Add up those things needed to compute
 *        -- mean,
 *      (both for all data and separately for data in a
 *      range defined in alloc_moments().
 *
 *    @param  mom    Pointer to previously allocated MOMENTS structure.
 *    @param  value  One measurement value
 *    @param  weight Weighting factor of this value
 *
 */

void fill_real_mean (MOMENTS *mom, HISTVALUE_REAL value, double weight)
{
   HISTSUM_REAL x1;
   
   if ( mom == (MOMENTS *) NULL )
      return; /* Histogram has not been allocated */

   x1 = (HISTSUM_REAL) (value * weight);
   
   mom->sum  += x1;
   mom->entries++;

   if ( value >= mom->lower_limit && value < mom->upper_limit )
   {
      mom->tsum += x1;
      mom->tentries++;
   }
   mom->level = 1;
}

/* ------------------------- stat_moments --------------------- */
/**
 *    Calculate moments (mean, rms, skewness, kurtosis) from
 *    the sums of powers of data values.
 *
 *    @param mom 'moments' structure with the sums of the powers of
 *                    data values (only 1st power if only mean to be
 *                    calculated, also 2nd power if r.m.s. to be calculated,
 *                    and also 3rd and 4th if skewness and kurtosis wanted.
 *    @param  stmom Pointer to structure for computed moments
 *
 *    @return  0 (o.k.), -1 and -2 (invalid data)
 *
 */

int stat_moments (MOMENTS *mom, struct momstat *stmom)
{
   double m3, m4, s2, sig2;

   stmom->mean = stmom->sigma = stmom->skewness = stmom->kurtosis = 0.;
   stmom->tmean = stmom->tsigma = stmom->tskewness = stmom->tkurtosis = 0.;

   if ( mom == (MOMENTS *) NULL )
      return(-1); /* Histogram has not been allocated */
   if ( mom->entries == 0 )
      return(-1); /* No sums available */
   stmom->mean = mom->sum / mom->entries;
   if ( mom->entries > 1 && mom->level > 1 )
   {
      s2 = stmom->mean*stmom->mean;
      sig2 = (mom->sum2-mom->sum*stmom->mean)/(mom->entries-1);
      stmom->sigma = sqrt(sig2);
      if ( sig2 != 0. && mom->level > 2 )
      {
         m3 = (mom->sum3 - stmom->mean*3.*mom->sum2) / mom->entries +
              2.*s2*stmom->mean;
         m4 = (mom->sum4 - stmom->mean * (4.*mom->sum3 -
              6.*stmom->mean*mom->sum2)) / mom->entries - 3.*s2*s2;
         stmom->skewness = m3/(sig2*stmom->sigma);
         stmom->kurtosis = m4/(sig2*sig2) - 3.;
      }
   }
   if ( mom->tentries == 0 )
      return(-2); /* No sums available */
   stmom->tmean = mom->tsum / mom->tentries;
   if ( mom->tentries > 1 && mom->level > 1)
   {
      s2 = stmom->tmean*stmom->tmean;
      sig2 = (mom->tsum2-mom->tsum*stmom->tmean) / (mom->tentries-1);
      stmom->tsigma = sqrt(sig2);
      if ( sig2 != 0. && mom->level > 2 )
      {
         m3 = (mom->tsum3 - stmom->tmean*3.*mom->tsum2) / mom->tentries +
             2.*s2*stmom->tmean;
         m4 = (mom->tsum4 - stmom->tmean * (4.*mom->tsum3 -
             6.*stmom->tmean*mom->tsum2)) / mom->tentries - 3.*s2*s2;
         stmom->tskewness = m3/(sig2*stmom->tsigma);
         stmom->tkurtosis = m4/(sig2*sig2) - 3.;
      }
   }
   return(0);
}
