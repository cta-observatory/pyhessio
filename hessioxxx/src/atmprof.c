/* ============================================================================

Copyright (C) 1990, 1997, 1998, 2008, 2009  Konrad Bernloehr

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

============================================================================ */

/* ================================================================== */
/**
 *  @file atmprof.c
 *  @short A stripped-down version of the interpolation of atmospheric
 *         profiles from the atmo.c file of the CORSIKA IACT/ATMO package.
 *
 *  The main differences are
 *    a) parameters are passed by value instead of FORTRAN by-reference way,
 *    b) the height is measured in meters.
 *
 *  The CORSIKA built-in profiles are not handled here.
 *
 *  @author  Konrad Bernloehr 
 *  @date    @verbatim CVS $Date: 2010/07/20 13:37:47 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.6 $ @endverbatim
 *
 */
/* ==================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "atmprof.h"
#include "fileopen.h"

/* Variables used for atmospheric profiles */
static int current_atmosphere;

#define MAX_PROFILE 50
static int num_prof;
static double p_alt[MAX_PROFILE], p_log_alt[MAX_PROFILE];
static double p_log_rho[MAX_PROFILE], p_rho[MAX_PROFILE];
static double p_log_thick[MAX_PROFILE];
static double p_log_n1[MAX_PROFILE];

static double top_of_atmosphere = 112.83e3;
static double bottom_of_atmosphere = 0.;

/* ================================================================== */
/*
   Linear interpolation functions from an older program of K.B.
   A binary search algorithm is used for fast interpolation.
*/

static void interp(double x, double *v, int n, int *ipl, double *rpl);

/* --------------------------- interp ------------------------------- */
/**
 *  @short Linear interpolation with binary search algorithm.
 *
 *  Linear interpolation between data point in sorted (i.e. monotonic
 *  ascending or descending) order. This function determines between
 *  which two data points the requested coordinate is and where between
 *  them. If the given coordinate is outside the covered range, the
 *  value for the corresponding edge is returned.
 *
 *  A binary search algorithm is used for fast interpolation.
 *
 *  @param  x Input: the requested coordinate
 *  @param  v Input: tabulated coordinates at data points
 *  @param  n Input: number of data points
 *  @param  ipl Output: the number of the data point following the requested
 *	    coordinate in the given sorting (1 <= ipl <= n-1)
 *  @param  rpl Output: the fraction (x-v[ipl-1])/(v[ipl]-v[ipl-1])
 *	    with 0 <= rpl <= 1
*/      

static void interp ( double x, double *v, int n, int *ipl, double *rpl )
{
   int i, l, m, j, lm;

#ifdef DEBUG_TEST_ALL
   if ( v == NULL || n <= 2 )
   {
      fprintf(stderr,"Invalid parameters for interpolation.\n");
      *ipl = 1;
      *rpl = 0.;
      return;
   }
#endif

   if ( v[0] < v[n-1] )
   {
      if (x <= v[0])
      {
         *ipl = 1;
         *rpl = 0.;
         return;
      }
      else if (x >= v[n-1])
      {
         *ipl = n-1;
         *rpl = 1.;
         return;
      }
      lm = 0;
   }
   else
   {
      if (x >= v[0])
      {
         *ipl = 1;
         *rpl = 0.;
         return;
      }
      else if (x <= v[n-1])
      {
         *ipl = n-1;
         *rpl = 1.;
         return;
      }
      lm = 1;
   }

   l = (n+1)/2-1;
   m = (n+1)/2;
   for (i=1; i<=30; i++ )
   {
      j = l;
      if (j < 1) j=1;
      if (j > n-1) j=n-1;
      if (x >= v[j+lm-1] && x <= v[j-lm])
      {
         *ipl = j;
         if ( v[j] != v[j-1] )
            *rpl = (x-v[j-1])/(v[j]-v[j-1]);
         else
            *rpl = 0.5;
         return;
      }
      m = (m+1)/2;
      if (x > v[j-1])
         l = l + (1-2*lm)*m;
      else
         l = l - (1-2*lm)*m;
   }
   fprintf(stderr,"Interpolation error.\n");
}

/* ----------------------------- rpol ------------------------------- */

static double rpol ( double *x, double *y, int n, double xp );

/**
 *  @short Linear interpolation with binary search algorithm.
 *
 *  Linear interpolation between data point in sorted (i.e. monotonic
 *  ascending or descending) order. The resulting interpolated value
 *  is returned as a return value.
 *
 *  This function calls interp() to find out where to interpolate.
 *  
 *  @param   x  Input: Coordinates for data table
 *  @param   y  Input: Corresponding values for data table
 *  @param   n  Input: Number of data points
 *  @param   xp Input: Coordinate of requested value
 *
 *  @return  Interpolated value
 *
*/

static double rpol ( double *x, double *y, int n, double xp )
{
   int ipl = 1;
   double rpl = 0.;

   interp ( xp, x, n, &ipl, &rpl );
   return y[ipl-1]*(1.-rpl) + y[ipl]*rpl;
}

/* ======================================================================= */

/* ------------------------ find_elsewhere ------------------------ */
/** 
 *  Find the atmospheric profiles elsewhere (in the sim_telarray configuration).
 */

static char *find_elsewhere (const char *fname, char *bf, size_t sz);

static char *find_elsewhere (const char *fname, char *bf, size_t sz)
{
   const char *sub_dirs[] = { "common", "hess", "hess2", "hess3", 
         "hess5000", "CTA" };
   const char *base_dirs[] = { "$SIM_TELARRAY_PATH", 
      "$HESSROOT/sim_telarray", "$CTA_PATH/sim_telarray",
      "../sim_telarray", ".", ".." };
   size_t ib, is, len;

   for ( ib=0; ib<sizeof(base_dirs)/sizeof(base_dirs[0]); ib++ )
   {
      if ( base_dirs[ib][0] == '$' )
      {
         const char *e = NULL, *f=base_dirs[ib]+1;
         const char *s = strchr(f,'/');
         if ( s != NULL && (s-f) < (int) sz )
         {
            strncpy(bf,f,s-f);
            bf[s-f] = '\0';
            e = getenv(bf);
         }
         else
            e = getenv(f);
         if ( e == NULL )
            continue;
         if ( (len=strlen(e)) < sz )
            strcpy(bf,e);
         if ( s != NULL )
         {
            len += strlen(s);
            if ( len >= sz )
               continue;
            strcat(bf,s);
         }
      }
      else
      {
         len = strlen(base_dirs[ib]);
         if ( len >= sz )
            continue;
         strcpy(bf,base_dirs[ib]);
      }
      for ( is=0; is<sizeof(sub_dirs)/sizeof(sub_dirs[0]); is++ )
      {
         if ( len + 5 + strlen(sub_dirs[is]+1+strlen(fname)) >= sz )
            continue;
         strcpy(bf+len,"/cfg/");
         strcpy(bf+len+5,sub_dirs[is]);
         strcat(bf+len+5,"/");
         strcat(bf+len+5,fname);
         if ( strpbrk(bf,";&<>\\") != NULL )
         {
            fprintf(stderr,"Invalid character in '%s'.\n", bf);
            continue;
         }
         if ( access(bf,R_OK) == 0 )
            return bf;
      }
   }

   return NULL;
}

/* ----------------------- init_atmprof ------------------------ */
/**
 *  @short Initialize atmospheric profiles.
 *
 *  Atmospheric models are read in from text-format tables.
 *  For the interpolation of relevant parameters (density, thickness,
 *  index of refraction, ...) all parameters are transformed such
 *  that linear interpolation can be easily used.
 *
 *  @param atmosphere Atmosphere number, to be expanded to the table file name.
 *
 *  @return 0 (OK) or -1 (error, e.g. table available)
 *
*/

int init_atmprof (int atmosphere)
{
   char fname[128];
   FILE *f;
   char line[1024];
   int count;
   double alt,rho,thick,n_1;
#ifdef LONG_ATMPROF
   double p,t,N,O3,H2O;
#endif
   char *other_fname;

   fflush(stdout);

   /* CORSIKA built-in atmospheres have atmosphere numbers <= 0 */
   if ( atmosphere <=0 )
   {
      fprintf(stderr,"init_atmprof: Cannot handle CORSIKA build-in atmospheric profiles.\n");
      return -1;
   }

   /* There are two different versions of data files. */
#ifndef LONG_ATMPROF
   sprintf(fname,"atmprof%d.dat",atmosphere);
#else
   sprintf(fname,"atm_profile_model_%d.dat",atmosphere);
#endif
   if ( (f=fileopen(fname,"r")) == NULL )
   {
      char bf[1024];
      perror(fname);

      if ( (other_fname = find_elsewhere(fname,bf,sizeof(bf)-1)) != NULL )
      {
         if ( (f=fileopen(other_fname,"r")) == NULL )
            perror(fname);
         else
            fprintf(stderr,"Found the file under '%s'.\n", other_fname);
      }
   }
   if ( f == NULL )
   {
#ifdef LONG_ATMPROF /* Try the other variant before giving up. */
      sprintf(fname,"atmprof%d.dat",atmosphere);
#else
      sprintf(fname,"atm_profile_model_%d.dat",atmosphere);
#endif
      fprintf(stderr,"Trying file %s instead.\n", fname);
      if ( (f=fileopen(fname,"r")) == NULL )
      {
         char bf[1024];
         perror(fname);
         if ( (other_fname = find_elsewhere(fname,bf,sizeof(bf)-1)) != NULL )
         {
            if ( (f=fileopen(other_fname,"r")) == NULL )
               perror(fname);
            else
               fprintf(stderr,"Found the file under '%s'.\n", other_fname);
         }
      }
   }
   if ( f == NULL )
      return -1;

   count = num_prof = 0;
   while ( fgets(line,sizeof(line)-1,f) != NULL && num_prof < MAX_PROFILE )
   {
      char *s;
      count++;

      for (s=line;*s==' ';s++)
         ;
      if ( *s=='#' ) /* Comment line */
         continue;
#ifndef LONG_ATMPROF
      /* The short files contain only data relevant for CORSIKA. */
      if ( sscanf(s,"%lf %lf %lf %lf",
        &alt,&rho,&thick,&n_1) != 4 )
#else
      /* The long files contain other data as well. */
      if ( sscanf(s,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
        &alt,&p,&t,&N,&rho,&thick,&O3,&H2O,&n_1) != 9 )
#endif
      {
         fprintf(stderr,"Syntax error in %s line %d.\n",fname,count);
         exit(1);
      }
      
      p_alt[num_prof] = alt*1e3; /* Altitude in file was in km */
      p_log_alt[num_prof] = (alt>0.)?log(alt*1e3):0.;
      p_log_rho[num_prof] = (rho>0.)?log(rho):-1000.;
      p_rho[num_prof] = rho;
      p_log_thick[num_prof] = (thick>0.)?log(thick):-1000.;
      p_log_n1[num_prof] = (n_1>0.)?log(n_1):-1000.;
      num_prof++;
   }
   
   fclose(f);
   fflush(stdout);
   printf("\n Atmospheric profile %d with %d levels read from file %s\n\n",
      atmosphere,num_prof,fname);

   if ( num_prof < 5 )
   {
      fprintf(stderr,
         "There are definitely too few atmospheric levels in this file.\n");
      fprintf(stderr,
         "Normally this kind of file should have 50 levels.\n");
      return -1;
   }

   bottom_of_atmosphere = p_alt[0];
   top_of_atmosphere    = p_alt[num_prof-1];
   
   current_atmosphere = atmosphere;
   
   return 0;
}

/* ---------------------------- rhofx ----------------------------- */
/**
 *
 *  @short Density of the atmosphere as a function of altitude.
 *
 *  @param  height altitude [m]
 *
 *  @return density [g/cm**3]
*/

double rhofx (double height)
{
   if ( current_atmosphere <= 0 )
      return 0.;
   return exp(rpol(p_alt,p_log_rho,num_prof,height));
}

/* ---------------------------- thickx ----------------------------- */
/**
 *
 *  @short Atmospheric thickness [g/cm**2] as a function of altitude.
 *
 *  @param  height altitude [m]
 *
 *  @return thickness [g/cm**2]
 *
*/

double thickx (double height)
{
   if ( current_atmosphere <= 0 )
      return 0.;
   return exp(rpol(p_alt,p_log_thick,num_prof,height));
}

/* ---------------------------- refidx_ ----------------------------- */
/**
 *
 *  @short Index of refraction as a function of altitude [cm].
 *
 *  @param height altitude [m]
 *
 *  @return index of refraction
 *
*/

double refidx (double height)
{
   if ( current_atmosphere <= 0 )
      return 1.;
   return 1.+exp(rpol(p_alt,p_log_n1,num_prof,height));
}

/* ---------------------------- heighx ----------------------------- */
/**
 *
 *  @short Altitude [m] as a function of atmospheric thickness [g/cm**2].
 *
 *  @param   thick atmospheric thickness [g/cm**2]
 *
 *  @return  altitude [m]
*/

double heighx (double thick)
{
   double h;
   if ( current_atmosphere <= 0 )
      return 0.;
   if ( (thick) <= 0. )
      return top_of_atmosphere;
   h = rpol(p_log_thick,p_alt,num_prof,log(thick));
   if ( h < top_of_atmosphere )
      return h;
   else
      return top_of_atmosphere;
}
