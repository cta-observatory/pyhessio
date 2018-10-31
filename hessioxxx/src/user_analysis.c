/* ============================================================================

Copyright (C) 2006, 2009, 2010, 2013  Konrad Bernloehr

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

/* ================================================================ */
/** @file user_analysis.c
 *  @short Code for analysis of simulated (and reconstructed)
 *         showers within the framework of the read_hess program.
 *
 *  Users wanting to make use of such analysis should modify the
 *  user_* functions provided here or the do_user_ana() function.
 *  Except for the do_user_ana() function and the user_set_...() functions, 
 *  all functions are declared as static to emphasize that their 
 *  interfaces can be changed here to the user's desires.
 *
 *  @author  Konrad Bernloehr 
 *  @date    initial version: August 2006
 *  @date    @verbatim CVS $Date: 2017/10/14 17:51:29 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.79 $ @endverbatim
 */

#include <limits.h>
#include "initial.h"
#include "io_basic.h"
#include "mc_tel.h"
#include "io_hess.h"
#include "io_histogram.h"
#include "fileopen.h"
#include "rec_tools.h"
#include "reconstruct.h"
#include "user_analysis.h"
#include "atmprof.h"
#include "straux.h"
#include "basic_ntuple.h"

#define MAX_TEL_TYPES 10

#ifndef PATH_MAX
# define PATH_MAX 4096
#endif

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

static int verbosity = 0;

struct tel_type_param
{
   int min_tel_id;
   int max_tel_id;
   double mirror_area;
   double flen;
   int num_pixels;
};

static int user_init_done = 0;
static int current_tel_type = 0;
static struct tel_type_param def_tel_type_param[MAX_TEL_TYPES];
static int saved_tel_type[H_MAX_TEL];
static char user_lookup_fname[1024]; 
static char hist_fname[1024];

/**
 *   Override the automatic naming for lookup files.
 */

void user_set_lookup_file (const char *fname)
{
   if ( fname != NULL && strlen(fname) < sizeof(user_lookup_fname) )
      strncpy(user_lookup_fname,fname,sizeof(user_lookup_fname)-1);
   else
   {
      fflush(stdout);
      fprintf(stderr,"Invalid lookup file name\n");
   }
}

/**
 *   Override the automatic naming for histogram files.
 */

void user_set_histogram_file (const char *fname)
{
   if ( fname != NULL && strlen(fname) < sizeof(hist_fname) )
      strncpy(hist_fname,fname,sizeof(hist_fname)-1);
   else
   {
      fflush(stdout);
      fprintf(stderr,"Invalid histogram file name\n");
   }
}

/**
 *  Select a specific telescope type for setting user parameters.
 */

void user_set_telescope_type (int itype)
{
   if ( itype > 0 && itype <= MAX_TEL_TYPES )
      current_tel_type = itype;
   else
      current_tel_type = 0;
}

/** 
 * Set telescope type parameters from a string (e.g. on the command line).
 *
 * Can be used to set all relevant parameters (others set to 0)
 * or just to switch the active type (no parameters other than the type number).
 */

int user_set_tel_type_param_by_str(const char *str)
{
   char word[6][32];
   int ipos = 0;
   int i;
   int itype = 0;
   struct tel_type_param *ttp = NULL;

   for (i=0; i<6; i++)
      word[i][0] = '\0';
   for (i=0; i<6; i++)
   {
      /* Take care that multiple commas or blanks will be treated like one. */
      /* So you should not provide something like '3,,,100.' to leave */
      /* some parameters at defaults but rather '3,0,0,100.' */
      if ( getword(str,&ipos,word[i],sizeof(word[i])-1,',','\0') <= 0 )
         break;
   }
   if ( word[0][0] == '\0' )
   {
      fflush(stdout);
      fprintf(stderr,"Bad telescope type parameters.\n");
      return -1;
   }
   itype = atoi(word[0]);
   if ( i == 1 && itype == 0 )
   {
      if ( current_tel_type != 0 )
         printf("Resetting to global parameter range.\n");
      current_tel_type = itype;
      return 0;
   }
   if ( itype < 1 || itype > MAX_TEL_TYPES )
   {
      fflush(stdout);
      fprintf(stderr,"Telescope type out of range.\n");
      return -1;
   }
   if ( i== 1 ) /* Just type number given but no parameters */
   {
      current_tel_type = itype;
      return 0;
   }

   ttp = &def_tel_type_param[itype-1];
   if ( word[1][0] != '\0' )
      ttp->min_tel_id = atoi(word[1]);
   else
      ttp->min_tel_id = 0;
   if ( word[2][0] != '\0' )
      ttp->max_tel_id = atoi(word[2]);
   else
      ttp->max_tel_id = 0;
   if ( word[3][0] != '\0' )
      ttp->mirror_area = atof(word[3]);
   else
      ttp->mirror_area = 0.;
   if ( word[4][0] != '\0' )
      ttp->flen = atof(word[4]);
   else
      ttp->flen = 0.;
   if ( word[5][0] != '\0' )
      ttp->num_pixels = atoi(word[5]);
   else
      ttp->num_pixels = 0;

   current_tel_type = itype;
   printf("Telescope type %d expected to match with:\n", itype);
   if ( ttp->min_tel_id != 0 || ttp->max_tel_id != 0 )
      printf("   telescope ID from %d to %d\n", ttp->min_tel_id, ttp->max_tel_id);
   if ( ttp->mirror_area != 0. )
      printf("   mirror area of about %3.1f m^2\n", ttp->mirror_area);
   if ( ttp->flen != 0. )
      printf("   focal length of about %4.2f m\n", ttp->flen);
   if ( ttp->num_pixels != 0 )
      printf("   number of pixels of about %d\n", ttp->num_pixels);
   if ( ttp->min_tel_id == 0 && ttp->max_tel_id == 0 && 
        ttp->mirror_area == 0. && ttp->flen == 0. && ttp->num_pixels == 0 )
   {
      printf("   this is a fallback type for any telescope.\n");
   }
   return 1;
}

/**
 *  Find out to which telescope type a telescope belongs, by best matching
 *  in the required parameters.
 */

int which_telescope_type (const struct hess_camera_settings_struct *cam_set);

int which_telescope_type (const struct hess_camera_settings_struct *cam_set)
{
   int best_type = 0;
   double best_match = -100.;
   int i;
   int idx = 0;

   for (i=0; i<MAX_TEL_TYPES; i++)
   {
      double m = 0.;
      struct tel_type_param *ttp = &def_tel_type_param[i];
      /* Matches in the telescope ID would give a huge boost. */
      if ( ttp->min_tel_id > 0 && ttp->max_tel_id > 0 )
      {
         if ( cam_set->tel_id >= ttp->min_tel_id && 
              cam_set->tel_id <= ttp->max_tel_id )
            m += 10.;
      }
      /* For the other parameters you can only get an improvement */
      /* if within a factor of 2 (up or down) from the desired value. */
      /* Missing parameters are ignored. */
      if ( ttp->mirror_area > 0. && cam_set->mirror_area > 0. )
      {
         m += 1. - fabs(log(ttp->mirror_area/cam_set->mirror_area)/log(2.)); /* Within factor of two */
      }
      if ( ttp->flen > 0. && cam_set->flen > 0. )
      {
         m += 1. - fabs(log(ttp->flen/cam_set->flen)/log(1.3)); /* Within 30% */
      }
      if ( ttp->num_pixels > 0 && cam_set->num_pixels > 0 )
      {
         m += 1. - fabs(log(ttp->num_pixels/(double)cam_set->num_pixels)/log(1.5)); /* Within 50% */
      }
      if ( m > best_match )
      {
         best_match = m;
         best_type = i+1;
      }
   }

   idx = find_tel_idx(cam_set->tel_id);
   if ( verbosity > 0 )
      printf("CT%d (#%d) is of type %d\n", cam_set->tel_id, idx, best_type);
   if ( idx >= 0 && idx < H_MAX_TEL )
      saved_tel_type[idx] = best_type;

   return best_type;
}

struct telescope_list
{
   size_t min_tel;
   size_t ntel;
   int *tel_id;
};

static struct telescope_list *alt_list = NULL;
static size_t n_list = 0;

static double max_theta = 0.2 * (M_PI/180.);
static double min_theta = 0.2 * (M_PI/180.);

static struct user_parameters up[MAX_TEL_TYPES+2];

/** Number of parameters, including: the gamma-ray source offset plus
 *    d_sp_idx, min_amp, tailcut_low, tailcut_high, min_pix, reco_flag,
 *    min_tel_img, max_tel_img, max_theta, theta_scale
 */
static int nparams, nparams_i, nparams_d;
static double *params;

struct user_parameters* user_get_parameters(int tp)
{
   if ( tp >= 0 && tp <= MAX_TEL_TYPES )
      return &up[tp];
   else
      return &up[0];
}

/** Get the best matching telescope type for a given telescope index.
 *  If user analysis is not activated, this will always be type 0.
 */

int user_get_type (int itel)
{
   if ( ! user_init_done )
      return 0;
   if ( itel >= 0 && itel < H_MAX_TEL )
      return saved_tel_type[itel];
   else
      return 0;
}

/** Evaluate energy-dependent cut parameters with
 *
 *   @param cut[0] the cut parameter at 1 TeV (lgE=0),
 *   @param cut[1] the slope of the cut parameters versus lgE,
 *   @param cut[2] the minimum cut parameter,
 *   @param cut[3] the maximum cut parameter.
 */

static double eval_cut_param(double *cut, double lgE);

static double eval_cut_param(double *cut, double lgE)
{
   double c = cut[0];
   if ( cut[1] == 0. )
      return c;
   c += cut[1]*lgE;
   if ( c < cut[2] )
      c = cut[2];
   if ( c > cut[3] )
      c = cut[3];
   return c;
}

void __attribute__ ((constructor)) user_init_parameters (void)
{
   static int user_init_param_done = 0;
   if ( ! user_init_param_done )
   {
     int i;
     size_t j, k;

     nparams_i = sizeof(up[i].i)/sizeof(int);
     nparams_d = sizeof(up[i].d)/sizeof(double);
     nparams = (nparams_i + nparams_d) * (MAX_TEL_TYPES+2) + 4;
     params = (double *) calloc(nparams,sizeof(double));

     /* We initialize more configurations than we expect to
        handle. The first one (i=0) will contain parameters
        only used globally but never used telescope-specific.
        It would also apply as a last resort if no reasonably
        matching type could be found. The actual telescope
        types should be numbers from 1 to (including) MAX_TEL_TYPES.
        The last one (n=MAX_TEL_TYPES+1) will just keep the 
        hardcoded defaults filled in here. It should not get 
        referenced later. */
     for ( i=0; i<= MAX_TEL_TYPES+1; i++ )
     {
      up[i].i.user_flags = 0;
      up[i].i.min_pix = 2.;
      up[i].i.reco_flag = 0;
      up[i].i.min_tel_img = 2;
      up[i].i.max_tel_img = 100;
      up[i].i.lref = 0;

      up[i].i.integrator = 0;
      up[i].i.integ_param[0] = up[i].i.integ_param[1] = 0;
      up[i].i.integ_no_rescale = 0;

      up[i].d.source_offset_deg = 0.0;
      up[i].d.minfrac = 0.;
      up[i].d.d_sp_idx = 0.0;
      up[i].d.min_amp = 80.;
      up[i].d.tailcut_low = 5.;
      up[i].d.tailcut_high = 10.;

      up[i].d.max_theta_deg = 0.2;
      up[i].d.theta_scale = 1.0;

      up[i].d.de2_cut_param[0] = 0.5;
      up[i].d.de2_cut_param[1] = 0.0;
      up[i].d.de2_cut_param[2] = 0.5;
      up[i].d.de2_cut_param[3] = 0.5;

      max_theta = up[i].d.max_theta_deg*(M_PI/180.);

      up[i].d.mscrw_min[0] = -2.0;
      up[i].d.mscrw_min[1] =  0.0;
      up[i].d.mscrw_min[2] = -2.0;
      up[i].d.mscrw_min[3] = -2.0;
      up[i].d.mscrw_max[0] = 0.6;
      up[i].d.mscrw_max[1] = 0.0;
      up[i].d.mscrw_max[2] = 0.6;
      up[i].d.mscrw_max[3] = 0.6;
      up[i].d.mscrl_min[0] = -2.0;
      up[i].d.mscrl_min[1] =  0.0;
      up[i].d.mscrl_min[2] = -2.0;
      up[i].d.mscrl_min[3] = -2.0;
      up[i].d.mscrl_max[0] = 1.2;
      up[i].d.mscrl_max[1] = 0.0;
      up[i].d.mscrl_max[2] = 1.2;
      up[i].d.mscrl_max[3] = 1.2;

      up[i].d.eres_cut_param[0] = 1.0;
      up[i].d.eres_cut_param[1] = 0.0;
      up[i].d.eres_cut_param[2] = 1.0;
      up[i].d.eres_cut_param[3] = 1.0;
      up[i].d.hmax_cut_param = 1.0;

      up[i].d.min_theta_deg = 0.2;
      min_theta = up[i].d.min_theta_deg*(M_PI/180.);
      up[i].d.camera_clipping_deg = 0.0;

      up[i].d.theta_escale[0] = 1.0;
      up[i].d.theta_escale[1] = 0.0;
      up[i].d.theta_escale[2] = 1.0;
      up[i].d.theta_escale[3] = 1.0;

      up[i].d.clip_amp = 0.0;

      for (k=0; k<2; k++)
         for (j=0; j<4; j++)
            up[i].d.d_integ_param[k][j] = 0.0;

      up[i].d.calib_scale = 0.;
      up[i].d.focal_length = 0.;
     }
     user_init_param_done = 1;
   }
}

/** Set user-defined flags: used to active HESS-style analysis 
 *
 *  @param uf 0: not exactly HESS-style analysis;
 *            1: HESS-style standard cuts;
 *            2: HESS-style hard cuts;
 *            3: HESS-style loose cuts.
 *            >=4: HESS-style (no re-scaling) but user-defined cut parameters.
 */

void user_set_flags (int uf)
{
   /* This particular parameter has global effect and not */
   /* a telescope-specific one. So it is applied to all setups */
   /* (except the one keeping the compiled-in defaults). */
   int i;
   printf("User flags: %d (%s).\n", uf,
      uf==0?"With rescaling of mscrw, mscrl":
      uf==1?"HESS-style standard cuts":
      uf==2?"HESS-style hard cuts":
      uf==3?"HESS-style loose cuts":
      "No rescaling of mscrw, mscrl");
   for (i=0; i<=MAX_TEL_TYPES; i++ )
   {
      up[i].i.user_flags = uf;
      if ( uf== 1 || uf == 2 || uf == 3 )
      {
         up[i].d.mscrw_min[0] = up[i].d.mscrw_min[1] = 
         up[i].d.mscrw_min[2] = up[i].d.mscrw_min[3] = 0.;
         up[i].d.mscrw_max[0] = up[i].d.mscrw_max[1] = 
         up[i].d.mscrw_max[2] = up[i].d.mscrw_max[3] = 0.;
         up[i].d.mscrl_min[0] = up[i].d.mscrl_min[1] = 
         up[i].d.mscrl_min[2] = up[i].d.mscrl_min[3] = 0.;
         up[i].d.mscrl_max[0] = up[i].d.mscrl_max[1] = 
         up[i].d.mscrl_max[2] = up[i].d.mscrl_max[3] = 0.;
      }
   }
}

/**
 *  @short Set the difference between generated MC spectrum and
 *        the assumed source spectrum.
 */

void user_set_spectrum(double di)
{
   /* This one cannot be telescope-type specific. */
   int i;
   for (i=0; i<=MAX_TEL_TYPES; i++ )
   {
      up[i].d.d_sp_idx = di;
   }
   printf("Events will be weighted by E^%5.3f\n", up[0].d.d_sp_idx);
}

/**
 *  @short Set the acceptable ranges for reconstructed impact positions.
 */

void user_set_impact_range (double *impact_range)
{
   /* This one cannot be telescope-type specific. */
   int i;
   for (i=0; i<=MAX_TEL_TYPES; i++ )
   {
      up[i].d.impact_range[0] = impact_range[0];
      up[i].d.impact_range[1] = impact_range[1];
      up[i].d.impact_range[2] = impact_range[2];
   }
   if ( impact_range[0] > 0. )
      printf("Array center must be at least %4.2f m from reconstructed shower axis.\n", 
         impact_range[0]);
   if ( impact_range[1] > 0. )
      printf("Reconstructed core position must have %4.2f <= x <= %4.2f\n",
         -impact_range[1], impact_range[1]);
   if ( impact_range[2] > 0. )
      printf("Reconstructed core position must have %4.2f <= y <= %4.2f\n",
         -impact_range[2], impact_range[2]);
}

/**
 *  @short Set the acceptable ranges for true impact positions.
 */

void user_set_true_impact_range (double *true_impact_range)
{
   /* This one cannot be telescope-type specific. */
   int i;
   for (i=0; i<=MAX_TEL_TYPES; i++ )
   {
      up[i].d.true_impact_range[0] = true_impact_range[0];
      up[i].d.true_impact_range[1] = true_impact_range[1];
      up[i].d.true_impact_range[2] = true_impact_range[2];
   }
   if ( true_impact_range[0] > 0. )
      printf("Array center must be at least %4.2f m from true shower axis.\n", 
         true_impact_range[0]);
   if ( true_impact_range[1] > 0. )
      printf("True core position must have %4.2f <= x <= %4.2f\n",
         -true_impact_range[1], true_impact_range[1]);
   if ( true_impact_range[2] > 0. )
      printf("True core position must have %4.2f <= y <= %4.2f\n",
         -true_impact_range[2], true_impact_range[2]);
}

/**
 *  Set the maximum core distance for telescopes if their images should
 *  be used beyond geometrical reconstruction.
 */

void user_set_max_core_distance (double rt)
{
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   up[c].d.max_core_distance = rt;
   if ( c > 0 )
      printf("Maximum core distance for telescope type %d: %f\n", c, rt);
   if ( c == 0 ) /* Setting unspecified type means setting all types. */
   {
      int i;
      printf("Maximum core distance for all telescope types: %f\n", rt);
      for (i=1; i<=MAX_TEL_TYPES; i++ )
      {
         up[i].d.max_core_distance = rt;
      }
   }
}

/**
 *  @short Set the minimum amplitude of images usable for the analysis.
 */

void user_set_min_amp(double a)
{
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   up[c].d.min_amp = a;
   if ( c > 0 )
      printf("Minimum image amplitude for telescope type %d: %f\n", c, a);
   if ( c == 0 ) /* Setting unspecified type means setting all types. */
   {
      int i;
      printf("Minimum image amplitude for all telescope types: %f\n", a);
      for (i=1; i<=MAX_TEL_TYPES; i++ )
      {
         up[i].d.min_amp = a;
      }
   }
}

/**
 *  @short Set the lower and upper tail cuts for the 
 *         standard two-level tail-cut scheme. 
 */

void user_set_tail_cuts(double tcl, double tch, int lref, double minfrac)
{
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   if ( tcl > tch ) /* wrong order: swap parameters */
   {
      double tt = tcl;
      tcl = tch;
      tch = tt;
   }
   up[c].d.tailcut_low = tcl;
   up[c].d.tailcut_high = tch;
   up[c].i.lref = lref;
   up[c].d.minfrac = minfrac;
   if ( c > 0 )
      printf("Tail-cut parameters for telescope type %d: %f, %f\n", c, tcl, tch);
   if ( c == 0 ) /* Setting unspecified type means setting all types. */
   {
      int i;
      printf("Tail-cut parameters for all telescope types: %f, %f\n", tcl, tch);
      for (i=1; i<=MAX_TEL_TYPES; i++ )
      {
         up[i].d.tailcut_low = tcl;
         up[i].d.tailcut_high = tch;
      }
   }
}

/**
 *  @short Set the minimum number of significant pixels in usable images.
 */

void user_set_min_pix(int mpx)
{
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   up[c].i.min_pix = mpx;
   if ( c > 0 )
      printf("Minimum number of pixels in image for tel. type %d: %d\n", c, mpx);
   if ( c == 0 ) /* Setting unspecified type means setting all types. */
   {
      int i;
      printf("Minimum number of pixels in image for all tel. types: %d\n", mpx);
      for (i=1; i<=MAX_TEL_TYPES; i++ )
      {
         up[i].i.min_pix = mpx;
      }
   }
}

/**
 *  @short Set the reconstruction level flag ('-r' option in read_hess).
 */

void user_set_reco_flag(int rf)
{
#if 0
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
#else
   int c = 0; /* This option is not used in a telescope-specific way so far. */
#endif

   up[c].i.reco_flag = rf;

   if ( c == 0 ) /* Setting unspecified type means setting all types. */
   {
      int i;
      for (i=1; i<=MAX_TEL_TYPES; i++ )
      {
         up[i].i.reco_flag = rf;
      }
   }
}

/**
 *  @short Set the minimum and maximum number of usable images for events used in analysis.
 */

void user_set_tel_img(int tmn, int tmx)
{
#if 0
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
#else
   int c = 0; /* This option is not used in a telescope-specific way so far. */
#endif

   up[c].i.min_tel_img = tmn;
   up[c].i.max_tel_img = tmx;

   if ( c > 0 )
      printf("Required telescope type %d images: %d to %d\n", c, tmn, tmx);
   if ( c == 0 ) /* Setting unspecified type means setting all types. */
   {
      int i;
      printf("Required telescope images: %d to %d\n", tmn, tmx);
      for (i=1; i<=MAX_TEL_TYPES; i++ )
      {
         up[i].i.min_tel_img = tmn;
         up[i].i.max_tel_img = tmx;
      }
   }
}

/**
 *   @short You may have alternative selections of (fewer) telescopes.
 */

void user_set_tel_list (size_t min_tel, size_t ntel, int *tel_id)
{
   size_t i;
   
   if ( ntel < 1 || tel_id == 0 )
      return;
   printf("Alternate telescope selection: any %zu of CTs", min_tel);
   for (i=0; i<ntel; i++)
      printf(" %d", tel_id[i]);
   printf("\n");

   if ( alt_list == NULL )
   {
      alt_list = (struct telescope_list *) 
                  malloc(sizeof(struct telescope_list));
      n_list = 0;
   }
   else
   {
      struct telescope_list *oa = alt_list;
      alt_list = (struct telescope_list *) 
                  realloc(alt_list,(n_list+1)*sizeof(struct telescope_list));
      /* Avoid cppcheck false positive (see exit(1) below!) */
      if ( alt_list == NULL )
         free(oa);
   }
   if ( alt_list == NULL )
   {
      perror("Allocation of telescope list failed");
      exit(1);
      n_list = 0;
   }
      
   alt_list[n_list].min_tel = min_tel;
   alt_list[n_list].ntel = ntel;
   alt_list[n_list].tel_id = (int *) malloc(ntel*sizeof(int));
   
   for (i=0; i<ntel; i++)
      alt_list[n_list].tel_id[i] = tel_id[i];
   
   n_list++;
}

/** Angular cut limit is multiplicity dependent. */
static double opt_theta_cut[7][H_MAX_TEL];

/**
 *  @short Set the maximum angle between source and reconstructed shower direction.
 */

void user_set_max_theta(double thmax, double thscale, double thmin)
{
   /* This is another global parameter, set for all types. */
   int i;
   for (i=0; i<=MAX_TEL_TYPES; i++)
   {
      up[i].d.min_theta_deg = thmin;
      if ( thmax > 0. )
         up[i].d.max_theta_deg = thmax;
      else if ( up[i].i.user_flags == 1 ) /* HESS standard cuts */
         up[i].d.max_theta_deg = up[i].d.min_theta_deg = sqrt(0.0125);
      else if ( up[i].i.user_flags == 2 ) /* HESS hard cuts */
         up[i].d.max_theta_deg = up[i].d.min_theta_deg = 0.1;
      else if ( up[i].i.user_flags == 3 ) /* HESS loose cuts */
         up[i].d.max_theta_deg = up[i].d.min_theta_deg = 0.2;
      up[i].d.theta_scale = thscale; /* Applies to both fixed and optimized theta cut. */
   }
   max_theta = up[0].d.max_theta_deg * (M_PI/180.);
   min_theta = up[0].d.min_theta_deg * (M_PI/180.);
   printf("Maximum theta value: %f deg (scaled by %f), minimum: %f deg\n",
      up[0].d.max_theta_deg, up[0].d.theta_scale, up[0].d.min_theta_deg);
}

/** By default the angular acceptance is the 80% containment radius.
 *  Performance may improve by using a smaller radius at low energies
 *  (stricter cut) and a larger radius at high energies (looser cut).
 *  This sets an additional lg(E) dependent scaling factor.
 */

void user_set_theta_escale (double *thes)
{
   /* This is another global parameter, set for all types. */
   int i, j;
   for (i=0; i<=MAX_TEL_TYPES; i++)
   {
      for ( j=0; j<4; j++ )
         up[i].d.theta_escale[j] = thes[j];
   }
   if ( up[0].d.theta_escale[2] != 1. || up[0].d.theta_escale[3] != 1. )
   {
      printf("Theta limit made energy-dependent (in addition to multiplicity-dependence).\n");
      printf("Scaling parameters: %f,%f,%f,%f\n",
         up[0].d.theta_escale[0], up[0].d.theta_escale[1],
         up[0].d.theta_escale[2], up[0].d.theta_escale[3]);
   }
}

/** The dE cut can be made more or less strict by a scale
 *  parameter which should be 1.0 by default and is
 *  below 1 for a stricter cut and above 1 for a looser cut.
 */

void user_set_de_cut (double *dec)
{
   /* This is another global parameter, set for all types. */
   if ( dec[0] > 0. )
   {
      int i;
      for (i=0; i<=MAX_TEL_TYPES; i++)
      {
         up[i].d.eres_cut_param[0] = dec[0];
         up[i].d.eres_cut_param[1] = dec[1];
         up[i].d.eres_cut_param[2] = dec[2];
         up[i].d.eres_cut_param[3] = dec[3];
      }
   }
   if ( dec[1] != 0. )
      printf("dE cut scaling parameters: %f,%f,%f,%f\n",dec[0],dec[1],dec[2],dec[3]);
   else
      printf("dE cut scaling parameter: %f\n",dec[0]);
}

/** Since the dE2 cut is not always of any help with default cut parameters,
 *  you can change the parameter to your needs.
 */

void user_set_de2_cut (double *de2c)
{
   /* This is another global parameter, set for all types. */
   if ( de2c[0] > 0. )
   {
      int i;
      for (i=0; i<=MAX_TEL_TYPES; i++)
      {
         up[i].d.de2_cut_param[0] = de2c[0];
         up[i].d.de2_cut_param[1] = de2c[1];
         up[i].d.de2_cut_param[2] = de2c[2];
         up[i].d.de2_cut_param[3] = de2c[3];
      }
   }
   if ( de2c[1] != 0. )
      printf("dE2 cut parameters: %f,%f,%f,%f\n",de2c[0],de2c[1],de2c[2],de2c[3]);
   else
      printf("dE2 cut parameter: %f\n",de2c[0]);
}

/** The hmax cut can be made or or less strict by a scale
 *  parameter which should be 1.0 by default and is
 *  below 1 for a stricter cut and above 1 for a looser cut.
 */

void user_set_hmax_cut (double hmaxc)
{
   /* This is another global parameter, set for all types. */
   if ( hmaxc > 0. )
   {
      int i;
      for (i=0; i<=MAX_TEL_TYPES; i++)
      {
         up[i].d.hmax_cut_param = hmaxc;
      }
   }
   printf("Hmax cut scaling parameter: %f\n",hmaxc);
}

/** Set shape cut parameters */

void user_set_shape_cuts (double wmin, double wmax, double lmin, double lmax)
{
#if 0
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
#else
   int c = 0; /* This option is not used in a telescope-specific way so far. */
#endif
   up[c].d.mscrw_min[0] = wmin; up[c].d.mscrw_min[1] = 0.;
   up[c].d.mscrw_max[0] = wmax; up[c].d.mscrw_max[1] = 0.;
   up[c].d.mscrl_min[0] = lmin; up[c].d.mscrl_min[1] = 0.;
   up[c].d.mscrl_max[0] = lmax; up[c].d.mscrl_max[1] = 0.;
   if ( c > 0 )
      printf("Shape cuts for tel. type %d: %f,%f,%f,%f\n",c,wmin,wmax,lmin,lmax);
   if ( c == 0 )
   {
      int i;
      printf("Shape cuts: %f,%f,%f,%f\n",wmin,wmax,lmin,lmax);
      for (i=1; i<=MAX_TEL_TYPES; i++)
      {
         up[i].d.mscrw_min[0] = wmin; up[i].d.mscrw_min[1] = 0.;
         up[i].d.mscrw_max[0] = wmax; up[i].d.mscrw_max[1] = 0.;
         up[i].d.mscrl_min[0] = lmin; up[i].d.mscrl_min[1] = 0.;
         up[i].d.mscrl_max[0] = lmax; up[i].d.mscrl_max[1] = 0.;
      }
   }
}

/** Set energy dependent scaled width limit. */

void user_set_width_max_cut (double *wmax)
{
#if 0
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
#else
   int c = 0; /* This option is not used in a telescope-specific way so far. */
#endif
   up[c].d.mscrw_max[0] = wmax[0];
   up[c].d.mscrw_max[1] = wmax[1];
   up[c].d.mscrw_max[2] = wmax[2];
   up[c].d.mscrw_max[3] = wmax[3];
   if ( c > 0 )
   {
      if ( wmax[1] != 0. )
         printf("Maximum width cut parameters for tel. type %d: %f,%f,%f,%f\n", 
            c,wmax[0],wmax[1],wmax[2],wmax[3]);
      else
         printf("Maximum width cut parameter for tel. type %d: %f\n", c, wmax[0]);
   }
   if ( c == 0 )
   {
      int i;
      if ( wmax[1] != 0. )
         printf("Maximum width cut parameters: %f,%f,%f,%f\n", 
            wmax[0],wmax[1],wmax[2],wmax[3]);
      else
         printf("Maximum width cut parameter: %f\n", wmax[0]);
      for (i=1; i<=MAX_TEL_TYPES; i++)
      {
         up[i].d.mscrw_max[0] = wmax[0];
         up[i].d.mscrw_max[1] = wmax[1];
         up[i].d.mscrw_max[2] = wmax[2];
         up[i].d.mscrw_max[3] = wmax[3];
      }
   }
}

/** Set energy dependent scaled length limit. */

void user_set_length_max_cut (double *lmax)
{
#if 0
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
#else
   int c = 0; /* This option is not used in a telescope-specific way so far. */
#endif
   up[c].d.mscrl_max[0] = lmax[0];
   up[c].d.mscrl_max[1] = lmax[1];
   up[c].d.mscrl_max[2] = lmax[2];
   up[c].d.mscrl_max[3] = lmax[3];
   if ( c > 0 )
   {
      if ( lmax[1] != 0. )
         printf("Maximum length cut parameters for tel. type %d: %f,%f,%f,%f\n", 
            c,lmax[0],lmax[1],lmax[2],lmax[3]);
      else
         printf("Maximum length cut parameter for tel. type %d: %f\n", c, lmax[0]);
   }
   if ( c == 0 )
   {
      int i;
      if ( lmax[1] != 0. )
         printf("Maximum length cut parameters: %f,%f,%f,%f\n", 
            lmax[0],lmax[1],lmax[2],lmax[3]);
      else
         printf("Maximum length cut parameter: %f\n", lmax[0]);
      for (i=1; i<=MAX_TEL_TYPES; i++)
      {
         up[i].d.mscrl_max[0] = lmax[0];
         up[i].d.mscrl_max[1] = lmax[1];
         up[i].d.mscrl_max[2] = lmax[2];
         up[i].d.mscrl_max[3] = lmax[3];
      }
   }
}

/** Set the telescope effective focal length. */

void user_set_focal_length (double f)
{
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   up[c].d.focal_length = f;
   if ( c == 0 )
   {
      int i;
      for ( i=1; i<=MAX_TEL_TYPES; i++ )
         up[i].d.focal_length = f;
   }
}

/** Set the maximum radius to be used of a camera. */

void user_set_clipping (double dc)
{
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   up[c].d.camera_clipping_deg = dc;
   if ( c == 0 )
   {
      int i;
      for ( i=1; i<=MAX_TEL_TYPES; i++ )
         up[i].d.camera_clipping_deg = dc;
   }
}

/** Set the maximum amplitude in a pixel */

void user_set_clipamp (double cpa)
{
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   up[c].d.clip_amp = cpa;
   if ( c == 0 )
   {
      int i;
      for ( i=1; i<=MAX_TEL_TYPES; i++ )
         up[i].d.clip_amp = cpa;
   }
}

/** Set the required trigger type(s) as a bit pattern */

void user_set_trg_req (int trg_req)
{
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   up[c].i.trg_req = trg_req;
   if ( c == 0 )
   {
      int i;
      for ( i=1; i<=MAX_TEL_TYPES; i++ )
         up[i].i.trg_req = trg_req;
   }
}

static int diffuse_mode = 0;
static double diffuse_off_axis_min = 0.;
static double diffuse_off_axis_max = M_PI/2.;

void user_set_diffuse_mode(int dm, double oar[])
{
   diffuse_mode = dm;
   diffuse_off_axis_min = oar[0] * (M_PI/180.);
   diffuse_off_axis_max = oar[1] * (M_PI/180.);
}

void user_set_verbosity(int v)
{
   verbosity = v;
}

static int event_selected = 0;

int user_selected_event ()
{
   return event_selected;
}

static int auto_lookup = 0;

void user_set_auto_lookup (int al)
{
   auto_lookup = al;
}

void user_set_integrator(int scheme)
{
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   up[c].i.integrator = scheme;
   if ( c == 0 )
   {
      int i;
      for ( i=1; i<=MAX_TEL_TYPES; i++ )
         up[i].i.integrator = scheme;
   }
}

void user_set_integ_window(int nsum, int noff, int ps_opt)
{
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   up[c].i.integ_param[0] = nsum;
   up[c].i.integ_param[1] = noff;
   up[c].i.integ_param[2] = ps_opt;
   if ( c == 0 )
   {
      int i;
      for ( i=1; i<=MAX_TEL_TYPES; i++ )
      {
         up[i].i.integ_param[0] = nsum;
         up[i].i.integ_param[1] = noff;
         up[i].i.integ_param[2] = ps_opt;
      }
   }
}

void user_set_integ_threshold(int ithg, int itlg)
{
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   up[c].i.integ_thresh[0] = ithg;
   up[c].i.integ_thresh[1] = itlg;
   if ( c == 0 )
   {
      int i;
      for ( i=1; i<=MAX_TEL_TYPES; i++ )
      {
         up[i].i.integ_thresh[0] = ithg;
         up[i].i.integ_thresh[1] = itlg;
      }
   }
}

void user_set_integ_no_rescale (int no)
{
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   up[c].i.integ_no_rescale = no;
   if ( c == 0 )
   {
      int i;
      for ( i=1; i<=MAX_TEL_TYPES; i++ )
      {
         up[i].i.integ_no_rescale = no;
      }
   }
}

void user_set_calib_scale (double s)
{
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   up[c].d.calib_scale = s;
   if ( c == 0 )
   {
      int i;
      for ( i=1; i<=MAX_TEL_TYPES; i++ )
      {
         up[i].d.calib_scale = s;
      }
   }
}

void user_set_nb_radius (double *r)
{
   int ib;
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   for ( ib=0; ib<3; ib++ )
      up[c].d.r_nb[ib] = r[ib];
   if ( c == 0 )
   {
      int i;
      for ( i=1; i<=MAX_TEL_TYPES; i++ )
      {
         for ( ib=0; ib<3; ib++ )
            up[i].d.r_nb[ib] = r[ib];
      }
   }
}

void user_set_nxt_radius (double r)
{
   int c = current_tel_type;
   if ( c < 0 || c > MAX_TEL_TYPES )
      return;
   up[c].d.r_ne = r;
   if ( c == 0 )
   {
      int i;
      for ( i=1; i<=MAX_TEL_TYPES; i++ )
      {
         
         up[i].d.r_ne = r;
      }
   }
}

/* ----------------------- expected_max_height ----------------------- */
/**
 *  @short Expected height of the shower maximum above the detector
 *         for gamma rays, based on simple analytical
 *         formula and exponential atmospheric profile.
 *
 *  @param E      The energy of the shower [TeV].
 *  @param theta  Then zenith angle of the shower [radians].
 *  @param height The height above sea level of the experiment [m].
 *  @return Height of shower maximum above detector [m]
 */

static double expected_max_height(double E, double theta, double height);

static double expected_max_height(double E, double theta, double height)
{
   double rho_0 = 1.2; /* (extrapolated) density at sea level, [kg/m^3] */
   double s = 8300.;   /* scale height, [m] */
   double X0 = 366.;   /* radiation length of air, [kg/m^2] */
   double Ec = 81.e-6; /* critical energy in air, [TeV] */
   double tmax;

   if ( E <= 0. )
      return -1.;

   tmax = X0*log(E/Ec+0.5); /* This is the PDG rather than the Rossi formula. */

   return s * log( (rho_0*s) / (cos(theta)*tmax) ) - height;
}

/* ----------------------- expected_max_distance ----------------------- */
/**
 *  @short Expected distance of the shower maximum from the detector
 *         for gamma rays, based on simple analytical
 *         formula and exponential atmospheric profile.
 *
 *  @param E      The energy of the shower [TeV].
 *  @param theta  Then zenith angle of the shower [radians].
 *  @param height The height above sea level of the experiment [m].
 *  @return Distance of shower maximum from detector [m]
 */

static double expected_max_distance(double E, double theta, double height);

static double expected_max_distance(double E, double theta, double height)
{
   return expected_max_height(E, theta, height) / cos(theta);
}

/* ------------------------  img_norm  -------------------------- */
/**
 *  @short Get scaled + reduced scaled image parameters (both HEGRA
 *         and HESS type scaling) as well as energy scaling from the
 *         lookups.
 *
 *  All variables for the results are optional. For variables which are
 *  of no interest, pass a NULL pointer.
 *
 *  @param w      Image width [rad].
 *  @param l      Image length [rad].
 *  @param A      Image amplitude [ peak p.e. ].
 *  @param lgA    log10(A)
 *  @param rc     Reconstructed core distance.
 *  @param tel_type Telescope type (for multiple lookups).
 *  @param scrw   Variable getting the scaled reduced width (HESS style).
 *  @param scrl   Variable getting the scaled reduced length (HESS style).
 *  @param scw    Variable getting the scaled width (HEGRA style).
 *  @param scl    Variable getting the scaled length (HEGRA style).
 *  @param sce    Variable getting the expected energy [TeV] for
 *                the given amplitude at the given core distance.
 *  @param scer   Variable getting the relative fluctuation of energy/amplitude at this point.
 *  @param rco    Variable getting the expected core distance based on width/length and amplitude.
 *  @param rcor   Variable getting the relative error in the core distance estimate.
 *  @param dimgo  Variable getting the expected distance in the image (as for rco).
 *  @param dimgor Variable getting the relative error in the image distance estimate.
 */

static int img_norm (double w, double l, double A, double lgA, double rc,
        int tel_type,
	double *scrw, double *scrl, double *scw, double *scl,
        double *sce, double *scer, 
        double *rco, double *rcor, double *dimgo, double *dimgor);

static int img_norm (double w, double l, double A, double lgA, double rc,
        int tel_type,
	double *scrw, double *scrl, double *scw, double *scl,
        double *sce, double *scer, 
        double *rco, double *rcor, double *dimgo, double *dimgor)
{
   static HISTOGRAM *h18113[MAX_TEL_TYPES+1], *h18114[MAX_TEL_TYPES+1], 
      *h18123[MAX_TEL_TYPES+1], *h18124[MAX_TEL_TYPES+1], 
      *h18153[MAX_TEL_TYPES+1], *h18154[MAX_TEL_TYPES+1];
   static HISTOGRAM *h18173[MAX_TEL_TYPES+1], *h18174[MAX_TEL_TYPES+1], 
      *h18183[MAX_TEL_TYPES+1], *h18184[MAX_TEL_TYPES+1];
   static double rl, rh, al, ah;
   static double sx, sy;
   static int nx=0, ny=0, nxo=0;
   int ix, iy, ixo;
   double wm, lm, em, ws, ls, es;
   double dom = 0., dos = 0., rom = 0., ros = 0., wol = 0.;

   if ( l > 0. )
      wol = w / l;

   if ( h18113[0] == 0 && h18113[1] == 0 ) /* Initialization after booking */
   {
      int tt = 0;
      for ( tt = 0; tt <= MAX_TEL_TYPES; tt++ )
      {
         int ht = 100000 * tt;
         if ( (h18113[tt] = get_histogram_by_ident(ht+18113)) == NULL ||
              (h18114[tt] = get_histogram_by_ident(ht+18114)) == NULL ||
              (h18123[tt] = get_histogram_by_ident(ht+18123)) == NULL ||
              (h18124[tt] = get_histogram_by_ident(ht+18124)) == NULL ||
              (h18153[tt] = get_histogram_by_ident(ht+18153)) == NULL ||
              (h18154[tt] = get_histogram_by_ident(ht+18154)) == NULL )
            continue;
         h18173[tt] = get_histogram_by_ident(ht+18173);
         h18174[tt] = get_histogram_by_ident(ht+18174);
         h18183[tt] = get_histogram_by_ident(ht+18183);
         h18184[tt] = get_histogram_by_ident(ht+18184);
         if ( h18173[tt] != NULL )
         {
            nxo = h18173[tt]->nbins;
         }

         rl = h18113[tt]->specific.real.lower_limit;
         rh = h18113[tt]->specific.real.upper_limit;
         al = h18113[tt]->specific_2d.real.lower_limit;
         ah = h18113[tt]->specific_2d.real.upper_limit;
         nx = h18113[tt]->nbins;
         ny = h18113[tt]->nbins_2d;
         sx = (rh-rl) / nx;
         sy = (ah-al) / ny;
      }
   }

   if ( tel_type < 0 || tel_type > MAX_TEL_TYPES || h18113[tel_type] == 0 )
      return -1;

   if ( lgA < al || lgA >= ah || rc < rl || rc >= rh )
      return 0;
   ix = (rc-rl) / sx;
   iy = (lgA-al) / sy;
   ixo= wol * nxo;
   if ( ix < 0 || ix >= nx || iy < 0 || iy >= ny )
      return 0;
   wm = h18113[tel_type]->extension->ddata[ix+nx*iy];
   ws = h18114[tel_type]->extension->ddata[ix+nx*iy];
   lm = h18123[tel_type]->extension->ddata[ix+nx*iy];
   ls = h18124[tel_type]->extension->ddata[ix+nx*iy];
   em = h18153[tel_type]->extension->ddata[ix+nx*iy];
   es = h18154[tel_type]->extension->ddata[ix+nx*iy];
   if ( h18173[tel_type] != NULL && ixo >= 0 && ixo < nxo )
   {
      rom = h18173[tel_type]->extension->ddata[ixo+nx*iy];
      ros = h18174[tel_type]->extension->ddata[ixo+nx*iy];
      dom = h18183[tel_type]->extension->ddata[ixo+nx*iy];
      dos = h18184[tel_type]->extension->ddata[ixo+nx*iy];
   }
   if ( scw != NULL )
   {
      if ( wm > 0. )
         *scw = w / wm;
      else
         *scw = 999.;
   }
   if ( scl != NULL )
   {
      if ( lm > 0. )
         *scl = l / lm;
      else
         *scl = 999.;
   }
   if ( scrw != NULL )
   {
      if ( wm > 0. && ws > 0. )
         *scrw = (w-wm) / ws;
      else
         *scrw = 999.;
   }
   if ( scrl != NULL )
   {
      if ( lm > 0. && ls > 0. )
         *scrl = (l-lm) / ls;
      else
         *scrl = 999.;
   }
   if ( sce != NULL )
   {
      if ( A > 0. && em > 0. )
         *sce = A/em;
      else
         *sce = 0.;
   }
   if ( scer != NULL )
   {
      if ( em > 0. && es > 0. )
         *scer = es/em;
      else
         *scer = 999.;
   }
   if ( rco != NULL )
   {
      *rco = rom;
   }
   if ( rcor != NULL && ros >= rom*rom )
   {
      *rcor = sqrt(ros-rom*rom);
   }
   if ( dimgo != NULL )
   {
      *dimgo = dom;
   }
   if ( dimgor != NULL && dos >= dom*dom )
   {
      *dimgor = sqrt(dos-dom*dom);
   }
   return 1;
}

/** Declare local (static) data here ... */

static int telescope_type[H_MAX_TEL];
static char lookup_fname[1024];

static double Az_src, Alt_src, Az_nom, Alt_nom, source_offset;

static MOMENTS *pixmom = NULL;

struct ebias_cor_data
{
   int ndat;
   double *lgE;
   double *lgDE;
};
static struct ebias_cor_data ebias;

/* ------------------------ ebias_correction ----------------------- */
/**
 *  @short Ask for a correction to log10(reconstructed energy), if available.
 *
 *  @return Bias in log10(energy), to be subtracted from log10(energy), or 0.
 */

double ebias_correction (double lgE);

double ebias_correction (double lgE)
{
   if ( ebias.ndat <= 0 || ebias.lgE == NULL || ebias.lgDE == NULL )
      return 0.; /* No correction available */
   return rpol(ebias.lgE, ebias.lgDE, ebias.ndat, lgE);
}


/* ----------------------- set_ebias_correction ---------------------- */
/**
 *  @short Set correction to log10(reconstructed energy), if available.
 *
 */

void set_ebias_correction (HISTOGRAM *h);

void set_ebias_correction (HISTOGRAM *h)
{
   ebias.ndat = 0;
   if ( ebias.lgE != NULL ) { free(ebias.lgE); ebias.lgE = NULL; }
   if ( ebias.lgDE != NULL ) { free(ebias.lgDE); ebias.lgDE = NULL; }

   if ( h == NULL ) /* No histogram available -> no correction later-on. */
   {
      printf("No energy-bias correction available.\n");
      return;
   }

   ebias.ndat = h->nbins;
   ebias.lgE = (double *) calloc(ebias.ndat, sizeof(double));
   ebias.lgDE = (double *) calloc(ebias.ndat, sizeof(double));

   if ( ebias.lgE != NULL && ebias.lgDE != NULL && ebias.ndat > 0 )
   {
      double e1 = h->specific.real.lower_limit;
      double e2 = h->specific.real.upper_limit;
      double de = (e2-e1)/h->nbins;
      int ibin;
      for ( ibin=0; ibin<h->nbins; ibin++ )
      {
         ebias.lgE[ibin] = e1 + (ibin+0.5)*de;
         if ( h->type == 'D' )
            ebias.lgDE[ibin] = h->extension->ddata[ibin];
         else if ( h->type == 'F' )
            ebias.lgDE[ibin] = h->extension->fdata[ibin];
      }
   }
}

static void init_telescope_types (AllHessData *hsdata);
static int tel_types_change = 0;
static int stat_type[MAX_TEL_TYPES+2];

/* --------------- init_telescope_types ----------------- */
/**
 *  @short Initialize what of type each telescope is. In normal 
 *         simulation data this is only needed once but in complex
 *         merged (via merge_simtel) data the necessary info may not be
 *         available for all of them when types for the first of
 *         them is needed.
 */

static void init_telescope_types (AllHessData *hsdata)
{
   int itel, tel_type, num_types = 0;
   int nas = 0;
   
   printf("\nInitialize (or re-initialize) available telescope types.\n");

   for ( itel=0; itel<hsdata->run_header.ntel; itel++ )
   {
      if ( hsdata->camera_set[itel].num_mirrors > 0 )
      {
         if ( telescope_type[itel] == 0 ) /* Report just once */
         {
            telescope_type[itel] = which_telescope_type(&hsdata->camera_set[itel]);
            if ( verbosity >= 0 )
             printf("CT%d of type %d at (%3.1f, %3.1f, %3.1f) m has %3.1f m^2 in %d mirrors (f=%3.1f m).\n",
              hsdata->run_header.tel_id[itel], 
              telescope_type[itel],
              hsdata->run_header.tel_pos[itel][0],
              hsdata->run_header.tel_pos[itel][1],
              hsdata->run_header.tel_pos[itel][2],
              hsdata->camera_set[itel].mirror_area,
              hsdata->camera_set[itel].num_mirrors,
              hsdata->camera_set[itel].flen);
            store_camera_radius(&hsdata->camera_set[itel],itel);
         }
         else if ( tel_types_change )
            telescope_type[itel] = which_telescope_type(&hsdata->camera_set[itel]);
      }
      else /* May report more than once in complex merged data, until info available. */
      {
         if ( verbosity >= 0 )
          printf("CT%d at (%3.1f, %3.1f, %3.1f) m with %3.1f m^2 mirror area not seen active (yet).\n",
            hsdata->run_header.tel_id[itel], 
            hsdata->run_header.tel_pos[itel][0],
            hsdata->run_header.tel_pos[itel][1],
            hsdata->run_header.tel_pos[itel][2],
            hsdata->camera_set[itel].mirror_area);
      }
   }

   /* Let's see for which type numbers we have any telescopes */
   for ( tel_type=0; tel_type<=MAX_TEL_TYPES+1; tel_type++ )
      stat_type[tel_type] = 0;
   for ( itel=0; itel<hsdata->run_header.ntel; itel++ )
   {
      tel_type = telescope_type[itel];
      if ( tel_type >= 0 && tel_type <= MAX_TEL_TYPES )
         stat_type[tel_type]++;
      if ( tel_type > 0 )
         nas++;
   }
   /* If there are gaps in the type presence, keep types actually used. */
   for ( tel_type=0; tel_type<=MAX_TEL_TYPES+1; tel_type++ )
   {
      if ( stat_type[tel_type] > 0 )
         num_types++;
   }
   if ( num_types > 1 )
   {
      printf("%d different types of telescopes were found.\n", num_types);
   }
   if ( nas == hsdata->run_header.ntel )
   {
      printf("All %d telescopes have a telescope type assigned.\n", nas);
   }
   else if ( nas == 0 )
   {
      printf("No telescope has a telescope type assigned.\n");
   }
   else
   {
      printf("%d out of %d telescopes have a telescope type assigned.\n", 
         nas, hsdata->run_header.ntel);
   }

   tel_types_change = 0;
}

static int init_hist_for_type[MAX_TEL_TYPES+2];
static int init_hist_global = 0;

static void book_hist_global(AllHessData *hsdata);
static void book_hist_for_type(AllHessData *hsdata, int itype);

static void book_hist_global(AllHessData *hsdata)
{
   double xylow[2], xyhigh[2];
   int nbins[2];
   int i, itel;

   if ( init_hist_global )
      return;

   printf("Booking global histograms.\n");
   
   xylow[0] = xylow[1] = -2000.;
   xyhigh[0] = xyhigh[1] = 2000.;
   nbins[0] = nbins[1] = 200;
   book_histogram(9100,
      "Core position of event triggered with >= 1 tel.",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9101,
      "Core position of event triggered with >= 1 tel. (10-20 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9102,
      "Core position of event triggered with >= 1 tel. (20-40 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9103,
      "Core position of event triggered with >= 1 tel. (40-80 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9104,
      "Core position of event triggered with >= 1 tel. (80-160 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9105,
      "Core position of event triggered with >= 1 tel. (160-316 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9106,
      "Core position of event triggered with >= 1 tel. (316-631 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9107,
      "Core position of event triggered with >= 1 tel. (0.631-1.26 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9108,
      "Core position of event triggered with >= 1 tel. (1.26-2.5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9109,
      "Core position of event triggered with >= 1 tel. (2.5-5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9110,
      "Core position of event triggered with >= 1 tel. (5-10 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9111,
      "Core position of event triggered with >= 1 tel. (10-20 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9112,
      "Core position of event triggered with >= 1 tel. (20-40 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9113,
      "Core position of event triggered with >= 1 tel. (40-80 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9114,
      "Core position of event triggered with >= 1 tel. (80-160 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9115,
      "Core position of event triggered with >= 1 tel. (>160 TeV)",
      "D", 2, xylow, xyhigh, nbins);

   book_histogram(9150,
      "Core position of event seen with >= 1 image",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9151,
      "Core position of event seen with >= 1 image (10-20 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9152,
      "Core position of event seen with >= 1 image (20-40 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9153,
      "Core position of event seen with >= 1 image (40-80 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9154,
      "Core position of event seen with >= 1 image (80-160 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9155,
      "Core position of event seen with >= 1 image (160-316 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9156,
      "Core position of event seen with >= 1 image (316-631 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9157,
      "Core position of event seen with >= 1 image (0.631-1.26 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9158,
      "Core position of event seen with >= 1 image (1.26-2.5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9159,
      "Core position of event seen with >= 1 image (2.5-5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9160,
      "Core position of event seen with >= 1 image (5-10 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9161,
      "Core position of event seen with >= 1 image (10-20 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9162,
      "Core position of event seen with >= 1 image (20-40 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9163,
      "Core position of event seen with >= 1 image (40-80 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9164,
      "Core position of event seen with >= 1 image (80-160 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9165,
      "Core position of event seen with >= 1 image (>160 TeV)",
      "D", 2, xylow, xyhigh, nbins);

   book_histogram(9200,
      "Core position of events triggered with >= 2 tel.",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9201,
      "Core position of event triggered with >= 2 tel. (10-20 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9202,
      "Core position of event triggered with >= 2 tel. (20-40 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9203,
      "Core position of event triggered with >= 2 tel. (40-80 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9204,
      "Core position of event triggered with >= 2 tel. (80-160 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9205,
      "Core position of event triggered with >= 2 tel. (160-316 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9206,
      "Core position of event triggered with >= 2 tel. (316-631 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9207,
      "Core position of event triggered with >= 2 tel. (0.631-1.26 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9208,
      "Core position of event triggered with >= 2 tel. (1.26-2.5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9209,
      "Core position of event triggered with >= 2 tel. (2.5-5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9210,
      "Core position of event triggered with >= 2 tel. (5-10 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9211,
      "Core position of event triggered with >= 2 tel. (10-20 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9212,
      "Core position of event triggered with >= 2 tel. (20-40 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9213,
      "Core position of event triggered with >= 2 tel. (40-80 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9214,
      "Core position of event triggered with >= 2 tel. (80-160 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9215,
      "Core position of event triggered with >= 2 tel. (>160 TeV)",
      "D", 2, xylow, xyhigh, nbins);

   book_histogram(9250,
      "Core position of event seen with >= 2 images",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9251,
      "Core position of event seen with >= 2 images (10-20 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9252,
      "Core position of event seen with >= 2 images (20-40 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9253,
      "Core position of event seen with >= 2 images (40-80 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9254,
      "Core position of event seen with >= 2 images (80-160 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9255,
      "Core position of event seen with >= 2 images (160-316 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9256,
      "Core position of event seen with >= 2 images (316-631 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9257,
      "Core position of event seen with >= 2 images (0.631-1.26 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9258,
      "Core position of event seen with >= 2 images (1.26-2.5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9259,
      "Core position of event seen with >= 2 images (2.5-5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9260,
      "Core position of event seen with >= 2 images (5-10 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9261,
      "Core position of event seen with >= 2 images (10-20 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9262,
      "Core position of event seen with >= 2 images (20-40 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9263,
      "Core position of event seen with >= 2 images (40-80 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9264,
      "Core position of event seen with >= 2 images (80-160 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9265,
      "Core position of event seen with >= 2 images (>160 TeV)",
      "D", 2, xylow, xyhigh, nbins);

   book_histogram(9300,
      "Core position of events triggered with >= 3 tel.",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9301,
      "Core position of event triggered with >= 3 tel. (10-20 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9302,
      "Core position of event triggered with >= 3 tel. (20-40 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9303,
      "Core position of event triggered with >= 3 tel. (40-80 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9304,
      "Core position of event triggered with >= 3 tel. (80-160 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9305,
      "Core position of event triggered with >= 3 tel. (160-316 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9306,
      "Core position of event triggered with >= 3 tel. (316-631 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9307,
      "Core position of event triggered with >= 3 tel. (0.631-1.26 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9308,
      "Core position of event triggered with >= 3 tel. (1.26-2.5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9309,
      "Core position of event triggered with >= 3 tel. (2.5-5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9310,
      "Core position of event triggered with >= 3 tel. (5-10 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9311,
      "Core position of event triggered with >= 3 tel. (10-20 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9312,
      "Core position of event triggered with >= 3 tel. (20-40 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9313,
      "Core position of event triggered with >= 3 tel. (40-80 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9314,
      "Core position of event triggered with >= 3 tel. (80-160 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9315,
      "Core position of event triggered with >= 3 tel. (>160 TeV)",
      "D", 2, xylow, xyhigh, nbins);

   book_histogram(9350,
      "Core position of event seen with >= 3 images",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9351,
      "Core position of event seen with >= 3 images (10-20 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9352,
      "Core position of event seen with >= 3 images (20-40 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9353,
      "Core position of event seen with >= 3 images (40-80 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9354,
      "Core position of event seen with >= 3 images (80-160 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9355,
      "Core position of event seen with >= 3 images (160-316 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9356,
      "Core position of event seen with >= 3 images (316-631 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9357,
      "Core position of event seen with >= 3 images (0.631-1.26 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9358,
      "Core position of event seen with >= 3 images (1.26-2.5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9359,
      "Core position of event seen with >= 3 images (2.5-5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9360,
      "Core position of event seen with >= 3 images (5-10 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9361,
      "Core position of event seen with >= 3 images (10-20 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9362,
      "Core position of event seen with >= 3 images (20-40 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9363,
      "Core position of event seen with >= 3 images (40-80 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9364,
      "Core position of event seen with >= 3 images (80-160 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9365,
      "Core position of event seen with >= 3 images (>160 TeV)",
      "D", 2, xylow, xyhigh, nbins);


   book_histogram(9400,
      "Core position of events triggered with >= 4 tel.",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9401,
      "Core position of event triggered with >= 4 tel. (10-20 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9402,
      "Core position of event triggered with >= 4 tel. (20-40 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9403,
      "Core position of event triggered with >= 4 tel. (40-80 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9404,
      "Core position of event triggered with >= 4 tel. (80-160 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9405,
      "Core position of event triggered with >= 4 tel. (160-316 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9406,
      "Core position of event triggered with >= 4 tel. (316-631 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9407,
      "Core position of event triggered with >= 4 tel. (0.631-1.26 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9408,
      "Core position of event triggered with >= 4 tel. (1.26-2.5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9409,
      "Core position of event triggered with >= 4 tel. (2.5-5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9410,
      "Core position of event triggered with >= 4 tel. (5-10 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9411,
      "Core position of event triggered with >= 4 tel. (10-20 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9412,
      "Core position of event triggered with >= 4 tel. (20-40 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9413,
      "Core position of event triggered with >= 4 tel. (40-80 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9414,
      "Core position of event triggered with >= 4 tel. (80-160 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9415,
      "Core position of event triggered with >= 4 tel. (>160 TeV)",
      "D", 2, xylow, xyhigh, nbins);

   book_histogram(9450,
      "Core position of event seen with >= 4 images",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9451,
      "Core position of event seen with >= 4 images (10-20 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9452,
      "Core position of event seen with >= 4 images (20-40 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9453,
      "Core position of event seen with >= 4 images (40-80 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9454,
      "Core position of event seen with >= 4 images (80-160 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9455,
      "Core position of event seen with >= 4 images (160-316 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9456,
      "Core position of event seen with >= 4 images (316-631 GeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9457,
      "Core position of event seen with >= 4 images (0.631-1.26 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9458,
      "Core position of event seen with >= 4 images (1.25-2.5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9459,
      "Core position of event seen with >= 4 images (2.5-5 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9460,
      "Core position of event seen with >= 4 images (5-10 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9461,
      "Core position of event seen with >= 4 images (10-20 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9462,
      "Core position of event seen with >= 4 images (20-40 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9463,
      "Core position of event seen with >= 4 images (40-80 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9464,
      "Core position of event seen with >= 4 images (80-160 TeV)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(9465,
      "Core position of event seen with >= 4 images (>160 TeV)",
      "D", 2, xylow, xyhigh, nbins);



   xylow[0] = -0.5; xyhigh[0] = hsdata->run_header.ntel+0.5; nbins[0] = hsdata->run_header.ntel+1;
   book_histogram(10001,
      "No. of triggered telescopes in triggered events", 
      "D", 1, xylow, xyhigh, nbins);
   book_histogram(10004,
      "No. of triggered telescopes in triggered events - no weights", 
      "D", 1, xylow, xyhigh, nbins);
   xylow[1]  = -3.; xyhigh[1] = 4.; nbins[1]  = 140;
   book_histogram(10003,
      "lg(true E) versus no. of triggered telescopes in triggered events", 
      "D", 2, xylow, xyhigh, nbins);
   xylow[1] = -0.5; xyhigh[1] = hsdata->run_header.ntel+0.5; nbins[1] = hsdata->run_header.ntel+1;
   book_histogram(10002,
      "No. of triggered telescopes versus participating Tel. ID in triggered events", 
      "D", 2, xylow, xyhigh, nbins);
   xylow[1]  = -0.5; xyhigh[1] = 10.5; nbins[1]  = 11;
   book_histogram(10005,
      "Telescope type versus no. of triggered telescopes per type in triggered events", 
      "D", 2, xylow, xyhigh, nbins);

   xylow[0]  = 0.;    xylow[1]  = -3.; /* 1 GeV */
   xyhigh[0] = 4000.; xyhigh[1] = 4.;  /* 10 PeV */
   nbins[0]  = 400;   nbins[1]  = 140;
   book_histogram(12001, 
      "lg(true E) versus array core distance, all", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22001,  "lg(true E), all - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22101,  "lg(true E), all - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12002, 
      "lg(true E) versus array core distance, triggered", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22002, "lg(true E), triggered - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22102, "lg(true E),triggered  - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12003, 
      "lg(true E) versus array core distance, passing amplitude cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22003, "lg(true E), passing amplitude cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22103, "lg(true E), passing amplitude cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12004, 
      "lg(true E) versus array core distance, passing ampl.+edge cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22004, "lg(true E), passing ampl.+edge cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22104, "lg(true E), passing ampl.+edge cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12005, 
      "lg(true E) versus array core distance, passing shape cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22005, "lg(true E), passing shape cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22105, "lg(true E), passing shape cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12006, 
      "lg(true E) versus array core distance, passing shape+angle cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22006, "lg(true E), passing shape+angle cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22106, "lg(true E), passing shape+angle cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12007, 
      "lg(true E) versus array core distance, passing shape+angle+dE cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22007, "lg(true E), passing shape+angle+dE cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22107, "lg(true E), passing shape+angle+dE cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12008, 
      "lg(true E) versus array core distance, passing shape+angle+dE+dE2 cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22008, "lg(true E), passing shape+angle+dE+dE2 cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22108, "lg(true E), passing shape+angle+dE+dE2 cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12009, 
      "lg(true E) versus array core distance, passing shape+angle+dE+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22009, "lg(true E), passing shape+angle+dE+dE2+hmax cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22109, "lg(true E), passing shape+angle+dE+dE2+hmax cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);


   book_histogram(12010, 
      "lg(true E) versus array core distance, passing angle+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22010, "lg(true E), passing angle+hmax cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22110, "lg(true E), passing angle+hmax cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12011, 
      "lg(true E) versus array core distance, passing shape+angle+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22011, "lg(true E), passing shape+angle+hmax cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22111, "lg(true E), passing shape+angle+hmax cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12012, 
      "lg(true E) versus array core distance, passing shape+angle+dE+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22012, "lg(true E), passing shape+angle+dE+hmax cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22112, "lg(true E), passing shape+angle+dE+hmax cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12013, 
      "lg(true E) versus array core distance, passing shape+angle+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22013, "lg(true E), passing shape+angle+dE2+hmax cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22113, "lg(true E), passing shape+angle+dE2+hmax cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12014, 
      "lg(true E) versus array core distance, passing angle+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22014, "lg(true E), passing angle+dE2+hmax cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22114, "lg(true E), passing angle+dE2+hmax cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);


   book_histogram(12015, 
      "lg(true E) versus array core distance, passing shape cuts, central 1 deg",
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22015, "lg(true E), passing shape cuts, central 1 deg - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22115, "lg(true E), passing shape cuts, central 1 deg - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);


   book_histogram(12053, 
      "lg(rec. E) versus array core distance, passing angle cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22053, "lg(rec. E), passing angle cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22153, "lg(rec. E), passing angle cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12054,
      "lg(rec. E) versus array core distance, any reconstructed energy", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22054, "lg(rec. E), any reconstructed energy - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22154, "lg(rec. E), any reconstructed energy - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12055, 
      "lg(rec. E) versus array core distance, passing shape cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22055, "lg(rec. E), passing shape cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22155, "lg(rec. E), passing shape cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12056, 
      "lg(rec. E) versus array core distance, passing shape+angle cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22056, "lg(rec. E), passing shape+angle cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22156, "lg(rec. E), passing shape+angle cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12057, 
      "lg(rec. E) versus array core distance, passing shape+angle+dE cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22057, "lg(rec. E), passing shape+angle+dE cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22157, "lg(rec. E), passing shape+angle+dE cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12058, 
      "lg(rec. E) versus array core distance, passing shape+angle+dE+dE2 cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22058, "lg(rec. E), passing shape+angle+dE+dE2 cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22158, "lg(rec. E), passing shape+angle+dE+dE2 cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12059, 
      "lg(rec. E) versus array core distance, passing shape+angle+dE+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22059, "lg(rec. E), passing shape+angle+dE+dE2+hmax cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22159, "lg(rec. E), passing shape+angle+dE+dE2+hmax cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);


   book_histogram(12060, 
      "lg(rec. E) versus array core distance, passing angle+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22060, "lg(rec. E), passing angle+hmax cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22160, "lg(rec. E), passing angle+hmax cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12061, 
      "lg(rec. E) versus array core distance, passing shape+angle+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22061, "lg(rec. E), passing shape+angle+hmax cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22161, "lg(rec. E), passing shape+angle+hmax cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12062, 
      "lg(rec. E) versus array core distance, passing shape+angle+dE+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22062, "lg(rec. E), passing shape+angle+dE+hmax cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22162, "lg(rec. E), passing shape+angle+dE+hmax cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12063, 
      "lg(rec. E) versus array core distance, passing shape+angle+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22063, "lg(rec. E), passing shape+angle+dE2+hmax cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22163, "lg(rec. E), passing shape+angle+dE2+hmax cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12064, 
      "lg(rec. E) versus array core distance, passing angle+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22064, "lg(rec. E), passing angle+dE2+hmax cuts - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22164, "lg(rec. E), passing angle+dE2+hmax cuts - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);


   book_histogram(12065, 
      "lg(rec. E) versus array core distance, passing shape cuts, central 1 deg",
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22065, "lg(rec. E), passing shape cuts, central 1 deg - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22165, "lg(rec. E), passing shape cuts, central 1 deg - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12067, 
      "lg(rec. E) versus array core distance, passing shape+dE cuts, central 1 deg",
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22067, "lg(rec. E), passing shape+dE cuts, central 1 deg - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22167, "lg(rec. E), passing shape+dE cuts, central 1 deg - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12068, 
      "lg(rec. E) versus array core distance, passing shape+dE+dE2 cuts, central 1 deg",
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22068, "lg(rec. E), passing shape+dE+dE2 cuts, central 1 deg - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22168, "lg(rec. E), passing shape+dE+dE2 cuts, central 1 deg - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);

   book_histogram(12069, 
      "lg(rec. E) versus array core distance, passing shape+dE+dE2+hmax cuts, central 1 deg",
      "D", 2, xylow, xyhigh, nbins);
   book_1d_histogram(22069, "lg(rec. E), passing shape+dE+dE2+hmax cuts, central 1 deg - no weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);
   book_1d_histogram(22169, "lg(rec. E), passing shape+dE+dE2+hmax cuts, central 1 deg - with weights",
      "D", xylow[1], xyhigh[1], nbins[1]*5);


   for ( itel=1; itel<9; itel++ )
   {
      int icut;
      const char *tcut[] = { "shape", "shape+angle", "shape+angle+dE+dE2+hmax", "shape+angle+hmax" };
      char title[1024];
      for ( icut=0; (size_t)icut<sizeof(tcut)/sizeof(tcut[0]); icut++ )
      {
         snprintf(title,sizeof(title)-1,
            "lg(E) versus array core distance, passing %s cuts, =%d good images",
            tcut[icut], itel);
         book_histogram(12100+100*icut+itel, title,
            "D", 2, xylow, xyhigh, nbins);
      }
   }

   xylow[0]  = -0.5;
   xyhigh[0] = (double)(nparams)-0.5;
   nbins[0]  = nparams;

   {
      int ip=0, itp;
      params[ip++] = 201.; /* Version */
      params[ip++] = MAX_TEL_TYPES+2; /* How many sets of each type */
      params[ip++] = nparams_i; /* Integers per set */
      for ( itp=0; itp <= MAX_TEL_TYPES+1; itp++ )
         for ( i=0; i<nparams_i; i++ )
            params[ip++] = *(((int *)&up[itp].i)+i);
      params[ip++] = nparams_d; /* Doubles per set */
      for ( itp=0; itp <= MAX_TEL_TYPES+1; itp++ )
         for ( i=0; i<nparams_d; i++ )
            params[ip++] = *(((double *)&up[itp].d)+i);
   }

   {  HISTOGRAM *h;
      if ( (h = get_histogram_by_ident(12099)) != NULL )
         free_histogram(h);
      h = book_histogram(12099,"Analysis parameters", 
            "D", 1, xylow, xyhigh, nbins);
      for (i=0; i<nparams; i++)
         fill_histogram(h,i,0.,params[i]);
   }

   xylow[0]  = 0.;    xylow[1]  = -3.; /* 1 GeV */
   xyhigh[0] = 30e3;  xyhigh[1] = 4.;  /* 10 PeV */
   nbins[0]  = 150.;  nbins[1]  = 140.;
   book_histogram(15001, "Hmax, lg(true E): ew, amp.+edge cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(15002, "Hmax, lg(true E): ew, passing shape cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(15003, "Hmax, lg(true E): ew, passing shape+angular cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(15004, "Hmax, lg(true E): ew, passing shape+angular+dE cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(15005, "Hmax, lg(true E): ew, passing shape+angular+dE+dE2 cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(15006, "Hmax, lg(true E): ew, passing shape+angular+dE+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);

   book_histogram(15101, "Hmax, lg(rec. E): ew, any reconstructed energy",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(15102, "Hmax, lg(rec. E): ew, passing shape cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(15103, "Hmax, lg(rec. E): ew, passing shape+angular cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(15104, "Hmax, lg(rec. E): ew, passing shape+angular+dE cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(15105, "Hmax, lg(rec. E): ew, passing shape+angular+dE+dE2 cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(15106, "Hmax, lg(rec. E): ew, passing shape+angular+dE+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);

   xylow[0]  = -5.;   xylow[1]  = -5.;
   xyhigh[0] = 20.;    xyhigh[1] = 20.;
   nbins[0]  = 250;   nbins[1]  = 250;
   book_histogram(17000, "MSRW, MSRL: ew", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(17001, "MSRW, MSRL: ew, within 1 deg", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(17002, "MSRW, MSRL: ew, passing angular cut", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(17003, "MSRW, MSRL: ew, passing dE cut", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(17004, "MSRW, MSRL: ew, passing angular+dE cut",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(17005, "MSRW, MSRL: ew, passing angular+dE+dE2 cut",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(17006, "MSRW, MSRL: ew, passing angular+dE+dE2+hmax cut",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(17500, "SRW, SRL: ew", 
      "D", 2, xylow, xyhigh, nbins);
   for (i=0;i<10;i++)
      book_histogram(17100+i, "MSRW, MSRL: ew for radial bins", 
         "D", 2, xylow, xyhigh, nbins);
   for (i=0;i<12;i++)
      book_histogram(17200+i, "MSRW, MSRL: ew for energy bins", 
         "D", 2, xylow, xyhigh, nbins);
   for (i=0;i<10;i++)
      book_histogram(17300+i, "MSRW, MSRL: ew for tel. multiplicity", 
         "D", 2, xylow, xyhigh, nbins);
   for (i=0;i<12;i++)
      book_histogram(17400+i, "MSRW, MSRL: ew for Xmax/50", 
         "D", 2, xylow, xyhigh, nbins);

   xylow[0]  = 0.;    xylow[1]  = 0.;
   xyhigh[0] = 100.;  xyhigh[1] = 1.;
   nbins[0]  = 100;   nbins[1]  = 200;
   book_histogram(17810, "hottest pixel 3(2) / pixel 1", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(17800, "hottest pixel / image amp.", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(17802, "pixel amp. / image amp. sigma", 
      "D", 2, xylow, xyhigh, nbins);
   xylow[1] = 0.; xyhigh[1] = 2.;
   book_histogram(17801, "log10(mean pixel amp.)", 
      "D", 2, xylow, xyhigh, nbins);
   xylow[1] = -5.; xyhigh[1] = 5.;
   book_histogram(17803, "pixel amp. / image amp. skewness", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(17804, "pixel amp. / image amp. kurtosis", 
      "D", 2, xylow, xyhigh, nbins);


   xylow[0]  = -3.;  xylow[1]  = 0.;
   xyhigh[0] = 4.;   xyhigh[1] = 800.;
   nbins[0]  = 70;   nbins[1]  = 20;
   book_histogram(18200,
      "Xmax versus Erec: no. of entries",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(18201,
      "Xmax versus Erec: event weights",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(18211,
      "Xmax versus Erec: e.w.*(lg Erec - lg Etrue)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(18212,
      "Xmax versus Erec: e.w.*(lg Erec - lg Etrue)**2",
      "D", 2, xylow, xyhigh, nbins);

   xylow[0]  = 0.;    xylow[1]  = 0.;
   xyhigh[0] = 100.;  xyhigh[1] = 2.;
   nbins[0]  = 100;   nbins[1]  = 200;
   book_histogram(19001, 
      "direction error versus no. of telescopes",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19101, 
      "direction error versus no. of telescopes, passing shape cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19201, 
      "direction error versus no. of telescopes, passing shape+dE cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19301, 
      "direction error versus no. of telescopes, passing shape+dE+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19401, 
      "direction error versus no. of telescopes, passing shape+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19501, 
      "direction error versus no. of telescopes, passing shape+dE+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19601, 
      "direction error versus no. of telescopes, passing shape+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19701, 
      "direction error versus no. of telescopes, passing shape+dE+dE2 cuts",
      "D", 2, xylow, xyhigh, nbins);

   xylow[0]  = -3.;    xylow[1]  = 0.;
   xyhigh[0] = 4.;  xyhigh[1] = 2.;
   nbins[0]  = 140;   nbins[1]  = 200;
   book_histogram(19002, 
      "direction error versus log10(energy)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19102, 
      "direction error versus log10(energy), passing shape cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19202, 
      "direction error versus log10(energy), passing shape+dE cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19302, 
      "direction error versus log10(energy), passing shape+dE+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19402, 
      "direction error versus log10(energy), passing shape+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19502, 
      "direction error versus log10(energy), passing shape+dE+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19602, 
      "direction error versus log10(energy), passing shape+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19702, 
      "direction error versus log10(energy), passing shape+dE+dE2 cuts",
      "D", 2, xylow, xyhigh, nbins);

   xylow[0]  = 0.;     xylow[1]  = 0.;
   xyhigh[0] = 2000.;  xyhigh[1] = 2.;
   nbins[0]  = 200;    nbins[1]  = 200;
   book_histogram(19003, 
      "Direction error versus core distance",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19103, 
      "Direction error versus core distance, passing shape cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19203, 
      "Direction error versus core distance, passing shape+dE cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19303, 
      "Direction error versus core distance, passing shape+dE+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19403, 
      "Direction error versus core distance, passing shape+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19503, 
      "Direction error versus core distance, passing shape+dE+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19603, 
      "Direction error versus core distance, passing shape+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19703, 
      "Direction error versus core distance, passing shape+dE+dE2 cuts",
      "D", 2, xylow, xyhigh, nbins);

   xylow[0]  = -3.;    xylow[1]  = -2.;
   xyhigh[0] = 4.;     xyhigh[1] = 1.;
   nbins[0]  = 140;    nbins[1]  = 150;
   book_histogram(19012, 
      "log10(rec. E/true E) versus log10(true E)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19013, 
      "log10(rec. E/true E) versus log10(rec. E)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19014, 
      "log10(rec. E0/true E) versus log10(rec. E0)",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19112, 
      "log10(rec. E/true E) versus log10(true E), passing shape cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19113, 
      "log10(rec. E/true E) versus log10(rec. E), passing shape cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19114,
      "log10(rec. E0/true E) versus log10(rec. E0), passing shape cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19212, 
      "log10(rec. E/true E) versus log10(true E), passing shape+dE cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19213, 
      "log10(rec. E/true E) versus log10(rec. E), passing shape+dE cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19214,
      "log10(rec. E0/true E) versus log10(rec. E0), passing shape+dE cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19312, 
      "log10(rec. E/true E) versus log10(true E), passing shape+dE+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19313, 
      "log10(rec. E/true E) versus log10(rec. E), passing shape+dE+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19314, 
      "log10(rec. E0/true E) versus log10(rec. E0), passing shape+dE+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19412, 
      "log10(rec. E/true E) versus log10(true E), passing shape+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19413, 
      "log10(rec. E/true E) versus log10(rec. E), passing shape+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19414, 
      "log10(rec. E0/true E) versus log10(rec. E0), passing shape+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19512, 
      "log10(rec. E/true E) versus log10(true E), passing shape+dE+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19513, 
      "log10(rec. E/true E) versus log10(rec. E), passing shape+dE+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19514, 
      "log10(rec. E0/true E) versus log10(rec. E0), passing shape+dE+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19612, 
      "log10(rec. E/true E) versus log10(true E), passing shape+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19613, 
      "log10(rec. E/true E) versus log10(rec. E), passing shape+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19614, 
      "log10(rec. E0/true E) versus log10(rec. E0), passing shape+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19712, 
      "log10(rec. E/true E) versus log10(true E), passing shape+dE+dE2 cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19713, 
      "log10(rec. E/true E) versus log10(rec. E), passing shape+dE+dE2 cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19714, 
      "log10(rec. E0/true E) versus log10(rec. E0), passing shape+dE+dE2 cuts",
      "D", 2, xylow, xyhigh, nbins);

   xylow[0]  = -5.;    xylow[1]  = -5.;
   xyhigh[0] = 5.;     xyhigh[1] = 5.;
   nbins[0]  = 200;    nbins[1]  = 200;
   book_histogram(19530, 
      "True direction (deg) in nominal plane",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19531, 
      "Rec. direction (deg) in nominal plane",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19532, 
      "Rec. direction (deg) in nominal plane, passing amplitude cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19533, 
      "Rec. direction (deg) in nominal plane, passing ampl.+edge cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19534, 
      "Rec. direction (deg) in nominal plane, passing shape cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19535, 
      "Rec. direction (deg) in nominal plane, passing shape+dE cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19536, 
      "Rec. direction (deg) in nominal plane, passing shape+dE+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19537, 
      "Rec. direction (deg) in nominal plane, passing shape+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19538, 
      "Rec. direction (deg) in nominal plane, passing shape+dE+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19539, 
      "Rec. direction (deg) in nominal plane, passing shape+dE2+hmax cuts",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(19540, 
      "Rec. direction (deg) in nominal plane, passing shape+dE+dE2 cuts",
      "D", 2, xylow, xyhigh, nbins);

   init_hist_global = 1;
}

static void book_hist_for_type(AllHessData *hsdata, int tel_type)
{
   double xylow[2], xyhigh[2];
   int nbins[2];
   int ht = 100000 * tel_type;

   if ( tel_type < 0 || tel_type > MAX_TEL_TYPES+1 )
      return;
   if ( init_hist_for_type[tel_type] )
      return;

   printf("Booking histograms for telescope type %d.\n", tel_type);

   xylow[0]  = 0.;    xylow[1]  = 0.;
   xyhigh[0] = 2000.; xyhigh[1] = 6.;
   nbins[0]  = 200;   nbins[1]  = 120;
   book_histogram(ht+18000, 
      "Log10(image ampl.) versus core distance: no. of entries",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18001, 
      "Log10(Image ampl.) versus core distance: event weights",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18011, 
      "Log10(Image ampl.) versus core distance: e.w.*width",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18012, 
      "Log10(Image ampl.) versus core distance: e.w.*width*width",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18021, 
      "Log10(Image ampl.) versus core distance: e.w.*length",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18022, 
      "Log10(Image ampl.) versus core distance: e.w.*length*length",
      "D", 2, xylow, xyhigh, nbins);

   book_histogram(ht+18051, 
      "Log10(Image ampl.) versus core distance: e.w.*Amp/E",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18052, 
      "Log10(Image ampl.) versus core distance: e.w.*(Amp/E)**2",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18061, 
      "Log10(Image ampl.) versus core distance, shape cuts: e.w.*Amp/E",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18062, 
      "Log10(Image ampl.) versus core distance, shape cuts: e.w.*(Amp/E)**2",
      "D", 2, xylow, xyhigh, nbins);

   xylow[0] = 0.; xyhigh[0] = 1.0; nbins[0] = 200.;
   book_histogram(ht+18005,
      "Log10(image ampl.) versus w/l: no. of entries",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18006,
      "Log10(image ampl.) versus w/l: event weights",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18071,
      "Log10(image ampl.) versus w/l: e.w.*rcore",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18072,
      "Log10(image ampl.) versus w/l: e.w.*rcore*rcore",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18081,
      "Log10(image ampl.) versus w/l: e.w.*rimg",
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18082,
      "Log10(image ampl.) versus w/l: e.w.*rimg*rimg",
      "D", 2, xylow, xyhigh, nbins);

   xylow[0]  = 0.;    xylow[1]  = -3;
   xyhigh[0] = 2000.; xyhigh[1] = 3.;
   nbins[0]  = 200;   nbins[1]  = 120;
   book_histogram(ht+18301,
      "lg(rec. E) vs. rec. core distance, any tel. in rec. event", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18302,
      "lg(rec. E) vs. rec. core distance, any tel. after shape cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18311,
      "lg(rec. E) vs. rec. core distance, trig. tel. in rec. event", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18312,
      "lg(rec. E) vs. rec. core distance, trig. tel. after shape cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18321,
      "lg(rec. E) vs. rec. core distance, min. amp. in rec. event", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18322,
      "lg(rec. E) vs. rec. core distance, min. amp. after shape cuts", 
      "D", 2, xylow, xyhigh, nbins);

   xylow[0]  = 0.;    xylow[1]  = -20;
   xyhigh[0] = 1000.; xyhigh[1] = 20.;
   nbins[0]  = 100;   nbins[1]  = 200;
   book_histogram(ht+18401,
      "time slope vs. true core distance, min amp.+edge cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18402,
      "time slope vs. true core distance, shape cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18403,
      "time slope vs. true core distance, shape+angle cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18404,
      "time slope vs. true core distance, shape+angle+dE+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);

   book_histogram(ht+18503,
      "time slope vs. rec. core distance, shape+angle cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18504,
      "time slope vs. rec. core distance, shape+angle+dE+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);

   xylow[0]  = 0.;    xylow[1]  = 0;
   xyhigh[0] = 1000.; xyhigh[1] = 7.5;
   nbins[0]  = 100;   nbins[1]  = 75;
   book_histogram(ht+18411,
      "time residual vs. true core distance, min amp.+edge cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18412,
      "time residual vs. true core distance, shape cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18413,
      "time residual vs. true core distance, shape+angle cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18414,
      "time residual vs. true core distance, shape+angle+dE+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);

   book_histogram(ht+18513,
      "time residual vs. rec. core distance, shape+angle cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18514,
      "time residual vs. rec. core distance, shape+angle+dE+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);

   xylow[0]  = 0.;    xylow[1]  = 0;
   xyhigh[0] = 1000.; xyhigh[1] = 15.;
   nbins[0]  = 100;   nbins[1]  = 75;
   book_histogram(ht+18421,
      "time width 1 vs. true core distance, min amp.+edge cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18422,
      "time width 1 vs. true core distance, shape cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18423,
      "time width 1 vs. true core distance, shape+angle cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18424,
      "time width 1 vs. true core distance, shape+angle+dE+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);

   book_histogram(ht+18523,
      "time width 1 vs. rec. core distance, shape+angle cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18524,
      "time width 1 vs. rec. core distance, shape+angle+dE+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);

   xylow[0]  = 0.;    xylow[1]  = 0;
   xyhigh[0] = 1000.; xyhigh[1] = 15.;
   nbins[0]  = 100;   nbins[1]  = 75;
   book_histogram(ht+18431,
      "time width 2 vs. true core distance, amp.+edge cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18432,
      "time width 2 vs. true core distance, shape cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18433,
      "time width 2 vs. true core distance, shape+angle cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18434,
      "time width 2 vs. true core distance, shape+angle+dE+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);

   book_histogram(ht+18533,
      "time width 2 vs. rec. core distance, shape+angle cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18534,
      "time width 2 vs. rec. core distance, shape+angle+dE+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);

   xylow[0]  = 0.;    xylow[1]  = 0;
   xyhigh[0] = 1000.; xyhigh[1] = 7.5;
   nbins[0]  = 100;   nbins[1]  = 75;
   book_histogram(ht+18441,
      "rise time vs. true core distance, min amp.+edge cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18442,
      "rise time vs. true core distance, shape cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18443,
      "rise time vs. true core distance, shape+angle cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18444,
      "rise time vs. true core distance, shape+angle+dE+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);

   book_histogram(ht+18543,
      "rise time vs. rec. core distance, shape+angle cuts", 
      "D", 2, xylow, xyhigh, nbins);
   book_histogram(ht+18544,
      "rise time vs. rec. core distance, shape+angle+dE+dE2+hmax cuts", 
      "D", 2, xylow, xyhigh, nbins);
 
   init_hist_for_type[tel_type] = 1;
}

/* -------------------------  user_init  ---------------------------- */
/**
 *  @short Initialisation of user analysis, booking of histograms etc.
 */

static void user_init (AllHessData *hsdata)
{
   /* One time initialisations like booking of histograms ... */

   int itel;
   int i;
   int tel_type;

   if ( user_init_done )
   {
      printf("\nUser analysis re-initialisation check\n");
      if ( tel_types_change )
      {
         printf("Telescope type definitions may have changed.\n");
         init_telescope_types(hsdata);
      }
      /* Adding telescope type specific histograms? */
      for ( tel_type=0; tel_type <= MAX_TEL_TYPES+1; tel_type++ )
      {
         if ( stat_type[tel_type] != 0 && init_hist_for_type[tel_type] == 0 )
            book_hist_for_type(hsdata, tel_type);
      }
      /* Nothing else to be done */
      return;
   }
   user_init_done = 1;

   printf("\n");
   printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
   printf("+ User analysis one-time initialisation called. +\n");
   printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");

   user_init_parameters();

   if ( hsdata != NULL && hsdata->mc_shower.primary_id == 0 && 
      up[0].d.theta_scale != 1.0 && hsdata->mc_run_header.viewcone[1] == 0. )
   {
      printf("Theta-cut scale factor forced to 1.0 for point-source gammas!\n");
      up[0].d.theta_scale = 1.0;
   }
   if ( diffuse_mode )
   {
      printf("Using diffuse mode, accepting showers between "
             "%4.2f and %4.2f deg from nominal viewing direction.\n",
         diffuse_off_axis_min*(180./M_PI),
         diffuse_off_axis_max*(180./M_PI));
   }

   Az_src  = (fabs(hsdata->mc_run_header.az_range[1] -
                     hsdata->mc_run_header.az_range[0]) < M_PI) ?
                0.5*(hsdata->mc_run_header.az_range[0] +
                     hsdata->mc_run_header.az_range[1]) :
                (hsdata->mc_run_header.az_range[0] +
                 hsdata->mc_run_header.az_range[1] < M_PI) ?
                0.5*(hsdata->mc_run_header.az_range[0] +
                     hsdata->mc_run_header.az_range[1] + M_PI) :
                0.5*(hsdata->mc_run_header.az_range[0] +
                     hsdata->mc_run_header.az_range[1] - M_PI);
   Alt_src = 0.5*(hsdata->mc_run_header.alt_range[0] + 
                      hsdata->mc_run_header.alt_range[1]);
   Az_nom  = hsdata->run_header.direction[0];
   Alt_nom = hsdata->run_header.direction[1];

   source_offset =angle_between(Az_src, Alt_src, Az_nom, Alt_nom);
   up[0].d.source_offset_deg = source_offset*(180./M_PI);

   printf("Simulated source direction is at Az=%4.2f deg, Alt %4.2f deg.\n",
      Az_src*(180./M_PI), Alt_src*(180./M_PI));
   printf("Nominal viewing direction is at Az=%4.2f deg, Alt %4.2f deg.\n",
      Az_nom*(180./M_PI), Alt_nom*(180./M_PI));
   if ( source_offset > 0. )
      printf("Source offset to nominal viewing direction is %5.3f deg.\n", up[0].d.source_offset_deg);

   init_telescope_types(hsdata);

   histogram_hashing(10000);
   if ( up[0].d.source_offset_deg < 0.01 ) 
      snprintf(lookup_fname, sizeof(lookup_fname),
         "lookups%d_%d-%d_%1.0f_%d_%1.0f_%1.0f.hdata.gz",
         up[0].i.reco_flag, up[0].i.min_tel_img, up[0].i.max_tel_img, 
         up[0].d.min_amp, up[0].i.min_pix,
         up[0].d.tailcut_low, up[0].d.tailcut_high);
   else
      snprintf(lookup_fname, sizeof(lookup_fname),
         "lookups%d_%d-%d_%1.0f_%d_%1.0f_%1.0f_off%4.2f.hdata.gz",
         up[0].i.reco_flag, up[0].i.min_tel_img, up[0].i.max_tel_img, 
         up[0].d.min_amp, up[0].i.min_pix,
         up[0].d.tailcut_low, up[0].d.tailcut_high, up[0].d.source_offset_deg);
   if ( user_lookup_fname[0] != '\0' )
      strncpy(lookup_fname,user_lookup_fname,sizeof(lookup_fname)-1);

   if (  read_histogram_file(lookup_fname,0) < 0 )
   {
      fflush(stdout);
      fprintf(stderr,"\nCannot read lookups file '%s'.\n\n", lookup_fname);
      sleep(2);
   }
   else
   {
      printf("Lookups read from '%s'.\n", lookup_fname);
      if ( up[0].i.user_flags == 0 ) /* not HESS-style */
      {
         int ires;
         HISTOGRAM *h;
         for ( i=0; i<H_MAX_TEL; i++ )
            opt_theta_cut[0][i] = 0.;
         for ( ires=0; ires<7 /* 3 */; ires++ )
         {
            if ( (h = get_histogram_by_ident(18911+ires)) != NULL )
            {
               for ( i=0; i<H_MAX_TEL && i<h->nbins; i++ )
	       {
                  opt_theta_cut[ires][i] = 
                     h->extension->ddata[i] * (M_PI/180.) * up[0].d.theta_scale;
		  if ( opt_theta_cut[ires][i] > max_theta * up[0].d.theta_scale && max_theta > 0. )
		     opt_theta_cut[ires][i] = max_theta * up[0].d.theta_scale;
		  if ( opt_theta_cut[ires][i] < min_theta * up[0].d.theta_scale )
		     opt_theta_cut[ires][i] = min_theta * up[0].d.theta_scale;
	       }
            }
            else if ( ires > 0 )
            {
               for ( i=0; i<H_MAX_TEL && i<h->nbins; i++ )
                  opt_theta_cut[ires][i] = opt_theta_cut[ires-1][i];
            }
         }
      }
   }
   {
      int ires;
      for ( ires=0; ires<7 /* 3 */; ires++ )
         for ( itel=0; itel<H_MAX_TEL; itel++ )
            if ( opt_theta_cut[ires][itel] <= 0. )
               opt_theta_cut[ires][itel] = max_theta * up[0].d.theta_scale;
   }
   set_ebias_correction(get_histogram_by_ident(19119)); /* If available */

#if 0
   for ( itel=0; itel<20; itel++ )
     printf("+++ ntel=%d: opt_theta_cut = %f, %f, %f deg after scaling by %f.\n",  itel,
      opt_theta_cut[0][itel]*(180./M_PI),  
      opt_theta_cut[1][itel]*(180./M_PI), 
      opt_theta_cut[2][itel]*(180./M_PI),
      up[0].d.theta_scale);
#endif

   book_hist_global(hsdata);

   /* Telescope type specific histograms */

   for ( tel_type=0; tel_type <= MAX_TEL_TYPES+1; tel_type++ )
   {

      if ( stat_type[tel_type] == 0 )
      {
         if ( tel_type <= MAX_TEL_TYPES )
            printf("No telescopes of type %d.\n", tel_type);
         continue; /* No telescopes of this type, so no need for histograms. */
      }
      book_hist_for_type(hsdata, tel_type);
   } /* End of telescope type specific histograms */

   pixmom = alloc_moments(0.,3000.);
}

/* ------------------------ user_mc_shower_fill ------------------------ */
/** Work to be done once per generated shower. */

static void user_mc_shower_fill (AllHessData *hsdata /* currently unused */)
{
   // double ewt = pow(hsdata->mc_shower.energy, d_sp_idx);
}

/* ---------------------- user_mc_event_fill ------------------------- */
/** 
 *  @short Work to be done once per shower usage.
 *
 *  Depending on sim_hessarray flags this might be called only 
 *    for triggered events or also for non-triggered events (default). 
 */

static void user_mc_event_fill (AllHessData *hsdata)
{
   double ewt = pow(hsdata->mc_shower.energy, up[0].d.d_sp_idx);

   /* Example code: */
      double rs;
      rs = line_point_distance(hsdata->mc_event.xcore, 
                     hsdata->mc_event.ycore, 0.,
                     cos(hsdata->mc_shower.altitude) *
                        cos(hsdata->mc_shower.azimuth),
                     cos(hsdata->mc_shower.altitude) *
                        sin(-hsdata->mc_shower.azimuth),
                     sin(hsdata->mc_shower.altitude),
                     0., 0., 0.);
      fill_histogram_by_ident(12001,rs,log10(hsdata->mc_shower.energy),ewt);
      fill_histogram_by_ident(22001,log10(hsdata->mc_shower.energy),0.,1.0);
      fill_histogram_by_ident(22101,log10(hsdata->mc_shower.energy),0.,ewt);
   /* - */
}

extern struct basic_ntuple bnt;

/* ------------------------- user_event_fill ------------------------- */
/** Fill (triggered) event specific histograms etc. */

static void user_event_fill (AllHessData *hsdata, int stage)
{
   double E_true = hsdata->mc_shower.energy;    /**< true energy [TeV] */
   if ( E_true <= 0. )
   {
      fflush(stdout);
      fprintf(stderr,"\nBad MC energy: %f\n", E_true);
      return;
   }
   double lg_E_true = log10(E_true);
   double ewt = pow(E_true, up[0].d.d_sp_idx);          /**< Event for desired spectral slope */
   double Az_true = hsdata->mc_shower.azimuth;
   double Alt_true= hsdata->mc_shower.altitude;
   double xc_true = hsdata->mc_event.xcore;
   double yc_true = hsdata->mc_event.ycore;
   double v_ts[H_MAX_TEL];                 /**< true core distance [m] */

   double energy = 0., lg_energy = -10., lg_energy0 = -10.;
   double Az = 0., Alt = 0.;
   double xc = 9999., yc = 9999.;
   double v_tr[H_MAX_TEL];                 /**< reconstructed core distance [m] */
   // double v_tel_az[H_MAX_TEL];             /**< telescope azimuth */
   // double v_tel_alt[H_MAX_TEL];            /**< telescope altitude */
   double v_amp[H_MAX_TEL];                /**< image amplitude [peak p.e.] */
   double v_w[H_MAX_TEL];                  /**< image width [rad] */
   double v_l[H_MAX_TEL];                  /**< image length [rad] */
   double v_rcog[H_MAX_TEL];               /**< radius of image c.o.g. in camera plane */
   double v_dimg[H_MAX_TEL];               /**< distance of image c.o.g. to source [rad] */
   int img_ok[H_MAX_TEL];                  /**< Amplitude and edge distance are ok */
   int sc_ok[H_MAX_TEL];
   int alt_ok = 0;
   
   /* Example code: */
      double rs, rr=9999., da;
      int itel, ntel=0;
      double w_sc = 0., s_scrw = 0., s_scrl = 0., s_scrw2 = 0., s_scrl2 = 0.;
      double mscrw = 999., mscrl = 999., mscrw_s = 999., mscrl_s = 999.;
      double s_sce = 0., w_sce = 0., eresol = 999.;
      double s_sce1 = 0., s_sce2 = 0., w_sce1 = 0., echi2 = 999.;
      double stheta, ctheta, x_nom = 999., y_nom = 999.; 
      int n_sc = 0, n_img = 0, n_img2 = 0, n_img3 = 0, npix_tot = 0;
      int shape_cuts_ok = 0, angle_cut_ok = 0, eres_cut_ok = 0;
      int icut, itp;
      int angle_cutx_ok[7] = {0, 0, 0, 0, 0, 0, 0};
      int eres2_cut_ok = 0, hmax_cut_ok = 0;
      double w_hmax = 0., s_hmax = 0., s_hmax2 = 0., hmax = 99999., hmax_err = 9999.;
      // double xsrc_nom = 0., ysrc_nom = 0.;
      double xrec = 999., yrec = 999.;
      double tr2s = 0., mdisp = 0., tss = 0., tsw = 0.;
      double r_lc = 125., r_lc_max = -1.;
      int num_trg_type[11];
      int ntp = sizeof(num_trg_type)/sizeof(num_trg_type[0]);

      static int have_atm_set = -1; /* For height to Xmax lookup */
      int num_trg = 0;

      if ( stage != 1 ) /* Only processing after our own shower reconstruction */
         return;
      
      if ( diffuse_mode )
      {
         double dn = angle_between(Az_nom, Alt_nom, Az_true, Alt_true);
         if ( dn < diffuse_off_axis_min || dn > diffuse_off_axis_max )
            return;
      }

      memset(&bnt,'\0',sizeof(bnt));
      bnt.primary = hsdata->mc_shower.primary_id;
      bnt.run = hsdata->run_header.run;
      bnt.event = hsdata->mc_event.event;
      bnt.lg_e_true = lg_E_true;
      bnt.weight = ewt;
      bnt.xfirst_true = thickx(hsdata->mc_shower.h_first_int) / sin(Alt_true);
      bnt.xmax_true = hsdata->mc_shower.xmax;
      bnt.xc_true = xc_true;
      bnt.yc_true = yc_true;
      bnt.az_true = Az_true;
      bnt.alt_true = Alt_true;
      
      for ( itp=0; itp<ntp; itp++ )
         num_trg_type[itp] = 0;

      for ( itel=0; itel<hsdata->run_header.ntel; itel++)
         if ( hsdata->event.teldata[itel].known )
	 {
	    num_trg++; // Actually rather the number of telescopes with data
            fill_histogram_by_ident(10002, hsdata->event.teldata[itel].tel_id, 0., ewt);
            itp = telescope_type[itel];
            if ( itp >= 0 && itp < ntp )
               num_trg_type[itp]++;
	 }

      if ( num_trg > 0 )
      {
         int itrg, ie = -1;
         for ( itel=0; itel<hsdata->run_header.ntel; itel++)
            fill_histogram_by_ident(10002, hsdata->event.teldata[itel].tel_id, num_trg, ewt);
         if ( E_true > 0.01 )
            ie = (int)(10.*log10(E_true/0.01)/3.) + 1; /* Ten bins per three decades in energy */
         if ( ie > 15 )
            ie = 15;
         for ( itrg=1; itrg<=num_trg && itrg<=4; itrg++ )
         {
            fill_histogram_by_ident(9000+100*itrg, xc_true, yc_true, ewt);
            if ( ie > 0 )
               fill_histogram_by_ident(9000+100*itrg+ie, xc_true, yc_true, ewt);
         }
      }

      fill_histogram_by_ident(10001, num_trg, 0., ewt);
      fill_histogram_by_ident(10004, num_trg, 0., 1.0);
      fill_histogram_by_ident(10003, num_trg, lg_E_true, ewt);
      for ( itp=0; itp<ntp; itp++ )
         fill_histogram_by_ident(10005, num_trg_type[itp], itp, ewt);

      /* True core offset w.r.t. array centre */
      rs = line_point_distance(xc_true, yc_true, 0.,
                     cos(Alt_true) * cos(Az_true),
                     cos(Alt_true) * sin(-Az_true),
                     sin(Alt_true),
                     0., 0., 0.);
      fill_histogram_by_ident(12002, rs, lg_E_true, ewt);
      fill_histogram_by_ident(22002, lg_E_true, 0., 1.0);
      fill_histogram_by_ident(22102, lg_E_true, 0., ewt);

      angles_to_offset(Az_true, Alt_true,
         Az_nom, Alt_nom, 180./M_PI, &x_nom, &y_nom); // In degrees!
      fill_histogram_by_ident(19530, x_nom, y_nom, ewt);

      /* Any limits in true impact position? */
      if ( up[0].d.true_impact_range[0] > 0. && rs > up[0].d.true_impact_range[0] )
         return;
      if ( up[0].d.true_impact_range[1] > 0. && fabs(xc_true) > up[0].d.true_impact_range[1] )
         return;
      if ( up[0].d.true_impact_range[2] > 0. && fabs(yc_true) > up[0].d.true_impact_range[2] )
         return;

      if ( hsdata->event.shower.known )
      {
         if ( (hsdata->event.shower.result_bits & 0x01) == 0x01 )
         {
            Az = hsdata->event.shower.Az;
            Alt = hsdata->event.shower.Alt;
         }
         else
            Az = Alt = -1.;
         if ( (hsdata->event.shower.result_bits & 0x04) == 0x04 )
         {
            xc = hsdata->event.shower.xc;
            yc = hsdata->event.shower.yc;
         }
         else
            xc = yc = 99999.;
      }

      bnt.xc = xc;
      bnt.yc = yc;
      bnt.az = Az;
      bnt.alt = Alt;

      for (itel=0; itel<hsdata->run_header.ntel; itel++)
      {
         double tel_az = hsdata->event.trackdata[itel].cor_known ?
                         hsdata->event.trackdata[itel].azimuth_cor :
                         hsdata->event.trackdata[itel].azimuth_raw;
         double tel_alt= hsdata->event.trackdata[itel].cor_known ?
                         hsdata->event.trackdata[itel].altitude_cor :
                         hsdata->event.trackdata[itel].altitude_raw;
         int j_img = 0; // For now always using the first set of image parameters.
         int ttyp = telescope_type[itel]; // Telescope type (1...n, or 0 as fallback)
         // v_tel_az[itel] = tel_az;
         // v_tel_alt[itel] = tel_alt;
         v_amp[itel]= 0.;
         v_rcog[itel] = 999.;
         img_ok[itel] = sc_ok[itel] = 0;
         v_dimg[itel] = 999.;

         if ( !hsdata->event.teldata[itel].known )
            continue;

         if ( hsdata->event.teldata[itel].num_image_sets > j_img &&
              hsdata->event.teldata[itel].img != NULL &&
              hsdata->event.teldata[itel].img[j_img].known &&
              hsdata->event.teldata[itel].img[j_img].amplitude > 0.)
         {
            /* Position of image c.o.g. in camera plane */
            double ximg = hsdata->event.teldata[itel].img[j_img].x;
            double yimg = hsdata->event.teldata[itel].img[j_img].y;
            ntel++;
            /* Distance of image c.o.g. from camera centre */
            v_rcog[itel] = sqrt(ximg*ximg+yimg*yimg);
            /* Image amplitude ('size') [p.e.] */
            v_amp[itel]  = hsdata->event.teldata[itel].img[j_img].amplitude;
            /* Image width */
            v_w[itel]    = hsdata->event.teldata[itel].img[j_img].w;
            /* Image length */
            v_l[itel]    = hsdata->event.teldata[itel].img[j_img].l;
            /* Distance of telescope from true core position */
            v_ts[itel]   = line_point_distance(xc_true, yc_true, 0.,
                     cos(Alt_true) * cos(Az_true),
                     cos(Alt_true) * sin(-Az_true),
                     sin(Alt_true),
                     hsdata->run_header.tel_pos[itel][0],
                     hsdata->run_header.tel_pos[itel][1],
                     hsdata->run_header.tel_pos[itel][2]);
            v_tr[itel]   = 9999.;
            if ( hsdata->event.shower.known && 
                 (hsdata->event.shower.result_bits & 0x05) == 0x05 )
               /* Distance of telescope from reconstructed core position */
               v_tr[itel]= line_point_distance(xc, yc, 0.,
                        cos(Alt) * cos(Az),
                        cos(Alt) * sin(-Az),
                        sin(Alt),
                        hsdata->run_header.tel_pos[itel][0],
                        hsdata->run_header.tel_pos[itel][1],
                        hsdata->run_header.tel_pos[itel][2]);
            if ( v_amp[itel] >= up[ttyp].d.min_amp && 
                 hsdata->event.teldata[itel].img[j_img].pixels >= up[ttyp].i.min_pix )
            {
               n_img++; /* Images with enough pixels and big enough amplitude */
               //if ( v_rcog[itel]+v_l[itel] <= 0.85*get_camera_radius(itel,0) )
               if ( ( up[ttyp].i.user_flags==0 && 
                      v_rcog[itel]+0.35*v_l[itel] <= 0.82*get_camera_radius(itel,0) ) ||
                    ( up[ttyp].i.user_flags > 0 && /* HESS-style */
                      v_rcog[itel] <= 0.764*get_camera_radius(itel,0) ) )
               {
                  if ( up[ttyp].d.max_core_distance > 0. && v_tr[itel] > up[ttyp].d.max_core_distance )
                     continue;
                  n_img2++; /* Images far enough from camera edge */
                  img_ok[itel]= 1;
                  if ( hsdata->event.shower.known && 
                       (hsdata->event.shower.result_bits & 0x05) == 0x05 )
                  {
                     double dimg;
                     double ww;
                     angles_to_offset(Az, Alt, tel_az, tel_alt, 1., &xrec, &yrec);
                     /* Distance between image c.o.g. and reconstructed shower direction
                        in the camera plane. */
                     dimg = sqrt((xrec-ximg)*(xrec-ximg)+(yrec-yimg)*(yrec-yimg));
                     ww = v_amp[itel] / (v_l[itel]*v_l[itel]) * dimg*dimg/(0.01+dimg)/(0.01+dimg) * v_tr[itel]/(125.+v_tr[itel]);
                     v_dimg[itel] = dimg;
                     s_hmax += ww * v_tr[itel] / dimg;
                     s_hmax2+= ww * (v_tr[itel] / dimg)*(v_tr[itel] / dimg);
                     w_hmax += ww;
                     n_img3++; /* Usable image in reconstructed showers */
                  }
               }
            }
         }
      }

      if ( n_img > 0 )
      {
         int iimg, ie = -1;
         if ( E_true > 0.01 )
            ie = (int)(10.*log10(E_true/0.01)/3.) + 1; /* Ten bins per three decades in energy */
         if ( ie > 15 )
            ie = 15;
         for ( iimg=1; iimg<=n_img && iimg<=4; iimg++ )
         {
            fill_histogram_by_ident(9050+100*iimg, xc_true, yc_true, ewt);
            if ( ie > 0 )
               fill_histogram_by_ident(9050+100*iimg+ie, xc_true, yc_true, ewt);
         }
      }

      hsdata->event.shower.result_bits &= 0x0f; /* We want to (re-)fill mean shape, energy, and Xmax parameters */
                                                /* but will keep direction and core position. */

      if ( n_img3 >= 2 && w_hmax > 0. )
      {
         /* Beware of rounding errors in sqrt argument */
         double s1 = s_hmax/w_hmax;
         double s2 = s_hmax2/w_hmax;
         double sd = s2-s1*s1;
         hmax = s1 * sin(Alt) + hsdata->mc_run_header.obsheight; /* Assume plane-parallel atmosphere */
         if ( sd > 0. )
            hmax_err = sqrt(sd) * sqrt(1./(n_img3-1.0)) * sin(Alt);
         else
            hmax_err = 99999.;
         if ( have_atm_set != hsdata->mc_run_header.atmosphere )
         {
            init_atmprof(hsdata->mc_run_header.atmosphere);
            have_atm_set = hsdata->mc_run_header.atmosphere;
         }
         hsdata->event.shower.xmax = thickx(hmax); /* Result is in g/cm^2 */
         hsdata->event.shower.err_xmax = 0.5*(thickx(hmax-hmax_err)-thickx(hmax+hmax_err));
         if ( verbosity > 0 )
            printf("Xmax = %5.2f +- %4.2f g/cm^2 (Hmax = %5.1f m) from %d images\n", 
               hsdata->event.shower.xmax, hsdata->event.shower.err_xmax, hmax, n_img3);
         if ( hsdata->event.shower.xmax > 0. )
            hsdata->event.shower.result_bits |= 0x300; /* We have Xmax and its error estimate */
         r_lc = (s_hmax/w_hmax) * sqrt(1.-1./pow(refidx(hmax),2.));
         if ( verbosity > 0 )
         {
            if ( r_lc_max < 0. )
            {
               int ih;
               for ( ih=2; ih<30; ih++ )
               {
                  double h = hsdata->mc_run_header.obsheight + ih*1000.;
                  double n = refidx(h);
                  double r = h/sin(Alt) * sqrt(1.-1./(n*n));
                  if ( r > r_lc_max )
                     r_lc_max = r;
               }
            }
            printf("Light cone has a radius of %4.2f m (max: %4.2f m)\n", r_lc, r_lc_max);
         }
      }

      bnt.xmax = hsdata->event.shower.xmax;
      bnt.sig_xmax = hsdata->event.shower.err_xmax * ((n_img3>1)?sqrt(n_img3-1.0):1.0);
      bnt.n_tsl0 = 0;
      bnt.acceptance = 0;

      for (itel=0; itel<hsdata->run_header.ntel; itel++)
      {
         int ht = 100000 * telescope_type[itel]; // Used for histograms
         int ttyp = telescope_type[itel];  // For lookup of user parameters

         int j_img = 0; // Like before, using only first image set.
         if ( v_amp[itel] > 0.)
         /* using only first image set */
         {
            double amp = v_amp[itel];
            double lg10_amp = (amp>0.) ? log10(amp) : -999.;
            double w = v_w[itel];
            double l = v_l[itel];
            double wol = (l>0.) ? w/l : -1.;
            double rcog = v_rcog[itel];
            double ts = v_ts[itel];
            double tr = v_tr[itel];
            double dimg = v_dimg[itel];
            double scrw = 999., scrl = 999., scw = 999., scl = 999.;
            double sce = 0., scer = 999.;
            double rwol = 0., rwol_err = 9999.; /* core distance estimated from w/l */
            double dwol = 0., dwol_err = 999.;  /* image distance estimated from w/l */

            if ( up[ttyp].i.user_flags > 0 ) /* HESS-style analysis */
            {
               if ( l <= 0.|| w < 0. || rcog > 0.764*get_camera_radius(itel,0) )
                  continue;
            }
            else
            {
               //if ( l <= 0.|| w < 0. || rcog+l > 0.85*get_camera_radius(itel,0) )
               if ( l <= 0.|| w < 0. || rcog+0.35*l > 0.82*get_camera_radius(itel,0) )
                  continue;
            }
            if ( up[ttyp].d.max_core_distance > 0. && v_tr[itel] > up[ttyp].d.max_core_distance )
               continue;

            /* Fill histograms for lookup table with true core distance. */
            if ( img_ok[itel] )
            {
               ImgData *img = &hsdata->event.teldata[itel].img[j_img];

               fill_histogram_by_ident(ht+18000,ts,lg10_amp,1.);
               fill_histogram_by_ident(ht+18001,ts,lg10_amp,ewt);
               fill_histogram_by_ident(ht+18011,ts,lg10_amp,ewt*w);
               fill_histogram_by_ident(ht+18012,ts,lg10_amp,ewt*w*w);
               fill_histogram_by_ident(ht+18021,ts,lg10_amp,ewt*l);
               fill_histogram_by_ident(ht+18022,ts,lg10_amp,ewt*l*l);
               fill_histogram_by_ident(ht+18051,ts,lg10_amp,ewt*amp/E_true);
               fill_histogram_by_ident(ht+18052,ts,lg10_amp,ewt*amp/E_true*amp/E_true);

               fill_histogram_by_ident(ht+18005,wol,lg10_amp,1.);
               fill_histogram_by_ident(ht+18006,wol,lg10_amp,ewt);
               fill_histogram_by_ident(ht+18071,wol,lg10_amp,ewt*ts);
               fill_histogram_by_ident(ht+18072,wol,lg10_amp,ewt*ts*ts);
               fill_histogram_by_ident(ht+18081,wol,lg10_amp,ewt*dimg);
               fill_histogram_by_ident(ht+18082,wol,lg10_amp,ewt*dimg*dimg);

               if ( img->tm_slope != 0. || img->tm_residual != 0. || img->tm_width1 == 0. )
               {
                  double xsh, ysh;
                  double tm_slope_deg = img->tm_slope * (M_PI/180.);
                  angles_to_offset(hsdata->event.shower.Az, hsdata->event.shower.Alt,
                     hsdata->event.trackdata[itel].azimuth_raw,
                     hsdata->event.trackdata[itel].altitude_raw, 1.0,
                     &xsh, &ysh);
                  if ( (xsh - hsdata->event.teldata[itel].img[j_img].x)*
                          cos(hsdata->event.teldata[itel].img[j_img].phi) +
                       (ysh - hsdata->event.teldata[itel].img[j_img].y)*
                          sin(hsdata->event.teldata[itel].img[j_img].phi) > 0. )
                     tm_slope_deg *= -1.;
                  /* FIXME: has to be generalized for different zenith angles. */
                  if ( tr > 200. && fabs(tm_slope_deg) < 1.5 )
                     bnt.n_tsl0++;
                  if ( tr > 50. )
                  {
                     double wt = pow(tr/(tr+100.),2.);
                     tss += pow(fabs(tm_slope_deg+4.82),1.21)*(125./tr) * wt;
                     tsw += wt;
                  }

                  fill_histogram_by_ident(ht+18401,ts,tm_slope_deg,ewt);
                  fill_histogram_by_ident(ht+18411,ts,img->tm_residual,ewt);
                  fill_histogram_by_ident(ht+18421,ts,img->tm_width1,ewt);
                  fill_histogram_by_ident(ht+18431,ts,img->tm_width2,ewt);
                  fill_histogram_by_ident(ht+18441,ts,img->tm_rise,ewt);
               }
            }

            /* All cut-related parameters use reconstructed core distance. */
            if ( img_norm(w, l, amp, lg10_amp, tr, telescope_type[itel],
                          &scrw, &scrl, &scw, &scl, &sce, &scer,
                          &rwol, &rwol_err, &dwol, &dwol_err) == 1 )
            {
               if ( amp > up[ttyp].d.min_amp && scrw < 990 && scrl < 990 )
               {
                  double ww = 1.;//sqrt(amp); // Not intuitive but 'sqrt(amp)' performs better than 'amp'.
                  sc_ok[itel] = 1;
                  /* Don't take individual odd images (e.g. with afterpulse
                     pixels somewhere in the camera) too seriously. */
                  if ( scrw > 10. )
                     scrw = 10.;
                  if ( scrl > 10. )
                     scrl = 10.;
                  n_sc++;
                  w_sc += ww;
                  s_scrw += scrw * ww;
                  s_scrw2+= scrw*scrw * ww;
                  s_scrl += scrl * ww;
                  s_scrl2+= scrl*scrl * ww;
                  if ( sce > 0. )
                  {
                     /* Weighting tries to account for to some extent for */
                     /* telescope-to-telescope fluctuations as well as for */
                     /* correlated fluctuations not fully represented in the */
                     /* look-ups. */
                     double wm, wchi2;

                     /* For the mean log of energy we weight according to */
                     /* relative error on energy, with systematics term added. */
                     wm = 1. / (0.01+(scer*scer));
                     s_sce  += log(sce) * wm;
                     w_sce  += wm;

#ifdef OLD_DE2
                     /* For the chi^2 estimate, we used to do this on */
                     /* absolute energy, with larger systematics term added. */
                     wchi2 = 1. / ((0.04+scer*scer) * (sce*sce));
                     s_sce1 += sce * wchi2;
                     s_sce2 += sce*sce * wchi2;
                     w_sce1 += wchi2;
#else
                     wchi2 = wm; // Use the same weighting as for energy estimate now
                     s_sce1 += log(sce) * wchi2;;
                     s_sce2 += log(sce)*log(sce) * wchi2;
                     w_sce1 += wchi2;
#endif
                  }
                  if ( img_ok[itel] )
                     npix_tot += hsdata->event.teldata[itel].img[j_img].pixels;
                  tr2s += tr * ww;
                  mdisp += (1.-v_w[itel]/v_l[itel]) * ww;
               }
               fill_histogram_by_ident(17500,scrw,scrl,ewt*sqrt(amp));

               if ( verbosity > 0 )
                  printf("Image with amplitude %f p.e. at core distance %3.1f (%3.1f (%3.1f+-%3.1f)) m:\n"
                         "    image c.o.g. distance to source = %f (%f+-%f)\n"
                         "    width = %f deg, scaled width = %f, scaled reduced width = %f\n"
                         "    length = %f deg, scaled length = %f, scaled reduced length = %f\n"
                         "    energy = %f +- %f TeV (instead of %f TeV)\n",
                         amp, ts, tr, rwol, rwol_err, 
                         dimg, dwol, dwol_err,
                         w*(180./M_PI), scw, scrw, 
                         l*(180./M_PI), scl, scrl,
                         sce, sce*scer, E_true); 
            }
         }
      }

      if ( n_sc > 0 && n_sc < up[0].i.min_tel_img )
      {
         /* Any alternate selection lists defined? */
         if ( alt_list != NULL )
         {
            size_t ialt, jtel, ktel;
            for (ialt=0; ialt<n_list; ialt++)
            {
               if ( (size_t)n_sc < alt_list[ialt].min_tel )
                  continue;
               ktel = 0;
               for (itel=0; itel<hsdata->run_header.ntel; itel++)
               {
                  if ( img_ok[itel] && sc_ok[itel] )
                  {
                     for (jtel=0; jtel<alt_list[ialt].ntel; jtel++)
                        if ( hsdata->event.teldata[itel].tel_id ==
                             alt_list[ialt].tel_id[jtel] )
                        {
                           ktel++;
                           break;
                        }
                  }
               }
               if ( ktel >= alt_list[ialt].min_tel )
               {
                  alt_ok = 1;
                  if ( verbosity > 0 )
                     printf("Alternate telescope selection %zu satisfied (%zu of %zu telescopes, with at least %zu required).\n", 
                        ialt, ktel, alt_list[ialt].ntel, alt_list[ialt].min_tel );
                  break;
               }
            }
         }
         /* Now try if lower multiplicities are required for any specific type. */
         if ( ! alt_ok )
         {
            int itype;
            size_t ctype[MAX_TEL_TYPES+1];
            int try_by_type = 0;
            for ( itype=1; itype <= MAX_TEL_TYPES; itype++ )
            {
               ctype[itype] = 0;
               if ( up[itype].i.min_tel_img < up[0].i.min_tel_img )
                  try_by_type = 1;
            }
            if ( try_by_type ) /* Is it worth accumulating the stats? */
            {
               for (itel=0; itel<hsdata->run_header.ntel; itel++)
               {
                  if ( img_ok[itel] && sc_ok[itel] && 
                       saved_tel_type[itel] >= 1 && 
                       saved_tel_type[itel] <= MAX_TEL_TYPES )
                     ctype[saved_tel_type[itel]]++;
               }
               for ( itype=1; itype <= MAX_TEL_TYPES; itype++ )
               {
                  if ( ctype[itype] >= (size_t) up[itype].i.min_tel_img )
                  {
                     alt_ok = 1;
                     if ( verbosity > 0 )
                        printf("Telescope selection satisfied for telescope type %d (%zu telescopes, with at least %d required).\n", 
                           itype, ctype[itype], up[itype].i.min_tel_img );
                    break;
                  }
               }
            }
         }
      }

      if ( (n_sc >= up[0].i.min_tel_img || alt_ok) && n_sc <= up[0].i.max_tel_img )
      {
         int ir, ie;
         
         mscrw = s_scrw / w_sc;
         mscrw_s = sqrt(s_scrw2/w_sc - (s_scrw/w_sc)*(s_scrw/w_sc)) * sqrt(n_sc/(n_sc-0.9999));
         mscrl = s_scrl / w_sc;
         mscrl_s = sqrt(s_scrl2/w_sc - (s_scrl/w_sc)*(s_scrl/w_sc)) * sqrt(n_sc/(n_sc-0.9999));

         if ( up[0].i.user_flags > 0 ) /* HESS-style */
         {
            /* No re-scaling used in HESS standard analysis */
         }
         else
         {
            /* Multiplicity dependent re-scaling results in almost 
               energy-independent shape cut efficiencies.*/
            double scw = sqrt(0.45+0.55*n_sc), scl = sqrt(0.57+0.43*n_sc);
            mscrw *= scw;
            mscrw_s *= scw;
            mscrl *= scl;
            mscrl_s *= scl;
         }

         hsdata->event.shower.mscl = mscrl; /* Use mean reduced scaled length as shape parameter. */
         hsdata->event.shower.mscw = mscrw; /* Use mean reduced scaled width as shape parameter. */
         hsdata->event.shower.result_bits |= 0x10; /* Mean (reduced) scaled width & length are known. */

         fill_histogram_by_ident(17000,mscrw,mscrl,ewt);
         ir = (int)(rs/50.);
         if ( ir >= 0 && ir < 10 )
            fill_histogram_by_ident(17100+ir,mscrw,mscrl,ewt);
         ie = (int)(floor(log10(E_true/0.01)*3.)+0.001);
         if ( ie >= 0 && ie < 12 )
            fill_histogram_by_ident(17200+ie,mscrw,mscrl,ewt);
         fill_histogram_by_ident(17300+ntel,mscrw,mscrl,ewt);

         if ( w_sce > 0. )
         {
            energy = exp(s_sce / w_sce);
            if ( energy > 0. )
            {
               lg_energy0 = log10(energy);
               lg_energy = lg_energy0 - ebias_correction(lg_energy0);
            }
            else
               lg_energy = lg_energy0 = -10.;
            
            eresol = 1./sqrt(w_sce);
            echi2  = (s_sce2/w_sce1 - (s_sce1/w_sce1)*(s_sce1/w_sce1)) 
                     * (double)n_sc / (n_sc - 0.9999);
            if ( echi2 < eval_cut_param(up[0].d.de2_cut_param,lg_energy) )
               eres2_cut_ok = 1;
         }
         else
         {
            energy = -10.;
            lg_energy = lg_energy0 = -10.;
            eresol = 999.;
            echi2  = 999.;
         }
      }
      if ( w_sc > 0. )
         mdisp /= w_sc;

      bnt.weight = ewt;
      bnt.rcm = (w_sc > 0.) ? (tr2s/w_sc) : 9999.;
      bnt.mdisp = mdisp;
      bnt.mscrw = mscrw;
      bnt.sig_mscrw = mscrw_s;
      bnt.mscrl = mscrl;
      bnt.sig_mscrl = mscrl_s;
      bnt.lg_e = lg_energy;
      bnt.sig_e = eresol;
      bnt.chi2_e = echi2;
      bnt.n_img = n_sc;
      bnt.n_trg = hsdata->event.central.num_teltrg;
      bnt.n_fail = 0; /* FIXME */
      bnt.n_pix = npix_tot;
      if ( tsw > 0. )
         bnt.tslope = tss/tsw;
      else
         bnt.tslope = 0.;
      bnt.tsphere = 0.;

      hsdata->event.shower.energy = energy;
      hsdata->event.shower.err_energy = energy*eresol;
      if ( energy > 0. )
         hsdata->event.shower.result_bits |= 0xc0; /* Both energy and error defined. */

      if ( energy > 0. && 
           eresol < eval_cut_param(up[0].d.eres_cut_param,lg_energy) * 
#if 1
              ( 0.20 - 0.09 * ( (lg_energy<0.0) ? lg_energy : tanh(lg_energy) ) ) )
#else
              ( ((energy>0.1) ? 0.25 - 0.05*lg_energy : 0.20 - 0.10*lg_energy) ) )
#endif
         eres_cut_ok = 1;

      if ( verbosity > 0 )
         printf("E=%f: Expected dE=%f < %f (%s); dE2=%f < %f (%s)\n",
	    energy, eresol, eval_cut_param(up[0].d.eres_cut_param,lg_energy) * 
              ( ((energy>0.1) ? 0.25 - 0.05*lg_energy :
                                0.20 - 0.10*lg_energy) ),
	    eres_cut_ok ? "yes" : "no",
	    echi2, eval_cut_param(up[0].d.de2_cut_param,lg_energy),
	    eres2_cut_ok ? "yes" : "no");

      if ( n_img3 >= 2 && w_hmax > 0. )
      {
         double hmax_exp = 
            expected_max_distance(energy, M_PI/2.-Alt,
               hsdata->mc_run_header.obsheight);
         hmax = s_hmax / w_hmax; /* It is a distance this time and not a height */
         if ( hmax > (1.0-0.1*up[0].d.hmax_cut_param)*hmax_exp-(1500.+500.*lg_energy+(lg_energy>0.?300.*lg_energy:0.)) && 
              hmax < (1.1+0.2*up[0].d.hmax_cut_param)*hmax_exp+(1500.+500.*lg_energy+(lg_energy>0.?300.*lg_energy:0.)) )
            hmax_cut_ok = 1;
      }
      else
      {
         hmax = -1.;
         hmax_cut_ok = 0;
      }

      if ( n_img >= up[0].i.min_tel_img && n_img <= up[0].i.max_tel_img )
      {
         fill_histogram_by_ident(12003,rs,lg_E_true,ewt);
         fill_histogram_by_ident(22003, lg_E_true, 0., 1.0);
         fill_histogram_by_ident(22103, lg_E_true, 0., ewt);
      }
      if ( n_img2 >= up[0].i.min_tel_img && n_img2 <= up[0].i.max_tel_img )
      {
         fill_histogram_by_ident(12004,rs,lg_E_true,ewt);
         fill_histogram_by_ident(22004, lg_E_true, 0., 1.0);
         fill_histogram_by_ident(22104, lg_E_true, 0., ewt);
         fill_histogram_by_ident(15001,hmax,lg_E_true,ewt);
      }

      /* In the following we need reconstructed shower parameters. */
      if ( !hsdata->event.shower.known || 
           (hsdata->event.shower.result_bits & 0x05) != 0x05 )
         return;

      /* Reconstructed  core offset w.r.t. array centre */
      rr = line_point_distance(xc, yc, 0.,
              cos(Alt)*cos(Az), cos(Alt)*sin(-Az), sin(Alt),
              0., 0., 0.);

      /* Any limits in reconstructed impact position? */
      if ( up[0].d.impact_range[0] > 0. && rr > up[0].d.impact_range[0] )
         return;
      if ( up[0].d.impact_range[1] > 0. && fabs(xc) > up[0].d.impact_range[1] )
         return;
      if ( up[0].d.impact_range[2] > 0. && fabs(yc) > up[0].d.impact_range[2] )
         return;

      angles_to_offset(Az, Alt,
         Az_nom, Alt_nom, 180./M_PI, &x_nom, &y_nom); // In degrees!
      fill_histogram_by_ident(19531,x_nom,y_nom,ewt);
      if ( n_img >= up[0].i.min_tel_img && n_img <= up[0].i.max_tel_img )
         fill_histogram_by_ident(19532,x_nom,y_nom,ewt);
      if ( n_img2 >= up[0].i.min_tel_img && n_img2 <= up[0].i.max_tel_img )
         fill_histogram_by_ident(19533,x_nom,y_nom,ewt);

      /* Angle between source and reconstructed direction */
      if ( diffuse_mode )
      {
         /* Diffuse data treated as a point source 
            (individual direction for every shower is the true shower direction). */
         stheta = angle_between(Az_true, Alt_true, Az, Alt);
      }
      else
      {
         /* Normal source where simulated with CORSIKA */
         stheta = angle_between(Az_src, Alt_src, Az, Alt);
      }

      if ( hsdata->event.shower.known )
      {
         double scaled = eval_cut_param(up[0].d.theta_escale,lg_energy);
         if ( stheta < max_theta * up[0].d.theta_scale * scaled )
            angle_cut_ok = 1;
         for ( icut=0; icut<7; icut++ )
         {
            if ( stheta < opt_theta_cut[icut][n_img2] * scaled )
               angle_cutx_ok[icut] = 1;
         }

         bnt.theta = stheta;
         bnt.sig_theta = sqrt((hsdata->event.shower.err_dir1*hsdata->event.shower.err_dir1+
            hsdata->event.shower.err_dir2*hsdata->event.shower.err_dir2) *
            ((hsdata->event.shower.num_img>=2)?(hsdata->event.shower.num_img-1.9999):0.) );
      }
      fill_histogram_by_ident(12054,rr,lg_energy,ewt);
      fill_histogram_by_ident(22054, lg_energy, 0., 1.0);
      fill_histogram_by_ident(22154, lg_energy, 0., ewt);

      /* Angle between viewing direction and reconstructed direction */
      ctheta = angle_between(Az_nom, Alt_nom, Az, Alt);

      if ( ctheta*(180./M_PI) < 1.0 )
         fill_histogram_by_ident(17001,mscrw,mscrl,ewt);
      if ( angle_cut_ok )
         fill_histogram_by_ident(17002,mscrw,mscrl,ewt);
      if ( eres_cut_ok )
         fill_histogram_by_ident(17003,mscrw,mscrl,ewt);
      if ( angle_cut_ok && eres_cut_ok )
         fill_histogram_by_ident(17004,mscrw,mscrl,ewt);
      if ( angle_cut_ok && eres_cut_ok && eres2_cut_ok )
         fill_histogram_by_ident(17005,mscrw,mscrl,ewt);
      if ( angle_cut_ok && eres_cut_ok && eres2_cut_ok && hmax_cut_ok )
         fill_histogram_by_ident(17006,mscrw,mscrl,ewt);

      fill_histogram_by_ident(15101,hmax,lg_energy,ewt);

      switch ( up[0].i.user_flags )
      {
         case 1:
            if ( mscrw < 0.9 && mscrl < 2.0 && 
                 mscrw > -2.0 && mscrl > -2.0 ) /* "standard cuts" */
               shape_cuts_ok = 1;
            break;
         case 2:
            if ( mscrw < 0.7 && mscrl < 2.0 && 
                 mscrw > -2.0 && mscrl > -2.0 ) /* "hard cuts" */
               shape_cuts_ok = 1;
            break;
         case 3:
            if ( mscrw < 1.2 && mscrl < 2.0 && 
                 mscrw > -2.0 && mscrl > -2.0 ) /* "loose cuts" */
               shape_cuts_ok = 1;
            break;
         default:
            if ( verbosity > 0 )
	       printf("Applying shape cuts %f < mscrw=%f < %f; %f < mscrl=%f < %f\n",
	          eval_cut_param(up[0].d.mscrw_min,lg_energy), mscrw,
	          eval_cut_param(up[0].d.mscrw_max,lg_energy),
	          eval_cut_param(up[0].d.mscrl_min,lg_energy), mscrl,
	          eval_cut_param(up[0].d.mscrl_max,lg_energy));
            if ( mscrw < eval_cut_param(up[0].d.mscrw_max,lg_energy) /*0.6*/ && 
                 mscrl < eval_cut_param(up[0].d.mscrl_max,lg_energy) /*1.2*/ && 
                 mscrw > eval_cut_param(up[0].d.mscrw_min,lg_energy) /*-2.0*/ && 
                 mscrl > eval_cut_param(up[0].d.mscrl_min,lg_energy) /*-2.0*/ )
               shape_cuts_ok = 1;
      }

      if ( verbosity > 0 )
         printf("n_img=%d, n_img2=%d, n_img3=%d, n_sc=%d, shape cut=%d, angle_ok=%d,%d,%d,%d,%d,%d,%d,%d, dE cut=%d, dE2 cut=%d, hmax cut=%d\n",
           n_img, n_img2, n_img3, n_sc, shape_cuts_ok, angle_cut_ok, 
           angle_cutx_ok[0], angle_cutx_ok[1], angle_cutx_ok[2],
           angle_cutx_ok[3], angle_cutx_ok[4], angle_cutx_ok[5], angle_cutx_ok[6],
           eres_cut_ok, eres2_cut_ok, hmax_cut_ok);

#ifdef FILL_ALL_TEL
      if ( hsdata->event.shower.known && 
           (hsdata->event.shower.result_bits & 0x05) == 0x05 )
      {
         for (itel=0; itel<hsdata->run_header.ntel; itel++)
         {
            int ht = 100000 * telescope_type[itel]; // Used for histograms
            double tr = v_tr[itel] = line_point_distance(xc, yc, 0.,
                     cos(Alt) * cos(Az),
                     cos(Alt) * sin(-Az),
                     sin(Alt),
                     hsdata->run_header.tel_pos[itel][0],
                     hsdata->run_header.tel_pos[itel][1],
                     hsdata->run_header.tel_pos[itel][2]);
            fill_histogram_by_ident(ht+18301,tr,lg_energy,ewt);
            if ( shape_cuts_ok )
               fill_histogram_by_ident(ht+18302,tr,lg_energy,ewt);
            if ( hsdata->event.teldata[itel].known )
            {
               fill_histogram_by_ident(ht+18311,tr,lg_energy,ewt);
               if ( shape_cuts_ok )
                  fill_histogram_by_ident(ht+18312,tr,lg_energy,ewt);
            }
            if ( img_ok[itel] && hsdata->event.teldata[itel].img != NULL )
            {
               fill_histogram_by_ident(ht+18321,tr,lg_energy,ewt);
               if ( shape_cuts_ok )
                  fill_histogram_by_ident(ht+18322,tr,lg_energy,ewt);
            }
         }
      }
#endif

      /* As a check for too strict shape cuts, fill histos without shape cuts applied. */
      if ( angle_cutx_ok[0] )
      {
         fill_histogram_by_ident(12053,rr,lg_energy,ewt);
         fill_histogram_by_ident(22053, lg_energy, 0., 1.0);
         fill_histogram_by_ident(22153, lg_energy, 0., ewt);
         if ( hmax_cut_ok )
         {
            fill_histogram_by_ident(12010,rs,lg_E_true,ewt);
            fill_histogram_by_ident(22010, lg_E_true, 0., 1.0);
            fill_histogram_by_ident(22110, lg_E_true, 0., ewt);
            fill_histogram_by_ident(12060,rr,lg_energy,ewt);
            fill_histogram_by_ident(22060, lg_energy, 0., 1.0);
            fill_histogram_by_ident(22160, lg_energy, 0., ewt);
            if ( eres2_cut_ok )
            {
               fill_histogram_by_ident(12014,rs,lg_E_true,ewt);
               fill_histogram_by_ident(22014, lg_E_true, 0., 1.0);
               fill_histogram_by_ident(22114, lg_E_true, 0., ewt);
               fill_histogram_by_ident(12064,rr,lg_energy,ewt);
               fill_histogram_by_ident(22064, lg_energy, 0., 1.0);
               fill_histogram_by_ident(22164, lg_energy, 0., ewt);
            }
         }
      }
      /* Normally, the shape cuts would be a minimum requirement for gamma candidates. */
      if ( shape_cuts_ok )
      {
         if ( verbosity > 0 && hsdata->event.shower.known )
            printf("Event passed shape cuts: mscrw=%5.3f, mscrl=%5.3f, E=%f+-%f (%f), c2/n=%f\n", 
               mscrw, mscrl, energy, energy/sqrt(w_sce), E_true, echi2);

         bnt.acceptance = 1;
         if ( angle_cutx_ok[0] )
         {
            bnt.acceptance = 2;
            if ( eres_cut_ok && angle_cutx_ok[1] )
            {
               bnt.acceptance = 3;
               if ( eres2_cut_ok && angle_cutx_ok[2] )
               {
                  bnt.acceptance = 4;
                  if ( hmax_cut_ok )
                     bnt.acceptance = 5;
               }
            }
         }

         fill_histogram_by_ident(12005,rs,lg_E_true,ewt);
         fill_histogram_by_ident(22005, lg_E_true, 0., 1.0);
         fill_histogram_by_ident(22105, lg_E_true, 0., ewt);
         fill_histogram_by_ident(12055,rr,lg_energy,ewt);
         fill_histogram_by_ident(22055, lg_energy, 0., 1.0);
         fill_histogram_by_ident(22155, lg_energy, 0., ewt);

         fill_histogram_by_ident(12100+n_img2,rs,lg_E_true,ewt);

         for (itel=0; itel<hsdata->run_header.ntel; itel++)
         {
            if ( img_ok[itel] && v_amp[itel] > 0. && hsdata->event.teldata[itel].img != NULL )
            {
               int j_img = 0;
               ImgData *img = &hsdata->event.teldata[itel].img[j_img];
               double ts = v_ts[itel];
               double tr = v_tr[itel];
               int ht = 100000 * telescope_type[itel]; // Used for histograms

               if ( img->tm_slope != 0. || img->tm_residual != 0. || img->tm_width1 == 0. )
               {
                  double xsh, ysh;
                  double tm_slope_deg = img->tm_slope * (M_PI/180.);
                  angles_to_offset(hsdata->event.shower.Az, hsdata->event.shower.Alt,
                     hsdata->event.trackdata[itel].azimuth_raw,
                     hsdata->event.trackdata[itel].altitude_raw, 1.0,
                     &xsh, &ysh);
                  if ( (xsh - hsdata->event.teldata[itel].img[j_img].x)*
                          cos(hsdata->event.teldata[itel].img[j_img].phi) +
                       (ysh - hsdata->event.teldata[itel].img[j_img].y)*
                          sin(hsdata->event.teldata[itel].img[j_img].phi) > 0. )
                     tm_slope_deg *= -1.;

                  fill_histogram_by_ident(ht+18402,ts,tm_slope_deg,ewt);
                  fill_histogram_by_ident(ht+18412,ts,img->tm_residual,ewt);
                  fill_histogram_by_ident(ht+18422,ts,img->tm_width1,ewt);
                  fill_histogram_by_ident(ht+18432,ts,img->tm_width2,ewt);
                  fill_histogram_by_ident(ht+18442,ts,img->tm_rise,ewt);

                  if ( angle_cutx_ok[0] )
                  {
                     fill_histogram_by_ident(ht+18403,ts,tm_slope_deg,ewt);
                     fill_histogram_by_ident(ht+18413,ts,img->tm_residual,ewt);
                     fill_histogram_by_ident(ht+18423,ts,img->tm_width1,ewt);
                     fill_histogram_by_ident(ht+18433,ts,img->tm_width2,ewt);
                     fill_histogram_by_ident(ht+18443,ts,img->tm_rise,ewt);

                     if ( hsdata->event.shower.known )
                     {
                        fill_histogram_by_ident(ht+18503,tr,tm_slope_deg,ewt);
                        fill_histogram_by_ident(ht+18513,tr,img->tm_residual,ewt);
                        fill_histogram_by_ident(ht+18523,tr,img->tm_width1,ewt);
                        fill_histogram_by_ident(ht+18533,tr,img->tm_width2,ewt);
                        fill_histogram_by_ident(ht+18543,tr,img->tm_rise,ewt);
                     }
                  }

                  if ( angle_cutx_ok[2] && eres_cut_ok && eres2_cut_ok && hmax_cut_ok )
                  {
                     fill_histogram_by_ident(ht+18404,ts,tm_slope_deg,ewt);
                     fill_histogram_by_ident(ht+18414,ts,img->tm_residual,ewt);
                     fill_histogram_by_ident(ht+18424,ts,img->tm_width1,ewt);
                     fill_histogram_by_ident(ht+18434,ts,img->tm_width2,ewt);
                     fill_histogram_by_ident(ht+18444,ts,img->tm_rise,ewt);

                     if ( hsdata->event.shower.known )
                     {
                        fill_histogram_by_ident(ht+18504,tr,tm_slope_deg,ewt);
                        fill_histogram_by_ident(ht+18514,tr,img->tm_residual,ewt);
                        fill_histogram_by_ident(ht+18524,tr,img->tm_width1,ewt);
                        fill_histogram_by_ident(ht+18534,tr,img->tm_width2,ewt);
                        fill_histogram_by_ident(ht+18544,tr,img->tm_rise,ewt);
                     }
                  }
               }
            }
         }

         if ( angle_cutx_ok[0] )
         {
            for (itel=0; itel<hsdata->run_header.ntel; itel++)
            {
               if ( img_ok[itel] && v_amp[itel] > 0. &&
                    hsdata->event.teldata[itel].raw != NULL &&
                    hsdata->event.teldata[itel].raw->known )
               {
                  struct user_parameters *up_tel = user_get_parameters(telescope_type[itel]);
                  double clip_amp = up_tel->d.clip_amp;
                  int jpix;
                  struct momstat stmom;
                  clear_moments(pixmom);
                  for ( jpix=0; jpix<hsdata->event.teldata[itel].image_pixels.pixels; jpix++)
                  {
                     int ipix = hsdata->event.teldata[itel].image_pixels.pixel_list[jpix];
                     double pixamp = calibrate_pixel_amplitude(hsdata,itel,ipix,0,-1,clip_amp);
                     //printf("Telescope %d pixel %d: %f peak p.e.\n", itel, ipix, pixamp);
                     fill_moments(pixmom,pixamp/v_amp[itel]);
                  }
                  stat_moments(pixmom,&stmom);
                  fill_histogram_by_ident(17801,
                     hsdata->event.teldata[itel].image_pixels.pixels,
                     log10(stmom.mean*v_amp[itel]),ewt);
                  fill_histogram_by_ident(17802,
                     hsdata->event.teldata[itel].image_pixels.pixels,
                     stmom.sigma,ewt);
                  fill_histogram_by_ident(17803,
                     hsdata->event.teldata[itel].image_pixels.pixels,
                     stmom.skewness,ewt);
                  fill_histogram_by_ident(17804,
                     hsdata->event.teldata[itel].image_pixels.pixels,
                     stmom.kurtosis,ewt);
               }
               if ( img_ok[itel] && hsdata->event.teldata[itel].img != NULL )
               {
                  if ( hsdata->event.teldata[itel].img[0].num_hot > 0 )
                  {
                     fill_histogram_by_ident(17800,
                        hsdata->event.teldata[itel].image_pixels.pixels,
                        hsdata->event.teldata[itel].img[0].hot_amp[0] / v_amp[itel], ewt);
                     if ( hsdata->event.teldata[itel].img[0].num_hot >= 3 )
                        fill_histogram_by_ident(17810,
                           hsdata->event.teldata[itel].image_pixels.pixels,
                           hsdata->event.teldata[itel].img[0].hot_amp[2] /
                           hsdata->event.teldata[itel].img[0].hot_amp[0], ewt);
                     else if ( hsdata->event.teldata[itel].img[0].num_hot >= 2 )
                        fill_histogram_by_ident(17810,
                           hsdata->event.teldata[itel].image_pixels.pixels,
                           hsdata->event.teldata[itel].img[0].hot_amp[1] /
                           hsdata->event.teldata[itel].img[0].hot_amp[0], ewt);
                  }
               }
            }

            fill_histogram_by_ident(12006,rs,lg_E_true,ewt);
            fill_histogram_by_ident(22006, lg_E_true, 0., 1.0);
            fill_histogram_by_ident(22106, lg_E_true, 0., ewt);
            fill_histogram_by_ident(12056,rr,lg_energy,ewt);
            fill_histogram_by_ident(22056, lg_energy, 0., 1.0);
            fill_histogram_by_ident(22156, lg_energy, 0., ewt);

            fill_histogram_by_ident(12200+n_img2,rs,lg_E_true,ewt);
            if ( eres_cut_ok )
            {
               if ( angle_cutx_ok[1] ) /* angle cut for shape+dE */
               {
                  fill_histogram_by_ident(12007,rs,lg_E_true,ewt);
                  fill_histogram_by_ident(22007, lg_E_true, 0., 1.0);
                  fill_histogram_by_ident(22107, lg_E_true, 0., ewt);
                  fill_histogram_by_ident(12057,rr,lg_energy,ewt);
                  fill_histogram_by_ident(22057, lg_energy, 0., 1.0);
                  fill_histogram_by_ident(22157, lg_energy, 0., ewt);
               }
               if ( eres2_cut_ok )
               {
                  if ( angle_cutx_ok[6] ) /* angle cut for shape+dE+dE2 */
                  {
                     fill_histogram_by_ident(12008,rs,lg_E_true,ewt);
                     fill_histogram_by_ident(22008, lg_E_true, 0., 1.0);
                     fill_histogram_by_ident(22108, lg_E_true, 0., ewt);
                     fill_histogram_by_ident(12058,rr,lg_energy,ewt);
                     fill_histogram_by_ident(22058, lg_energy, 0., 1.0);
                     fill_histogram_by_ident(22158, lg_energy, 0., ewt);
                     if ( hsdata->event.shower.xmax > 0 )
                     {
                        /* Look for Hmax induced biases */
                        int ihmax = (int)(hsdata->event.shower.xmax/50.);
                        if ( ihmax >= 0 && ihmax < 12 )
                           fill_histogram_by_ident(17400+ihmax,mscrw,mscrl,ewt);
                        fill_histogram_by_ident(18200,lg_energy,hsdata->event.shower.xmax,1.);
                        fill_histogram_by_ident(18201,lg_energy,hsdata->event.shower.xmax,ewt);
                        fill_histogram_by_ident(18211,lg_energy,hsdata->event.shower.xmax,
                           ewt*(lg_energy-lg_E_true));
                        fill_histogram_by_ident(18212,lg_energy,hsdata->event.shower.xmax,
                           ewt*(lg_energy-lg_E_true)*(lg_energy-lg_E_true));
                     }
                  }
                  if ( hmax_cut_ok && angle_cutx_ok[3]  ) /* angle cut for shape+dE+dE2+hmax */
                  {
                     event_selected = 1; /* Select for DST level 1x extraction */
                     fill_histogram_by_ident(12009,rs,lg_E_true,ewt);
                     fill_histogram_by_ident(22009, lg_E_true, 0., 1.0);
                     fill_histogram_by_ident(22109, lg_E_true, 0., ewt);
                     fill_histogram_by_ident(12059,rr,lg_energy,ewt);
                     fill_histogram_by_ident(22059, lg_energy, 0., 1.0);
                     fill_histogram_by_ident(22159, lg_energy, 0., ewt);
                     fill_histogram_by_ident(12300+n_img2,rs,lg_E_true,ewt);
                  }
               }
            }
            /* Similar for other final cut order */
            if ( hmax_cut_ok )
            {
               if ( angle_cutx_ok[3] ) /* angle cut for shape+hmax */
               {
                  fill_histogram_by_ident(12011,rs,lg_E_true,ewt);
                  fill_histogram_by_ident(22011, lg_E_true, 0., 1.0);
                  fill_histogram_by_ident(22111, lg_E_true, 0., ewt);
                  fill_histogram_by_ident(12061,rr,lg_energy,ewt);
                  fill_histogram_by_ident(22061, lg_energy, 0., 1.0);
                  fill_histogram_by_ident(22161, lg_energy, 0., ewt);
                  fill_histogram_by_ident(12400+n_img2,rs,lg_E_true,ewt);
               }
               if ( eres_cut_ok && angle_cutx_ok[4] ) /* angle cut for shape+dE+hmax */
               {
                  fill_histogram_by_ident(12012,rs,lg_E_true,ewt);
                  fill_histogram_by_ident(22012, lg_E_true, 0., 1.0);
                  fill_histogram_by_ident(22112, lg_E_true, 0., ewt);
                  fill_histogram_by_ident(12062,rr,lg_energy,ewt);
                  fill_histogram_by_ident(22062, lg_energy, 0., 1.0);
                  fill_histogram_by_ident(22162, lg_energy, 0., ewt);
               }
               if ( eres2_cut_ok && angle_cutx_ok[5] ) /* angle cut for shape+dE2+hmax */
               {
                  fill_histogram_by_ident(12013,rs,lg_E_true,ewt);
                  fill_histogram_by_ident(22013, lg_E_true, 0., 1.0);
                  fill_histogram_by_ident(22113, lg_E_true, 0., ewt);
                  fill_histogram_by_ident(12063,rr,lg_energy,ewt);
                  fill_histogram_by_ident(22063, lg_energy, 0., 1.0);
                  fill_histogram_by_ident(22163, lg_energy, 0., ewt);
               }
            }
         }
         if ( ctheta*(180./M_PI) < 1.0 )
         {
            fill_histogram_by_ident(12015,rs,lg_E_true,ewt);
            fill_histogram_by_ident(22015, lg_E_true, 0., 1.0);
            fill_histogram_by_ident(22115, lg_E_true, 0., ewt);
            fill_histogram_by_ident(12065,rr,lg_energy,ewt);
            fill_histogram_by_ident(22065, lg_energy, 0., 1.0);
            fill_histogram_by_ident(22165, lg_energy, 0., ewt);
            if ( eres_cut_ok )
            {
               fill_histogram_by_ident(12067,rr,lg_energy,ewt);
               fill_histogram_by_ident(22067, lg_energy, 0., 1.0);
               fill_histogram_by_ident(22167, lg_energy, 0., ewt);
               if ( eres2_cut_ok )
               {
                  fill_histogram_by_ident(12068,rr,lg_energy,ewt);
                  fill_histogram_by_ident(22068, lg_energy, 0., 1.0);
                  fill_histogram_by_ident(22168, lg_energy, 0., ewt);
                  if ( hmax_cut_ok )
                  {
                     fill_histogram_by_ident(12069,rr,lg_energy,ewt);
                     fill_histogram_by_ident(22069, lg_energy, 0., 1.0);
                     fill_histogram_by_ident(22169, lg_energy, 0., ewt);
                  }
               }
            }
         }

         for (itel=0; itel<hsdata->run_header.ntel; itel++)
         {
            if ( v_amp[itel]> 0. )
            {
               int ht = telescope_type[itel] * 100000;
               fill_histogram_by_ident(ht+18061,v_ts[itel],log10(v_amp[itel]),
                  ewt*(v_amp[itel]/E_true));
               fill_histogram_by_ident(ht+18062,v_ts[itel],log10(v_amp[itel]),
                  ewt*(v_amp[itel]/E_true)*(v_amp[itel]/E_true));
            }
         }
         fill_histogram_by_ident(15002,hmax,lg_E_true,ewt);
         fill_histogram_by_ident(15102,hmax,lg_energy,ewt);
         if ( angle_cut_ok )
         {
            fill_histogram_by_ident(15003,hmax,lg_E_true,ewt);
            fill_histogram_by_ident(15103,hmax,lg_energy,ewt);
            if ( eres_cut_ok )
            {
               fill_histogram_by_ident(15004,hmax,lg_E_true,ewt);
               fill_histogram_by_ident(15104,hmax,lg_energy,ewt);
               if ( eres2_cut_ok )
               {
                  fill_histogram_by_ident(15005,hmax,lg_E_true,ewt);
                  fill_histogram_by_ident(15105,hmax,lg_energy,ewt);
                  if ( hmax_cut_ok )
                  {
                     fill_histogram_by_ident(15006,hmax,lg_E_true,ewt);
                     fill_histogram_by_ident(15106,hmax,lg_energy,ewt);
                  }
               }
            }
         }
         fill_histogram_by_ident(19534,x_nom,y_nom,ewt);
         if ( hmax_cut_ok )
         {
            fill_histogram_by_ident(19537,x_nom,y_nom,ewt);
            if ( eres2_cut_ok )
               fill_histogram_by_ident(19539,x_nom,y_nom,ewt);
         }
         if ( eres_cut_ok )
         {
            fill_histogram_by_ident(19535,x_nom,y_nom,ewt);
            if ( eres2_cut_ok )
            {
               fill_histogram_by_ident(19540,x_nom,y_nom,ewt);
               if ( hmax_cut_ok )
                  fill_histogram_by_ident(19506,x_nom,y_nom,ewt);
            }
            if ( hmax_cut_ok )
               fill_histogram_by_ident(19538,x_nom,y_nom,ewt);
         }
      }
      else if ( verbosity > 0 )
         printf("Event failed shape cuts: mscrw=%5.3f, mscrl=%5.3f, E=%f+-%f (%f), c2/n=%f\n",
             mscrw, mscrl, energy,energy/sqrt(w_sce+1e-10), E_true, echi2);

      if ( hsdata->event.shower.known )
      {
         da = angle_between(hsdata->event.shower.Az, hsdata->event.shower.Alt,
                            hsdata->mc_shower.azimuth, hsdata->mc_shower.altitude) 
              * (180./M_PI);
         /* any reconstructed shower */
         fill_histogram_by_ident(19001,hsdata->event.shower.num_img,da,ewt);
         fill_histogram_by_ident(19002,lg_E_true,da,ewt);
         fill_histogram_by_ident(19003,rs,da,ewt);
         fill_histogram_by_ident(19012,lg_E_true,lg_energy-lg_E_true,ewt);
         fill_histogram_by_ident(19013,lg_energy,lg_energy-lg_E_true,ewt);
         fill_histogram_by_ident(19014,lg_energy0,lg_energy0-lg_E_true,ewt);
         if ( shape_cuts_ok )
         {  /* shape */
            fill_histogram_by_ident(19101,hsdata->event.shower.num_img,da,ewt);
            fill_histogram_by_ident(19102,lg_E_true,da,ewt);
            fill_histogram_by_ident(19103,rs,da,ewt);
            fill_histogram_by_ident(19112,lg_E_true,lg_energy-lg_E_true,ewt);
            fill_histogram_by_ident(19113,lg_energy,lg_energy-lg_E_true,ewt);
            fill_histogram_by_ident(19114,lg_energy0,lg_energy0-lg_E_true,ewt);
            if ( hmax_cut_ok )
            {  /* shape+hmax */
               fill_histogram_by_ident(19401,hsdata->event.shower.num_img,da,ewt);
               fill_histogram_by_ident(19402,log10(hsdata->mc_shower.energy),da,ewt);
               fill_histogram_by_ident(19403,rs,da,ewt);
               fill_histogram_by_ident(19412,lg_E_true,lg_energy-lg_E_true,ewt);
               fill_histogram_by_ident(19413,lg_energy,lg_energy-lg_E_true,ewt);
               fill_histogram_by_ident(19414,lg_energy0,lg_energy0-lg_E_true,ewt);
               if ( eres2_cut_ok )
               {  /* shape+dE2+hmax */
                  fill_histogram_by_ident(19601,hsdata->event.shower.num_img,da,ewt);
                  fill_histogram_by_ident(19602,lg_E_true,da,ewt);
                  fill_histogram_by_ident(19603,rs,da,ewt);
                  fill_histogram_by_ident(19612,lg_E_true,lg_energy-lg_E_true,ewt);
                  fill_histogram_by_ident(19613,lg_energy,lg_energy-lg_E_true,ewt);
                  fill_histogram_by_ident(19614,lg_energy0,lg_energy0-lg_E_true,ewt);
               }
            }
            if ( eres_cut_ok )
            {  /* shape+dE */
               fill_histogram_by_ident(19201,hsdata->event.shower.num_img,da,ewt);
               fill_histogram_by_ident(19202,log10(hsdata->mc_shower.energy),da,ewt);
               fill_histogram_by_ident(19203,rs,da,ewt);
               fill_histogram_by_ident(19212,lg_E_true,lg_energy-lg_E_true,ewt);
               fill_histogram_by_ident(19213,lg_energy,lg_energy-lg_E_true,ewt);
               fill_histogram_by_ident(19214,lg_energy0,lg_energy0-lg_E_true,ewt);
               if ( hmax_cut_ok )
               {  /* shape+dE+hmax */
                  fill_histogram_by_ident(19501,hsdata->event.shower.num_img,da,ewt);
                  fill_histogram_by_ident(19502,lg_E_true,da,ewt);
                  fill_histogram_by_ident(19503,rs,da,ewt);
                  fill_histogram_by_ident(19512,lg_E_true,lg_energy-lg_E_true,ewt);
                  fill_histogram_by_ident(19513,lg_energy,lg_energy-lg_E_true,ewt);
                  fill_histogram_by_ident(19514,lg_energy0,lg_energy0-lg_E_true,ewt);
               }
               if ( eres2_cut_ok )
               {  /* shape+dE+dE2 */
                  fill_histogram_by_ident(19701,hsdata->event.shower.num_img,da,ewt);
                  fill_histogram_by_ident(19702,lg_E_true,da,ewt);
                  fill_histogram_by_ident(19703,rs,da,ewt);
                  fill_histogram_by_ident(19712,lg_E_true,lg_energy-lg_E_true,ewt);
                  fill_histogram_by_ident(19713,lg_energy,lg_energy-lg_E_true,ewt);
                  fill_histogram_by_ident(19714,lg_energy0,lg_energy0-lg_E_true,ewt);
                  if ( hmax_cut_ok )
                  {  /* shape+dE+dE2+hmax */
                     fill_histogram_by_ident(19301,hsdata->event.shower.num_img,da,ewt);
                     fill_histogram_by_ident(19302,lg_E_true,da,ewt);
                     fill_histogram_by_ident(19303,rs,da,ewt);
                     fill_histogram_by_ident(19312,lg_E_true,lg_energy-lg_E_true,ewt);
                     fill_histogram_by_ident(19313,lg_energy,lg_energy-lg_E_true,ewt);
                     fill_histogram_by_ident(19314,lg_energy0,lg_energy0-lg_E_true,ewt);
                  }
               }
            }
         }
      }
   /* - */
}

/* ------------------------- user_done ------------------------------ */
/** After all data for a file (usually one run) was processed. */

static void user_done (AllHessData *hsdata)
{
}

/* ------------------------- prog_path ------------------------- */
/** Find the path from which the current program was started. */

static char *prog_path (void);

static char *prog_path ()
{
   /* If the following two variables are provided by the caller,
      this function can be made reentrant. */
   static char bf[PATH_MAX+2];
   size_t bs = PATH_MAX+1;

   char cmdl_fname[128], program[PATH_MAX+2], path_elem[PATH_MAX+2];
   pid_t i = getpid();
   int n;
   size_t nu;
   FILE *f;
   char *s, *path;

   /* First action is to get argv[0] without help by main(). */
   snprintf(cmdl_fname, sizeof(cmdl_fname)-1,"/proc/%d/cmdline", i); /* Linux-specific */
   if ( (f = fopen(cmdl_fname,"r")) == NULL )
   {
      strncpy(bf,".",bs); /* Probably not Linux. Need another trick. */
      return bf;
   }
   for (nu=0; nu+1<bs; nu++)
      if ( (bf[nu] = fgetc(f)) == '\0' )
         break;
   fclose(f);

   /* Now lets see if the program was called with an explicit path name on the command line. */
   if ( (s = strrchr(bf,'/')) != NULL ) /* Explicit path in cmdline? */
   {
      if ( s != bf )
         *s = '\0';
      else
         s[1] = '\0'; /* A program in the root directory */
      return bf;
   }

   /* We have to scan the PATH environment variables for the location. */
   if ( (path = getenv("PATH")) == NULL )
   {
      strncpy(bf,".",bs); /* Perhaps the program itself destroyed PATH. */
      return bf;
   }
   strncpy(program,bf,PATH_MAX);
   n = 0;
   while ( getword(path,&n,path_elem,PATH_MAX,':','\0') > 0 )
   {
      size_t len = strlen(path_elem);
      if ( len+1+strlen(program) < PATH_MAX )
      {
         path_elem[len] = '/';
         strcpy(path_elem+len+1,program);
         if ( access(path_elem,X_OK) == 0 )
         {
            path_elem[len] = '\0';
            strncpy(bf,path_elem,bs);
            return bf;
         }
      }
   }

   strncpy(bf,".",bs); /* Probably wrong but we shouldn't come along here actually. */
   return bf;
}

/* ------------------------- user_finish --------------------------- */
/** Final call before program terminates. */

static void user_finish (AllHessData *hsdata)
{
   if ( hist_fname[0] == '\0' )
   {
      strcpy(hist_fname,"user_histograms.hdata.gz");
      if ( getenv("USER_ANALYSIS") != NULL )
      {
         if ( strlen(getenv("USER_ANALYSIS"))+16 < sizeof(hist_fname) )
         {
            snprintf(hist_fname,sizeof(hist_fname),"%s.hdata.gz",getenv("USER_ANALYSIS"));
         }
      }
   }
   /* Example code (here auto-generating lookup tables): */
   fflush(NULL);
   fprintf(stderr,"Writing histograms to '%s'\n", hist_fname);
   write_all_histograms(hist_fname);
   if ( hsdata != NULL && hsdata->mc_shower.primary_id == 0 )
   {
      if ( auto_lookup )
      {
         char cmd[4096];
         char *bsdir = prog_path(); /* gen_lookup should be in the same place as read_hess */
         char *hessroot = getenv("HESSROOT");
         if ( hessroot == NULL )
            hessroot = getenv("CTA_PATH");
         snprintf(cmd,sizeof(cmd),"%s/gen_lookup",bsdir);
         if ( access(cmd,X_OK) == 0 )
         {
            snprintf(cmd,sizeof(cmd),"%s/gen_lookup %s %s", bsdir, hist_fname, lookup_fname);
         }
         else if ( strcmp(bsdir,".") == 0 && hessroot != NULL )
         {
            snprintf(cmd,sizeof(cmd),"%s/bin/gen_lookup %s %s", hessroot, hist_fname, lookup_fname);
         }
         else
            snprintf(cmd,sizeof(cmd),"gen_lookup %s %s", hist_fname, lookup_fname);
         fprintf(stderr,"Auto-generating lookups with command '%s'\n",cmd);
         system(cmd);
      }
      else
         fprintf(stderr,"Generate a matching lookup file (if these were gamma showers) with:\n"
          "    gen_lookup %s %s\n", hist_fname, lookup_fname);
      
      if ( pixmom != NULL )
      {
         free_moments(pixmom);
         pixmom = NULL;
      }
   }
   /* - */
}

/* ======================= do_user_ana ==================== */
/* @short All-purpose user function, delegates work to static
 *        functions where we can easily change interfaces.
 *
 * @param hsdata Pointer to the all-including data structure (see io_hess.h).
 * @param item_type The item type of the data block that was read before.
 * @param stage Only relevent for the IO_TYPE_HESS_EVENT item type
 *              (stage 0 is with shower parameters reconstructed in
 *              sim_hessarray and stage 1 is with parameters reconstructed
 *              in read_hess; see '-r' option there) and for the
 *              0 (=end of data) item type (stage 0 after each file,
 *              stage 1 just before read_hess terminates).
 */

int do_user_ana (AllHessData *hsdata, unsigned long item_type, int stage)
{
   static int init_done = 0;
   event_selected = 0;

   switch ( item_type )
   {
      case IO_TYPE_HESS_RUNHEADER:
         /* New run. This may be too early for initialisation. */
         /* Relevant new data: hsdata->run_header */
         break;
      case IO_TYPE_HESS_MCRUNHEADER:
         /* MC run header following run header. */
         /* Relevant new data: hsdata->mc_run_header */
         break;
      case IO_TYPE_HESS_MC_SHOWER:
         if ( !init_done || tel_types_change )
         {
            user_init(hsdata);
            init_done = 1;
         }
         if ( tel_types_change )
            init_telescope_types(hsdata);
         /* Relevant new data: hsdata->mc_shower */
         user_mc_shower_fill(hsdata);
         break;
      case IO_TYPE_HESS_MC_EVENT:
         /* Relevant new data: hsdata->mc_event */
         user_mc_event_fill(hsdata);
         break;
      case IO_TYPE_HESS_CAMSETTINGS:
         tel_types_change = 1;
         break;
      case IO_TYPE_HESS_EVENT:
         if ( !init_done || tel_types_change )
         {
            user_init(hsdata);
            init_done = 1;
         }
         if ( tel_types_change )
            init_telescope_types(hsdata);
         /* Here comes the real beef. */
         /*   Stage 0: with reconstruction from sim_hessarray */
         /*   Stage 1: optional, with new reconstruction in read_hess. */
         user_event_fill(hsdata,stage);
         break;
      case IO_TYPE_HESS_RUNSTAT:
         /* End-of-run stats for triggered events */
         /* Relevant new data: hsdata->run_stat */
         break;
      case IO_TYPE_HESS_MC_RUNSTAT:
         /* End-of-run stats for MC showers */
         /* Relevant new data: hsdata->mc_run_stat */
         break;
      case 0:
         if ( stage == 0 )
            user_done(hsdata);
         else
            user_finish(hsdata);
         break;
      default:
         break; /* Everything else is ignored. */
   }

   return 0;
}
