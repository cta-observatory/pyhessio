/* ============================================================================

Copyright (C) 2003, 2009, 2011, 2015  Konrad Bernloehr

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

/** @file reconstruct.c
 *  @short Second moments type image analysis.
 *
 *  @date    @verbatim CVS $Revision: 1.70 $ @endverbatim
 *  @version @verbatim CVS $Date: 2017/06/07 14:33:27 $ @endverbatim
 */

#include "initial.h"      /* This file includes others as required. */
#include "io_hess.h"
#include "rec_tools.h"
#include "reconstruct.h"
#include "user_analysis.h"
#ifdef WITH_RANDFLAT
#include "rndm2.h"
#endif

/** The factor needed to transform from mean p.e. units to units of the single-p.e. peak:
    Depends on the collection efficiency, the asymmetry of the single p.e. amplitude 
    distribution and the electronic noise added to the signals. Default value is for HESS. */
#define CALIB_SCALE 0.92

static int px_shape_type[H_MAX_TEL];

#define H_MAX_NB 50

struct camera_nb_list
{
   int npix;           ///< Number of pixels in camera.
   int nbsize;         ///< Number of neighbours in list (elements in nblist).
   int *pix_num_nb;    ///< Number of neighbours for each pixel.
   int *pix_first_nb;  ///< Where in list is the first of the neighbours for each pixel.
   int *nblist;        ///< The actual packed list of all neighbours for all pixels.
};

static struct camera_nb_list nb_lists[H_MAX_TEL][3]; ///< To be filled with up to 3 neighbour lists for each telescope.
static struct camera_nb_list ext_list[H_MAX_TEL];    ///< Optional extension lists beyond image cleaning.

int allocate_nb_list(int itel, int npix, int shape_type, int nnbs, int *nbs);
int deallocate_nb_list(int itel);

int deallocate_nb_list(int itel)
{
   int k;
   if ( itel < 0 || itel >= H_MAX_TEL )
      return -1;
   if ( ext_list[itel].nblist != NULL )
   {
      free(ext_list[itel].nblist);
      ext_list[itel].nblist = NULL;
   }
   if ( ext_list[itel].pix_num_nb != NULL )
   {
      free(ext_list[itel].pix_num_nb);
      ext_list[itel].pix_num_nb = NULL;
   }
   if ( ext_list[itel].pix_first_nb != NULL )
   {
      free(ext_list[itel].pix_first_nb);
      ext_list[itel].pix_first_nb = NULL;
   }
   ext_list[itel].nbsize = ext_list[itel].npix = 0;
   for ( k=0; k<3; k++ )
   {
      if ( nb_lists[itel][k].nblist != NULL )
      {
         free(nb_lists[itel][k].nblist);
         nb_lists[itel][k].nblist = NULL;
      }
      if ( nb_lists[itel][k].pix_num_nb != NULL )
      {
         free(nb_lists[itel][k].pix_num_nb);
         nb_lists[itel][k].pix_num_nb = NULL;
      }
      if ( nb_lists[itel][k].pix_first_nb != NULL )
      {
         free(nb_lists[itel][k].pix_first_nb);
         nb_lists[itel][k].pix_first_nb = NULL;
      }
      nb_lists[itel][k].nbsize = nb_lists[itel][k].npix = 0;
   }
   
   return 0;
}

static int image_list[H_MAX_TEL][H_MAX_PIX];
static int image_numpix[H_MAX_TEL];

static double pixel_amp[H_MAX_TEL][H_MAX_PIX];
static int show_total_amp = 0;
static int pixel_sat[H_MAX_TEL];

static char pixel_disabled[H_MAX_TEL][H_MAX_PIX];
static int any_disabled[H_MAX_TEL];

static double camera_radius_eff[H_MAX_TEL];
static double camera_radius_max[H_MAX_TEL];

static double integration_correction[H_MAX_TEL][H_MAX_GAINS];

static int verbosity = 0;

/* ---------------------- set_disabled_pixels ------------------------ */

/** Set up pixels to be ignored (regarded as zero amplitude) in the analysis
    if they either have HV disabled or the camera active radius is clipped. 
    
    @param hsdata Pointer to all available data and configurations.
    @param itel Telescope index where we set new values.
    @param broken_pixels_fraction Optional fraction of additional pixels to be
              set like dead pixels (not usable for analysis).
    
    Disabled pixels are ignored in the evaluation of the camera radius. 
 */

int set_disabled_pixels(AllHessData *hsdata, int itel, double broken_pixels_fraction)
{
   int npix, ipix;
   CameraSettings *camset;
   PixelDisabled *pixdis;
   double clip_radius = 0;
   int ttype;
   UserParameters *up;
   any_disabled[itel] = 0;

   if ( hsdata == NULL || itel < 0 || itel >= H_MAX_TEL )
      return -1;

   camset = &hsdata->camera_set[itel];
   pixdis = &hsdata->pixel_disabled[itel];
   ttype = which_telescope_type(camset);
   up = user_get_parameters(ttype);
   clip_radius = tan(up->d.camera_clipping_deg*(M_PI/180.)) * camset->flen;
   npix = camset->num_pixels;

   for ( ipix=0; ipix<npix; ipix++ )
   {
      int nvd, ivd;
      pixel_disabled[itel][ipix] = 0;

      /* Disable due to clipping the camera size from simulation to analysis */
      if ( clip_radius > 0. )
      {
         double xp = camset->xpix[ipix];
         double yp = camset->ypix[ipix];
         double r = sqrt(xp*xp+yp*yp);
         if ( r+0.5*camset->size[ipix] >= clip_radius )
         {
            pixel_disabled[itel][ipix] = 1;
            any_disabled[itel] = 1;
         }
      }
      
      /* Disable randomly broken pixels */
      if ( broken_pixels_fraction > 0. )
      {
         /* Would have preferred to use the same random number generator */
	 /* as in sim_telarray but using generator from stdlib to avoid */
	 /* reverse dependence, at least for now. */
#ifdef WITH_RANDFLAT
      	 if ( RandFlat() < broken_pixels_fraction )
#else
      	 if ( drand48() < broken_pixels_fraction )
#endif
         {
            pixel_disabled[itel][ipix] = 1;
            any_disabled[itel] = 1;
         }
      }

      /* Disable pixels that were simulated with HV off (or disabled in an earlier processing stage) */
      nvd = pixdis->num_HV_disabled;
      for ( ivd=0; ivd<nvd; ivd++ )
      {
         if ( pixdis->HV_disabled[ivd] == ipix )
         {
            pixel_disabled[itel][ipix] = 1;
            any_disabled[itel] = 1;
         }
      }
   }

   /* The combined set of disabled pixels is stored as HV disabled */
   if ( any_disabled[itel] )
   {
      int nd = 0;
      for ( ipix=0; ipix<npix; ipix++ )
      {
         if ( pixel_disabled[itel][ipix] )
         {
            pixdis->HV_disabled[nd] = ipix;
            nd++;
         }
      }
      pixdis->num_HV_disabled = nd;
   }

   return 0;
}

/* ------------------------ guess_pixel_shape --------------------------- */
/**
 *  Guess the common pixel shape type from relative positions of neighbours.
 */

static int guess_pixel_shape(CameraSettings *camset, int itel);

static int guess_pixel_shape(CameraSettings *camset, int itel)
{
   int npix = camset->num_pixels;
   int i, j;
   int stat_st[6] = {0, 0, 0, 0, 0, 0 };
   double asum = 0., dsum = 0., aod2 = 0.;
   int px_shape = camset->pixel_shape[0];

   /* If we know the shape and know that all pixels have the same shape, there is no guessing */
   if ( camset->common_pixel_shape && px_shape >= 0 )
      return px_shape;

   for (i=1; i<npix; i++)
   {
      if ( camset->pixel_shape[i] != px_shape )
      {
         px_shape = -2;
         break;
      }
   }
   if ( px_shape >= 0 )
      return px_shape;

   if ( px_shape == -2 )
   {
      printf("Different pixel shapes in telescope ID  %d, need to guess.\n", camset->tel_id);
   }

   /* If shape type is unknown or differs, we apply some heuristics */

   for (i=0; i<npix; i++)
   {
      if ( pixel_disabled[itel][i] )
         continue;
      asum += camset->area[i];
      dsum += camset->size[i];

      for (j=0; j<i; j++)
      {
         double ds, dx, dy, d2, a;
         int ia;
         if ( pixel_disabled[itel][j] )
            continue;
         ds = camset->size[i] + camset->size[j];
         dx = camset->xpix[i] - camset->xpix[j];
         dy = camset->ypix[i] - camset->ypix[j];
         d2 = dx*dx + dy*dy;
         if ( d2 < 0.5*ds*ds ) 
         {
            a  = (180./M_PI) * atan2(dy,dx);
            if ( a < -1. )
               a += 180.;
            ia = ((int) ((a+0.5)/5.)) * 5;
            if ( ia == 0 )
               stat_st[0]++;
            else if ( ia == 60 )
               stat_st[1]++;
            else if ( ia == 90 )
               stat_st[2]++;
            else if ( ia == 120 )
               stat_st[3]++;
            else if ( ia == 30 )
               stat_st[4]++;
            else if ( ia == 150 )
               stat_st[5]++;
         }
      }
   }

   asum /= (((double) npix)+1e-10);
   dsum /= (((double) npix)+1e-10);
   aod2 = asum/(dsum*dsum);

   if ( stat_st[0] > 0 && stat_st[2] > 0 && 
        stat_st[1] == 0 && stat_st[3] == 0 )
   {
      px_shape = 2; /* probably square pixels */
      if ( aod2 < 0.99 || aod2 > 1.01 )
      {
         fprintf(stderr,
            "Pixel positions in telescope %d indicate square pixels but area/size^2 does not match.\n",
            camset->tel_id);
      }
   }
   else /* most likely hexagonal, orientation? */
   {
      if ( 4*stat_st[2] < (stat_st[1] + stat_st[3]) )
         px_shape = 1;
      else if ( stat_st[2] > 0 && stat_st[0] == 0 )
         px_shape = 3;
      else
      {
         px_shape = 0; /* Fall back to circular */
         if ( aod2 < 0.99*M_PI/4. || aod2 > 1.01*M_PI/4. )
         {
            fprintf(stderr,
               "Pixel positions in telescope %d indicate round pixels but area/size^2 does not match.\n",
               camset->tel_id);
         }
      }
      if ( px_shape != 0 )
      {
         if ( aod2 < 0.99*sqrt(3.)/2. || aod2 > 1.01*sqrt(3.)/2. )
         {
            if ( aod2 >= 0.99*M_PI/4. && aod2 <= 1.01*M_PI/4. )
            {
               px_shape = 0; // Round pixels on hexagonal pattern
            }
            else
               fprintf(stderr,
                  "Pixel positions in telescope %d indicate hexagonal pixels but area/size^2 does not match.\n"
                  "Expected %f (or %f for round pixels) but found %f\n",
                  camset->tel_id, sqrt(3.)/2., M_PI/4., aod2);
         }
      }
   }
   for (i=0; i<npix; i++)
   {
      if ( camset->pixel_shape[i] == -1 )
         camset->pixel_shape[i] = px_shape;
   }
#ifdef DEBUG_PIXEL_NB
   fprintf(stderr,"Pixel shape type of telescope #%d (ID %d) seems to be %d (stat = %d, %d, %d, %d)\n",
      itel, camset->tel_id, px_shape, stat_st[0], stat_st[1], stat_st[2], stat_st[3]);
#endif

   return px_shape;
}

/* --------------------------- find_neighbours ---------------------------- */

/** Find the list of neighbours for each pixel. */

static int find_neighbours(CameraSettings *camset, int itel);

static int find_neighbours(CameraSettings *camset, int itel)
{
   int nb[3][H_MAX_PIX][H_MAX_NB];  ///< Temporary neighbour lists for one telescope.
   int xt[H_MAX_PIX][H_MAX_NB];
   int nnb[3][H_MAX_PIX], nxt[H_MAX_PIX];
   int ntot[3], ntxt;
   double r2[3] = { 1.0, 0., 0. }, rxt2 = 0.;
   int ttype = which_telescope_type(camset);
   int i, j, k, n;
   int npix = camset->num_pixels;
   UserParameters *up = user_get_parameters(ttype);
   px_shape_type[itel] = guess_pixel_shape(camset, itel);

#ifndef DEBUG_PIXEL_NB
 if ( verbosity > 0 )
#endif
      printf("Start of neighbour pixel finding for telescope ID %d\n",camset->tel_id);

   for ( k=0; k<3; k++ )
   {
      if ( nb_lists[itel][k].nblist != NULL || nb_lists[itel][k].pix_num_nb != NULL || 
           nb_lists[itel][k].pix_first_nb != NULL )
      {
         fprintf(stderr,"Invalid pixel neighbour list initialization for telescope ID %d\n", 
            camset->tel_id);
         return -1;
      }
   }
   if ( ext_list[itel].nblist != NULL || ext_list[itel].pix_num_nb != NULL ||
        ext_list[itel].pix_first_nb != NULL )
   {
      fprintf(stderr,"Invalid pixel neighbour extension list initialization for telescope ID %d\n", 
         camset->tel_id);
      return -1;
   }

   for ( k=0; k<3; k++ )
   {
      ntot[k] = 0;
      for ( i=0; i<npix; i++ )
         nnb[k][i] = 0;
   }
   for ( i=0; i<npix; i++ )
      nxt[i] = 0;

   for ( k=0; k<3; k++ )
      r2[k] = up->d.r_nb[k] * up->d.r_nb[k];
   rxt2 = up->d.r_ne * up->d.r_ne;
   if ( r2[0] <= 0. ) /* Apply shape-type dependent default neighbour limit (first neighbour only). */
   {
      /* The default radius used to be sqrt(2.) times pixel diameter without extra margin for gaps. */
      /* New defaults include diagonal neighbours for square pixels */
      /* but do normally not extend across module gaps (also did not before). */
      if ( px_shape_type[itel] == 2 )
         r2[0] = 1.6*1.6; /* Sqrt(2.) plus some margin for gaps, includes diagonal neighbours. */
      else
         r2[0] = 1.2*1.2; /* With 20% margin for gaps: effectively as old defaults. */
   }

   /* Detect neighbours, up to some maximum, not packed yet. */
   for (i=0; i<npix; i++)
   {
      if ( pixel_disabled[itel][i] )
         continue;
      for (j=0; j<i; j++)
      {
         double ds, dx, dy, d2;
         if ( pixel_disabled[itel][j] )
            continue;
         ds = 0.5*(camset->size[i] + camset->size[j]);
         dx = camset->xpix[i] - camset->xpix[j];
         dy = camset->ypix[i] - camset->ypix[j];
         d2 = dx*dx + dy*dy;
         /* Immediate neighbours */
         if ( d2 < r2[0]*(ds*ds) )
         {
            if ( nnb[0][i] < H_MAX_NB )
               nb[0][i][nnb[0][i]++] = j;
            if ( nnb[0][j] < H_MAX_NB )
               nb[0][j][nnb[0][j]++] = i;
         }
         /* Further neighbours only on request (not used with classical 2-level cleaning) */
         else if ( d2 < r2[1]*(ds*ds) )
         {
            if ( nnb[1][i] < H_MAX_NB )
               nb[1][i][nnb[1][i]++] = j;
            if ( nnb[1][j] < H_MAX_NB )
               nb[1][j][nnb[1][j]++] = i;
         }
         /* There can be a third set of even more distant neighbours. */
         else if ( d2 < r2[2]*(ds*ds) )
         {
            if ( nnb[2][i] < H_MAX_NB )
               nb[2][i][nnb[2][i]++] = j;
            if ( nnb[2][j] < H_MAX_NB )
               nb[2][j][nnb[2][j]++] = i;
         }
         /* The extension list beyond image cleaning is independent of the neighbours for the cleaning itself. */
         if ( d2 < rxt2*(ds*ds) )
         {
            if ( nxt[i] < H_MAX_NB )
               xt[i][nxt[i]++] = j;
            if ( nxt[j] < H_MAX_NB )
               xt[j][nxt[j]++] = i;
         }
      }
   }

#ifdef DEBUG_PIXEL_NB
   for (i=0; i<npix && i<30; i++)
   {
      printf("Pixel %d has %d direct neighbours, %d in second set, %d in third set, %d in extension list.\n",
         i, nnb[0][i], nnb[1][i], nnb[2][i], nxt[i]);
#ifdef DEBUG_PIXEL_NB2
      { int j;
        if ( nnb[0][i] > 0 )
        { printf("   direct:"); for (j=0; j<nnb[0][i]; j++) printf(" %d", nb[0][i][j]); printf("\n"); }
        if ( nnb[1][i] > 0 )
        { printf("   second:"); for (j=0; j<nnb[1][i]; j++) printf(" %d", nb[1][i][j]); printf("\n"); }
        if ( nnb[2][i] > 0 )
        { printf("   third:"); for (j=0; j<nnb[2][i]; j++) printf(" %d", nb[2][i][j]); printf("\n"); }
        if ( nxt[i] > 0 )
        { printf("   extension:"); for (j=0; j<nxt[i]; j++) printf(" %d", xt[i][j]); printf("\n"); }
      }
#endif
   }
#endif

   /* Set up packed neighbour list */
   for ( k=0; k<3; k++ )
   {
      ntot[k] = 0;
      for ( i=0; i<npix; i++ )
         ntot[k] += nnb[k][i];
      if ( ntot[k] > 0 )
      {
#ifndef DEBUG_PIXEL_NB
         if ( verbosity > 1 )
#endif
         printf("Neighbour pixel list %d in telescope ID %d has size %d from %d pixels.\n",
            k, camset->tel_id, ntot[k], npix);
         if ( (nb_lists[itel][k].nblist = (int *) calloc(ntot[k],sizeof(int))) == NULL ||
              (nb_lists[itel][k].pix_num_nb = (int *) calloc(npix,sizeof(int))) == NULL ||
              (nb_lists[itel][k].pix_first_nb = (int *) calloc(npix,sizeof(int))) == NULL )
         {
            fprintf(stderr,"Allocation of neighbour list %d failed for telescope ID %d\n", k, camset->tel_id);
            return -1;
         }
         nb_lists[itel][k].npix = npix;
         nb_lists[itel][k].nbsize = ntot[k];
         n = 0;
         for ( i=0; i<npix; i++ )
         {
            nb_lists[itel][k].pix_first_nb[i] = n;
            nb_lists[itel][k].pix_num_nb[i] = nnb[k][i];
            /* Fill packed list of neighbours */
            for ( j=0; j<nnb[k][i]; j++ )
               nb_lists[itel][k].nblist[n+j] = nb[k][i][j];
            n += nnb[k][i];
         }
      }
      else
      {
         nb_lists[itel][k].npix = npix;
         nb_lists[itel][k].nbsize = 0;
      }
   }
   /* Set up packed extension list */
   ntxt = 0;
   for ( i=0; i<npix; i++ )
      ntxt += nxt[i];
   if ( ntxt > 0 )
   {
#ifndef DEBUG_PIXEL_NB
      if ( verbosity > 1 )
#endif
      printf("Extension pixel list in telescope ID %d has size %d from %d pixels.\n",
         camset->tel_id, ntxt, npix);
      if ( (ext_list[itel].nblist = (int *) calloc(ntxt,sizeof(int))) == NULL ||
           (ext_list[itel].pix_num_nb = (int *) calloc(npix,sizeof(int))) == NULL ||
           (ext_list[itel].pix_first_nb = (int *) calloc(npix,sizeof(int))) == NULL )
      {
         fprintf(stderr,"Allocation of neighbour extension list failed for telescope ID %d\n", camset->tel_id);
         return -1;
      }
      ext_list[itel].npix = npix;
      ext_list[itel].nbsize = ntxt;
      n = 0;
      for ( i=0; i<npix; i++ )
      {
         ext_list[itel].pix_first_nb[i] = n;
         ext_list[itel].pix_num_nb[i] = nxt[i];
         /* Fill packed list of neighbours */
         for ( j=0; j<nxt[i]; j++ )
            ext_list[itel].nblist[n+j] = xt[i][j];
         n += nxt[i];
      }
   }
   else
   {
      ext_list[itel].npix = npix;
      ext_list[itel].nbsize = 0;
   }

#ifdef DEBUG_PIXEL_NB
   for (i=0; i<npix && i<30; i++)
   {
      printf("Pixel %d has packed %d direct neighbours, %d in second set, %d in third set, %d in extension list.\n",
         i, nb_lists[itel][0].nbsize>0 ? nb_lists[itel][0].pix_num_nb[i] : 0, 
         nb_lists[itel][1].nbsize>0 ? nb_lists[itel][1].pix_num_nb[i] : 0, 
         nb_lists[itel][2].nbsize>0 ? nb_lists[itel][2].pix_num_nb[i] : 0, 
         ext_list[itel].nbsize>0 ? ext_list[itel].pix_num_nb[i] : 0);
#ifdef DEBUG_PIXEL_NB2
      { int j;
        if ( nb_lists[itel][0].nbsize>0 && nb_lists[itel][0].pix_num_nb[i] > 0 )
        { printf("   direct:"); for (j=0; j<nb_lists[itel][0].pix_num_nb[i]; j++) 
            printf(" %d", nb_lists[itel][0].nblist[nb_lists[itel][0].pix_first_nb[i]+j]); printf("\n"); }
        if ( nb_lists[itel][1].nbsize>0 && nb_lists[itel][1].pix_num_nb[i] > 0 )
        { printf("   second:"); for (j=0; j<nb_lists[itel][1].pix_num_nb[i]; j++) 
            printf(" %d", nb_lists[itel][1].nblist[nb_lists[itel][1].pix_first_nb[i]+j]); printf("\n"); }
        if ( nb_lists[itel][2].nbsize>0 && nb_lists[itel][2].pix_num_nb[i] > 0 )
        { printf("   third:"); for (j=0; j<nb_lists[itel][2].pix_num_nb[i]; j++) 
            printf(" %d", nb_lists[itel][2].nblist[nb_lists[itel][2].pix_first_nb[i]+j]); printf("\n"); }
        if ( ext_list[itel].nbsize>0 && ext_list[itel].pix_num_nb[i] > 0 )
        { printf("   extension:"); for (j=0; j<ext_list[itel].pix_num_nb[i]; j++) 
            printf(" %d", ext_list[itel].nblist[ext_list[itel].pix_first_nb[i]+j]); printf("\n"); }
      }
#endif
   }
#endif

   return 0;
}

/* ----------------------------- store_camera_radius ---------------------- */

/* Determine the size of the camera in units of radians (i.e. for unit focal length). */

int store_camera_radius(CameraSettings *camset, int itel);

int store_camera_radius(CameraSettings *camset, int itel)
{
   int npix = camset->num_pixels, actpix = 0;
   int i;
   double sr = 0., rmx = 0.;
   
   if ( itel < 0 || itel >= H_MAX_TEL || npix < 2 )
      return -1;

   for (i=0; i<npix; i++)
   {
      if ( pixel_disabled[itel][i] )
         continue;
      actpix++;
      double x = camset->xpix[i];
      double y = camset->ypix[i];
      double r = sqrt(x*x+y*y);
      sr += r;
      if ( r > rmx )
         rmx = r; 
   }
   camera_radius_max[itel] = rmx / camset->flen;
   camera_radius_eff[itel] = 1.5 * sr / (double) actpix / camset->flen;
   
   if ( verbosity > 0 )
      printf("CT%d with %d pixels (%d active) has effective radius %5.3f deg, max. radius %5.3f deg.\n",
   	camset->tel_id, npix, actpix,
   	camera_radius_eff[itel]*(180./M_PI),
   	camera_radius_max[itel]*(180./M_PI) );

   return 0;
}

/* ----------------------- get_camera_radius ------------------------------ */

double get_camera_radius(int itel, int maxflag)
{
   if ( itel < 0 || itel >= H_MAX_TEL )
      return -1.;
   if ( maxflag )
      return camera_radius_max[itel];
   else
      return camera_radius_eff[itel];
}


/* -------------------------- select_calibration_channel ------------------ */

static int no_low_gain = 0, no_high_gain = 0;

/** Control if only low-gain or high-gain should get used instead of both.
 *
 *  @param chn 0 (both channels), 1 (only high gain), 2 (only low gain)
 */

void select_calibration_channel (int chn)
{
   if ( chn == 0 )
      no_low_gain = no_high_gain = 0;
   else if ( chn == 1 )
   {
      no_low_gain = 1;
      no_high_gain = 0;
   }
#if ( H_MAX_GAINS >= 2 )
   else if ( chn == 2 )
   {
      no_low_gain = 0;
      no_high_gain = 1;
   }
#endif
}

/* ----------------------- calibrate_amplitude ---------------------------- */

/** @short Calibrate amplitudes in all pixels of a camera. 
 *
 *  This function is operating only on pulse sums, either from normal raw data
 *  or from timing/pulse shape analysis. Use calibrate_pixel_amplitude() for
 *  calibration of individual samples.
 *
 *  @param hsdata Pointer to all available data and configurations.
 *  @param itel   Index of telescope in the relevant arrays (not the ID).
 *  @param flag_amp_tm 0: Use normal integrated amplitude.
 *                     1: Use integration around global peak position from
 *                        pulse shape analysis. May include all pixels or only selected.
 *                     2: Use integration around local peak position from
 *                        pulse shape analysis. Return 0 for pixels without
 *                        a fairly significant peak.
 *  @param clip_amp: if >0, any calibrated amplitude is clipped not to exceed this value [mean p.e.].
 */

int calibrate_amplitude(AllHessData *hsdata, int itel, 
   int flag_amp_tm, double clip_amp)
{
   int npix = hsdata->camera_set[itel].num_pixels;
   int i, j;
   AdcData *raw;
   TelEvent *te;
   LasCalData *lcal;
   TelMoniData *moni;
   int ns = 0;
   PixelTiming *pixtm = NULL;
   PixelList *tpl = NULL;
   int with_tm = 0, glob_only_selected = 1, check_triggered = 0;
   int tel_type;
   struct user_parameters *up;
   double calib_scale;
   int no_raw = 0;

   if ( itel < 0 || itel >= H_MAX_TEL )
      return -1;
   
   /* Internal calibration is to units of mean single-p.e. amplitude but for */
   /* better comparison with experiment typically use units where the most likely */
   /* single-p.e. amplitude got fitted ('single p.e. peak'). */
   tel_type = user_get_type(itel);
   up = user_get_parameters(tel_type);
   calib_scale = (up->d.calib_scale > 0. ? up->d.calib_scale : CALIB_SCALE);

   for (i=0; i<npix; i++)
      pixel_amp[itel][i] = 0.;
   pixel_sat[itel] = 0;

   te = &hsdata->event.teldata[itel];
   tpl = &te->trigger_pixels;
   raw = te->raw;
   if ( raw == NULL )
      no_raw = 1;
   else if ( !raw->known )
      no_raw = 1;
   if ( no_raw ) /* If no raw data is available but calibrated data we use that */
   {
      if ( te->pixcal != NULL && te->pixcal->known )
      {
         for ( i=0; i<npix && i<te->pixcal->num_pixels; i++ )
         {
            if ( te->pixcal->significant[i] )
               pixel_amp[itel][i] = te->pixcal->pixel_pe[i];
         }
         return 0;
      }
      else
         return -1; /* Neither raw nor calibrated data available */
   }

   /* Continue with calibrating the raw data. This is the preferred way since we are more flexible. */
   lcal = &hsdata->tel_lascal[itel];
   moni = &hsdata->tel_moni[itel];

   if ( flag_amp_tm && te->pixtm != NULL && te->pixtm->known && 
      te->pixtm->before_peak >= 0 && te->pixtm->after_peak >= 0 )
   {
      pixtm = te->pixtm;
      glob_only_selected = (pixtm->threshold < 0 ? 1 : 0);
      with_tm = 1;
   }
   
   for (i=0; i<npix; i++)
   {
      double npe, npe_hg=0., sig_hg;
#if ( H_MAX_GAINS >= 2 )
      double npe_lg=0., sig_lg;
#endif

      int significant = raw->significant[i];
      int hg_known = (no_high_gain && raw->num_gains > 1) ? 0 : raw->adc_known[HI_GAIN][i];
#if ( H_MAX_GAINS >= 2 )
      int lg_known = ((raw->num_gains>=2 && !no_low_gain) ? raw->adc_known[LO_GAIN][i] : 0);
#endif

      if ( with_tm )
      {
         if ( pixtm->timval[i][0] >= 0. /* We have a significant peak */
              || ((flag_amp_tm == 1) && (!glob_only_selected)) /* All global sums are stored. */ )
         {
            if ( flag_amp_tm == 1 ) /* Sum around global peak position */
            {
               sig_hg = pixtm->pulse_sum_glob[HI_GAIN][i] + 1e-12;
#if ( H_MAX_GAINS >= 2 )
               sig_lg = pixtm->pulse_sum_glob[LO_GAIN][i];
#endif
            }
            else /* Sum around local peak position (expected to be biased) */
            {
               sig_hg = pixtm->pulse_sum_loc[HI_GAIN][i];
#if ( H_MAX_GAINS >= 2 )
               sig_lg = pixtm->pulse_sum_loc[LO_GAIN][i];
#endif
            }
         }
         else
#if ( H_MAX_GAINS >= 2 )
            sig_lg = sig_hg = 0.;
#else
            sig_hg = 0.;
#endif

         npe_hg = sig_hg * lcal->calib[HI_GAIN][i];
#if ( H_MAX_GAINS >= 2 )
         npe_lg = sig_lg * lcal->calib[LO_GAIN][i];
#endif
      }
      else
      {
         sig_hg = hg_known ? (raw->adc_sum[HI_GAIN][i] -
	      moni->pedestal[HI_GAIN][i]) : 0.;
         npe_hg = sig_hg * lcal->calib[HI_GAIN][i];

#if ( H_MAX_GAINS >= 2 )
         sig_lg = lg_known ? (raw->adc_sum[LO_GAIN][i] -
	      moni->pedestal[LO_GAIN][i]) : 0.;
         npe_lg = sig_lg * lcal->calib[LO_GAIN][i];
#endif
      }

      if ( pixel_disabled[itel][i] )
      {
         raw->significant[i] = 0;
         // raw->adc_known[HI_GAIN][i] = 0;
         // raw->adc_known[LO_GAIN][i] = 0;
         pixel_amp[itel][i] = 0.;
         continue;
      }

      if ( !significant ) 
      	 npe = 0.;
      else if ( hg_known && sig_hg < 10000 && sig_hg > -1000 )
      	 npe = npe_hg;
#if ( H_MAX_GAINS >= 2 )
      else if ( raw->num_gains >= 2 )
      	 npe = npe_lg;
#endif
      else
         npe = npe_hg;

      if ( significant )
         ns++;

      if ( clip_amp > 0. )
         if ( npe > clip_amp )
         {
            npe = clip_amp;
            pixel_sat[itel]++;
         }

      /* npe is in units of 'mean photo-electrons' (unit = mean p.e. signal). */
      /* We convert to experimentalist's 'peak photo-electrons' */
      /* now (unit = most probable p.e. signal after experimental resolution). */
      /* Keep in mind: peak(10 p.e.) != 10*peak(1 p.e.) */
      pixel_amp[itel][i] = calib_scale * npe;
   }

   /* In case we want to keep the calibrated data for further use or storage */
   if ( te->pixcal != NULL )
   {
      PixelCalibrated *pc = te->pixcal;
      pc->tel_id = raw->tel_id;
      pc->num_pixels = npix;
      for (i=0; i<npix; i++)
      {
         pc->pixel_pe[i] = pixel_amp[itel][i];
         pc->significant[i] = raw->significant[i];
      }
      if ( raw->list_known )
      {
         pc->list_known = raw->list_known;
         pc->list_size = raw->list_size;
         for ( j=0; j<raw->list_size; j++ )
            pc->pixel_list[j] = raw->adc_list[j];
      }
      else
      {
         pc->list_known = 0;
         pc->list_size = 0;
      }
      pc->int_method = 0;
      if ( flag_amp_tm > 0 )
         pc->int_method = -flag_amp_tm;
      pc->known = raw->known;
   }

   /* Any of the previously triggered pixels may be disabled now. */
   if ( any_disabled[itel] )
      for ( i=0, j=0; i<tpl->pixels; i++ )
      {
         if ( pixel_disabled[itel][tpl->pixel_list[i]] )
         {
            // printf("Triggered pixel %d in tel. %d, event %d is disabled.\n", 
            //    tpl->pixel_list[i], te->tel_id, te->loc_count);
            tpl->pixel_list[i] = -1;
            check_triggered = 1;
         }
      }
   /* The triggered status of the telescope may have changed. */
   if ( check_triggered )
   {
      for ( i=0, j=0; i<tpl->pixels; i++ )
      {
         if ( tpl->pixel_list[i] >= 0 )
         {
            if ( i>j )
               tpl->pixel_list[j] = tpl->pixel_list[i];
            j++;
         }
      }
      tpl->pixels = j;

      /* Less pixels surviving than needed for telescope trigger. */
      if ( tpl->pixels < hsdata->pixel_set[itel].min_pixel_mult )
      {
         CentralEvent *ce = &hsdata->event.central;
         // printf("Telescope %d no longer triggered in event %d (%d) because of disabled pixels.\n", 
         //    te->tel_id, ce->glob_count, hsdata->mc_event.event);
         te->known = 0;
         /* If we have a list of triggered telescopes, remove this one. */
         for ( i=0, j=0; i<ce->num_teltrg; i++ )
         {
            if ( ce->teltrg_list[i] != te->tel_id )
            {
               if ( i>j )
                  ce->teltrg_list[j] = ce->teltrg_list[i];
               j++;
            }
         }
         ce->num_teltrg = j;
         /* If we have a list of telescopes with data, remove this one. */
         for ( i=0, j=0; i<ce->num_teldata; i++ )
         {
            if ( ce->teldata_list[i] != te->tel_id )
            {
               if ( i>j )
                  ce->teldata_list[j] = ce->teldata_list[i];
               j++;
            }
         }
         ce->num_teldata = j;
         /* Update old-style, small array bit patterns. */
         if ( itel < 32 )
         {
            int mask = ~(1 << itel);
            ce->teltrg_pattern &= mask;
            ce->teldata_pattern &= mask;
         }
      }
   }

   return 0;
}

/* ---------------------- calibrate_pixel_amplitude ----------------------- */

/** Calibrate a single pixel amplitude.
 *
 *  @param hsdata Pointer to all available data and configurations.
 *  @param itel   Index of telescope in the relevant arrays (not the ID).
 *  @param ipix   The pixel number (0 ... npix-1).
 *  @param flag_amp_tm 0: Use normal integrated amplitude.
 *                     1: Use integration around global peak position from
 *                        pulse shape analysis. May include all pixels or only selected.
 *                     2: Use integration around local peak position from
 *                        pulse shape analysis. Return 0 for pixels without
 *                        a fairly significant peak.
 *  @param itime -1: sum of samples of type as given in flag_amp_tm
 *               0...(nsamples-1): sample data (if available) for one time slice
 *  @param clip_amp: if >0, any calibrated amplitude is clipped not to exceed this value [mean p.e.].
 *
 * @return Pixel amplitude in peak p.e. units (based on conversion factor from H.E.S.S.).
 */

double calibrate_pixel_amplitude(AllHessData *hsdata, int itel, 
   int ipix, int flag_amp_tm, int itime, double clip_amp)
{
   int i = ipix, npix, significant, hg_known;
   double npe, sig_hg, npe_hg;
#if (H_MAX_GAINS >= 2)
   double sig_lg, npe_lg;
   int lg_known;
#endif
   AdcData *raw;
   TelEvent *te;
   LasCalData *lcal;
   TelMoniData *moni;
   int tel_type;
   struct user_parameters *up;
   double calib_scale;
   int no_raw = 0;

   if ( hsdata == NULL || itel < 0 || itel >= H_MAX_TEL )
      return 0.;

   /* Internal calibration is to units of mean single-p.e. amplitude but for */
   /* better comparison with experiment typically use units where the most likely */
   /* single-p.e. amplitude got fitted ('single p.e. peak'). */
   tel_type = user_get_type(itel);
   up = user_get_parameters(tel_type);
   calib_scale = (up->d.calib_scale > 0. ? up->d.calib_scale : CALIB_SCALE);

   if ( pixel_disabled[itel][ipix] )
      return 0.;
   npix = hsdata->camera_set[itel].num_pixels;
   if ( ipix < 0 || ipix >= npix )
      return 0.;
   te = &hsdata->event.teldata[itel];
   raw = te->raw;
   if ( raw == NULL )
      no_raw = 1;
   else if ( !raw->known )
      no_raw = 1;
   if ( no_raw ) /* If no raw data is available but calibrated data we use that */
   {
      if ( te->pixcal != NULL && te->pixcal->known )
      {
         if ( ipix >= 0 && ipix < npix && ipix < te->pixcal->num_pixels )
            if ( te->pixcal->significant[ipix] )
               return te->pixcal->pixel_pe[ipix];
      }
      return 0.; /* Neither raw nor calibrated data available for this pixel */
   }

   lcal = &hsdata->tel_lascal[itel];
   moni = &hsdata->tel_moni[itel];
 
   significant = raw->significant[i];

   hg_known = (no_high_gain && raw->num_gains>1) ? 0 : raw->adc_known[HI_GAIN][i];
#if (H_MAX_GAINS >= 2)
   if ( raw->num_gains >= 2 && !no_low_gain )
      lg_known = raw->adc_known[LO_GAIN][i];
   else
      lg_known = 0;
#endif

   if ( itime >= 0 && flag_amp_tm == 0 )
   {
      if ( raw->num_samples <= 1 || itime >= raw->num_samples || !significant )
         return 0.;
      /* For zero-suppressed sample mode data check relevant bit */
      if ( (raw->zero_sup_mode & 0x20) != 0 && (raw->significant[i] & 0x020) == 0 )
         return 0.;
      sig_hg = hg_known ? (raw->adc_sample[HI_GAIN][i][itime] -
           moni->pedestal[HI_GAIN][i]/(double)raw->num_samples) : 0;
      npe = npe_hg = sig_hg * lcal->calib[HI_GAIN][i];
#if (H_MAX_GAINS >= 2 )
      if ( lg_known )
      {
         sig_lg = raw->adc_sample[LO_GAIN][i][itime] -
            moni->pedestal[LO_GAIN][i]/(double)raw->num_samples;
         npe_lg = sig_lg * lcal->calib[LO_GAIN][i];
         /* FIXME: need to make the high/low switch-over point flexible */
         if ( hg_known && sig_hg < 2000 && sig_hg > -300 ) /* Lower thresholds than for sums, here assuming 12-bit ADC */
            npe = npe_hg;
         else if ( raw->num_gains >= 2 )
            npe = npe_lg;
         else
            npe = npe_hg;
      }
#endif
      return calib_scale * npe;
   }

   if ( flag_amp_tm && te->pixtm != NULL && te->pixtm->known && 
      te->pixtm->before_peak >= 0 && te->pixtm->after_peak >= 0 )
   {
      PixelTiming *pixtm = te->pixtm;
      int glob_only_selected = (pixtm->threshold < 0 ? 1 : 0);

      if ( pixtm->timval[i][0] >= 0. /* We have a significant peak */
           || ((flag_amp_tm == 1) && (!glob_only_selected)) /* All global sums are stored. */ )
      {
         if ( flag_amp_tm == 1 ) /* Sum around global peak position */
         {
            sig_hg = pixtm->pulse_sum_glob[HI_GAIN][i] + 1e-12 /* fuzz needed for camera_image*/ ;
#if (H_MAX_GAINS >= 2 )
            sig_lg = pixtm->pulse_sum_glob[LO_GAIN][i];
#endif
         }
         else /* Sum around local peak position (expected to be biased) */
         {
            sig_hg = pixtm->pulse_sum_loc[HI_GAIN][i];
#if (H_MAX_GAINS >= 2 )
            sig_lg = pixtm->pulse_sum_loc[LO_GAIN][i];
#endif
         }
      }
      else
#if (H_MAX_GAINS >= 2 )
         sig_lg = sig_hg = 0.;
#else
         sig_hg = 0.;
#endif

      npe_hg = sig_hg * lcal->calib[HI_GAIN][i];
#if (H_MAX_GAINS >= 2 )
      npe_lg = sig_lg * lcal->calib[LO_GAIN][i];
#endif
   }
   else
   {
      sig_hg = hg_known ? (raw->adc_sum[HI_GAIN][i] -
           moni->pedestal[HI_GAIN][i]) : 0.;
      npe_hg = sig_hg * lcal->calib[HI_GAIN][i];

#if (H_MAX_GAINS >= 2 )
      sig_lg = lg_known ? (raw->adc_sum[LO_GAIN][i] -
           moni->pedestal[LO_GAIN][i]) : 0.;
      npe_lg = sig_lg * lcal->calib[LO_GAIN][i];
#endif
   }

   if ( !significant ) 
      npe = 0.;
#if (H_MAX_GAINS >= 2 )
   /* FIXME: need to make the high/low switch-over point flexible */
   else if ( hg_known && sig_hg < 10000 && sig_hg > -1000 )
      npe = npe_hg;
   else if ( raw->num_gains >= 2 )
      npe = npe_lg;
#endif
   else
      npe = npe_hg;

   if ( clip_amp > 0. )
   {
      if ( npe > clip_amp )
         npe = clip_amp;
   }

   /* npe is in units of 'mean photo-electrons' (unit = mean p.e. signal.). */
   /* We convert to experimentalist's 'peak photo-electrons' */
   /* now (unit = most probable p.e. signal after experimental resolution). */
   /* Keep in mind: peak(10 p.e.) != 10*peak(1 p.e.) */
   return calib_scale * npe;
}

/* --------------------------- simple_integration -------------------------- */
/**
 *  @short Integrate sample-mode data (traces) over a common and fixed interval.
 *
 *  The integration window can be anywhere in the available length of the traces.
 *  Since the calibration function subtracts a pedestal that corresponds to the
 *  total length of the traces we may also have to add a pedestal contribution
 *  for the samples not summed up.
 *  No weighting of individual samples is applied.
 *
 *  @param hsdata Pointer to all available data and configurations.
 *  @param itel   Sequence number of the telescope being processed.
 *  @param nsum   Number of samples to sum up (is reduced if exceeding available length).
 *  @param nskip  Number of initial samples skipped (adapted such that interval fits into what is available).
 *                  Note: for multiple gains, this results in identical integration regions.
 */ 

static int simple_integration(AllHessData *hsdata, int itel, int nsum, int nskip);

static int simple_integration(AllHessData *hsdata, int itel, int nsum, int nskip)
{
   int isamp, ipix, igain;
   TelEvent *teldata = NULL;
   AdcData *raw;
   TelMoniData *moni;

   if ( hsdata == NULL || itel < 0 || itel >= H_MAX_TEL )
      return -1;

   teldata = &hsdata->event.teldata[itel];
   moni = &hsdata->tel_moni[itel];
   if ( !teldata->known || teldata->raw == NULL )
      return -1;
   raw = teldata->raw;
   if ( !raw->known )
      return -1;
   if ( raw->num_samples <= 1 || !(raw->known&2) )
      return 0;
   if ( nsum + nskip > raw->num_samples )
   {
      if ( nsum >= raw->num_samples )
      {
         nskip = 0;
         nsum = raw->num_samples;
      }
      else
         nskip = raw->num_samples - nsum;
   }
   for (igain=0; igain<raw->num_gains; igain++)
   {
      for (ipix=0; ipix<raw->num_pixels; ipix++)
      {
         /* For zero-suppressed sample mode data check relevant bit */
         if ( (raw->zero_sup_mode & 0x20) != 0 && (raw->significant[ipix] & 0x020) == 0 )
            raw->adc_sum[igain][ipix] = 0;
         else if ( raw->significant[ipix] && raw->adc_known[igain][ipix] )
         {
            int sum = 0;
            for ( isamp=0; isamp<nsum; isamp++ )
               sum += raw->adc_sample[igain][ipix][isamp+nskip];
            if ( nsum != raw->num_samples )
            {
               /* Keep in mind that the calibration functions subtract a sum pedestal */
               /* corresponding to num_samples bins. Add remaining pedestal. */
               sum += (int) ((raw->num_samples-nsum)*moni->pedestal[igain][ipix]/(double)raw->num_samples+0.5);
            }
            if ( integration_correction[itel][igain] > 0. )
               sum = (int)((sum-moni->pedestal[igain][ipix]) * 
                     integration_correction[itel][igain] + 
                     moni->pedestal[igain][ipix] + 0.5);
            raw->adc_sum[igain][ipix] = sum;
         }
         else
            raw->adc_sum[igain][ipix] = 0;
      }
   }
   return 0;
}

/* --------------------------- global_peak_integration -------------------------- */
/**
 *  @short Integrate sample-mode data (traces) over a common interval around a global signal peak.
 *
 *  The integration window can be anywhere in the available length of the traces.
 *  Since the calibration function subtracts a pedestal that corresponds to the
 *  total length of the traces we may also have to add a pedestal contribution
 *  for the samples not summed up.
 *  No weighting of individual samples is applied.
 *
 *  @param hsdata Pointer to all available data and configurations.
 *  @param itel   Sequence number of the telescope being processed.
 *  @param nsum   Number of samples to sum up (is reduced if exceeding available length).
 *  @param nbefore  Start the integration a number of samples before the peak,
 *                  as long as it fits into the available data range.
 *                  Note: for multiple gains, this results in identical integration regions.
 *  @param sigamp Amplitude in ADC counts above pedestal at which a signal is
 *                  considered as significant (separate for high gain/low gain).
 */ 

static int global_peak_integration(AllHessData *hsdata, int itel, int nsum, int nbefore, int *sigamp);

static int global_peak_integration(AllHessData *hsdata, int itel, int nsum, int nbefore, int *sigamp)
{
   int isamp, ipix, igain;
   TelEvent *teldata = NULL;
   AdcData *raw;
   TelMoniData *moni;
   int jpeak[H_MAX_PIX], ppeak[H_MAX_PIX], npeaks=0, peakpos_hg=-1;

   if ( hsdata == NULL || itel < 0 || itel >= H_MAX_TEL )
      return -1;

   teldata = &hsdata->event.teldata[itel];
   moni = &hsdata->tel_moni[itel];
   if ( !teldata->known || teldata->raw == NULL )
      return -1;
   raw = teldata->raw;
   if ( !raw->known )
      return -1;
   if ( raw->num_samples <= 1 || !(raw->known&2) )
      return 0;
   if ( nsum > raw->num_samples )
   {
      nsum = raw->num_samples;
   }
   for (igain=0; igain<raw->num_gains; igain++)
   {
      int peakpos = -1, start = 0;
      npeaks = 0;
      for (ipix=0; ipix<raw->num_pixels; ipix++)
      {
         raw->adc_sum[igain][ipix] = 0;
         /* For zero-suppressed sample mode data check relevant bit */
         if ( (raw->zero_sup_mode & 0x20) != 0 && (raw->significant[ipix] & 0x020) == 0 )
            continue;
         else if ( raw->significant[ipix] && raw->adc_known[igain][ipix] )
         {
            int pedsamp = (int) (moni->pedestal[igain][ipix]/(double)raw->num_samples+0.5);
            int significant = 0, ipeak = -1, p=0;
            for ( isamp=0; isamp<raw->num_samples; isamp++ )
            {
               if ( raw->adc_sample[igain][ipix][isamp] - pedsamp >= sigamp[igain] )
               {
                  int isamp2;
                  significant = 1;
                  ipeak = isamp;
                  p = raw->adc_sample[igain][ipix][isamp];
                  for ( isamp2=isamp+1; isamp2<raw->num_samples; isamp2++ )
                  {
                     if ( raw->adc_sample[igain][ipix][isamp2] > p )
                     {
                        ipeak = isamp2;
                        p = raw->adc_sample[igain][ipix][isamp2];
                     }
                  }
                  break;
               }
            }
            if ( significant )
            {
               jpeak[npeaks] = ipeak;
               ppeak[npeaks] = p - pedsamp;
               npeaks++;
            }
         }
      }
      if ( npeaks > 0 )
      {
         double ps=0., pjs=0.;
         int ipeak;
         for (ipeak=0; ipeak<npeaks; ipeak++)
         {
            ps += ppeak[ipeak];
            pjs += ppeak[ipeak] * jpeak[ipeak];
         }
         if ( ps > 0. )
            peakpos = (int) (pjs/ps+0.5);
         else
            peakpos = 0;
         start = peakpos - nbefore;
         if ( start < 0 )
            start = 0;
         if ( start + nsum > raw->num_samples )
            start = raw->num_samples - nsum;
         if ( igain == 0 )
            peakpos_hg = peakpos;
      }
      /* If low-gain does not have a significant peak of its own, use the high-gain peak */
      else if ( igain > 0 && peakpos_hg >= 0 )
         peakpos = peakpos_hg;
      else /* No idea where to sum up. We could leave it all empty or use an ad-hoc value. */
         continue; /* Leaving it empty make it obvious that the threshold may be too low */
         // peakpos = (2*raw->num_samples)/3; /* This ad-hoc value may still work in many cases */

      start = peakpos - nbefore;
      if ( start < 0 )
         start = 0;
      if ( start + nsum > raw->num_samples )
         start = raw->num_samples - nsum;
#if 0
printf("** Event %d, telescope %d: npeaks=%d, peakpos = %d, start = %d, correction = %f\n",
hsdata->event.central.glob_count,
teldata->tel_id, npeaks, peakpos, start, integration_correction[itel][0]);
#endif
      for (ipix=0; ipix<raw->num_pixels; ipix++)
      {
         /* For zero-suppressed sample mode data check relevant bit */
         if ( (raw->zero_sup_mode & 0x20) != 0 && (raw->significant[ipix] & 0x020) == 0 )
            raw->adc_sum[igain][ipix] = 0;
         else if ( raw->significant[ipix] && raw->adc_known[igain][ipix] )
         {
            int sum = 0;
            for ( isamp=0; isamp<nsum; isamp++ )
               sum += raw->adc_sample[igain][ipix][isamp+start];
            if ( nsum != raw->num_samples )
            {
               /* Keep in mind that the calibration functions subtract a sum pedestal */
               /* corresponding to num_samples bins. Add remaining pedestal. */
               sum += (int) ((raw->num_samples-nsum)*moni->pedestal[igain][ipix]/(double)raw->num_samples+0.5);
            }
            if ( integration_correction[itel][igain] > 0. )
               sum = (int)((sum-moni->pedestal[igain][ipix]) * 
                     integration_correction[itel][igain] + 
                     moni->pedestal[igain][ipix] + 0.5);
            raw->adc_sum[igain][ipix] = sum;
         }
         else
            raw->adc_sum[igain][ipix] = 0;
      }
   }
   return 0;
}

/* --------------------------- local_peak_integration -------------------------- */
/**
 *  @short Integrate sample-mode data (traces) around a pixel-local signal peak.
 *
 *  The integration window can be anywhere in the available length of the traces.
 *  Since the calibration function subtracts a pedestal that corresponds to the
 *  total length of the traces we may also have to add a pedestal contribution
 *  for the samples not summed up.
 *  No weighting of individual samples is applied.
 *
 *  @param hsdata Pointer to all available data and configurations.
 *  @param itel   Sequence number of the telescope being processed.
 *  @param nsum   Number of samples to sum up (is reduced if exceeding available length).
 *  @param nbefore  Start the integration a number of samples before the peak,
 *                  as long as it fits into the available data range.
 *                  Note: for multiple gains, this may result in identical integration regions (depending on signal).
 *  @param sigamp Amplitude in ADC counts above pedestal at which a signal is
 *                  considered as significant (separate for high gain/low gain).
 */ 

static int local_peak_integration(AllHessData *hsdata, int itel, int nsum, int nbefore, int *sigamp);

static int local_peak_integration(AllHessData *hsdata, int itel, int nsum, int nbefore, int *sigamp)
{
   int isamp, ipix, igain;
   TelEvent *teldata = NULL;
   AdcData *raw;
   TelMoniData *moni;
   int peakpos = -1, start = 0, peakpos_hg=-1;

   if ( hsdata == NULL || itel < 0 || itel >= H_MAX_TEL )
      return -1;

   teldata = &hsdata->event.teldata[itel];
   moni = &hsdata->tel_moni[itel];
   if ( !teldata->known || teldata->raw == NULL )
      return -1;
   raw = teldata->raw;
   if ( !raw->known )
      return -1;
   if ( raw->num_samples <= 1 || !(raw->known&2) )
      return 0;
   if ( nsum > raw->num_samples )
   {
      nsum = raw->num_samples;
   }

   for (ipix=0; ipix<raw->num_pixels; ipix++)
   {
      for (igain=0; igain<raw->num_gains; igain++)
         raw->adc_sum[igain][ipix] = 0;
      /* For zero-suppressed sample mode data check relevant bit */
      if ( (raw->zero_sup_mode & 0x20) != 0 && (raw->significant[ipix] & 0x020) == 0 )
         continue;
      if ( raw->significant[ipix] && raw->adc_known[HI_GAIN][ipix] )
      {
         int pedsamp = (int) (moni->pedestal[HI_GAIN][ipix]/(double)raw->num_samples+0.5);
         int significant = 0, ipeak = -1, p=0;
         for ( isamp=0; isamp<raw->num_samples; isamp++ )
         {
            if ( raw->adc_sample[HI_GAIN][ipix][isamp] - pedsamp >= sigamp[HI_GAIN] )
            {
               int isamp2;
               significant = 1;
               ipeak = isamp;
               p = raw->adc_sample[HI_GAIN][ipix][isamp];
               for ( isamp2=isamp+1; isamp2<raw->num_samples; isamp2++ )
               {
                  if ( raw->adc_sample[HI_GAIN][ipix][isamp2] > p )
                  {
                     ipeak = isamp2;
                     p = raw->adc_sample[HI_GAIN][ipix][isamp2];
                  }
               }
               break;
            }
         }
         peakpos = peakpos_hg = ipeak;
         if ( significant  && peakpos >= 0 )
         {
            int sum = 0;
            start = peakpos - nbefore;
            if ( start < 0 )
               start = 0;
            if ( start + nsum > raw->num_samples )
               start = raw->num_samples - nsum;
            for ( isamp=0; isamp<nsum; isamp++ )
               sum += raw->adc_sample[HI_GAIN][ipix][isamp+start];
            if ( nsum != raw->num_samples )
            {
               /* Keep in mind that the calibration functions subtract a sum pedestal */
               /* corresponding to sum_samples bins. Add remaining pedestal. */
               sum += (int) ((raw->num_samples-nsum)*moni->pedestal[HI_GAIN][ipix]/(double)raw->num_samples+0.5);
            }
            if ( integration_correction[itel][HI_GAIN] > 0. )
               sum = (int) ((sum-moni->pedestal[HI_GAIN][ipix]) * 
                     integration_correction[itel][HI_GAIN] + 
                     moni->pedestal[HI_GAIN][ipix] + 0.5);
            raw->adc_sum[HI_GAIN][ipix] = sum;
         }
      }
#if (H_MAX_GAINS > 1)
      if ( raw->num_gains > 1 && raw->significant[ipix] && raw->adc_known[LO_GAIN][ipix] )
      {
         /* Normally, low gain would be integrated over the same interval */
         /* but the high-gain signal may be missing or in complete saturation. */
         /* Thus we first try to see if the low-gain channel has a significant signal by itself. */
         int pedsamp = (int) (moni->pedestal[LO_GAIN][ipix]/(double)raw->num_samples+0.5);
         int significant = 0, ipeak = -1, p=0;
         for ( isamp=0; isamp<raw->num_samples; isamp++ )
         {
            if ( raw->adc_sample[LO_GAIN][ipix][isamp] - pedsamp >= sigamp[LO_GAIN] )
            {
               int isamp2;
               significant = 1;
               ipeak = isamp;
               p = raw->adc_sample[LO_GAIN][ipix][isamp];
               for ( isamp2=isamp+1; isamp2<raw->num_samples; isamp2++ )
               {
                  if ( raw->adc_sample[LO_GAIN][ipix][isamp2] > p )
                  {
                     ipeak = isamp2;
                     p = raw->adc_sample[LO_GAIN][ipix][isamp2];
                  }
               }
               break;
            }
         }
         if ( significant )
            peakpos = ipeak;
         else
            peakpos = peakpos_hg;
         if ( peakpos >= 0 )
         {
            int sum = 0;
            start = peakpos - nbefore;
            if ( start < 0 )
               start = 0;
            if ( start + nsum > raw->num_samples )
               start = raw->num_samples - nsum;
            for ( isamp=0; isamp<nsum; isamp++ )
               sum += raw->adc_sample[LO_GAIN][ipix][isamp+start];
            if ( nsum != raw->num_samples )
            {
               /* Keep in mind that the calibration functions subtract a sum pedestal */
               /* corresponding to num_samples bins. Add remaining pedestal. */
               sum += (int) ((raw->num_samples-nsum)*moni->pedestal[LO_GAIN][ipix]/(double)raw->num_samples+0.5);
            }
            if ( integration_correction[itel][LO_GAIN] > 0. )
               sum = (int)((sum-moni->pedestal[LO_GAIN][ipix]) * 
                     integration_correction[itel][LO_GAIN] + 
                     moni->pedestal[LO_GAIN][ipix] + 0.5);
            raw->adc_sum[LO_GAIN][ipix] = sum;
         }
      }
#endif
   }
   return 0;
}


/* --------------------------- nb_peak_integration -------------------------- */
/**
 *  @short Integrate sample-mode data (traces) around a peak in the signal sum of neighbouring pixels.
 *
 *  The integration window can be anywhere in the available length of the traces.
 *  Since the calibration function subtracts a pedestal that corresponds to the
 *  total length of the traces we may also have to add a pedestal contribution
 *  for the samples not summed up.
 *  No weighting of individual samples is applied.
 *
 *  @param hsdata Pointer to all available data and configurations.
 *  @param lwt    Weight of the local pixel (0: peak from neighbours only,
 *                  1: local pixel counts as much as any neighbour).
 *  @param itel   Sequence number of the telescope being processed.
 *  @param nsum   Number of samples to sum up (is reduced if exceeding available length).
 *  @param nbefore  Start the integration a number of samples before the peak,
 *                  as long as it fits into the available data range.
 *                  Note: for multiple gains, this results in identical integration regions.
 *  @param sigamp Amplitude in ADC counts above pedestal at which a signal is
 *                  considered as significant (separate for high gain/low gain).
 */ 

static int nb_peak_integration(AllHessData *hsdata, int lwt, int itel, int nsum, int nbefore, int *sigamp);

static int nb_peak_integration(AllHessData *hsdata, int lwt, int itel, int nsum, int nbefore, int *sigamp)
{
   int isamp, ipix, igain, ipeak, p;
   TelEvent *teldata = NULL;
   AdcData *raw;
   TelMoniData *moni;
   int peakpos = -1, start = 0, peakpos_hg=-1;
   struct camera_nb_list *nbl;

   if ( hsdata == NULL || itel < 0 || itel >= H_MAX_TEL )
      return -1;

   teldata = &hsdata->event.teldata[itel];
   moni = &hsdata->tel_moni[itel];
   if ( !teldata->known || teldata->raw == NULL )
      return -1;
   raw = teldata->raw;
   if ( !raw->known )
      return -1;
   if ( raw->num_samples <= 1 || !(raw->known&2) )
      return 0;
   if ( nsum > raw->num_samples )
   {
      nsum = raw->num_samples;
   }

   /* For this integration scheme we need the list of neighbours early on */
   if ( nb_lists[itel][0].nblist == NULL )
   {
      find_neighbours(&hsdata->camera_set[itel],itel);
   }
   
   if ( nb_lists[itel][0].nblist == NULL || 
        nb_lists[itel][0].nbsize <= 0 )
      return -1;
   nbl = &nb_lists[itel][0];

   for (ipix=0; ipix<raw->num_pixels; ipix++)
   {
      int nb_samples[H_MAX_SLICES], inb, knb=0;
      for (igain=0; igain<raw->num_gains; igain++)
         raw->adc_sum[igain][ipix] = 0;
      /* For zero-suppressed sample mode data check relevant bit of the current pixel */
      if ( (raw->zero_sup_mode & 0x20) != 0 && (raw->significant[ipix] & 0x020) == 0 )
         continue;
      for (isamp=0; isamp<raw->num_samples; isamp++ )
         nb_samples[isamp] = 0;
      for ( inb=0; inb<nbl->pix_num_nb[ipix]; inb++ )
      {
         int ipix_nb = nbl->nblist[nbl->pix_first_nb[ipix]+inb];
         /* For zero-suppressed sample mode data ALSO check relevant bit of the neighbour */
         if ( (raw->zero_sup_mode & 0x20) != 0 && (raw->significant[ipix_nb] & 0x020) == 0 )
            continue;
         if ( raw->significant[ipix_nb] && raw->adc_known[HI_GAIN][ipix_nb] )
         {
            /* No need for (flat) pedestal subtraction here since we just look for the peak position. */
            for (isamp=0; isamp<raw->num_samples; isamp++ )
               nb_samples[isamp] += raw->adc_sample[HI_GAIN][ipix_nb][isamp];
            knb++;
         }
      }
      if ( lwt > 0 )
      {
         if ( raw->significant[ipix] && raw->adc_known[HI_GAIN][ipix] )
         {
            /* This plain summation assumes pixels have roughly similar response */
            for (isamp=0; isamp<raw->num_samples; isamp++ )
               nb_samples[isamp] += raw->adc_sample[HI_GAIN][ipix][isamp] * lwt;
            knb++;
         }
      }
      if ( knb == 0 ) /* No integration window available for truely isolated pixels */
         continue;
      ipeak = 0;
      p = nb_samples[0];
      for ( isamp=1; isamp<raw->num_samples; isamp++ )
      {
         if ( nb_samples[isamp] > p )
         {
            p = nb_samples[isamp];
            ipeak = isamp;
         }
      }
      peakpos = peakpos_hg = ipeak;
      start = peakpos - nbefore;
      if ( start < 0 )
         start = 0;
      if ( start + nsum > raw->num_samples )
         start = raw->num_samples - nsum;
      if ( raw->significant[ipix] && raw->adc_known[HI_GAIN][ipix] )
      {
         // int pedsamp = (int) (moni->pedestal[HI_GAIN][ipix]/(double)raw->num_samples+0.5);
         {
            int sum = 0;
            for ( isamp=0; isamp<nsum; isamp++ )
               sum += raw->adc_sample[HI_GAIN][ipix][isamp+start];
            if ( nsum != raw->num_samples )
            {
               /* Keep in mind that the calibration functions subtract a sum pedestal */
               /* corresponding to num_samples bins. Add remaining pedestal. */
               sum += (int) ((raw->num_samples-nsum)*moni->pedestal[HI_GAIN][ipix]/(double)raw->num_samples+0.5);
            }
            if ( integration_correction[itel][HI_GAIN] > 0. )
               sum = (int)((sum-moni->pedestal[HI_GAIN][ipix]) * 
                     integration_correction[itel][HI_GAIN] + 
                     moni->pedestal[HI_GAIN][ipix] + 0.5);
            raw->adc_sum[HI_GAIN][ipix] = sum;
         }
      }
#if (H_MAX_GAINS > 1)
      if ( raw->num_gains > 1 && raw->significant[ipix] && raw->adc_known[LO_GAIN][ipix] )
      {
         /* Low gain is integrated over the same interval here. */
         // int pedsamp = (int) (moni->pedestal[LO_GAIN][ipix]/(double)raw->num_samples+0.5);
         {
            int sum = 0;
            for ( isamp=0; isamp<nsum; isamp++ )
               sum += raw->adc_sample[LO_GAIN][ipix][isamp+start];
            if ( nsum != raw->num_samples )
            {
               /* Keep in mind that the calibration functions subtract a sum pedestal */
               /* corresponding to sum_samples bins. Add remaining pedestal. */
               sum += (int) ((raw->num_samples-nsum)*moni->pedestal[LO_GAIN][ipix]/(double)raw->num_samples+0.5);
            }
            if ( integration_correction[itel][LO_GAIN] > 0. )
               sum = (int)((sum-moni->pedestal[LO_GAIN][ipix]) * 
                     integration_correction[itel][LO_GAIN] + 
                     moni->pedestal[LO_GAIN][ipix] + 0.5);
            raw->adc_sum[LO_GAIN][ipix] = sum;
         }
      }
#endif
   }
   return 0;
}

/* ------------------------------ gradient_integration -------------------------- */
/** @short Fit gradient of pixel pulse peak times along image and evaluate
 *         the fitted line for getting the time around which pulses get integrated.
 *
 *  There are basically three problems: a) bootstrap problem for finding significant pixels,
 *  b) robustness of the fit in case of pixels that don't follow the time gradient, and
 *  c) what to do with pixels that have a large enough signal at a time not consistent
 *  with the fitted line.
 */

static int gradient_integration(AllHessData *hsdata, int itel, int nsum, int nbefore, int *sigamp);

static int gradient_integration(AllHessData *hsdata, int itel, int nsum, int nbefore, int *sigamp)
{
#if 0
   int isamp, ipix, igain, ipeak, p;
   TelEvent *teldata = NULL;
   AdcData *raw;
   TelMoniData *moni;
   // int peakpos = -1, start = 0, peakpos_hg=-1;
   // struct camera_nb_list *nbl;

   if ( hsdata == NULL || itel < 0 || itel >= H_MAX_TEL )
      return -1;

   teldata = &hsdata->event.teldata[itel];
   raw = teldata->raw;
   moni = &hsdata->tel_moni[itel];
   if ( !teldata->known || teldata->raw == NULL )
      return -1;
   raw = teldata->raw;
   if ( !raw->known )
      return -1;
   if ( raw->num_samples <= 1 || !(raw->known&2) )
      return 0;
   if ( nsum > raw->num_samples )
   {
      nsum = raw->num_samples;
   }
#endif
   /* FIXME: Not really implemented yet ... */
   return -1;
}

#define WITH_PZPSA 1

#ifdef WITH_PZPSA

/* This is reformatted code contributed by Thomas Kihm for FlashCam */

/* ---------------- PzpsaSmoothUpsampleU16 ---------------- */
/**
   @short Upsample (expand the n input values to us samples each)
   Subtract baseline bl and correct for a single pole 
   decay with the decay time pz and smooth the resulting trace 
   with two moving averages with a width of us. 
   The output is placed in array op and returns the new number 
   of samples (n*us). 
   
   This function derived from code by T.Kihm, using uint16_t for
   input array element type and double for output.
   Example: PzpsaSmoothUpsampleU16(50,4,tti,0.,mpz,tto,&mxop,&imxop);
   
   @param n   Number of elements in input array ip
   @param us  Upsampling factor (use '4' to upsample from 250 MHZ to one GHz).
   @param ip  Pointer to input array of ADC raw data of type uint16_t
   @param bl  Baseline (pedestal) on input per sample
   @param pz  Pole-zero compensation factor in differencing (0<=pz<=1)
   @param op  Pointer to output array of type double
   @param max Maximum content in output array (only filled if not NULL)
   @param at  Position of maximum bin in output array (only filled if not NULL)
 */

static int PzpsaSmoothUpsampleU16 (int n, int us, uint16_t *ip, double bl, 
   double pz, double *op, double *max, int *at);

static int PzpsaSmoothUpsampleU16 (int n, int us, uint16_t *ip, double bl, 
   double pz, double *op, double *max, int *at)
{
   int   i, i1;          ///< running indices
   double v2, v1;        ///< the next and prev. input samples
   double sum1, sum2;    ///< the running sum of 1.st and 2.nd average
   double tmp;           ///< a temp var for intermediate copy 
   double pzc2, pzc1;    ///< the next and prev. pz corrected value
   double *out1 = op;    ///< the out pointer of the first runsum
   double *out2 = op;    ///< the out pointer of the second runsum
   double mult = 1./(us*us); ///< the multiplier to correct the two runsums
   double peakmax = -1e30; ///< peak maximum
   int   peakat  = 0;   ///< peak position

   v1 = v2 = (ip[0]-bl)*mult; 
   pzc2 = pzc1 = v2;
   sum1 = pzc2*us;
   sum2 = sum1*us;

   for ( i=0; i<us; i++) 
      *out1++ = sum1; 

   for (i=1; i<n; i++) 
   {
      v2 = (ip[i]-bl)*mult; 
      pzc2 = (v2-v1); 
      v1 = v2;
      for (i1=0; i1<us; i1++) 
      {
         sum1 += pzc2 - pzc1*pz; 
         *out1++ = sum1;
         tmp = *out2; 
         *out2++ = sum2; 
         if (sum2 > peakmax)
         {
            peakmax = sum2;
            peakat = out2-op-1;
         }
         sum2 += sum1 - tmp; 
      }
      pzc1 = pzc2;
   }  
   n*=us;

   for (v2=op[n-1]; out2<(op+n); ) 
   {
      tmp = *out2; 
      *out2++ = sum2;
      if (sum2 > peakmax)
      {
         peakmax = sum2;
         peakat = out2 - op - 1;
      }
      sum2 += v2 - tmp; 
   }

   if (max != NULL)
      *max = peakmax; 
   if (at != NULL)
      *at = peakat;

   return n; 
}


/* --------------------- PzpsaPeakProperty ---------------------- */
/**
    @short Calculates the peak property of the signal in (n samples) 
    at position pos. 

    The signal is integrated from sample pos-w to pos+w and the 
    result is stored in intsum. 

    The cog is the center of gravity calculated by the area 
    above the minumum of the signal from pos-w to pos+w

    Returns a quality value for the signal which is defined as
    in[pos]-(in[start]+in[stop])/2. Negativ values indicate that
    no positive signal was found.
 */

static double PzpsaPeakProperty (int n, double *in, int pos, int w, double *intsum, double *cog);

static double PzpsaPeakProperty (int n, double *in, int pos, int w, double *intsum, double *cog)
{
   int i;
   double v, sum, isum, min;
   int start = pos - w;
   int stop  = pos + w;
   if (start < 0) 
      start = 0;
   if (stop >= n) 
      stop = n-1;

   for (sum=0, i=start; i<=stop; i++)
      sum += in[i];

   if (intsum != NULL)
      *intsum = sum;

   for (min=in[start], i=start; i<=stop; i++)
      min = (in[i]<min) ? in[i] : min;

   for (isum=sum=0, i=start; i<=stop; i++)
   {
      v = in[i] - min;
      isum += v*i;
      sum += v;
   }
   if (sum != 0.)
      v = isum/sum; 
   else 
      v = 0.; 
   if (cog) 
      *cog = v; 
   return in[pos] - 0.5*(in[start]+in[stop]); 
} 
  

#endif

/* ------------------------------ nb_fc_shaped_peak_integration -------------------------- */

/** @short Pulse integration based on peaks in neighbour pixel signals after FlashCam-style
    pulse shaping.

    Basically like nb_peak_integration for lwt=0 but pulses are all
    upscaled in sampling frequency by a factor of four and one several variants
    for FlashCam-style pulse shaping is applied first. Signal extraction = integration
    also allows for different variants.
    There are actually way more variants available than necessary,
    intended for evaluation and testing.

    Note that the psopt parameter is specified with the '--integration-window'
    command line option as the third value. (Recommended values for the first two
    are 1,0 (=nsum,nbefore). Nsum=0 means nsum=1.)
    Interpret psopt as decimal MHTO (with M=psopt/1000, H=(psopt%1000)/100, T=(psopt%100)/10, O=psopt%10):
      O = -1 : Full pzpsa shaping and peak finding over full readout range, no neighbours involved.
               This results in a significant bias for positive NSB fluctuations.
           0 : Full pzpsa shaping but peak finding in signal of neighbours,
               avoiding the beginning (first 7) and end (last 3) of the upsampled
               signal because these are noisier and result in artifacts. (OK to use)
           9 : Like '0' but include the beginning and end for peak finding. (better use 0 or 1)
           1 : Like '9' but do own differencing to have smooth start and end. (recommended)
           2 : Like '1' but do differencing between second-to-next original samples.
           3 : Like '1' but do pulse shaping with own, more explicit code
               (differs in the beginning and the end but otherwise the same).
           4 : Like '2' but do pulse shaping with own, more explicit code.
      T =  0 : Use integration from nbefore the peak for nsum upsampled samples
               and determine the pixel timing as the peak position close to the
               peak times in the signal of neighbour pixels (except for O = -1).
        >  0 : Use the PzpsaPeakProperty code for summation and center-of-gravity
               determination of peaks, with T as width parameter. Nsum and nbefore are ignored.
      H =  0 : Summation region is entirely determined by the peak in the signal
               of neighbourin pixels, without any bias for NSB fluctuations.
           1 : Summation region allows for small adjustement in peak position
               in the signal of the pixel itself. Small NSB bias.
      M =  0 : Not touching pixel timing structure.
           1 : Re-evaluate and refill pixel timing from shaped signals and peaks.
               Unless a new 'integration threshold' is given, the old threshold
               for significant pixel timings gets re-used (but pixel is list still new).
 *
 *  @param hsdata Pointer to all available data and configurations.
 *  @param itel   Sequence number of the telescope being processed.
 *  @param nsum   Number of samples to sum up (is reduced if exceeding available length).
 *  @param nbefore  Start the integration a number of samples before the peak,
 *                  as long as it fits into the available data range.
 *                  Note: for multiple gains, this may result in identical integration regions (depending on signal).
 *  @param sigamp (not used)
 *  @param psopt  Pulse shaping option as described
 *  @param ithr   Integration threshold in ADC counts gets actually used for significance in pixel timing.
 *
 *  @return 0 (OK), -1 (error)
 */

static int nb_fc_shaped_peak_integration(AllHessData *hsdata, int itel, int nsum, int nbefore, int *sigamp, int psopt, int ithr);

static int nb_fc_shaped_peak_integration(AllHessData *hsdata, int itel, int nsum, int nbefore, int *sigamp, int psopt, int ithr)
{
   int isamp, ipix, igain, ipeak;
   TelEvent *teldata = NULL;
   AdcData *raw;
   TelMoniData *moni;
   int peakpos = -1, start = 0, nsamp4 = -1;
   struct camera_nb_list *nbl;
   static double *buffer = NULL;
   static size_t bfsize = 0;
   size_t bfreq;
   size_t off_gain, off_pix;
   double mpz1 = 0.758, mpz2 = 0.758*0.758; /* Hardcoded to match FlashCam pulse fall-off */
   /* double sc1 = 0.25/(1.-mpz1), sc2 = 0.25/(1.-mpz2); // for equal area */
   double sc1 = 1.0, sc2 = 1.0; /* Needs to match values in set_integration_correction */
   double pixtim_thr = 1e10;
   int pixtim_flag = psopt/1000;
   if ( pixtim_flag )
      psopt = psopt%1000;
   
   if ( hsdata == NULL || itel < 0 || itel >= H_MAX_TEL )
      return -1;

   teldata = &hsdata->event.teldata[itel];
   moni = &hsdata->tel_moni[itel];

   if ( !teldata->known || teldata->raw == NULL )
      return -1;

   if (  teldata->pixtm == NULL )
      pixtim_flag = 0;
   else if ( pixtim_flag )
   {
// printf("Prepare to refill pixel timing structure ...\n");
      pixtim_thr = fabs((double) teldata->pixtm->threshold);
      /* If any --integration-threshold was given, we use it here: */
      if ( ithr > 0 )
         pixtim_thr = ithr;
      teldata->pixtm->threshold = -pixtim_thr;
      teldata->pixtm->list_type = 0; /* Re-evaluate pixel list when writing */
      teldata->pixtm->list_size = 0;
      teldata->pixtm->before_peak = nbefore;
      teldata->pixtm->after_peak = nsum - nbefore - 1;
      teldata->pixtm->num_types = 7;
      teldata->pixtm->time_type[0] = 1;
      teldata->pixtm->time_type[1] = 2;
      teldata->pixtm->time_type[2] = 2;
      teldata->pixtm->time_type[3] = 2;
      teldata->pixtm->time_type[4] = 4;
      teldata->pixtm->time_type[5] = 4;
      teldata->pixtm->time_type[6] = 5;
      teldata->pixtm->time_level[0] = 1.0;
      teldata->pixtm->time_level[1] = 0.2;
      teldata->pixtm->time_level[2] = 0.5;
      teldata->pixtm->time_level[3] = 0.8;
      teldata->pixtm->time_level[4] = 0.5;
      teldata->pixtm->time_level[5] = 0.2;
      teldata->pixtm->time_level[6] = teldata->pixtm->threshold;
      teldata->pixtm->known = 0; /* Set to known when done */
   }

   raw = teldata->raw;
   if ( !raw->known )
      return -1;
   if ( raw->num_samples <= 1 || !(raw->known&2) )
      return 0;
   nsamp4 = 4*raw->num_samples;
   if ( nsum <= 0 )
      nsum = 1;
   if ( nsum > nsamp4 )
   {
      nsum = nsamp4;
   }
   /* Required buffer size (in bytes) for pulse shaping of this camera data */
   bfreq = nsamp4 * raw->num_pixels * raw->num_gains * sizeof(double);
   if ( bfreq > bfsize ) /* Need to (re-)allocate? */
   {
      double *bfb = buffer;
      if ( (buffer = (double *) realloc(buffer,bfreq)) == NULL )
      {
         buffer = bfb;
         return -1;
      }
      bfsize = bfreq;
   }
   /* Offsets per pixel and per gain (in doubles) */
   off_pix = nsamp4;
   off_gain = off_pix * raw->num_pixels;

   for ( igain=0; igain<raw->num_gains; igain++ )
   {
      double *bg = buffer + igain*off_gain, *bpx;
      double ped; ///< Pedestal in raw signal, per sample.
      double bpx0[4*H_MAX_SLICES], bpx1[4*H_MAX_SLICES];
      uint16_t *smp;

      for ( ipix=0; ipix<raw->num_pixels; ipix++ )
      {
         /* Initialize pulse sum to (sum) pedestal in case of any 'continue' statement ... */
         raw->adc_sum[igain][ipix] = (uint32_t) (moni->pedestal[igain][ipix]+0.5);

         if ( ! raw->significant[ipix] || ! raw->adc_known[igain][ipix] )
            continue;
         if ( (raw->zero_sup_mode & 0x20) != 0 && (raw->significant[ipix] & 0x020) == 0 )
            continue;

         bpx = bg + off_pix*ipix;
         ped = moni->pedestal[igain][ipix] / (double) raw->num_samples;
         smp = &raw->adc_sample[igain][ipix][0];

#ifdef WITH_PZPSA
         if ( psopt%10 == 0 || psopt%10 == 9 ) /* Pulse shaping completely with pzpsa code but no peak detection */
         {
            PzpsaSmoothUpsampleU16(raw->num_samples,4,smp,ped,mpz1,bpx,NULL,NULL);
#if 0
            /* Code for debugging and testing only ... */
            if ( teldata->pixtm->timval[ipix][0] >= 0 &&
                 teldata->pixtm->pulse_sum_glob[0][ipix] >= 200 )
            {
               for ( isamp=0; isamp<nsamp4; isamp++ )
                  printf("==> %4d %3d  %f\n", ipix, isamp, bpx[isamp]);
               exit(9);
            }
#endif
         }
         else if ( psopt%10 == 1 ) /* Pulse shaping with pzpsa code but own differentiation
                                      (with pole-zero) which differs from psopt==0 in up to
                                      seven slices at the beginning and up to three at the end */ 
         {
            PzpsaSmoothUpsampleU16(raw->num_samples,4,smp,ped,0.,bpx1,NULL,NULL);
            /* Smooth start, no sudden peak at start of window */
            bpx[0] = sc1 * ( bpx1[0] - mpz1*bpx1[0] ) * 0.0;
            bpx[1] = sc1 * ( bpx1[1] - mpz1*bpx1[0] ) * 0.25;
            bpx[2] = sc1 * ( bpx1[2] - mpz1*bpx1[0] ) * 0.5;
            bpx[3] = sc1 * ( bpx1[3] - mpz1*bpx1[0] ) * 0.75;
            /* These should be the same differences like in pzpsa code */
            for ( isamp=4; isamp<nsamp4; isamp++ )
               bpx[isamp] = sc1 * ( bpx1[isamp] - mpz1*bpx1[isamp-4] );
#if 0
            /* Code for debugging and testing only ... */
            if ( teldata->pixtm->timval[ipix][0] >= 0 &&
                 teldata->pixtm->pulse_sum_glob[0][ipix] >= 200 )
            { 
                double bpx2[4*H_MAX_SLICES+8];
                PzpsaSmoothUpsampleU16(raw->num_samples,4,smp,ped,mpz1,bpx2,NULL,NULL);
                for ( isamp=0; isamp<nsamp4; isamp++ )
                  printf("==> %4d %3d  %f  %f  %f\n", ipix, isamp, bpx1[isamp], bpx[isamp], bpx2[isamp]);
               exit(9);
            }
#endif
         }
         else if ( psopt%10 == 2 ) /* Pulse shaping with pzpsa code but instead of differences
                                      of subsequent slices (before upscaling, e.g. 4 ns apart) 
                                      it takes differences of second-next slices (e.g. 8 ns apart),
                                      with adapted pole-zero compensation. */
         {
            PzpsaSmoothUpsampleU16(raw->num_samples,4,smp,ped,0.,bpx1,NULL,NULL);
            /* Smooth start, no sudden peak at start of window */
            for ( isamp=0; isamp<4 && isamp<nsamp4; isamp++ )
               bpx[isamp] = sc2 * ( bpx1[isamp+4] - mpz2*bpx1[0] ) * isamp/4.;
            /* Note that the differences are for symmetric offsets */
            for ( isamp=4; isamp+4<nsamp4; isamp++ )
               bpx[isamp] = sc2 * ( bpx1[isamp+4] - mpz2*bpx1[isamp-4] );
            /* Smooth end to avoid random peaks at the end of the window */
            for ( ; isamp<nsamp4; isamp++ )
               bpx[isamp] = sc2 * ( bpx[nsamp4-1] - mpz2*bpx1[isamp-4] ) * (nsamp4-1-isamp)/4.;
#if 0
            /* Code for debugging and testing only ... */
            if ( teldata->pixtm->timval[ipix][0] >= 0 &&
                 teldata->pixtm->pulse_sum_glob[0][ipix] >= 200 )
            {
               for ( isamp=0; isamp<nsamp4; isamp++ )
                  printf("==> %4d %3d  %f\n", ipix, isamp, bpx[isamp]);
               exit(9);
            }
#endif
         }
         else if ( psopt == -1 ) /* Special case of pulse shaping completely with pzpsa code
                                    plus peak detection without involving neighbours.
                                    As a result there is a significant bias by positive
                                    NSB fluctuations, depending on the level of NSB. */
         {
            double bmax = 0., sum = 0.;
            int ipmx = 0;
            PzpsaSmoothUpsampleU16(raw->num_samples,4,smp,ped,mpz1,bpx,&bmax,&ipmx);
            sum = bmax;
            if ( integration_correction[itel][igain] > 0. )
               sum *= integration_correction[itel][igain];
            sum += moni->pedestal[igain][ipix];
            raw->adc_sum[igain][ipix] = (sum>0.) ? (int) (sum+0.5) : 0;
#if 0
            /* Code for debugging and testing only ... */
            if ( teldata->pixtm->timval[ipix][0] >= 0 &&
                 teldata->pixtm->pulse_sum_glob[0][ipix] >= 500 )
            {
               LasCalData *lcal = &hsdata->tel_lascal[itel];
               printf("xxx> ipix=%d, bmax=%f, ipmx=%d, sum=%f, p.e.=%f (glob=%f, loc=%f)\n", ipix, bmax, ipmx, sum,
                  (sum-moni->pedestal[igain][ipix])*lcal->calib[igain][ipix],
                  teldata->pixtm->pulse_sum_glob[0][ipix]*lcal->calib[0][ipix],
                  teldata->pixtm->pulse_sum_loc[0][ipix]*lcal->calib[0][ipix]);
             //   for ( isamp=0; isamp<nsamp4; isamp++ )
             //      printf("==> %4d %3d  %f\n", ipix, isamp, bpx[isamp]);
             //   exit(9);
            }
#endif
         }
         else
#endif
         {  /* Not using pzpsa code - slower but except at start and end identical to pzpsa. */
            /* This code is here for verification and better understanding what the
               pzpsa code does. Except for such tests there is no point in using it. */

            /* Upsampling with linear interpolation, subtract pedestal */
            double a0 = smp[0] - ped, a1=a0, a2=a0;
            for ( isamp=0; isamp<raw->num_samples; isamp++ )
            {
               a0 = a1;
               a1 = a2;
               a2 = smp[isamp] - ped;
               bpx0[4*isamp+0] = 0.5*(a0+a1);
               bpx0[4*isamp+1] = 0.25*a0 + 0.75*a1;
               bpx0[4*isamp+2] = a1;
               bpx0[4*isamp+3] = 0.75*a1 + 0.25*a2;
            }
            /* Running average of four upscaled samples */
            bpx1[0] = 0.25*(bpx0[0]+bpx0[0]+bpx0[1]+bpx0[2]);
            for ( isamp=1; isamp+3<nsamp4; isamp++ )
               bpx1[isamp] = 0.25*(bpx0[isamp-1]+bpx0[isamp]+bpx0[isamp+1]+bpx0[isamp+2]);
            bpx1[nsamp4-3] = 0.25*(bpx0[nsamp4-4]+bpx0[nsamp4-3]+bpx0[nsamp4-2]+bpx0[nsamp4-1]);
            bpx1[nsamp4-2] = 0.25*(bpx0[nsamp4-3]+bpx0[nsamp4-2]+2*bpx0[nsamp4-1]);
            bpx1[nsamp4-1] = 0.25*(bpx0[nsamp4-2]+3*bpx0[nsamp4-1]);

            /* Note that after interpolation and running average the
               center-of-gravity of any pulses (minus pedestal) is unchanged
               and the resulting signal of a delta-shaped raw pulse is symmetric. */

            /* Use one of two differencing options */
            bpx = bg + off_pix*ipix;
            if ( psopt%10 == 3 ) /* Difference of values one original sample apart (like 'pzpsa.c' code) */
            {
               /* Same smooth start as for psopt == 1 */
               for ( isamp=0; isamp<4 && isamp<nsamp4; isamp++ )
                  bpx[isamp] = sc1 * (bpx1[isamp]-mpz1*bpx1[0]) * isamp/4.;
               for ( isamp=4; isamp<nsamp4; isamp++ )
               {
                  bpx[isamp] = sc1*(bpx1[isamp]-mpz1*bpx1[isamp-4]);
               }
            }
            else if ( psopt%10 == 4 ) /* Difference of values two original samples apart (retain more signal but slightly longer shape) */
            {
               /* Same smooth start as for psopt == 2 */
               for ( isamp=0; isamp<4 && isamp<nsamp4; isamp++ )
                  bpx[isamp] = sc2 * (bpx1[isamp+4]-mpz2*bpx1[0]) * isamp/4.;
               /* Same differences as for psopt == 2 */
               for ( isamp=4; isamp+4<nsamp4; isamp++ )
               {
                  bpx[isamp] = sc2*(bpx1[isamp+4]-mpz2*bpx1[isamp-4]);
               }
               /* and same smooth lead-out. */
               for ( ; isamp<nsamp4; isamp++ )
                  bpx[isamp] = sc2*(bpx1[nsamp4-1]-mpz2*bpx1[isamp-4]) * (nsamp4-1-isamp)/4.;
            }
            else
            {
               fflush(stdout);
               fprintf(stderr,"Invalid pulse shaping option %d\n", psopt);
               return -1;
            }

#if 0
            /* Code for debugging and testing only ... */
            if ( teldata->pixtm->timval[ipix][0] >= 0 &&
                 teldata->pixtm->pulse_sum_glob[0][ipix] >= 200 )
            {
               bpx = bg + off_pix*ipix;
               for ( isamp=0; isamp<nsamp4; isamp++ )
                  printf("==> %4d %3d  %f\n", ipix, isamp, bpx[isamp]);
               exit(9);
            }
#endif
         }      

      } /* Not pzpsa */  
   }
   
   if ( psopt < 0 ) /* Not actually using neighbour pixels? */
      return 0;


   /* We need the list of neighbours now */
   if ( nb_lists[itel][0].nblist == NULL )
   {
      find_neighbours(&hsdata->camera_set[itel],itel);
   }
   
   if ( nb_lists[itel][0].nblist == NULL || 
        nb_lists[itel][0].nbsize <= 0 )
      return -1;
   nbl = &nb_lists[itel][0];

   /* For simplicity we do this separately and in the same way for
      each gain channel although FlashCam/DigiCam at least have only
      one channel. Inside this loop we have the same peak finding as
      in nb_peak_integration with zero local pixel weight (lwt=0),
      although now on the upscaled and shaped signal. */
   for (igain=0; igain<raw->num_gains; igain++)
   {
      double *bg = buffer + igain*off_gain, *bpx, *bpx_nb;
      for (ipix=0; ipix<raw->num_pixels; ipix++)
      {
         double nb_samples[H_MAX_SLICES], p = 0.;
         int inb, knb=0, ifirst, ilast;
         
         /* Keep in mind that adc_sum has already been initialized in first loop */

         /* Require available sample-mode data in this pixel, not zero-suppressed. */
         if ( ! raw->significant[ipix] || ! raw->adc_known[igain][ipix] ||
              ( (raw->zero_sup_mode & 0x20) != 0 && (raw->significant[ipix] & 0x020) == 0 ) )
         {
            if ( pixtim_flag )
            {
               teldata->pixtm->pulse_sum_glob[igain][ipix] = 
               teldata->pixtm->pulse_sum_loc[igain][ipix] = 0.;
               if ( igain == 0 )
                  teldata->pixtm->timval[ipix][0] = -1; /* Others not rewritten */
            }
            continue;
         }

         bpx = bg + off_pix*ipix;
         for (isamp=0; isamp<4*raw->num_samples; isamp++ )
            nb_samples[isamp] = 0.;

         /* Add up signals of suitable neighbour pixels */
         for ( inb=0; inb<nbl->pix_num_nb[ipix]; inb++ )
         {
            int ipix_nb = nbl->nblist[nbl->pix_first_nb[ipix]+inb];
            /* For zero-suppressed sample mode data ALSO check relevant bit of the neighbour */
            if ( (raw->zero_sup_mode & 0x20) != 0 && (raw->significant[ipix_nb] & 0x020) == 0 )
               continue;
            if ( raw->significant[ipix_nb] && raw->adc_known[igain][ipix_nb] )
            {
               bpx_nb = bg + off_pix*ipix_nb;
               for (isamp=0; isamp<4*raw->num_samples; isamp++ )
                  nb_samples[isamp] += bpx_nb[isamp];
               knb++;
            }
         }

         if ( knb == 0 ) /* No integration window available for truely isolated pixels */
         {
            if ( pixtim_flag )
            {
               teldata->pixtm->pulse_sum_glob[igain][ipix] = 
               teldata->pixtm->pulse_sum_loc[igain][ipix] = 0.;
               if ( igain == 0 )
                  teldata->pixtm->timval[ipix][0] = -1; /* Others not rewritten */
            }
            continue; /* Perhaps we should also reset the 0x01 raw->significant bit? */
         }

         if ( psopt%10 > 0 )
         {
            ifirst = 0;
            ilast = nsamp4;
         }
         else /* Avoid artifacts/noise at beginning and end of original pzpsa shaped pulses */
         {
            ifirst = 7;
            ilast = nsamp4 - 3;
         }
         ipeak = ifirst;
         p = nb_samples[ifirst];
         for ( isamp=ifirst+1; isamp<ilast; isamp++ )
         {
            if ( nb_samples[isamp] > p )
            {
               p = nb_samples[isamp];
               ipeak = isamp;
            }
         }
         peakpos = ipeak;

         if ( raw->significant[ipix] && raw->adc_known[igain][ipix] )
         {
            double sum = 0., cog = 0.;
            int icog = (peakpos < 1 ) ? 1 : (peakpos+1 >= nsamp4 ? nsamp4-2 : peakpos);
            if ( psopt%100 < 10 ) /* sum of shaped signal */
            {
               if ( psopt/100 != 1 )
               {
                  /* Summation range completely independent of signal in pixel, no bias at all. */
                  start = peakpos - nbefore;
                  if ( start < 0 )
                     start = 0;
                  if ( start + nsum > nsamp4 )
                     start = nsamp4 - nsum;
                  for ( isamp=0; isamp<nsum; isamp++ )
                     sum += bpx[isamp+start];
               }

               /* The peak position in this pixel should be close to neighbours but may differ a bit */
               if ( icog > ifirst+1 && bpx[icog-1] > bpx[icog] && bpx[icog-1] > bpx[icog+1] )
               {
                  icog--;
                  if ( icog > ifirst+1 && bpx[icog-1] > bpx[icog] && bpx[icog-1] > bpx[icog+1] )
                  {
                     icog--;
                     if ( icog > ifirst+1 && bpx[icog-1] > bpx[icog] && bpx[icog-1] > bpx[icog+1] )
                        icog--;
                  }
               }
               else if ( icog < ilast-2 && bpx[icog+1] > bpx[icog] && bpx[icog+1] > bpx[icog-1] )
               {
                  icog++;
                  if ( icog < ilast-2 && bpx[icog+1] > bpx[icog] && bpx[icog+1] > bpx[icog-1] )
                  {
                     icog++;
                     if ( icog < ilast-2 && bpx[icog+1] > bpx[icog] && bpx[icog+1] > bpx[icog-1] )
                        icog++;
                  }
               }
               ipeak = icog;
               /* Despite the name, this is not (unlike in pzpsa code) a center-of-gravity but
                  the position of the peak of an assumed parabola through the icog slice and 
                  the preceding and following slices. */
               cog = icog + (bpx[icog+1]-bpx[icog-1])/(4.*bpx[icog]-2.*(bpx[icog+1]+bpx[icog-1]));
               /* If we didn't start from a peak or negative signals were involved we
                  may actually get a result too far away from where we started. */
               if ( fabs(cog-icog) > 1.0 )
                  cog = icog;

               if ( psopt/100 == 1 )
               {
                  /* Summation range defined after small adjustments for peak position, small NSB bias */
                  start = icog - nbefore;
                  if ( start < 0 )
                     start = 0;
                  if ( start + nsum > nsamp4 )
                     start = nsamp4 - nsum;
                  for ( isamp=0; isamp<nsum; isamp++ )
                     sum += bpx[isamp+start];
               }
               
            }
            else /* Use the pzpsa summation and center-of-gravity scheme */
            {
               int w = (psopt%100)/10; /**< Extension of summation/cog region [peakpos-w : peakpos+w] */
               sum = 0.;
               (void) PzpsaPeakProperty (nsamp4, bpx, peakpos, w, &sum, &cog);
               sum /= (1.+2*w);
               icog = (int) (cog+0.5);
            }

            /* In any case we have to apply a correction factor, depending on the
               pulse-shaping option in use, such that later calibration returns
               a meaningful result -  although there are still order of 10% 
               systematic effects depending on the option used. */
            if ( integration_correction[itel][igain] > 0. )
               sum *= integration_correction[itel][igain];
            
            if ( pixtim_flag ) /* Asked to re-write the pixel timing structure */
            {
               if ( ipeak < 0 || ipeak >= nsamp4 ) /* No proper peak known */
               {
                  teldata->pixtm->pulse_sum_glob[igain][ipix] = 
                  teldata->pixtm->pulse_sum_loc[igain][ipix] = 0.;
                  if ( igain == 0 )
                     teldata->pixtm->timval[ipix][0] = -1; /* Others not rewritten */
               }
               else if ( bpx[ipeak] < pixtim_thr ) /* Pixel timing not considered significant */
               {
                  teldata->pixtm->pulse_sum_glob[igain][ipix] = 
                  teldata->pixtm->pulse_sum_loc[igain][ipix] = 0.;
                  if ( igain == 0 )
                     teldata->pixtm->timval[ipix][0] = -1; /* Others not rewritten */
               }
               else
               {
                  teldata->pixtm->pulse_sum_glob[igain][ipix] = 
                  teldata->pixtm->pulse_sum_loc[igain][ipix] = (sum>0.) ? (int) (sum+0.5) : 0;
                  if ( igain == 0 )
                  {
                     double ts20, ts50, ts80, te20, te50, tsthr, tethr;
                     double apeak = bpx[ipeak];
                     int i;
                     ts20 = ts50 = ts80 = tsthr = 0; 
                     te20 = te50 = tethr = nsamp4-1;
                     for ( i=ipeak-1; i>=ifirst; i-- )
                     {
                        if ( bpx[i] <= 0.2 * apeak )
                        {
                           ts20 = i + (0.2*apeak - bpx[i])/(bpx[i+1]-bpx[i]);
                           break;
                        }
                     }
                     for ( i=ipeak-1; i>=ifirst; i-- )
                     {
                        if ( bpx[i] <= 0.5 * apeak )
                        {
                           ts50 = i + (0.5*apeak - bpx[i])/(bpx[i+1]-bpx[i]);
                           break;
                        }
                     }
                     for ( i=ipeak-1; i>=ifirst; i-- )
                     {
                        if ( bpx[i] <= 0.8 * apeak )
                        {
                           ts80 = i + (0.8*apeak - bpx[i])/(bpx[i+1]-bpx[i]);
                           break;
                        }
                     }
                     for ( i=ipeak-1; i>=ifirst; i-- )
                     {
                        if ( bpx[i] <= pixtim_thr )
                        {
                           tsthr = i + (pixtim_thr - bpx[i])/(bpx[i+1]-bpx[i]);
                           break;
                        }
                     }
                     for ( i=ipeak+1; i<ilast; i++ )
                     {
                        if ( bpx[i] <= pixtim_thr )
                        {
                           tethr = i - (pixtim_thr - bpx[i])/(bpx[i-1]-bpx[i]);
                           break;
                        }
                     }
                     for ( i=ipeak+1; i<ilast; i++ )
                     {
                        if ( bpx[i] <= 0.5 * apeak )
                        {
                           te50 = i - (0.5*apeak - bpx[i])/(bpx[i-1]-bpx[i]);
                           break;
                        }
                     }
                     for ( i=ipeak+1; i<ilast; i++ )
                     {
                        if ( bpx[i] <= 0.2 * apeak )
                        {
                           te20 = i - (0.2*apeak - bpx[i])/(bpx[i-1]-bpx[i]);
                           break;
                        }
                     }
                     teldata->pixtm->timval[ipix][0] = 0.25*cog; /* In units of original sampling time step */
                     teldata->pixtm->timval[ipix][1] = 0.25*ts20;
                     teldata->pixtm->timval[ipix][2] = 0.25*ts50;
                     teldata->pixtm->timval[ipix][3] = 0.25*ts80;
                     teldata->pixtm->timval[ipix][4] = 0.25*(te50-ts50);
                     teldata->pixtm->timval[ipix][5] = 0.25*(te20-ts20);
                     teldata->pixtm->timval[ipix][6] = 0.25*(tethr-tsthr);
                  }
               }
            }

            /* We add the original pedestal, again to keep calibration as-is. */
            sum += moni->pedestal[igain][ipix];
            /* Round the result to the nearest unsigned integer */
            raw->adc_sum[igain][ipix] = (sum>0.) ? (uint32_t) (sum+0.5) : 0.;
            /* Alternative: (uint32_t)nearbyint(sum) */
#if 0
            /* Code for debugging and testing only ... */
            if ( teldata->pixtm->timval[ipix][0] >= 0 &&
                 teldata->pixtm->pulse_sum_glob[0][ipix] >= 500 )
            {
               LasCalData *lcal = &hsdata->tel_lascal[itel];
               printf("xxx> ipix=%d bpx[%d]=%f, cog=%f, sum=%f, p.e.=%f (glob=%f, loc=%f)\n", ipix, ipeak, bpx[ipeak], cog, sum,
                  (sum-moni->pedestal[igain][ipix])*lcal->calib[igain][ipix],
                  teldata->pixtm->pulse_sum_glob[0][ipix]*lcal->calib[0][ipix],
                  teldata->pixtm->pulse_sum_loc[0][ipix]*lcal->calib[0][ipix]);
            }
#endif
         }
      }
   }

   if ( pixtim_flag )
      teldata->pixtm->known = 1;

   return 0;
}

/* -------------------------- qpol ------------------------------------- */
/** Quick interpolation in array of points equidistant in x coordinate */

static double qpol(double x, int np, double *yval);

static double qpol(double x, int np, double *yval)
{
   int ix = (int) x;
   if ( x < 0 || x >= (double) np )
      return 0.;
   if ( ix+1 >= np )
      return 0.;
   return yval[ix]*(ix+1-x) + yval[ix+1]*(x-ix);
}

/* ------------------ set_integration_correction ---------------------- */

/** With partial pulse integration we extract a correction factor
    from partial to full pulse area from the reference pulse shape
    provided by MC. Since actual pulses may have an intrinsic width 
    (and as a result are wider than the reference pulse) this can 
    still lead to a bit underestimated p.e. values. But this is
    hard to fix without knowing the true width of light pulses. */

static int set_integration_correction(AllHessData *hsdata, int itel, int integrator, int *intpar);

static int set_integration_correction(AllHessData *hsdata, int itel, int integrator, int *intpar)
{
   int igain, ibin, iphase;
   int nbins = intpar[0];
   int noff  = intpar[1];
   int psopt = intpar[2];
   int nsamp;
   TelEvent *teldata = NULL;
   AdcData *raw = NULL;
   if ( hsdata == NULL || itel < 0 || itel >= H_MAX_TEL )
      return -1;
   teldata = &hsdata->event.teldata[itel];
   raw = teldata->raw;
   if ( raw == NULL ) /* No point in integration correction if there is nothing to integrate */
      return -1;
   nsamp = raw->num_samples;
   
   PixelSetting *ps = &hsdata->pixel_set[itel];
   CameraOrganisation *co = &hsdata->camera_org[itel];

// printf("###> # nsamp=%d, nbins=%d, noff=%d, psopt=%d, time_slice=%f, ref_step=%f\n", 
//   nsamp, nbins, noff, psopt,ps->time_slice, ps->ref_step);

   for (igain=0; igain<co->num_gains; igain++)
   {
      double sum = 0., asum = 0., speak = 0.;
      int ipeak = 0;
      double st = ps->time_slice / ps->ref_step;
      double sr = 1./st;

      integration_correction[itel][igain] = 1.0; /* Fall-back to avoid repeated attempts. */
      if ( ps->nrefshape <= igain || ps->time_slice == 0. || ps->ref_step == 0. )
         continue;

      /* Sum over all the pulse we have */
      for ( ibin=0; ibin<ps->lrefshape; ibin++ )
      {
         asum += ps->refshape[igain][ibin];
         if ( ps->refshape[igain][ibin] > speak )
         {
            speak = ps->refshape[igain][ibin];
            ipeak = ibin;
         }
      }
      /* Rescale to original time step */
      asum *= sr;

      if ( integrator == 7 ) /* Window is in upsampled and shaped pulse */
      {
         uint16_t adc[2*H_MAX_SLICES];
//         double traw[2*H_MAX_SLICES];
         double shaped1[8*H_MAX_SLICES], shaped[8*H_MAX_SLICES];
//         double shaped0[8*H_MAX_SLICES];
         int nsamp4 = 4 * nsamp, istart;
         double mpz1 = 0.758;
         double mpz2 = mpz1*mpz1;
         double sc1 = 1.0;  /* like pzpsa not scaled; equal area would be: 0.25/(1.-mpz1) */
         double sc2 = 1.0;  /* equal area: 0.25/(1.-mpz2) */
         double speakraw = speak;
         double asum1 = 0., speak1 = 0.;
         int ipeak1 = 0;
         if ( nbins <= 0 )
            nbins = 1;
         if ( nbins > nsamp4 )
            nbins = nsamp4;
         if ( noff > nbins )
            noff = nbins;

// printf("###> st=%f, sr=%f, ipeak=%d, mpz1=%f, sc1=%f, mpz2=%f, sc2=%f\n", st, sr, ipeak, mpz1, sc1, mpz2, sc2);
         for ( iphase=0; iphase<10; iphase++ )
         {
            double bmx = 0.;
            double ti = ((iphase*0.1-0.5)-noff) * st + ipeak;
            int ipmx = 0;
            for ( ibin=0; ibin<2*nsamp; ibin++ )
            {
//               traw[ibin] = ((ibin-nsamp/2)*st + ti) * ps->ref_step;
               adc[ibin] = (uint16_t) (qpol((ibin-nsamp/2)*st + ti, ps->lrefshape, ps->refshape[igain])/speakraw * 3000. + 1000);
            }

#ifdef WITH_PZPSA

// PzpsaSmoothUpsampleU16(2*nsamp,4,adc,1000.,0.,shaped0,NULL,NULL);
// PzpsaSmoothUpsampleU16(2*nsamp,4,adc,1000.,mpz1,shaped1,NULL,NULL);

            switch ( psopt%10 )
            {
               case -1: /* Full pzpsa shaping with full range peak search included */
                  PzpsaSmoothUpsampleU16(2*nsamp,4,adc,1000.,mpz1,shaped,&bmx,&ipmx);
// printf("xxx> bmx=%f, shaped(%d)=%f\n",bmx,ipmx,shaped[ipmx]);
                  sum += bmx;
                  speak1 = bmx;
                  ipeak1 = ipmx;
                  break;
               case 0: /* Full pzpsa shaping+diff with nb peak search excluding edges */
               case 9: /* Full pzpsa shaping+diff with nb peak search including edges */
                  PzpsaSmoothUpsampleU16(2*nsamp,4,adc,1000.,mpz1,shaped,NULL,NULL);
                  break;
               case 1: /* Pzpsa shaping with own differencing */
               case 3: /* For simplicity also using pzpsa code for psopt=3 */
                  PzpsaSmoothUpsampleU16(2*nsamp,4,adc,1000.,0.,shaped1,NULL,NULL);
                  shaped[0] = sc1 * shaped1[2];
                  shaped[1] = sc1 * shaped1[3];
                  for ( ibin=2; ibin+2<2*nsamp4; ibin++ )
                     shaped[ibin] = sc1 * ( shaped1[ibin+2] - mpz1*shaped1[ibin-2] );
                  for ( ; ibin<2*nsamp4; ibin++ )
                     shaped[ibin] = -sc1 * mpz1*shaped1[ibin-2];
                  break;
               case 2: /* Pzpsa shaping with own 2nd-next differencing */
               case 4: /* For simplicity also using pzpsa code for psopt=4 */
                  PzpsaSmoothUpsampleU16(2*nsamp,4,adc,1000.,0.,shaped1,NULL,NULL);
                  for ( ibin=0; ibin<4 && ibin<2*nsamp4; ibin++ )
                     shaped[ibin] = sc2 * shaped1[ibin+4];
                  for ( ibin=4; ibin+4<2*nsamp4; ibin++ )
                     shaped[ibin] = sc2 * ( shaped1[ibin+4] - mpz2*shaped1[ibin-4] );
                  for ( ; ibin<2*nsamp4; ibin++ )
                     shaped[ibin] = -sc2 * mpz2*shaped1[ibin-4];
                  break;
               default:
                  return -1;
            }
            if ( psopt >= 0 )
            {
               speak1 = 0.;
               ipeak1 = noff;
               for ( ibin=0; ibin<2*nsamp4; ibin++ )
               {
// if ( iphase == 0 )
// printf("###> %f:  %f  %f  %f  %f\n", traw[ibin/4]+ibin%4, (adc[ibin/4]-1000.), shaped0[ibin], shaped1[ibin], shaped[ibin]);
                  asum1 += shaped[ibin];
                  if ( shaped[ibin] > speak1 )
                  {
                     speak1 = shaped[ibin];
                     ipeak1 = ibin;
                  }
               }
               istart = ipeak1 - nbins;
               if ( istart < 0 )
                  istart = 0;
               if ( (psopt%100)/10 > 0 )
               {
                  double peaksum=0., cog=0.;
                  int w = (psopt%100)/10;
                  (void) PzpsaPeakProperty (2*nsamp4, shaped, ipeak1, w, &peaksum, &cog);
                  sum += peaksum/(1.+2*w);
               }
               else
               {
                  for ( ibin=istart; ibin<istart+nbins; ibin++ )
                     sum += shaped[ibin];
               }
            }

#else
            fflush(NULL);
            fprintf(stderr,"... Setting up integration correction factor for integration scheme 7 "
               "without pzpsa code (for psopt 3 and 4) not implemented ...\n");
            exit(1);
#endif

         }
         sum *= speakraw / 3000. * 0.1; /* Tried 10 phases, had amplitude scaled by factor 3000/speakraw */

printf("Integration correction factor: asum=%f, asum1=%f, sum=%f, speakraw=%f at %d, speak=%f at %d of %d\n",
   asum, asum1/3000., sum*0.1, speakraw, ipeak, speak1, ipeak1, 2*nsamp4);

      }

      else /* Window is in raw pulse shape (and original sampling) */
      {
         /* Now sum up given interval starting from peak, averaging over phase */
         for ( iphase=0; iphase<10; iphase++ )
         {
            double ti = ((iphase*0.1-0.5)-noff) * st + ipeak;
            for ( ibin=0; ibin<nbins; ibin++ )
               sum += qpol(ibin*st + ti, ps->lrefshape, ps->refshape[igain]);
         }
         sum *= 0.1; /* Tried 10 phases */
      }

      if ( sum > 0. && asum > 0. )
         integration_correction[itel][igain] = asum/sum;

printf("Integration correction factor for telescope #%d (ID=%d) gain %d (scheme %d, with %d samples, offset %d, option %d) is %f\n", 
   itel, raw->tel_id, igain, integrator, nbins, noff, psopt, integration_correction[itel][igain]);

   }
   return 0;
}

/* ------------------------------ pixel_integration -------------------------- */
/**
 *  @short Pixel integration steering function. Work is done in selected integration function.
 */

static int pixel_integration(AllHessData *hsdata, int itel, struct user_parameters *up);

static int pixel_integration(AllHessData *hsdata, int itel, struct user_parameters *up)
{
#ifdef DEBUG_PIXEL_NB
printf("Pixel integration for telescope #%d\n",itel);
#endif

   if ( hsdata == NULL || up == NULL )
      return -1;

   if ( integration_correction[itel][0] == 0. )
      set_integration_correction(hsdata, itel, up->i.integrator, up->i.integ_param);

   switch ( up->i.integrator )
   {
      case 1: /* Fixed integration region */
         return simple_integration(hsdata, itel, up->i.integ_param[0], up->i.integ_param[1]);
         break;
      case 2: /* Integration region by global peak of significant pixels */
         return global_peak_integration(hsdata, itel, up->i.integ_param[0], up->i.integ_param[1], up->i.integ_thresh);
         break;
      case 3: /* Peak in each pixel determined independently */
         return local_peak_integration(hsdata, itel, up->i.integ_param[0], up->i.integ_param[1], up->i.integ_thresh);
         break;
      case 4: /* Integration region determined by signal in neighbours only. */
         return nb_peak_integration(hsdata, 0, itel, up->i.integ_param[0], up->i.integ_param[1], up->i.integ_thresh);
         break;
      case 5: /* Mix of neighbours and local pixel */
         return nb_peak_integration(hsdata, 3, itel, up->i.integ_param[0], up->i.integ_param[1], up->i.integ_thresh);
         break;
      case 6: /* Time gradient fit used to define placement of integration window */
         return gradient_integration(hsdata, itel, up->i.integ_param[0], up->i.integ_param[1], up->i.integ_thresh);
         break;
      case 7: /* Integration of FlashCam-shaped signal in region defined by neighbours only. */
         return nb_fc_shaped_peak_integration(hsdata, itel, up->i.integ_param[0], up->i.integ_param[1], up->i.integ_thresh,
            up->i.integ_param[2], up->i.integ_thresh[0]);
         break;
      default:
         fprintf(stderr,"Invalid integration method %d.\n", up->i.integrator);
         return -1;
   }
   return 0;
}

/* --------------------------- clean_image_tailcut ------------------------ */
/** 
 *  @short Use dual-level tail-cut image cleaning procedure to get pixel list. 
 *
 *  In contrast to the classical dual-level tail-cuts this function
 *  has an optional restriction to only those pixels having an
 *  amplitude above a given fraction of the n-th hottest pixel.
 *  This should almost stop the increase of width and length
 *  with increasing intensity after some point.
 *
 *  @param hsdata Pointer to all available data and configurations.
 *  @param itel   Sequence number of the telescope being processed.
 *  @param al     The lower of the two tail-cut thresholds.
 *  @param ah     The higher of the two tail-cut thresholds.
 *  @param lref   Determines which pixel, after sorting by amplitude,
 *                will be used as providing the reference amplitude.
 *                Example: use 3 for the third hottest pixel.
 *                If this number is <= 0, the classical scheme is used.
 *  @param minfrac Which fraction of the reference amplitude is required
 *                for pixels to be included in the final image.
 *                If this number is <= 0.0, the classical scheme is used.
 */

static int clean_image_tailcut(AllHessData *hsdata, int itel, 
   double al, double ah, int lref, double minfrac);

static int clean_image_tailcut(AllHessData *hsdata, int itel, 
   double al, double ah, int lref, double minfrac)
{
   int pass_low[H_MAX_PIX], pass_high[H_MAX_PIX];
   int npix;
   int i;
   TelEvent *teldata = NULL;
   struct camera_nb_list *nbl = NULL;

   if ( hsdata == NULL || itel < 0 || itel >= H_MAX_TEL )
      return -1;
   npix = hsdata->camera_set[itel].num_pixels;
   
   teldata = &hsdata->event.teldata[itel];
   teldata->image_pixels.pixels = image_numpix[itel] = 0;
   if ( !teldata->known || teldata->raw == NULL )
      return -1;
   if ( !teldata->raw->known )
      return -1;

   if ( nb_lists[itel][0].nblist == NULL )
   {
      find_neighbours(&hsdata->camera_set[itel],itel);
   }
   if ( nb_lists[itel][0].nbsize <= 0 )
      return -1;
   nbl = &nb_lists[itel][0];

   for (i=0; i<npix; i++)
   {
      if ( pixel_amp[itel][i] < al )
         pass_low[i] = pass_high[i] = 0;
      else if ( pixel_amp[itel][i] < ah )
      {
         pass_low[i] = 1;
         pass_high[i] = 0;
      }
      else
         pass_low[i] = pass_high[i] = 1;
   }

   for (i=0; i<npix; i++)
   {
      if ( pass_high[i] )
      {
         int j;
         for (j=0; j<nbl->pix_num_nb[i]; j++)
            if ( pass_low[nbl->nblist[nbl->pix_first_nb[i]+j]] )
	    {
               image_list[itel][image_numpix[itel]] = i;
               teldata->image_pixels.pixel_list[image_numpix[itel]] = i;
               image_numpix[itel]++;
	       break;
	    }
      }
      else if ( pass_low[i] )
      {
         int j;
         for (j=0; j<nbl->pix_num_nb[i]; j++)
            if ( pass_high[nbl->nblist[nbl->pix_first_nb[i]+j]] )
	    {
               image_list[itel][image_numpix[itel]] = i;
               teldata->image_pixels.pixel_list[image_numpix[itel]] = i;
               image_numpix[itel]++;
	       break;
	    }
      }
   }

   /* If a minimum fraction of the amplitude of the n-th hottest pixel */
   /* is required, we sort the pixels by amplitude first. */
   if ( lref > 0 && lref < image_numpix[itel] && minfrac > 0. )
   {
      int j;
      double refamp;
      for (i=0; i<image_numpix[itel]; i++)
      {
         for (j=i+1; j<image_numpix[itel]; j++)
         {
            int ipix = teldata->image_pixels.pixel_list[i];
            int jpix = teldata->image_pixels.pixel_list[j];
            if ( pixel_amp[itel][jpix] > pixel_amp[itel][ipix] )
            {
               teldata->image_pixels.pixel_list[i] = jpix;
               teldata->image_pixels.pixel_list[j] = ipix;
            }
         }
      }
      refamp = pixel_amp[itel][teldata->image_pixels.pixel_list[lref-1]];
      for ( i=lref; i<image_numpix[itel]; i++ )
      {
         int ipix = teldata->image_pixels.pixel_list[i];
         if ( pixel_amp[itel][ipix] < minfrac*refamp )
         {
#if 0
int k;
printf("\nSorted pixel list with");
for (k=0; k<image_numpix[itel]; k++)
{
   int kpix = teldata->image_pixels.pixel_list[k];
   printf(" %f",pixel_amp[itel][kpix]);
}
printf(" gets truncated after %d hottest pixels.\n", i);
#endif
            image_numpix[itel] = i;
            break;
         }
      }
   }

   teldata->image_pixels.pixels = image_numpix[itel];

   return 0;
}

/* ---------------------------- second_moments ---------------------------- */

/** Reconstruction of second moments parameters from cleaned image. */

/* Note: Can only be used after clean_image_tailcut()
   because it is using external static storage filled there. */

static int second_moments(AllHessData *hsdata, int itel, int cut_id, int nimg, double clip_amp);

static int second_moments(AllHessData *hsdata, int itel, int cut_id, int nimg, double clip_amp)
{
   static int iprint = 0;
   CameraSettings *camset = &hsdata->camera_set[itel];
   int i, j;
   double dx = 0., dy = 0.; // Source offsets in camera
   double sx = 0.,sy = 0., sxx = 0., sxy = 0., syy = 0., sA = 0.;
   double img_scale = 180./M_PI/hsdata->camera_set[itel].flen;
   double a, b, sx3, sx4, beta, cb, sb;
   double alpha, distance, miss, width, length;
   double xmean, ymean, orientation, direction;
   double skewness, kurtosis;
   int hot_pixel[5];
   double hot_amp[5];
   TelEvent *teldata = NULL;
   double stot=0.;
   int npix = hsdata->camera_set[itel].num_pixels;
   ImgData *img = NULL;

   if ( itel < 0 || itel >= H_MAX_TEL )
      return -1;

   teldata = &hsdata->event.teldata[itel];
   if ( teldata->raw == NULL )
      return -1;
   if ( !teldata->raw->known )
      return -1;
   if ( teldata->img == NULL )
      return -1;
   if ( nimg < 0 || nimg >= teldata->max_image_sets )
      return -1;
   if ( nimg >= teldata->num_image_sets )
      teldata->num_image_sets = nimg+1;
   img = &teldata->img[nimg];
   img->known = 0;
   img->amplitude = 0.;
   img->pixels = 0;
   img->num_sat = pixel_sat[itel];
   img->clip_amp = clip_amp;

   if ( show_total_amp )
      for (i=0; i<npix; i++)
         stot += pixel_amp[itel][i];

   if ( image_numpix[itel] < 2 ) // Minimum 2 pixels
      return -1;

   for (j=0; j<5; j++)
   {
      hot_pixel[j] = -1;
      hot_amp[j] = -1;
   }

   for (j=0; j<image_numpix[itel]; j++)
   {
      i = image_list[itel][j];
      sA += pixel_amp[itel][i];
      if ( pixel_amp[itel][i] > hot_amp[4] )
      {
         int k, l;
         for (k=0; k<5; k++ )
            if ( pixel_amp[itel][i] > hot_amp[k] )
            {
               for ( l=4; l>k; l-- )
               {
                  hot_amp[l] = hot_amp[l-1];
                  hot_pixel[l] = hot_pixel[l-1];
               }
               hot_amp[k] = pixel_amp[itel][i];
               hot_pixel[k] = i;
               break;
            }
      }
   }
   
   if ( sA < 1. )
      return -1;

   for (j=0; j<image_numpix[itel]; j++)
   {
      int ipix = image_list[itel][j];
      double x = camset->xpix[ipix] - dx;
      double y = camset->ypix[ipix] - dy;
      double A = pixel_amp[itel][ipix];
      sx  += (A * x);
      sxx += (A * x) * x;
      sxy += (A * x) * y;
      sy  += (A * y);
      syy += (A * y) * y;
   }
   sx /= sA;
   sy /= sA;
   sxx = sxx/sA - sx*sx;
   sxy = sxy/sA - sx*sy; 
   syy = syy/sA - sy*sy;
   
   // if ( sxy != 0. )
   if ( fabs(sxy) > 1e-8*fabs(sxx) && fabs(sxy) > 1e-8*fabs(syy) )
   {
      double p1 = syy - sxx, p2 = sxy*sxy;
      double q, r1, r2;
      if ( p2 > 1e-8*(p1*p1) )
         q = p1 + sqrt(p1*p1+4.*p2);
      else
         q = 2.*p2;
      b = 0.5 * q/sxy;
      a = sy - b*sx;
      if ( (r1 = syy + 2.*p2/q) > 0. )
	 length = img_scale * sqrt(r1);
      else
	 length = 0.;
      if ( (r2 = sxx - 2.*p2/q) > 0. )
	 width = img_scale * sqrt(r2);
      else
	 width = 0.;
   }
   else
   {
      if ( fabs(syy) < 1e-8*fabs(sxx) )
         syy = 0.;
      else if ( fabs(sxx) < 1e-8*fabs(syy) )
         sxx = 0.;
      if ( sxx > syy && syy >= 0. )
      {
	 length = img_scale * sqrt(sxx);
	 width = img_scale * sqrt(syy);
	 b = 0.;
	 a = sy;
      }
      else if ( syy >= 0. && sxx >= 0. )
      {
	 length = img_scale * sqrt(syy);
	 width = img_scale * sqrt(sxx);
	 b = 100000.;
	 a = sy - b*sx;
      }
      else
      {
         a = b = 0.;
         length = 
         width = img_scale * 0.001;
      }
   }
   beta = atan(b);
   cb = cos(beta);
   sb = sin(beta);
   if ( (distance = img_scale * sqrt(sx*sx+sy*sy)) == 0. )
      distance = img_scale * 0.001;
   miss = img_scale * fabs(a)/sqrt(b*b+1.);
   if ( miss/distance <= 1. )
      alpha = 180./M_PI * asin(miss/distance);
   else
      alpha = 90.;
   xmean = img_scale * (sx + dx);
   ymean = img_scale * (sy + dy);
   direction = 180./M_PI * (beta + camset->cam_rot);
   orientation = 180./M_PI * (atan2(sy,sx) + camset->cam_rot);
   if ( camset->cam_rot != 0. )
   {
      double rmean = sqrt(xmean*xmean + ymean*ymean);
      double rphi = atan2(ymean,xmean) + camset->cam_rot;
      xmean = rmean * cos(rphi);
      ymean = rmean * sin(rphi);
   }

   sxx = sx3 = sx4 = 0.;
   for (j=0; j<image_numpix[itel]; j++)
   {
      int ipix = image_list[itel][j];
      double x = camset->xpix[ipix] - dx;
      double y = camset->ypix[ipix] - dy;
      double A = pixel_amp[itel][ipix];
      double xp;
      xp =  cb*(x-sx) + sb*(y-sy);
      /* yp = -sb*(x-sx) + cb*(y-sy); */ /* Not used */
      sxx += (A*xp) * xp;
      sx3 += ((A*xp) * xp) * xp;
      sx4 += (((A*xp) * xp) * xp) * xp;
   }

   if ( verbosity >= 0 )
   {
      if ( iprint++ == 0 )
         printf("#@* Lines starting with '@*' contain the following columns:\n"
                "#@*  (1): event\n"
                "#@*  (2): telescope\n"
                "#@*  (3): energy\n"
                "#@*  (4): core distance to telescope\n"
                "#@*  (5): image size (amplitude) [p.e.]\n"
                "#@*  (6): number of pixels in image\n"
                "#@*  (7): width [deg.]\n"
                "#@*  (8): length [deg.]\n"
                "#@*  (9): distance [deg.]\n"
                "#@* (10): miss [deg.]\n"
                "#@* (11): alpha [deg.]\n"
                "#@* (12): orientation [deg.]\n"
                "#@* (13): direction [deg.]\n"
                "#@* (14): image c.o.g. x [deg.]\n"
                "#@* (15): image c.o.g. y [deg.]\n"
                "#@* (16): Xmax [g/cm^2]\n"
                "#@* (17): Hmax [m]\n"
                "#@* (18): Size without tail-cuts (all-pixels sum)\n"
                "#@* (19-23): Hot pixels\n");

      printf("@* %d %d %6.3f %7.2f %7.1f %d %7.4f %7.4f %7.4f %7.4f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.2f  %7.2f  %3.1f %3.1f %3.1f %3.1f %3.1f  %d %d %d %d %d\n",
      hsdata->mc_event.event, 
      camset->tel_id,
      hsdata->mc_shower.energy, 
      line_point_distance(hsdata->mc_event.xcore, hsdata->mc_event.ycore,0.,
                     cos(hsdata->mc_shower.altitude)*cos(hsdata->mc_shower.azimuth),
                     cos(hsdata->mc_shower.altitude)*sin(-hsdata->mc_shower.azimuth),
                     sin(hsdata->mc_shower.altitude), 
                     hsdata->run_header.tel_pos[itel][0],
                     hsdata->run_header.tel_pos[itel][1], 
                     hsdata->run_header.tel_pos[itel][2]),
      sA, image_numpix[itel],
      width, length, distance, miss, alpha, orientation, direction, 
      xmean, ymean,
      hsdata->mc_shower.xmax, hsdata->mc_shower.hmax,
      stot,
      hot_amp[0], hot_amp[1], hot_amp[2], hot_amp[3], hot_amp[4],
      hot_pixel[0], hot_pixel[1], hot_pixel[2], hot_pixel[3], hot_pixel[4]);
   }

   skewness = sx3/pow(sxx,1.5);
   if ( skewness < 0. )
   {
      alpha = 180. - alpha;
      direction += 180.;
   }
   if ( sxx > 0. )
      kurtosis = sx4/(sxx*sxx) - 3.;
   else
      kurtosis = 0.;

   /* Just filling into the first image set. May overwrite existing image data. */
   img->pixels = image_numpix[itel];
   img->cut_id = cut_id;
   img->amplitude = sA;
   img->x = xmean * (M_PI/180.);
   img->y = ymean * (M_PI/180.);
   img->phi = direction * (M_PI/180.);
   img->w = width * (M_PI/180.);
   img->l = length * (M_PI/180.);
   img->skewness = skewness;
   img->kurtosis = kurtosis;
   img->known = 1;

   return 0;
}

/* -------------------------- pixel_timing_analysis ----------------------- */

/** Calculate summary results from pixel timing data. */

/* Note: Can only be used after clean_image_tailcut()
   because it is using external static storage filled there. */

static int pixel_timing_analysis(AllHessData *hsdata, int itel, int nimg);

static int pixel_timing_analysis(AllHessData *hsdata, int itel, int nimg)
{
   TelEvent *teldata = NULL;
   PixelTiming *pixtm = NULL;
   ImgData *img = NULL;
   CameraSettings *camset = NULL;
   PixelSetting *pixset = NULL;
   double sw=0., sx=0., sy=0., sxx=0., sxy=0., syy=0., swd1=0., swd2=0., srt=0.;
   double sdt=0., sdt2=0.;
   int n = 0, j, k, kpeak = -1, kwa = -1, kwr50 = -1, kwr20 = -1, ksr20 = -1, ksr80 = -1;
   double img_scale, b, m, d, time_slice = 1.0 /* A good guess for older data */;
   double res = 0.;

   if ( itel < 0 || itel >= H_MAX_TEL )
      return -1;
   teldata = &hsdata->event.teldata[itel];
   camset = &hsdata->camera_set[itel];
   img_scale = 1./hsdata->camera_set[itel].flen;
   pixset = &hsdata->pixel_set[itel];
   if ( pixset->time_slice > 0. )
      time_slice = pixset->time_slice;

   if ( teldata->img == NULL || nimg < 0 || nimg >= teldata->num_image_sets )
      return -1;
   img = &teldata->img[nimg];

   /* Reset resulting data */
   img->tm_slope = img->tm_residual = img->tm_width1 = img->tm_width2 = img->tm_rise = 0.;

   if ( (pixtm = teldata->pixtm) == NULL )
      return -1;
   if ( !pixtm->known )
      return -1;

   /* Identify where the things that we want actually are */
   for (k=0; k<pixtm->num_types && k<H_MAX_PIX_TIMES; k++)
   {
      if ( pixtm->time_type[k] == PIX_TIME_PEAKPOS_TYPE )
         kpeak = k;
      else if ( pixtm->time_type[k] == PIX_TIME_WIDTH_ABS_TYPE )
         kwa = k;
      else if ( pixtm->time_type[k] == PIX_TIME_WIDTH_REL_TYPE )
      {
         if ( pixtm->time_level[k] > 0.4 && pixtm->time_level[k] < 0.6 )
            kwr50 = k;
         else if ( pixtm->time_level[k] > 0.15 && pixtm->time_level[k] < 0.25 )
            kwr20 = k;
      }
      else if ( pixtm->time_type[k] == PIX_TIME_STARTPOS_REL_TYPE )
      {
         if ( pixtm->time_level[k] > 0.7 && pixtm->time_level[k] < 0.9 )
            ksr80 = k;
         else if ( pixtm->time_level[k] > 0.15 && pixtm->time_level[k] < 0.25 )
            ksr20 = k;
      }
   }
   if ( kpeak < 0 ) /* We need at least the peak position */
      return -1;

   for (j=0; j<image_numpix[itel]; j++)     // External data: image_numpix
   {
      int ipix = image_list[itel][j];       // External data: image_list
      double A, x, y, xr, t, wi, wd1=0., wd2=0., rt=0.;
      if ( ipix < 0 || ipix > pixtm->num_pixels )
         continue;
      if ( pixtm->timval[ipix][0] < 0. )
         continue;
      if ( pixel_disabled[itel][ipix] )
         continue;
      wi = 0.;
      if ( (A=pixel_amp[itel][ipix]) > 0. ) // External data: pixel_amp
         wi = A / (A+100.); /* Assume errors ~1/sqrt(A), levels off for large amplitudes. */
      else
         continue;
      x = camset->xpix[ipix]*img_scale - img->x;
      y = camset->ypix[ipix]*img_scale - img->y;
      xr =  cos(img->phi)*x + sin(img->phi)*y;
      t = pixtm->timval[ipix][kpeak];      /* Pixel peak time */
      if ( kwr50 >= 0 )
         wd1 = pixtm->timval[ipix][kwr50]; /* Pulse width at 50% of peak */
      else if ( kwa >= 0 )
         wd1 = pixtm->timval[ipix][kwa];   /* Pulse width above threshold */
      if ( kwr20 >= 0 )
         wd2 = pixtm->timval[ipix][kwr20]; /* Pulse width at 20% of peak */
      if ( ksr20 > 0 && ksr80 > 0 )
         rt = pixtm->timval[ipix][ksr80] - pixtm->timval[ipix][ksr20]; /* Rise time */
      n++;
      /* Sum up for weighted linear least squares fit */
      sw  += wi;
      sx  += wi * xr;
      sy  += wi * t;
      sxx += wi * xr*xr;
      sxy += wi * xr*t;
      syy += wi * t*t;
      /* Sum up the other stuff, also weighted */
      swd1+= wi * wd1;
      swd2+= wi * wd2;
      srt += wi * rt;
   }

   if ( sw == 0. || n < 2 ) /* Slope is definitely undetermined */
      return 0;
   d = sw*sxx-sx*sx;
   if ( d<1e-10 ) /* Undetermined within numerical accuracy */
      return 0;
   m = (sw*sxy-sx*sy) / d;
   b = (sxx*sy-sx*sxy) / d;

   img->tm_slope = m * time_slice;
   img->tm_width1 = swd1 / sw * time_slice;
   img->tm_width2 = swd2 / sw * time_slice;
   img->tm_rise = srt / sw * time_slice;

   /* Second round for r.m.s. residuals [Do we need that round?] */
   for (j=0; j<image_numpix[itel]; j++)
   {
      int ipix = image_list[itel][j];
      double A, x, y, xr, t, dt, wi;
      if ( ipix < 0 || ipix > pixtm->num_pixels )
         continue;
      if ( pixtm->timval[ipix][0] < 0. )
         continue;
      if ( pixel_disabled[itel][ipix] )
         continue;
      wi = 0.;
      if ( (A=pixel_amp[itel][ipix]) > 0. )
         wi = A / (A+100.); /* Same as above, no need to sum up. */
      else
         continue;
      x = camset->xpix[ipix]*img_scale - img->x;
      y = camset->ypix[ipix]*img_scale - img->y;
      xr =  cos(img->phi)*x + sin(img->phi)*y;
      t = pixtm->timval[ipix][kpeak]; /* Pixel peak time */
      dt = t - b - m*xr;
      sdt += wi * dt;
      sdt2+= wi * dt*dt;
   }
   res = sdt2/sw-(sdt/sw)*(sdt/sw);
   if ( res > 0. )
      img->tm_residual = sqrt(res*(double)n/(double)(n-1)) 
            * time_slice;

   return 0;
}

/* ------------------------------ image_reconstruct ----------------------- */

/** Calibrate and clean image pixels and reconstruct second moments parameters from images. */
/*
 *  @param hsdata Pointer to all available data and configurations.
 *  @param itel   Index of telescope in the relevant arrays (not the ID).
 *  @param cut_id Identifier of the method and levels of tail-cuts in image cleaning.
 *  @param tcl    The lower of the two tail-cut thresholds.
 *  @param tch    The higher of the two tail-cut thresholds.
 *  @param minfrac Which fraction of the reference amplitude is required
 *                for pixels to be included in the final image.
 *                If this number is <= 0.0, the classical scheme is used.
 *  @param nimg   -2: normal image as image 0, image from timing as image 1.
 *                -1: any free image (replacing existing image of matching cut-id if available,)
 *                >= 0: results go to specific image struct in telescope event data.
 *  @param flag_amp_tm 0: Use normal integrated amplitude.
 *                     1: Use integration around global peak position from
 *                        pulse shape analysis. May include all pixels or only selected.
 *                     2: Use integration around local peak position from
 *                        pulse shape analysis. Return 0 for pixels without
 *                        a fairly significant peak.
 *  @param clip_amp: if >0, any calibrated amplitude is clipped not to exceed this value [mean p.e.].
 */

static int image_reconstruct(AllHessData *hsdata, int itel, int cut_id, 
	double tcl, double tch, int lref, double minfrac, int nimg, 
        int flag_amp_tm, double clip_amp);

static int image_reconstruct(AllHessData *hsdata, int itel, int cut_id, 
	double tcl, double tch, int lref, double minfrac, int nimg, 
        int flag_amp_tm, double clip_amp)
{
   TelEvent *teldata = NULL;
   int second_image_from_timing = (nimg==-2 ? 1 : 0);

   if ( itel < 0 || itel >= H_MAX_TEL )
      return -1;

   if ( nb_lists[itel][0].nblist == NULL )
   {
      find_neighbours(&hsdata->camera_set[itel],itel);
   }

   teldata = &hsdata->event.teldata[itel];
   if ( teldata->num_image_sets < 2 )
   {
      nimg = 0; // Only place for one image set, no choice.
      second_image_from_timing = 0;
   }
   if ( nimg == -1 && teldata->img != NULL )  /* No specific image struct? */
   {
      int j, k=-1;
      for ( j=0; j<teldata->num_image_sets; j++ )
      {
         if ( !teldata->img[j].known && k == -1 )
            k = j;
         if ( teldata->img[j].cut_id == cut_id )
         {
            nimg = j;
            break;
         }
      }
      if ( nimg == -1 ) /* No image struct of matching cut id known? */
      {
         if ( k >= 0 )  /* Any image struct not in use? */
            nimg = k;
         else
            nimg = 0;   /* Last resort: use the first image struct */
      }
   }
   if ( nimg == -2 )
      nimg = 0;
   if ( teldata->raw != NULL && teldata->raw->known )
   {
      calibrate_amplitude(hsdata, itel, second_image_from_timing?0:flag_amp_tm, clip_amp);
      if ( clean_image_tailcut(hsdata, itel, tcl, tch, lref, minfrac) < 0 )
         return -1;
      if ( second_moments(hsdata, itel, cut_id, nimg, clip_amp) < 0 )
         return -1;
      pixel_timing_analysis(hsdata,itel,nimg); /* Optional; may fail */

      if ( second_image_from_timing && teldata->num_image_sets > 1) 
      /* We reconstruct both image types */
      {
         calibrate_amplitude(hsdata, itel, flag_amp_tm==0?2:flag_amp_tm, clip_amp);
         if ( clean_image_tailcut(hsdata, itel, tcl, tch, lref, minfrac) < 0 )
            return -1;
         if ( second_moments(hsdata, itel, 2, 1, clip_amp) < 0 )
            return -1;
      }

      return 0;
   }
   else
      return -1;
}


/* ------------------------- clean_raw_data ------------------------------- */

int clean_raw_data (AllHessData *hsdata, int itel, int clean_flag, 
   int tcl, int tch, struct user_parameters *up);
   
int clean_raw_data (AllHessData *hsdata, int itel, int clean_flag, 
   int tcl, int tch, struct user_parameters *up)
{
   int siglev[H_MAX_PIX];
   TelEvent *teldata = &hsdata->event.teldata[itel];
   CameraSettings *camset = &hsdata->camera_set[itel];
   AdcData *raw = teldata->raw;
   PixelTiming *pixtm = teldata->pixtm;
   PixelCalibrated *pixcal = teldata->pixcal;
   int npix = camset->num_pixels;
   int nimg = teldata->image_pixels.pixels;
   int ipix, j, k, kpix;
   int nsig = 0, nstr = 0;
//   int siglist[H_MAX_PIX];
   int have_nb = 1, have_ext = 1;
   
#ifdef CLEAN_DEBUG
printf("** Cleaning data in event %d telescope ID %d with method %d\n",
   teldata->glob_count, teldata->tel_id, clean_flag);
#endif

   if ( clean_flag == 9 ) /* Poor man's solution: rely on pixel timing significance */
   {
      teldata->readout_mode = 9; /* Not advertized and not competitive in most cases */
      return 0;
   }

   /* Was the neighbour finding done yet? */
   if ( nb_lists[itel][0].nblist == NULL )
   {
      find_neighbours(&hsdata->camera_set[itel],itel);
   }

   if ( nb_lists[itel][0].nblist == NULL ||
        nb_lists[itel][0].pix_num_nb == NULL || nb_lists[itel][0].nbsize == 0 )
   {
#ifdef CLEAN_DEBUG
      printf("No neighbour list for telescope ID %d", hsdata->camera_set[itel].tel_id);
      if ( nb_lists[itel][0].nblist != NULL )
         printf(" (npix=%d, nbsize=%d)\n", nb_lists[itel][0].npix, nb_lists[itel][0].nbsize);
      else
         printf(" (NULL)\n");
#endif
      have_nb = 0;
   }
   if ( ext_list[itel].nblist == NULL || 
        ext_list[itel].pix_num_nb == NULL || ext_list[itel].nbsize == 0 )
   {
#ifdef CLEAN_DEBUG
      printf("No extension neighbour list for telescope ID %d\n", hsdata->camera_set[itel].tel_id);
#endif
      have_ext = 0;
   }

   /* Default is not to be in cleaned up data */
   for ( ipix=0; ipix<npix; ipix++ )
   {
      siglev[ipix] = 0;
   }
#ifdef CLEAN_DEBUG
printf("** Pixels in image(s): %d\n",teldata->image_pixels.pixels);
#endif
   /* Pixels passing image cleaning form the core of the cleaned list, */
   /* adding all extension neighbours of those. */
   for ( j=0; j<nimg; j++ )
   {
      int nx, bx;
      ipix = teldata->image_pixels.pixel_list[j];
#ifdef CLEAN_DEBUG
      if ( ipix < 0 || ipix >= npix )
        printf("Pixel outside range\n");
#endif
      siglev[ipix] = 1; /* Pixel is in set passing image cleaning */

      if ( have_ext && ipix<ext_list[itel].npix )
      { 
         nx = ext_list[itel].pix_num_nb[ipix];
         bx = ext_list[itel].pix_first_nb[ipix];
         if ( bx+nx <= ext_list[itel].nbsize )
         for ( k=0; k<nx; k++ )
         {
            kpix = ext_list[itel].nblist[bx+k];
            if ( siglev[kpix] == 0 )
               siglev[kpix] = 4; /* Pixel is an extended neighbour of a cleaned image pixel */
         }
         else
           printf("Extension pixels outside range\n");
      }
      if ( (clean_flag%10) == 4 )
      {
         if ( have_nb && ipix<nb_lists[itel][0].npix )
         {
            nx = nb_lists[itel][0].pix_num_nb[ipix];
            bx = nb_lists[itel][0].pix_first_nb[ipix];
            if ( bx+nx <= nb_lists[itel][0].nbsize )
            for ( k=0; k<nx; k++ )
            {
               kpix = nb_lists[itel][0].nblist[bx+k];
               siglev[kpix] |= 2; /* Pixel is a normal neighbour of a cleaned image pixel */
            }
         else
           printf("Neighbour pixels outside range\n");
         }
      }
   }
   
   /* Disabled pixels (like HV off, no signal) remain off the list */
   if ( any_disabled[itel] )
   {
      for ( ipix=0; ipix<npix; ipix++ )
      {
         if ( pixel_disabled[itel][ipix] && siglev[ipix] )
            siglev[ipix] = 0;
      }
   }
   /* This should contain these pixels as well but make sure we get them all. */
   for ( j=0; j<hsdata->pixel_disabled[itel].num_HV_disabled; j++ )
   {
      ipix = hsdata->pixel_disabled[itel].HV_disabled[j];
      if ( ipix < 0 || ipix >= npix )
        printf("Pixel outside range\n");
      siglev[ipix] = 0;
   }   
   
   /* We may have pixel intensities in 'raw' (AdcData), 'pixtm' (PixelTiming),
      and/or 'pixcal' (PixelCalibrated). Reset the available but now 
      insignificant ones. */
   if ( raw != NULL )
      raw->list_known = raw->list_size = 0;
   for ( ipix=0; ipix<npix; ipix++ )
   {
      /* Now mark the selected pixels everywhere possible as significant */
      if ( siglev[ipix] )
      {
//         siglist[nsig] = ipix;
         if ( raw != NULL )
         {
            int known = 0;
            for ( j=0; j<raw->num_gains; j++ )
               if ( raw->adc_known[j][ipix] )
                  known = 1;
            if ( !known )
            {
               raw->significant[ipix] = 0;
               continue;
            }
            raw->adc_list[nsig] = ipix;
            if ( (clean_flag%10) == 1 )      /* Sums for clean+extension but no traces */
               raw->significant[ipix] = 1;
            else if ( (clean_flag%10) == 2 ) /* Sums for all pixels and traces only for clean+extension */
            {
               raw->significant[ipix] = 0x21;
               nstr++;
            }
            else if ( (clean_flag%10) == 3 ) /* Sums and traces only for clean+extension */
            {
               raw->significant[ipix] = 0x21;
               nstr++;
            }
            else if ( (clean_flag%10) == 4 ) /* Sums for clean+extension, traces only for clean+nb */
            {
               if ( (siglev[ipix] & 0x03) )
               {
                  raw->significant[ipix] = 0x21; /* both if part of cleaned image */
                  nstr++;
               }
               else
                  raw->significant[ipix] = 1; /* only sum if in extension */
            }
            else if ( (clean_flag%10) == 5 ) /* Sums for clean+extension, traces only for clean */
            {
               if ( siglev[ipix] == 1 )
               {
                  raw->significant[ipix] = 0x21; /* both if part of cleaned image */
                  nstr++;
               }
               else
                  raw->significant[ipix] = 1; /* only sum if in extension */
            }
            nsig++;
         }
      }
      else
      {
         if ( raw != NULL )
         {
            if ( (clean_flag%10) == 2 ) /* Sums for all pixels and traces only for clean+extension */
            {
               int known = 0;
               for ( j=0; j<raw->num_gains; j++ )
                  if ( raw->adc_known[j][ipix] )
                     known = 1;
               if ( known )
               {
                  raw->significant[ipix] = 1;
                  raw->adc_list[nsig] = ipix; /* only sum */
//                  siglist[nsig] = ipix;
                  nsig++;
               }
            }
            else /* Nothing stored outside of clean+extension */
               raw->significant[ipix] = 0;
         }
      }
   }
   raw->list_size = nsig;
   raw->list_known = 1;
#ifdef CLEAN_DEBUG
printf("** Telescope %d has %d significant pixels after cleaning (%d with traces)\n",
raw->tel_id, nsig, nstr);
#endif
   if ( (clean_flag%10) >= 2 && (teldata->readout_mode & 0xff) < 2 &&
         raw->num_samples > 1 )
      teldata->readout_mode = 2;
   if ( (clean_flag%10) == 1 && (teldata->readout_mode & 0xff) > 0 )
      teldata->readout_mode = 0;
   if ( nstr == 0 ) /* With no significant pixels for traces, turn off storing samples. */
      teldata->readout_mode = 0;

   if ( pixtm != NULL && pixtm->known )
   {
      int igain, ntm = 0;
      pixtm->list_type = pixtm->list_size = 0;
      for ( ipix=0; ipix<npix; ipix++ )
      {
         if ( siglev[ipix] == 0 )
         {
            pixtm->timval[ipix][0] = -1;
            for ( igain=0; igain < pixtm->num_gains; igain++ )
            {
               pixtm->pulse_sum_loc[igain][ipix] = 0;
               pixtm->pulse_sum_glob[igain][ipix] = 0;
            }
         }
         else if ( pixtm->timval[ipix][0] != -1 )
            ntm++;
      }
      if ( ntm == 0 ) /* If nothing left, no need to write a pixel timing block */
         pixtm->known = 0;
   }
   
   if ( pixcal != NULL && pixcal->known )
   {
      int ncal = 0;
      pixcal->list_known = 0;
      for ( ipix=0; ipix<npix; ipix++ )
      {
         if ( pixcal->significant[ipix] )
         {
            if ( siglev[ipix] )
            {
               pixcal->pixel_list[ncal++] = ipix;
            }
            else
               pixcal->significant[ipix] = 0;
         }
      }
      pixcal->list_size = ncal;
      pixcal->list_known = 1;
      if ( ncal == 0 )
         pixcal->known = 0;
   }

   return 0;
}


/* ----------------------------- shower_reconstruct ----------------------- */

/** Shower reconstruction (geometrical reconstruction only) */

static int shower_reconstruct (AllHessData *hsdata, const double *min_amp_tel, 
      const size_t *min_pix_tel, int cut_id);

static int shower_reconstruct (AllHessData *hsdata, const double *min_amp_tel, 
      const size_t *min_pix_tel, int cut_id)
{
   double amp[H_MAX_TEL];
   double ximg[H_MAX_TEL], yimg[H_MAX_TEL], phi[H_MAX_TEL], disp[H_MAX_TEL];
   double xtel[H_MAX_TEL], ytel[H_MAX_TEL], ztel[H_MAX_TEL];
   double az[H_MAX_TEL], alt[H_MAX_TEL], flen[H_MAX_TEL], cam_rot[H_MAX_TEL];
   double ref_az, ref_alt;
   double shower_az=0., shower_alt=0., xcore=0., ycore=0., var_dir=0., var_core=0.;
   int flag = 0, ntel = 0, itel, rc, ntrg=0;
   int list[H_MAX_TEL];
   static int iprint = 0;

   int pattern = 0;
   
   hsdata->event.shower.known = 0;

   for (itel=0; itel<hsdata->event.num_tel && itel<H_MAX_TEL; itel++)
   {
      TelEvent *teldata = &hsdata->event.teldata[itel];
      ImgData *img;
      double r_cog;
      double min_amp = min_amp_tel[itel];
      if ( min_amp <= 0. )
         min_amp = 80.;

      if ( !teldata->known || teldata->img == NULL )
         continue;
      ntrg++;
      int s_img = -1, j_img;
      /* Check if we have an image of exactly the requested type. */
      for ( j_img=0; j_img < teldata->num_image_sets; j_img++ )
         if ( cut_id == teldata->img[j_img].cut_id &&
              teldata->img[j_img].known )
            s_img = j_img;
      /* If not exactly that type available but the type asked for is
         rather general, just take the first image. */
      if ( s_img < 0 && (cut_id == 0 || cut_id == 1) ) 
         for ( j_img=0; j_img < teldata->num_image_sets; j_img++ )
            if ( teldata->img[j_img].known )
            {
               s_img = j_img;
               break;
            }
      if ( s_img < 0 )
         continue;
      img = &teldata->img[s_img];
      if ( (amp[ntel] = img->amplitude) < min_amp || 
           (size_t)img->pixels < min_pix_tel[itel] )
         continue;
      if ( camera_radius_eff[itel] == 0. )
         store_camera_radius(&hsdata->camera_set[itel],itel);
      
      r_cog = sqrt(img->x*img->x + img->y*img->y);
      if ( r_cog > 0.8 * camera_radius_eff[itel] )
         continue;
      ximg[ntel] = img->x;
      yimg[ntel] = img->y;
      phi[ntel]  = img->phi;
      disp[ntel] = (img->l>0. && img->l>img->w) ? 1.-img->w/img->l : 1e-3;
      flen[ntel] = 1.0; // not: hsdata->camera_set[itel].flen;
      xtel[ntel] = hsdata->run_header.tel_pos[itel][0];
      ytel[ntel] = hsdata->run_header.tel_pos[itel][1];
      ztel[ntel] = hsdata->run_header.tel_pos[itel][2];
      if ( hsdata->event.trackdata[itel].cor_known )
      {
         az[ntel]  = hsdata->event.trackdata[itel].azimuth_cor;
         alt[ntel] = hsdata->event.trackdata[itel].altitude_cor;
      }
      else
      {
         az[ntel]  = hsdata->event.trackdata[itel].azimuth_raw;
         alt[ntel] = hsdata->event.trackdata[itel].altitude_raw;
      }
      cam_rot[ntel] = 0.;
      list[ntel] = hsdata->event.teldata[itel].tel_id;
      ntel++;
      if ( itel < 16 )
         pattern |= (1<<itel);
   }
   
   if ( ntel < 2 )
      return 0;

   ref_az  = hsdata->run_header.direction[0];
   ref_alt = hsdata->run_header.direction[1];
   
   rc = shower_geometric_reconstruction(ntel, amp, ximg, yimg, phi, disp,
         xtel, ytel, ztel, az, alt, flen, cam_rot, ref_az, ref_alt, flag,
         &shower_az, &shower_alt, &var_dir, &xcore, &ycore, &var_core);
   if ( rc >= 1 )
   {
      ShowerParameters *shower = &hsdata->event.shower;
      shower->num_trg = shower->num_read = ntrg;
      shower->num_img = ntel;
      shower->img_pattern = pattern; // Cannot be represented this way with many telescopes.
      for ( itel=0; itel<ntel; itel++ )
         shower->img_list[itel] = list[itel];
      shower->known = 1;
      if ( rc == 1 )
      {
         if ( ntel > 2 )
            shower->result_bits = 3;
         else
            shower->result_bits = 1;
      }
      else if ( rc >= 2 )
      {
         if ( ntel > 2 )
            shower->result_bits = 15;
         else
            shower->result_bits = 5;
      }
      shower->Az  = shower_az;
      shower->Alt = shower_alt;
      shower->xc  = xcore;
      shower->yc  = ycore;
      shower->mscl = shower->mscw = -1.;
      shower->energy = -1.;
      shower->xmax = 0.;
      if ( ntel > 2 )
      {
         /* Error estimates based on variance, assuming equal contributions
            in both coordinates. */
         shower->err_dir1 = shower->err_dir2 = sqrt(var_dir/(ntel-2.)/2.);
         shower->err_dir3 = 0.; /* We don't have an error matrix, assume independent */
         shower->err_core1 = shower->err_core2 = sqrt(var_core/(ntel-2.)/2.);
         shower->err_core3 = 0.;
      }
      else
      {
         shower->err_dir1 = shower->err_dir2 = 0.;
         shower->err_dir3 = 0.;
         shower->err_core1 = shower->err_core2 = 0.,
         shower->err_core3 = 0.;
      }

if ( verbosity > 0 )
 printf("### Shower reconstruction with %d telescopes: Az=%lf deg, Alt=%lf deg, xc=%lf m, yc=%lf m\n",
   ntel,shower_az*(180./M_PI),shower_alt*(180./M_PI),xcore,ycore);

      if ( verbosity >= 0 )
      {
         if ( iprint++ == 0 )
            printf("#@; Lines starting with '@;' contain the following columns:\n"
                   "#@;  (1): event\n"
                   "#@;  (2): number of telescopes triggered\n"
                   "#@;  (3): number of images used in shower reconstruction\n"
                   "#@;  (4): results bit pattern\n"
                   "#@;  (5): shower azimuth [deg] (with bit 0)\n"
                   "#@;  (6): shower altitude [deg] (with bit 0)\n"
                   "#@;  (7): angle between reconstructed and true direction [deg]\n"
                   "#@;  (8): core position x [m] (with bit 2)\n"
                   "#@;  (9): core position y [m] (with bit 2)\n"
                   "#@; (10): horizontal displacement between reconstructed and true core [m]\n"
                   "#@; (11): MSCL [deg] (with bit 4)\n"
                   "#@; (12): MSCW [deg] (with bit 4)\n"
                   "#@; (13): energy [TeV] (with bit 6)\n"
                   "#@; (14): Xmax [g/cm^2] (with bit 8)\n");

         if ( hsdata->event.shower.known )
            printf("@; %d %d %d %d   %f %f %f   %f %f %f   %f %f  %f %f\n",
               hsdata->mc_event.event, 
               hsdata->event.shower.num_trg, hsdata->event.shower.num_img,
               hsdata->event.shower.result_bits,
               hsdata->event.shower.Az*(180./M_PI), hsdata->event.shower.Alt*(180./M_PI),
               angle_between(hsdata->event.shower.Az, hsdata->event.shower.Alt,
                  hsdata->mc_shower.azimuth, hsdata->mc_shower.altitude) * (180./M_PI),
               hsdata->event.shower.xc, hsdata->event.shower.yc,
               sqrt((hsdata->event.shower.xc-hsdata->mc_event.xcore)*
                    (hsdata->event.shower.xc-hsdata->mc_event.xcore) +
                    (hsdata->event.shower.yc-hsdata->mc_event.ycore)*
                    (hsdata->event.shower.yc-hsdata->mc_event.ycore)),
               hsdata->event.shower.mscl*(180./M_PI), 
               hsdata->event.shower.mscw*(180./M_PI),
               hsdata->event.shower.energy, hsdata->event.shower.xmax);
      }
   }

   return ntel;
}

/* --------------------------------- reconstruct -------------------------- */

/** Image/shower reconstruction function
 *
 *  @param hsdata Pointer to all available data and configurations.
 *  @param reco_flag If >= 3 then redo image cleaning before shower reconstruction.
 *                 If >= 4 then the total image intensities are re-determined
 *                 and that may change which images are used or not in
 *                 the shower reconstruction.
 *  @param min_amp The minimum amplitude required in images (telescope-specific,
 *                 that means requiring an array of at least size H_MAX_TEL).
 *  @param min_pix The minimum number of pixels required in images (telescope-specific).
 *  @param tcl    The lower of the two tail-cut thresholds (telescope-specific).
 *  @param tch    The higher of the two tail-cut thresholds (telescope-specific).
 *  @param lref   Determines which pixel, after sorting by amplitude,
 *                will be used as providing the reference amplitude (telescope-specific).
 *                Example: use 3 for the third hottest pixel.
 *                If this number is <= 0, the classical scheme is used.
 *  @param minfrac Which fraction of the reference amplitude is required
 *                for pixels to be included in the final image (telescope-specific).
 *                If this number is <= 0.0, the classical scheme is used.
 *  @param nimg   Which of (sometimes) several images should be filled?
 *                Use -1 to replace an existing image of the same cut id
 *                (if such an image exists) or add another image (if there
 *                is free space for it) or replace the first image (if
 *                all else fails).
 *                Use -2 to indicate that image analysis from normal
 *                integrated amplitude should go into first image and
 *                (if available) that from pixel timing (around local
 *                peak position or otherwise global peak position) should
 *                go into the second image.
 *  @param flag_amp_tm 0: Use normal integrated amplitude.
 *                     1: Use integration around global peak position from
 *                        pulse shape analysis. May include all pixels or only selected.
 *                     2: Use integration around local peak position from
 *                        pulse shape analysis. Return 0 for pixels without
 *                        a fairly significant peak.
 *
 */

int reconstruct (AllHessData *hsdata, int reco_flag, 
      const double *min_amp, const size_t *min_pix, const double *tcl, const double *tch, 
      const int *lref, const double *minfrac, int nimg, int flag_amp_tm, int clean_flag);

int reconstruct (AllHessData *hsdata, int reco_flag, 
      const double *min_amp, const size_t *min_pix, const double *tcl, const double *tch, 
      const int *lref, const double *minfrac, int nimg, int flag_amp_tm, int clean_flag)
{
   int itel, cut_id = 0;
   if ( min_amp == NULL || min_pix == NULL || 
        tcl == NULL || tch == NULL ||
        lref == NULL || minfrac == NULL )
      return -1;

#ifdef DEBUG_PIXEL_NB
printf("Reconstruct data\n");
#endif

   /* Note: cut_id must be in the range 0 to 255 */
   cut_id = 1; /* Unspecified classical two-level tailcut */
   /* Note: We used to assign cut IDs based on lower tailcut level but don't anymore. */ 

   if ( reco_flag >= 4 )
      show_total_amp = 1;
   else
      show_total_amp = 0;

   if ( reco_flag >= 3 )
   {
      int have_raw_data = 0;
      for (itel=0; itel<hsdata->run_header.ntel; itel++)
      {
         if ( hsdata->event.teldata[itel].known &&
              hsdata->event.teldata[itel].raw != NULL &&
              hsdata->event.teldata[itel].raw->known )
         {
            have_raw_data = 1;
            break;
         }
      }
      if ( have_raw_data )
      {
         for (itel=0; itel<hsdata->run_header.ntel; itel++)
         {
            if ( hsdata->event.teldata[itel].known &&
                 hsdata->event.teldata[itel].raw != NULL &&
                 hsdata->event.teldata[itel].raw->known )
            {
               int tel_type = user_get_type(itel);
               struct user_parameters *up = user_get_parameters(tel_type);

               /* (Re-) summing pixel intensities, if samples available and integrator given: */
               if ( up->i.integrator > 0 )
               {
                  pixel_integration(hsdata, itel, up);
               }

               /* (Re-) calculate Hillas parameters: */
               if ( hsdata->event.teldata[itel].img != NULL )
               {
                  double clip_amp = up->d.clip_amp;
                  if ( up->d.focal_length != 0. ) /* Override MC nominal focal length */
                  {
                     if ( verbosity >= 0 && hsdata->camera_set[itel].flen != up->d.focal_length )
                        printf("Effective focal length of telescope %d of type %d changed from %5.3f to %5.3f m.\n",
                           hsdata->event.teldata[itel].tel_id, tel_type, 
                           hsdata->camera_set[itel].flen, up->d.focal_length);
                     hsdata->camera_set[itel].flen = up->d.focal_length;
                  }
                  image_reconstruct(hsdata, itel, cut_id, 
                     tcl[itel], tch[itel], lref[itel], minfrac[itel], 
                     nimg, flag_amp_tm, clip_amp);
               }
               
               if ( clean_flag )
               {
                  clean_raw_data(hsdata, itel, clean_flag, tcl[itel], tch[itel], up);
               }
            }
         }
      }
   }

   shower_reconstruct(hsdata, min_amp, min_pix, cut_id);

   return 0;
}

/* ----------------------------- set_reco_verbosity ----------------------- */

void set_reco_verbosity(int v)
{
   verbosity = v;
}
