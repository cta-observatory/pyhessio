/* ============================================================================

Copyright (C) 2003, 2009, 2011  Konrad Bernloehr

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
 *  @date    @verbatim CVS $Revision: 1.58 $ @endverbatim
 *  @version @verbatim CVS $Date: 2015/05/27 11:36:48 $ @endverbatim
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

#define H_MAX_NB1 8
#define H_MAX_NB2 24
static int neighbours1[H_MAX_TEL][H_MAX_PIX][H_MAX_NB1]; 
//static int neighbours2[H_MAX_TEL][H_MAX_PIX][H_MAX_NB2];
static int nnb1[H_MAX_TEL][H_MAX_PIX];
//static int nnb2[H_MAX_TEL][H_MAX_PIX];
static int has_nblist[H_MAX_TEL];
static int px_shape_type[H_MAX_TEL];

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

      nvd = pixdis->num_HV_disabled;
      for ( ivd=0; ivd<nvd; ivd++ )
         if ( pixdis->HV_disabled[ivd] == ipix )
         {
            pixel_disabled[itel][ipix] = 1;
            any_disabled[itel] = 1;
         }
   }

   return 0;
}

/* --------------------------- find_neighbours ---------------------------- */

/** Find the list of neighbours for each pixel. */

static int find_neighbours(CameraSettings *camset, int itel);

static int find_neighbours(CameraSettings *camset, int itel)
{
   int npix = camset->num_pixels;
   int i, j;
   int stat_st[6] = {0, 0, 0, 0, 0, 0 };
   double asum = 0., dsum = 0., aod2 = 0.;

   for (i=0; i<npix; i++)
   {
      if ( pixel_disabled[itel][i] )
         continue;
      asum += camset->area[i];
      dsum += camset->size[i];
      for (j=0; j<H_MAX_NB1; j++)
         neighbours1[itel][i][j] = -1;
      nnb1[itel][i] = 0;
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
            if ( nnb1[itel][i] < H_MAX_NB1 )
               neighbours1[itel][i][nnb1[itel][i]++] = j;
            if ( nnb1[itel][j] < H_MAX_NB1 )
               neighbours1[itel][j][nnb1[itel][j]++] = i;
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
   has_nblist[itel] = 1;

   asum /= (((double) npix)+1e-10);
   dsum /= (((double) npix)+1e-10);
   aod2 = asum/(dsum*dsum);

   if ( stat_st[0] > 0 && stat_st[2] > 0 && 
        stat_st[1] == 0 && stat_st[3] == 0 )
   {
      px_shape_type[itel] = 2; /* probably square pixels */
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
         px_shape_type[itel] = 1;
      else if ( stat_st[2] > 0 && stat_st[0] == 0 )
         px_shape_type[itel] = 3;
      else
      {
         px_shape_type[itel] = 0; /* Fall back to circular */
         if ( aod2 < 0.99*M_PI/4. || aod2 > 1.01*M_PI/4. )
         {
            fprintf(stderr,
               "Pixel positions in telescope %d indicate round pixels but area/size^2 does not match.\n",
               camset->tel_id);
         }
      }
      if ( px_shape_type[itel] != 0 )
      {
         if ( aod2 < 0.99*sqrt(3.)/2. || aod2 > 1.01*sqrt(3.)/2. )
         {
            if ( aod2 >= 0.99*M_PI/4. && aod2 <= 1.01*M_PI/4. )
            {
               px_shape_type[itel] = 0; // Round pixels on hexagonal pattern
            }
            else
               fprintf(stderr,
                  "Pixel positions in telescope %d indicate hexagonal pixels but area/size^2 does not match.\n"
                  "Expected %f (or %f for round pixels) but found %f\n",
                  camset->tel_id, sqrt(3.)/2., M_PI/4., aod2);
         }
      }
   }
#if 0
   fprintf(stderr,"Pixel shape type of telescope #%d seems to be %d (stat = %d, %d, %d, %d)\n",
      itel, px_shape_type[itel], stat_st[0], stat_st[1], stat_st[2], stat_st[3]);
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
               /* corresponding to sum_samples bins. Add remaining pedestal. */
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
               /* corresponding to sum_samples bins. Add remaining pedestal. */
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
   if ( !has_nblist[itel] )
      find_neighbours(&hsdata->camera_set[itel],itel);

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
      for ( inb=0; inb<nnb1[itel][ipix]; inb++ )
      {
         int ipix_nb = neighbours1[itel][ipix][inb];
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
      if ( knb == 0 )
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
               /* corresponding to sum_samples bins. Add remaining pedestal. */
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

static int set_integration_correction(AllHessData *hsdata, int itel, int nbins, int noff);

static int set_integration_correction(AllHessData *hsdata, int itel, int nbins, int noff)
{
   int igain, ibin, iphase;
   if ( hsdata == NULL || itel < 0 || itel >= H_MAX_TEL )
      return -1;
   PixelSetting *ps = &hsdata->pixel_set[itel];
   CameraOrganisation *co = &hsdata->camera_org[itel];

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

      /* Now sum up given interval starting from peak, averaging over phase */
      for ( iphase=0; iphase<5; iphase++ )
      {
         double ti = ((iphase*0.2-0.4)-noff) * st + ipeak;
         for ( ibin=0; ibin<nbins; ibin++ )
            sum += qpol(ibin*st + ti, ps->lrefshape, ps->refshape[igain]);
      }
      sum *= 0.2;
      if ( sum > 0. && asum > 0. )
         integration_correction[itel][igain] = asum/sum;
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
   if ( hsdata == NULL || up == NULL )
      return -1;
   if ( integration_correction[itel][0] == 0. )
      set_integration_correction(hsdata, itel, up->i.integ_param[0], up->i.integ_param[1]);
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

   if ( hsdata == NULL || itel < 0 || itel >= H_MAX_TEL )
      return -1;
   npix = hsdata->camera_set[itel].num_pixels;
   
   teldata = &hsdata->event.teldata[itel];
   teldata->image_pixels.pixels = image_numpix[itel] = 0;
   if ( !teldata->known || teldata->raw == NULL )
      return -1;
   if ( !teldata->raw->known )
      return -1;

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
         for (j=0; j<nnb1[itel][i]; j++)
            if ( pass_low[neighbours1[itel][i][j]] )
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
         for (j=0; j<nnb1[itel][i]; j++)
            if ( pass_high[neighbours1[itel][i][j]] )
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
   if ( nimg < 0 || nimg >= teldata->num_image_sets )
      return -1;
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

   if ( !has_nblist[itel] )
      find_neighbours(&hsdata->camera_set[itel],itel);

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
   double shower_az, shower_alt, xcore, ycore, var_dir, var_core;
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
      const int *lref, const double *minfrac, int nimg, int flag_amp_tm);

int reconstruct (AllHessData *hsdata, int reco_flag, 
      const double *min_amp, const size_t *min_pix, const double *tcl, const double *tch, 
      const int *lref, const double *minfrac, int nimg, int flag_amp_tm)
{
   int itel, cut_id = 0;
   if ( min_amp == NULL || min_pix == NULL || 
        tcl == NULL || tch == NULL ||
        lref == NULL || minfrac == NULL )
      return -1;

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
                  image_reconstruct(hsdata, itel, cut_id, 
                     tcl[itel], tch[itel], lref[itel], minfrac[itel], 
                     nimg, flag_amp_tm, clip_amp);
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
