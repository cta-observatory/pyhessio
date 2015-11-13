/* ============================================================================

Copyright (C) 2006, 2009  Konrad Bernloehr

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

#ifdef __cplusplus
extern "C" {
#endif

struct user_parameters
{
   struct /* Only type int */
   {
      /** 1: HESS-style analysis standard cuts; 2: hard cuts; 3: loose cuts. */
      int user_flags;
      /** The minimum number of significant pixels in usable images. */
      int min_pix;
      /** Reconstruction level flag */
      int reco_flag;
      /** Minimum and maximum number of usable images for events used in analysis. */
      int min_tel_img;
      int max_tel_img;
      /** Which pixel's amplitude is used as reference */
      int lref;

      /** The type of pixel intensity integration scheme */
      /** 0: none (implicitly all samples), 1: simple, 2: around global peak, */
      /** 3: around local peak, 4: around peak in neighbour pixels. */
      int integrator;
      /** Integration-scheme-specific integer parameters, typically:  */
      /** number of bins to integrate and some offset value from start or back from detected peak. */
      int integ_param[2];
      /** Integer type thresholds for significance in ADC units (one per gain) */
      int integ_thresh[2];
      /** Set to 1 if integration over small window should not rescale for fraction of single p.e. trace. */
      int integ_no_rescale;
      
      /** Required trigger type (bit pattern: bit 0 = majo, 1=asum, 2=dsum) */
      int trg_req;
   } i;

   struct /* Only type double */
   {
      double source_offset_deg;
      /** Difference between generated MC spectrum (e.g. E^-2.0) and
       *  assumed source spectrum (e.g. E^-2.5), e.g. case d_sp_idx = -0.5 .
       */
      double d_sp_idx;
      /** The minimum amplitude [ peak p.e. ] of images usable for the analysis. */
      double min_amp;
      /** The lower and upper tail cuts for the standard two-level tail-cut scheme. */
      double tailcut_low;
      double tailcut_high;
      /** Minimum fraction of reference amplitude is needed */
      double minfrac;
      double max_theta_deg;
      double theta_scale;
      double de2_cut_param[4];

      double mscrw_min[4];
      double mscrw_max[4];
      double mscrl_min[4];
      double mscrl_max[4];

      double eres_cut_param[4];
      double hmax_cut_param;
      double min_theta_deg;
      /** Pixel outside this radius (if > 0) should be ignored in image reconstruction. */
      double camera_clipping_deg;

      double theta_escale[4]; /**< If the angular acceptance deviates from the 80% containment. */

      double clip_amp; /**< Pixel intensity clipped to this value after calibration, if this param is not zero. */

      /** Integration-scheme- and gain-specific floating-point parameters. */
      double d_integ_param[2][4];

      /** Calibration scale from mean-p.e. units to experimental units (0.0: like HESS). */
      double calib_scale;

      /** Radii for initial neighbour pixel search */
      double r_nb[3]; /**< Maximum search radii for neighbours [pixel diameter] */
      double r_ne; /**< Radius for extending significant pixels in image cleaning [pixel diameter] */
      
      /* Ranges for reconstructed impact position in array (global). */
      double impact_range[3]; /**< [0]: maximum distance of array center from shower axis, [1],[2]: max. |x|,|y| of core in ground plane. */
      /* Ranges for true impact position in array (global). */
      double true_impact_range[3]; /**< As for impact_ranhe */
      /* Maximum core distance of telescopes used in analysis after geometric reconstruction. */
      double max_core_distance;
   } d;
};
typedef struct user_parameters UserParameters;

void user_init_parameters(void);
struct user_parameters* user_get_parameters(int itype);

void user_set_lookup_file (const char *fname);
void user_set_histogram_file (const char *fname);

void user_set_telescope_type (int itype);
int user_set_tel_type_param_by_str(const char *str);
int which_telescope_type (const struct hess_camera_settings_struct *cam_set);

int user_get_type (int itel);

void user_set_spectrum(double di);
void user_set_impact_range (double *impact_range);
void user_set_true_impact_range (double *true_impact_range);
void user_set_max_core_distance (double rt);
void user_set_min_amp(double a);
void user_set_tail_cuts(double tcl, double tch, int lref, double minfrac);
void user_set_min_pix(int mpx);
void user_set_reco_flag(int rf);
void user_set_tel_img(int tmn,int tmx);
void user_set_tel_list (size_t min_tel, size_t ntel, int *tel_id);
void user_set_max_theta(double thmax, double thscale, double thmin);
void user_set_de_cut (double *dec);
void user_set_de2_cut (double *de2c);
void user_set_hmax_cut (double hmaxc);
void user_set_shape_cuts (double wmin, double wmax, double lmin, double lmax);
void user_set_width_max_cut (double *wmx);
void user_set_length_max_cut (double *lmx);
void user_set_clipping (double dc);
void user_set_clipamp (double cpa);
void user_set_verbosity(int v);
void user_set_flags (int uf);
void user_set_auto_lookup (int al);
void user_set_theta_escale (double *the);
void user_set_diffuse_mode(int dm, double oar[]);
void user_set_integrator(int scheme);
void user_set_integ_window(int nsum, int noff);
void user_set_integ_threshold(int ithg, int itlg);
void user_set_trg_req(int trg_req);
void user_set_integ_no_rescale (int no);
void user_set_calib_scale (double s);
void user_set_nb_radius (double *r);
void user_set_nxt_radius (double r);
int user_selected_event(void);
int do_user_ana (AllHessData *hsdata, unsigned long item_type, int stage);

#ifdef __cplusplus
}
#endif
