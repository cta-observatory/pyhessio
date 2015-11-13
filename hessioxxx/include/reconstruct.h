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

/* reconstruct.c */
double line_point_distance(double xp1, double yp1, double zp1, double cx, double cy, double cz, double x, double y, double z);
int reconstruct(AllHessData *hsdata, int reco_flag, 
   const double *min_amp, const size_t *min_pix, const double *tcl, const double *tch, 
   const int *lref, const double *minfrac, int nimg, int flag_amp_tm);
int store_camera_radius(CameraSettings *camset, int itel);
double get_camera_radius(int itel, int maxflag);
void select_calibration_channel (int chn);
int calibrate_amplitude(AllHessData *hsdata, int itel, 
   int flag_amp_tm, double clip_amp);
double calibrate_pixel_amplitude(AllHessData *hsdata, int itel, 
   int ipix, int flag_amp_tm, int itime, double clip_amp);
double calibrate_pixel_sample_amplitude(AllHessData *hsdata, int itel, 
   int ipix, int flag_amp_tm, int itime, double clip_sample_amp);
void set_reco_verbosity(int v);
int set_disabled_pixels(AllHessData *hsdata, int itel, double broken_pixels_fraction);

#ifdef __cplusplus
}
#endif
