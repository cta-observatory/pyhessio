/* ============================================================================

Copyright (C) 2009  Konrad Bernloehr

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

#ifndef CAMERA_IMAGE_H__LOADED
#define CAMERA_IMAGE_H__LOADED 1

#ifdef __cplusplus
extern "C" {
#endif

/* In camera_image.c: */
void hesscam_ps_plot(const char *image_fname, AllHessData *hsdata, 
   int itel, int type, int amp_tm, double clip_amp);

void hesscam_type_sum_plot(const char *image_fname, AllHessData *hsdata, int teltype);

#ifdef __cplusplus
}
#endif

#endif
