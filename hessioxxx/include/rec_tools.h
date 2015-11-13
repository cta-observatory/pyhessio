/* ============================================================================

Copyright (C) 2001, 2009  Konrad Bernloehr

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

#ifndef REC_TOOLS_H__LOADED
#define REC_TOOLS_H__LOADED 1

#ifdef __cplusplus
extern "C" {
#endif

void angles_to_offset(double obj_azimuth, double obj_altitude, 
   double azimuth, double altitude, double focal_length, 
   double *xoff, double *yoff);
void offset_to_angles(double xoff, double yoff, double azimuth, double altitude, 
   double focal_length, double *obj_azimuth, double *obj_altitude);
void get_shower_trans_matrix(double azimuth, double altitude, double trans[][3]);
void cam_to_ref(double ximg, double yimg, double phi, 
   double ref_azimuth, double ref_altitude, double cam_rot, 
   double azimuth, double altitude, double focal_length, 
   double *axref, double *ayref, double *phiref);
int intersect_lines(double xp1, double yp1, double phi1, 
   double xp2, double yp2, double phi2, double *xs, double *ys, double *sang);
int shower_geometric_reconstruction(int ntel, const double *amp, 
   const double *ximg, const double *yimg, 
   const double *phi, const double *disp,
   const double *xtel, const double *ytel, const double *ztel, 
   const double *az, const double *alt, 
   const double *flen, const double *cam_rot, 
   double ref_az, double ref_alt,  int flag,
   double *shower_az, double *shower_alt, double *var_dir,
   double *xc, double *yc, double *var_core);
double angle_between(double azimuth1, double altitude1, 
   double azimuth2, double altitude2);
double line_point_distance (double xp1, double yp1, double zp1, 
   double cx, double cy, double cz,
   double x, double y, double z);

#ifdef __cplusplus
}
#endif

#endif
