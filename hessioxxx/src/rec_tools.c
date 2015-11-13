/* ============================================================================

Copyright (C) 2000, 2009  Konrad Bernloehr

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

/** @file rec_tools.h
 *  @short Tools for shower geometric reconstruction.
 *
 *  Shower geometric reconstruction based on the major axes of the
 *  telescope images. The image parameters from each telescope are
 *  transformed to a common reference frame first before the
 *  average intersection point of all images is calculated in
 *  plane coordinates.
 *
 *  @author  Konrad Bernloehr 
 *  @date    2000, 2009
 *  @date    @verbatim CVS $Date: 2014/05/07 13:08:25 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.17 $ @endverbatim
 */
/* ================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "initial.h"
#include "rec_tools.h"
#include "io_hess.h"

/* ------------------- line_point_distance --------------------- */
/**
 *  Distance between a straight line and a point in space.
 *
 *  @param xp1, yp1, zp1:  reference point on the line
 *  @param cx, cy, cz:  direction cosines of the line
 *  @param x, y, z:	point in space
 *
 *  @return distance
*/

double line_point_distance (double xp1, double yp1, double zp1, 
   double cx, double cy, double cz,
   double x, double y, double z)
{
   double a, a1, a2, a3, b;
   
   a1 = (y-yp1)*cz - (z-zp1)*cy;
   a2 = (z-zp1)*cx - (x-xp1)*cz;
   a3 = (x-xp1)*cy - (y-yp1)*cx;
   a  = a1*a1 + a2*a2 + a3*a3;
   b = cx*cx + cy*cy + cz*cz;
   if ( a<0. || b<= 0. )
      return -1;
   return sqrt(a/b);
}

#if 0
static void rotate_xy (double cam_rot, double *x, double *y)
{
   double r = (*x)*(*x) + (*y)*(*y);
   double p = atan2(*y,*x) + cam_rot;
   
   *x = r * cos(p);
   *y = r * sin(p);
}
#endif

/* =================== angles_to_offset ====================== */
/**
 *  @short Transform telescope and object Alt/Az to offset in camera.
 *
 *  Transform from given telescope and object angles (Az/Alt) to
 *  the offset the object has in the camera plane.
 *
 *  This does not account for any rotation of the camera and its pixels.
*/

void angles_to_offset (double obj_azimuth, double obj_altitude,
   double azimuth, double altitude, double focal_length, 
   double *xoff, double *yoff)
{
   double daz = obj_azimuth - azimuth;
   double coa = cos(obj_altitude);

   double xp0 = -cos(daz) * coa;
   double yp0 = sin(daz) * coa;
   double zp0 = sin(obj_altitude);

   double cx = sin(altitude);
   double sx = cos(altitude);

   double xp1 = cx*xp0 + sx*zp0;
   double yp1 = yp0;
   double zp1 = -sx*xp0 + cx*zp0;
   
   if ( xp1 == 0 && yp1 == 0 ) /* On-axis ? */
   {
      *xoff = *yoff = 0.;
      return;
   }
#if 0 /* The more intuitive way but taking more CPU time */
   double q = acos(zp1); // off-axis angle
   double d = tan(q) * focal_length; // distance to camera center
   double alpha = atan2(yp1,xp1); // orientation

   *xoff = d * cos(alpha);
   *yoff = d * sin(alpha);
#else /* Actually we don't need any further function calls */
   double s = focal_length / zp1;

   *xoff = s * xp1;
   *yoff = s * yp1;
#endif
}

/* =================== offset_to_angles ====================== */
/**
 *  @short Transform from offset in camera to corresponding Az/Alt.
 *
 *  Transform from the offset an object or image has in the
 *  camera plane of a telescope to the corresponding Az/Alt.
 *
 *  This does not account for any rotation of the camera and its pixels.
 *  (xoff and yoff are assumed to be corrected for camera rotation).
*/

void offset_to_angles (double xoff, double yoff, 
   double azimuth, double altitude, double focal_length, 
   double *obj_azimuth, double *obj_altitude)
{
   if ( xoff == 0. && yoff == 0. ) /* Avoid division by zero */
   {
      *obj_azimuth = azimuth;
      *obj_altitude = altitude;
      return;
   }
   else
   {
      double d = sqrt(xoff*xoff+yoff*yoff);
      double q = atan(d/focal_length);

      double sq = sin(q);
      double xp1 = xoff * (sq/d);
      double yp1 = yoff * (sq/d);
      double zp1 = cos(q);

      double cx = sin(altitude);
      double sx = cos(altitude);

      double xp0 = cx*xp1 - sx*zp1;
      double yp0 = yp1;
      double zp0 = sx*xp1 + cx*zp1;

      *obj_altitude = asin(zp0);
      *obj_azimuth  = atan2(yp0,-xp0) + azimuth;
      if ( *obj_azimuth < 0. )
         *obj_azimuth += 2.*M_PI;
      else if ( *obj_azimuth >= (2.*M_PI ) )
         *obj_azimuth -= 2.*M_PI;
   }
}

/* ================= get_shower_trans_matrix ================= */
/**
 *  @short Calculate transformation matrix.
 *
 *  Calculate transformation matrix from horizontal reference
 *  frame to one z axis in the given Az/Alt direction and
 *  the x axis in the plane defined by Az/Alt and zenith.
*/

void get_shower_trans_matrix (double azimuth, double altitude,
   double trans[][3])
{
   double cos_z = sin(altitude);
   double sin_z = cos(altitude);
   double cos_az = cos(azimuth);
   double sin_az = sin(azimuth);
   
   trans[0][0] = cos_z*cos_az;
   trans[1][0] = sin_az;
   trans[2][0] = sin_z*cos_az;
   
   trans[0][1] = -cos_z*sin_az;
   trans[1][1] = cos_az;
   trans[2][1] = -sin_z*sin_az;
   
   trans[0][2] = -sin_z;
   trans[1][2] = 0.;
   trans[2][2] = cos_z;
}

/* =================== cam_to_ref ====================== */
/**
 *  @short Transform from one camera to common reference frame.
 *
 *  Transform from the camera plane coordinate system of
 *  a telescope looking to altitude/azimuth to a plane
 *  coordinate system of a potential telescope looking
 *  to a reference direction ref_azimuth,ref_altitude
 *  and having unit focal length. Rotation of image
 *  angles is accounted for but not imaging errors.
*/

void cam_to_ref (double ximg, double yimg, double phi,
   double ref_azimuth, double ref_altitude, double cam_rot,
   double azimuth, double altitude, double focal_length,
   double *axref, double *ayref, double *phiref)
{
   double s, c;
   double ximg_rot, yimg_rot;
   double azm_img, alt_img, dphi1, dphi2;
   
   if ( fabs(cam_rot) > 1e-14 )
   {
      c = cos(cam_rot);
      s = sin(cam_rot);
      ximg_rot = ximg*c + yimg*s;
      yimg_rot = yimg*c - ximg*s;
   }
   else
   {
      ximg_rot = ximg;
      yimg_rot = yimg;
   }

   offset_to_angles(ximg_rot, yimg_rot, azimuth, altitude, focal_length,
      &azm_img, &alt_img);
   dphi1 = -atan(tan(azm_img-azimuth)*sin(alt_img));
   angles_to_offset(azm_img, alt_img, ref_azimuth, ref_altitude, 1.0,
      axref, ayref);
   dphi2 = -atan(tan(azm_img-ref_azimuth)*sin(alt_img));

   *phiref = phi + cam_rot + (dphi2-dphi1); 
}

/* =================== intersect_lines ====================== */
/**
 *  @short Intersect pairs of lines.
 *
 *  Intersect a pair of straight lines in a plane and return
 *  the intersection point and the angle at which the lines intersect.
*/

int intersect_lines (double xp1, double yp1, double phi1,
   double xp2, double yp2, double phi2, double *xs, double *ys, double *sang)
{
   double A1, B1, C1;
   double A2, B2, C2;
   double detAB, detBC, detCA;
   double s1, c1, s2, c2;
   
   /* Hesse normal form for line 1 */
   s1 = sin(phi1);
   c1 = cos(phi1);
   A1 = s1;
   B1 = -c1;
   C1 = yp1*c1 - xp1*s1;

   /* Hesse normal form for line 2 */
   s2 = sin(phi2);
   c2 = cos(phi2);
   A2 = s2;
   B2 = -c2;
   C2 = yp2*c2 - xp2*s2;

   detAB = (A1*B2-A2*B1);
   detBC = (B1*C2-B2*C1);
   detCA = (C1*A2-C2*A1);
   
   if ( fabs(detAB) < 1e-14 ) /* parallel */
   {
      if ( sang )
      	 *sang = 0.;
      if ( fabs(detBC) < 1e-14 && fabs(detCA) < 1e-14 ) /* same lines */
      {
      	 /* We could take any point on the line but use the middle here. */
      	 *xs = 0.5*(xp1+xp2);
	 *ys = 0.5*(yp1+yp2);
	 return 2;
      }
      *xs = *ys = 0.;
      return 0;
   }
   
   *xs = detBC / detAB;
   *ys = detCA / detAB;
   
   if ( sang != NULL )
   {
      double dx1 = (*xs-xp1);
      double dx2 = (*xs-xp2);
      double dy1 = (*ys-yp1);
      double dy2 = (*ys-yp2);
      double dr1 = sqrt(dx1*dx1+dy1*dy1);
      double dr2 = sqrt(dx2*dx2+dy2*dy2);
      double cos_ang;
      if ( dr1*dr2 == 0. )
      	 *sang = 0.;
      else
      {
         cos_ang = (dx1*dx2+dy1*dy2) / (dr1*dr2);
         if ( cos_ang >= 1. )
            *sang = 0.;
         else if ( cos_ang <= -1. )
            *sang = M_PI;
         else
            *sang = acos(cos_ang);
      }
   }
   
   return 1;
}

#ifndef MAX_TEL
# ifdef H_MAX_TEL
#  define MAX_TEL H_MAX_TEL
# else
#  define MAX_TEL 100
# endif
#endif

static inline double square(double a) { return a*a; }

/* ================ shower_geometric_reconstruction ================ */
/**
 *  @short Simple reconstruction by intersecting pairs of lines.
 *
 *  Simple geometric shower reconstruction by intersecting pairs
 *  of straigh lines (from major axis of second moments ellipses
 *  after transformation to a common plane), first for the
 *  shower direction and then for the core position.
 *  No errors on reconstructed direction or core position are
 *  calculated.
 *  This should sooner or later be superceded by a fit procedure
 *  taking advantage of estimated errors on image positions and angles.
 *
 *  @param ntel   The number of telescopes with suitable images.
 *  @param amp    The image amplitudes in each suitable telescope [p.e.].
 *  @param ximg   The image c.o.g. x positions in the local camera coordinate systems.
 *  @param yimg   The image c.o.g. y positions in the local camera coordinate systems.
 *  @param phi    The image major axis direction [rad].
 *  @param disp   The DISP parameter (1.-width/length), used for giving preference
 *                to elongated images. Set all to 1.0 if unknown or no preference wanted.
 *                Can also be passed as a NULL pointer instead.
 *  @param xtel   The x coordinate of the telescope positions within array [m].
 *  @param ytel   The y coordinate of the telescope positions within array [m].
 *  @param ztel   The z coordinate of the telescope positions within array [m].
 *  @param az     The azimuth angles to which the telescopes are pointing
 *                (N->E->S->W) [rad].
 *  @param alt    The altitude angles to which the telescopes are pointing [rad].
 *  @param flen   The focal length to which ximg and yimg are scaled
 *                (1.0 if in units of radians, otherwise flen is in meters).
 *  @param cam_rot Camera rotation angle [rad].
 *  @param ref_az The reference azimuth angle (system nominal azimuth) [rad].
 *  @param ref_alt The reference altitude angle (system nominal altitude) [rad].
 *  @param flag   Use the reconstucted direction to derive the core position (0)
 *                or use the nominal direction for that (1 or any other non-zero). 
 *                The second version may sightly improve core distance and thus 
 *                energy accuracy for well-defined point sources.
 *  @param shower_az Return the reconstructed shower azimuth angle (N->E->S->W) [rad].
 *  @param shower_alt Return the reconstructed shower altitude angle [rad].
 *  @param var_dir Variance (dx**2+dy**2)/ntel of reconstructed direction for more than two images.
 *                Can be NULL if you are not interested in it.
 *  @param xc     Return the reconstructed core position x coordinate (at z=0) [m].
 *  @param yc     Return the reconstructed core position y coordinate (at z=0) [m].
 *  @param var_core Variance (dx**2+dy**2)/ntel of reconstructed core position for more than two images.
 *                Can be NULL if you are not interested in it.
*/

int shower_geometric_reconstruction (int ntel, 
   const double *amp, const double *ximg, const double *yimg, 
   const double *phi, const double *disp,
   const double *xtel, const double *ytel, const double *ztel,
   const double *az, const double *alt, 
   const double *flen, const double *cam_rot,
   double ref_az, double ref_alt, int flag,
   double *shower_az, double *shower_alt, double *var_dir,
   double *xc, double *yc, double *var_core)
{
   double xang[MAX_TEL], yang[MAX_TEL], aphi[MAX_TEL];
   double xt[MAX_TEL], yt[MAX_TEL];
   double xs, ys, sa, w, xh, yh, zh;
   double sum_xs = 0., sum_ys = 0., sum_w = 0., sum_xs2 = 0., sum_ys2 = 0.;
   int itel, jtel;
   double trans[3][3];
   
   if ( ntel < 2 )
   {
      fprintf(stderr,"Not enough images.\n");
      return 0;
   }
   else if ( ntel > MAX_TEL )
   {
      fprintf(stderr,"Too many images, current limit is %d.", MAX_TEL);
      return 0;
   }
   
   /* Convert positions of images to a common reference frame. */
   for ( itel=0; itel<ntel; itel++ )
   {
      xang[itel] = yang[itel] = aphi[itel] = 0.;
      if ( amp[itel] <= 10. )
      	 continue;
      cam_to_ref(ximg[itel], yimg[itel], phi[itel], ref_az, ref_alt,
      	 cam_rot[itel], az[itel], alt[itel], flen[itel],
      	 &xang[itel], &yang[itel], &aphi[itel]);

      /* Note: at this point the angles are only corrected for */
      /* camera rotation but not for imaging errors. */
   }
   sum_xs = sum_ys = sum_w = 0.;
   for ( itel=0; itel<ntel; itel++ )
   {
      double amp_red;
      if ( amp[itel] <= 10. )
      	 continue;
      for ( jtel=0; jtel<itel; jtel++ )
      {
	 if ( amp[jtel] <= 10. )
      	    continue;
         /* Do the lines have an intersection point? */
         if ( intersect_lines(xang[itel],yang[itel],aphi[itel],
               xang[jtel],yang[jtel],aphi[jtel],
               &xs, &ys, &sa) != 1 )
            continue;
         /* "Reduced amplitude" like reduced mass in celestial dynamics. */
         amp_red = (amp[itel]*amp[jtel])/(amp[itel]+amp[jtel]);
#define WT_DISP 1
#ifdef WT_DISP
         /* New: lower weight for round images, low amplitude on low
            amplitude, and again on small intersection angles. */
         /* FIXME: We might also want to use the NSB rates where different
            telescope types are involved. */
         if ( disp != NULL )
            w = square(amp_red * sin(sa) * disp[itel] * disp[jtel]);
         else
            w = square(amp_red * sin(sa));
#else
         w = (amp[itel]<amp[jtel]?amp[itel]:amp[jtel]) * sin(sa); /* Old style */
#endif
         sum_xs += xs * w;
         sum_xs2+= xs*xs * w;
         sum_ys += ys * w;
         sum_ys2+= ys*ys * w;
         sum_w  += w;
      }
   }

   if ( fabs(sum_w) < 1e-10 )
      return -1;

   /* Weighted average of intersection points. */
   sum_xs /= sum_w;
   sum_ys /= sum_w;

   /* Variance between the different intersection coordinates. */
   if ( var_dir != NULL )
   {
      if ( ntel > 2 )
         *var_dir = ( (sum_xs2/sum_w - sum_xs*sum_xs) +
                      (sum_ys2/sum_w - sum_ys*sum_ys) );
      else
         *var_dir = 0.;
   }

   /* Convert this planar position to az/alt angles. */
   offset_to_angles(sum_xs, sum_ys, ref_az, ref_alt, 1.0,
      shower_az, shower_alt);

   *shower_az -= (2.*M_PI)*floor(*shower_az/(2.*M_PI));

   /* Get transformation matrix from horizontal rectangular */
   /* coordinates to shower plane. */
   if ( flag == 0 ) /* assume reconstructed direction */
      get_shower_trans_matrix(*shower_az, *shower_alt, trans);
   else             /* assume reference direction */
      get_shower_trans_matrix(ref_az, ref_alt, trans);
      
   for ( itel=0; itel<ntel; itel++ )
   {
      xt[itel] = trans[0][0]*xtel[itel] + 
                 trans[0][1]*ytel[itel] +
      	         trans[0][2]*ztel[itel];
      yt[itel] = trans[1][0]*xtel[itel] + 
                 trans[1][1]*ytel[itel] +
      	         trans[1][2]*ztel[itel];
      /* Assume parallel projection onto the shower plane (zt[itel]=0) */
      /* Given some assumption on the distance to the shower maximum, */
      /* e.g. 10 km/cos(z), we could also try a central projection */
      /* by calculating the intersection (x,y) at z=0. */
   }
   
   sum_xs = sum_ys = sum_w = sum_xs2 = sum_ys2 = 0.;
   for ( itel=0; itel<ntel; itel++ )
      for ( jtel=0; jtel<itel; jtel++ )
      {
         double amp_red;
         /* Do the lines have an intersection point? */
         if ( intersect_lines(xt[itel],yt[itel],aphi[itel],
               xt[jtel],yt[jtel],aphi[jtel],
               &xs, &ys, &sa) != 1 )
            continue;

         /* Weighting of intersection point of pair of lines. */
         amp_red = (amp[itel]*amp[jtel])/(amp[itel]+amp[jtel]);
#ifdef WT_DISP
         if ( disp != NULL )
            w = square(amp_red * sin(sa) * disp[itel] * disp[jtel]);
         else
            w = square(amp_red * sin(sa));
#else
         w = (amp[itel]<amp[jtel]?amp[itel]:amp[jtel]) * sin(sa); /* Old style */
#endif
         sum_xs += xs * w;
         sum_xs2+= xs*xs * w;
         sum_ys += ys * w;
         sum_ys2+= ys*ys * w;
         sum_w  += w;
      }

   if ( sum_w == 0. )
      return 1;

   /* Reverse transformation matrix is just transposed but z!=0 */
   xs = sum_xs / sum_w;
   ys = sum_ys / sum_w;
   xh = trans[0][0] * xs +
        trans[1][0] * ys;
   yh = trans[0][1] * xs +
        trans[1][1] * ys;
   zh = trans[0][2] * xs +
        trans[1][2] * ys;

   /* Extrapolation to the ground (detection level) */
   *xc = xh - trans[2][0]*zh/trans[2][2];
   *yc = yh - trans[2][1]*zh/trans[2][2];

   if ( var_core != NULL )
   {
      if ( ntel > 2 )
         *var_core = ( (sum_xs2/sum_w - sum_xs*sum_xs) +
                       (sum_ys2/sum_w - sum_ys*sum_ys) );
      else
         *var_core = 0.;
   }

   return 2;
}

/* ==================== angle_between ====================== */
/**
 * @short Calculate the angle between two directions given in spherical coordinates.
 *
 * @return The angle between the two directions in units of radians.
 */

double angle_between (double azimuth1, double altitude1, double azimuth2, double altitude2)
{
   double ax1 = cos(azimuth1)*cos(altitude1);
   double ay1 = sin(-azimuth1)*cos(altitude1);
   double az1 = sin(altitude1);
   double ax2 = cos(azimuth2)*cos(altitude2);
   double ay2 = sin(-azimuth2)*cos(altitude2);
   double az2 = sin(altitude2);
   double cos_ang = ax1*ax2 + ay1*ay2 + az1*az2;
   /* Check for rounding errors pushing us outside the valid range. */
   if ( cos_ang <= -1. )
      return M_PI;
   else if ( cos_ang >= 1. )
      return 0.;
   else
      return acos(cos_ang);
}

