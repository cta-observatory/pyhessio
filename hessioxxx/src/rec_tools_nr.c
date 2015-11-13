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

/** @file rec_tools_nr.h
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
 *  @date    @verbatim CVS $Date: 2010/11/15 15:35:00 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.4 $ @endverbatim
 */
/* ================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rec_tools.h"

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
   double x, double y, double z);

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

/* =================== angles_to_offset ====================== */
/**
 *  @short Transform telescope and object Alt/Az to offset in camera.
 *
 *  Transform from given telescope and object angles (Az/Alt) to
 *  the offset the object has in the camera plane.
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

