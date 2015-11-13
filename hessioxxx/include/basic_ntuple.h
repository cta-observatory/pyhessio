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

/** @file basic_ntuple.h
 *  @short Desclaration of the basic_ntuple struct.
 */

#ifndef BASIC_NTUPLE_HEADER
#define BASIC_NTUPLE_HEADER 1

#ifdef __cplusplus
extern "C" {
#endif

/** A struct with basic per-shower parameters, to be used as an
 *  n-tuple in the event selection.
 */

struct basic_ntuple
{
   // True MC data, to be ignored for selection:
   int primary;     /**< Primary particle ID. */
   int run;         /**< Simulation run number. */
   int event;       /**< Event number (100*shower number + array number) */
   double weight;   /**< Event weight, not to be used for selection (based on true energy). */
   double lg_e_true;/**< log10(true energy of primary). */
   double xfirst_true;/**< Atmospheric depth of first interaction. */
   double xmax_true;/**< True shower maximum atmospheric depth (not well defined with few particles). */
   double xc_true;  /**< True core position at detection level (x coordinate). */
   double yc_true;  /**< True core position at detection level (y coordinate). */
   double az_true;  /**< True shower direction (Azimuth). */
   double alt_true; /**< True shower direction (Altitude). */

   // Reconstructed parameters not suitable for selection:
   double xc;       /**< Reconstructed core position at detection level (x coordinate). */
   double yc;       /**< Reconstructed core position at detection level (y coordinate). */
   double az;       /**< Reconstructed shower direction (Azimuth). */
   double alt;      /**< Reconstructed shower direction (Altitude). */

   // Geometrical parameters: 
   double rcm;      /**< Mean core distance of telescopes used in reconstruction. */
   double mdisp;    /**< Mean DISP (1.-width/length) of usable images. */
   double theta;    /**< Angle between source position and rec. shower direction. */
   double sig_theta;/**< R.m.s. of theta of telescopes pairs (if > 2 tel.). */
   double mscrw;    /**< Mean scaled reduced width. */
   double sig_mscrw;/**< R.m.s. of scaled reduced widths of individual images. */
   double mscrl;    /**< Mean scaled reduced length. */
   double sig_mscrl;/**< R.m.s. of scaled reduced lengths of indvidual images. */
   double xmax;     /**< Depth of shower maximum. */
   double sig_xmax; /**< R.m.s. of Xmax from individual telescopes/images. */
   // Energy parameters (lookups based):
   double lg_e;     /**< Log10 of reconstructed energy. */
   double sig_e;    /**< Relative error estimate on E (NOT the r.m.s. of individual estimates). */
   double chi2_e;   /**< Consistency of individual energy estimates as reduced chi**2 value. */
   // Image and trigger timing parameters:
   double tslope;   /**< Core distance corrected mean time slope (deg/ns/100 m). */
   double tsphere;  /**< R.m.s. of trigger times from spherical propagation from shower max. */

   // Integer-type parameters:
   size_t n_img;    /**< Number of used images. */
   size_t n_trg;    /**< Number of triggered telescopes. */
   size_t n_fail;   /**< Number of failed triggers (telescopes expected to trigger). */
   size_t n_tsl0;   /**< Number of images with zero time slope well outside light pool. */
   size_t n_pix;    /**< Total number of used pixels in all used images. */

   /// Event acceptance level by standard selection scheme
   /// (0: no; 1: shape cuts; 2: +angular cut; 3: +dE cut; 4: +dE2 cut; 5: +Hmax cut.
   size_t acceptance;
};


int list_ntuple(FILE *f, const struct basic_ntuple *b, int wtr);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif

