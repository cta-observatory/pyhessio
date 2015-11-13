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

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "basic_ntuple.h"

static int list_init = 0;
static int with_true = 0;

/**
 * List the parameters useful for event selection plus some
 * more parameters which should not be used for event selection.
 *
 * @param f Output file, to be opened beforehand.
 * @param b Pointer to the struct containing all the relevant numbers.
 * @param wtr Non-zero on first call to write also true MC parameters.
 */

int list_ntuple(FILE *f, const struct basic_ntuple *b, int wtr)
{
   if ( f == NULL || ferror(f) )
      return -1;
   if ( ! list_init )
   {
      int n = 0;
      with_true = wtr;
      fprintf(f,"# Basic selection n-tuple: version 2 (%s true MC data)\n",
         with_true?"with":"without");
      fprintf(f,"# Column %d: Log10(E/TeV)\n", ++n);
      fprintf(f,"# Column %d: Mean core distance [m]\n", ++n);
      fprintf(f,"# Column %d: Mean DISP parameter\n", ++n);
      fprintf(f,"# Column %d: Angle between shower direction and test position [deg]\n", ++n);
      fprintf(f,"# Column %d: Scatter of pair-wise intersections [deg]\n", ++n);
      fprintf(f,"# Column %d: Mean scaled reduced width (mscrw)\n", ++n);
      fprintf(f,"# Column %d: Sigma of scrw in different telescopes\n", ++n);
      fprintf(f,"# Column %d: Mean scaled reduced length (mscrl)\n", ++n);
      fprintf(f,"# Column %d: Sigma of scrl in different telescopes\n", ++n);
      fprintf(f,"# Column %d: Xmax (g/cm^2)\n", ++n);
      fprintf(f,"# Column %d: Sigma of Xmax between telescopes\n", ++n);
      fprintf(f,"# Column %d: Sigma(E)/E as used for dE cut\n", ++n);
      fprintf(f,"# Column %d: Chi2(E)/(ntel-1) as used for dE2 cut\n", ++n);
      fprintf(f,"# Column %d: Corrected mean time slope [ns/deg/100 m]\n", ++n);
      fprintf(f,"# Column %d: (central trigger scatter, not filled yet)????\n", ++n);
      fprintf(f,"# Column %d: No. if images used\n", ++n);
      fprintf(f,"# Column %d: No. of telescopes triggered.\n", ++n);
      fprintf(f,"# Column %d: No. of telescopes which should have triggered but did not\n", ++n);
      fprintf(f,"# Column %d: No. of telescopes with zero time slope well outside light pool\n", ++n);
      fprintf(f,"# Column %d: Total number of used pixels in all used telescopes\n", ++n);
      fprintf(f,"# Column %d: Acceptance level\n", ++n);
    if ( with_true )
    {
      fprintf(f,"# Separator '|', followed by true MC values, not to be used for selection:\n");
      fprintf(f,"# Column %d: Primary particle ID\n", ++n);
      fprintf(f,"# Column %d: Run\n", ++n);
      fprintf(f,"# Column %d: Event\n", ++n);
      fprintf(f,"# Column %d: Event weight\n", ++n);
      fprintf(f,"# Column %d: Log10(true energy/TeV)\n", ++n);
      fprintf(f,"# Column %d: Depth of first interaction\n", ++n);
      fprintf(f,"# Column %d: Depth of maximum\n", ++n);
      fprintf(f,"# Column %d: True core x position [m]\n", ++n);
      fprintf(f,"# Column %d: True core y position [m]\n", ++n);
      fprintf(f,"# Column %d: True shower azimuth [deg]\n", ++n);
      fprintf(f,"# Column %d: True shower altitude [deg]\n", ++n);
      fprintf(f,"# Separator '|', followed by auxilliary reconstructed data:\n");
      fprintf(f,"# Column %d: Core x position [m]\n", ++n);
      fprintf(f,"# Column %d: Core y position [m]\n", ++n);
      fprintf(f,"# Column %d: Shower azimuth [deg]\n", ++n);
      fprintf(f,"# Column %d: Shower altitude [deg]\n", ++n);
    }
      fprintf(f,"#\n");
      list_init = 1;
   }

   // fprintf(f,"%f\t", b->weight);
   fprintf(f,"%f", b->lg_e);
   fprintf(f,"\t%f", b->rcm);
   fprintf(f,"\t%f", b->mdisp);
   fprintf(f,"\t%f\t%f", b->theta*(180./M_PI), b->sig_theta*(180./M_PI));
   fprintf(f,"\t%f\t%f", b->mscrw, b->sig_mscrw);
   fprintf(f,"\t%f\t%f", b->mscrl, b->sig_mscrl);
   fprintf(f,"\t%f\t%f", b->xmax, b->sig_xmax);
   fprintf(f,"\t%f\t%f", b->sig_e, b->chi2_e);
   fprintf(f,"\t%f\t%f", b->tslope, b->tsphere);
   fprintf(f,"\t%zu\t%zu\t%zu", b->n_img, b->n_trg, b->n_fail);
   fprintf(f,"\t%zu", b->n_tsl0);
   fprintf(f,"\t%zu", b->n_pix);
   fprintf(f,"\t%zu", b->acceptance);
  if ( with_true )
  {
   fprintf(f,"|%d",  b->primary);
   fprintf(f,"\t%d", b->run);
   fprintf(f,"\t%d", b->event);
   fprintf(f,"\t%f", b->weight);
   fprintf(f,"\t%f", b->lg_e_true);
   fprintf(f,"\t%f", b->xfirst_true);
   fprintf(f,"\t%f", b->xmax_true);
   fprintf(f,"\t%f", b->xc_true);
   fprintf(f,"\t%f", b->yc_true);
   fprintf(f,"\t%f", b->az_true);
   fprintf(f,"\t%f", b->alt_true);
   fprintf(f,"|%f",  b->xc);
   fprintf(f,"\t%f", b->yc);
   fprintf(f,"\t%f", b->az);
   fprintf(f,"\t%f", b->alt);
  }
   fprintf(f,"\n");

   return 0;
}
