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

/* ================================================================ */
/** @file gen_lookup.c
 *  @short Generate image shape and energy lookups for user analysis 
 *         in read_hess.
 *
 *  Read_hess must be run with user analysis once and the generated
 *  histogram file is used by this program to generate the lookups.
 *  The lookup file is used in the next round of read_hess user analysis,
 *  if found under the desired name. Look at the last lines of output
 *  from read_hess (or at the beginning, right after the history) to
 *  see how the lookup file should be called (depends on tail cut
 *  parameters, and so on).
 *
 *  @date    @verbatim CVS $Revision: 1.21 $ @endverbatim
 *  @version @verbatim CVS $Date: 2012/05/11 13:18:48 $ @endverbatim
 */
/* ================================================================ */

#include "initial.h"
#include "io_basic.h"
#include "histogram.h"
#include "io_histogram.h"
#include "fileopen.h"

HISTOGRAM *h18000, *h18001, *h18011, *h18012, *h18021, *h18022, *h18051, *h18052;
HISTOGRAM *h18100, *h18101, *h18111, *h18112, *h18121, *h18122, *h18151, *h18152;
HISTOGRAM *h18113, *h18114, *h18123, *h18124, *h18140, *h18141, *h18153, *h18154;
HISTOGRAM *h18005, *h18006, *h18071, *h18072, *h18081, *h18082;
HISTOGRAM *h18105, *h18106, *h18171, *h18172, *h18181, *h18182;
HISTOGRAM *h18173, *h18174, *h18183, *h18184;
HISTOGRAM *h18200, *h18201, *h18211, *h18212;
HISTOGRAM *h18301, *h18311, *h18321, *h18322;

void fill_gaps(void);

/**
  * @short Fill gaps in those histograms used for generating the lookups.
  *
  * Depending on the physical quantities we have different
  * strategies for interpolation/extrapolation/smoothing.
  */

void fill_gaps()
{
   int min_count = 50;
   int nx = h18000->nbins;
   int ny = h18000->nbins_2d;
   int iy, ix, iy2;
   // double *ty, *ty2;

   if ( nx == 0 || ny == 0 )
      exit(1);
   // ty = (double *) calloc(ny,sizeof(double)); // ???
   // ty2 = (double *) calloc(ny,sizeof(double)); // ???

   for (iy=0; iy<ny; iy++)
   {
      double y = h18000->specific_2d.real.lower_limit + 
             (iy+0.5) * (h18000->specific_2d.real.upper_limit
                - h18000->specific_2d.real.lower_limit) /
                 h18000->nbins_2d;
      for (ix = 0; ix <nx; ix++ )
      {
         int ix1 = ix, ix2 = ix, kx = 0;
         double x = h18000->specific.real.lower_limit + 
                (ix+0.5) * (h18000->specific.real.upper_limit
                   - h18000->specific.real.lower_limit) /
                    h18000->nbins;
         double ns = h18000->extension->ddata[ix+nx*iy],
                nxs= h18000->extension->ddata[ix+nx*iy];
         double s  = h18001->extension->ddata[ix+nx*iy],
                sw = h18011->extension->ddata[ix+nx*iy],
                sww= h18012->extension->ddata[ix+nx*iy],
                sl = h18021->extension->ddata[ix+nx*iy],
                sll= h18022->extension->ddata[ix+nx*iy],
                xs = h18001->extension->ddata[ix+nx*iy],
                se = h18051->extension->ddata[ix+nx*iy],
                see= h18052->extension->ddata[ix+nx*iy];
         while ( ns < min_count && kx<=5 )
         {
            kx++;
            if ( ix >= kx )
            {
               ix1 = ix-kx;
               ns += h18000->extension->ddata[ix1+nx*iy];
               s  += h18001->extension->ddata[ix1+nx*iy];
               sw += h18011->extension->ddata[ix1+nx*iy];
               sww+= h18012->extension->ddata[ix1+nx*iy];
               sl += h18021->extension->ddata[ix1+nx*iy];
               sll+= h18022->extension->ddata[ix1+nx*iy];
            }
            if ( ix+kx < nx )
            {
               ix2 = ix+kx;
               ns += h18000->extension->ddata[ix2+nx*iy];
               s  += h18001->extension->ddata[ix2+nx*iy];
               sw += h18011->extension->ddata[ix2+nx*iy];
               sww+= h18012->extension->ddata[ix2+nx*iy];
               sl += h18021->extension->ddata[ix2+nx*iy];
               sll+= h18022->extension->ddata[ix2+nx*iy];
            }
         }
         if ( ns >= min_count-1. )
         {
            fill_histogram(h18100,x,y,ns);
            fill_histogram(h18101,x,y,s);
            fill_histogram(h18111,x,y,sw);
            fill_histogram(h18112,x,y,sww);
            fill_histogram(h18121,x,y,sl);
            fill_histogram(h18122,x,y,sll);
         }
         for ( iy2=iy-3; iy2<=iy+3; iy2++ )
         {
            double w = (iy2<iy) ? 1.-(iy-iy2)/4. : (iy2>iy) ? 1.-(iy2-iy)/4. : 0.;
            if ( iy2 >= 0 && iy2 < ny )
            {
               nxs += w * h18000->extension->ddata[ix+nx*iy2];
               xs  += w * h18001->extension->ddata[ix+nx*iy2];
               se  += w * h18051->extension->ddata[ix+nx*iy2],
               see += w * h18052->extension->ddata[ix+nx*iy2];
            }
         }
         if ( nxs >= min_count-1. )
         {
            fill_histogram(h18140,x,y,nxs);
            fill_histogram(h18141,x,y,xs);
            fill_histogram(h18151,x,y,se);
            fill_histogram(h18152,x,y,see);
         }
      }

   }


   if ( h18005 == NULL || h18171 == NULL )
      return;

   nx = h18005->nbins;
   ny = h18005->nbins_2d;
   
   for (iy=0; iy<ny; iy++)
   {
      double y = h18005->specific_2d.real.lower_limit + 
             (iy+0.5) * (h18005->specific_2d.real.upper_limit
                - h18005->specific_2d.real.lower_limit) /
                 h18005->nbins_2d;
      for (ix = 0; ix <nx; ix++ )
      {
         double x = h18005->specific.real.lower_limit + 
                (ix+0.5) * (h18005->specific.real.upper_limit
                   - h18005->specific.real.lower_limit) /
                    h18005->nbins;
         double ns = h18005->extension->ddata[ix+nx*iy];
         double s  = h18006->extension->ddata[ix+nx*iy],
                sr = h18071->extension->ddata[ix+nx*iy],
                srr= h18072->extension->ddata[ix+nx*iy],
                sd = h18081->extension->ddata[ix+nx*iy],
                sdd= h18082->extension->ddata[ix+nx*iy];
         
         /* No gap filling / smoothing here yet */
         
         fill_histogram(h18105,x,y,ns);
         fill_histogram(h18106,x,y,s);
         fill_histogram(h18171,x,y,sr);
         fill_histogram(h18172,x,y,srr);
         fill_histogram(h18181,x,y,sd);
         fill_histogram(h18182,x,y,sdd);
      }
   }
}

void gen_image_lookups(void);

/**
 *  @short Generate the lookups for image shape parameters and energy.
 */

void gen_image_lookups()
{
   int nx = h18000->nbins;
   int ny = h18000->nbins_2d;
   int iy, ix;
   if ( nx == 0 || ny == 0 )
      exit(1);

   clear_histogram(h18113);
   clear_histogram(h18114);
   clear_histogram(h18123);
   clear_histogram(h18124);
   clear_histogram(h18153);
   clear_histogram(h18154);
   clear_histogram(h18173);
   clear_histogram(h18174);
   clear_histogram(h18183);
   clear_histogram(h18184);

   for (iy=0; iy<ny; iy++)
   {
      double y = h18000->specific_2d.real.lower_limit + 
             (iy+0.5) * (h18000->specific_2d.real.upper_limit
                - h18000->specific_2d.real.lower_limit) /
                 h18000->nbins_2d;
      for (ix = 0; ix <nx; ix++ )
      {
         double x = h18000->specific.real.lower_limit + 
                (ix+0.5) * (h18000->specific.real.upper_limit
                   - h18000->specific.real.lower_limit) /
                    h18000->nbins;
         double ns = h18100->extension->ddata[ix+nx*iy];
         double s  = h18101->extension->ddata[ix+nx*iy];
         double sw = h18111->extension->ddata[ix+nx*iy];
         double sww= h18112->extension->ddata[ix+nx*iy];
         double sl = h18121->extension->ddata[ix+nx*iy];
         double sll= h18122->extension->ddata[ix+nx*iy];
         double nxs= h18140->extension->ddata[ix+nx*iy];
         double xs = h18141->extension->ddata[ix+nx*iy];
         double se = h18151->extension->ddata[ix+nx*iy];
         double see= h18152->extension->ddata[ix+nx*iy];

         double wm = 0., ws = 1e-10, lm = 0., ls = 1e-10, em = 0., es = 1e-10;
         if ( ns >= 2 && s > 0. )
         {
            wm = sw / s;
            ws = sqrt((sww/s - wm*wm) * ns/(ns-1));
            lm = sl / s;
            ls = sqrt((sll/s - lm*lm) * ns/(ns-1));
         }
         if ( nxs >= 2 && xs > 0. )
         {
	    em = se / xs;
	    es = sqrt((see/xs - em*em) * nxs/(nxs-1));
         }

         fill_histogram(h18113,x,y,wm);
         fill_histogram(h18114,x,y,ws);
         fill_histogram(h18123,x,y,lm);
         fill_histogram(h18124,x,y,ls);
         fill_histogram(h18153,x,y,em);
         fill_histogram(h18154,x,y,es);
      }
   }

   if ( h18005 == NULL || h18171 == NULL )
      return;

   nx = h18005->nbins;
   ny = h18005->nbins_2d;
   for (iy=0; iy<ny; iy++)
   {
      double y = h18005->specific_2d.real.lower_limit + 
             (iy+0.5) * (h18005->specific_2d.real.upper_limit
                - h18005->specific_2d.real.lower_limit) /
                 h18005->nbins_2d;
      for (ix = 0; ix <nx; ix++ )
      {
         double x = h18005->specific.real.lower_limit + 
                (ix+0.5) * (h18005->specific.real.upper_limit
                   - h18005->specific.real.lower_limit) /
                    h18005->nbins;
         double nwol = h18105->extension->ddata[ix+nx*iy];
         double swol = h18106->extension->ddata[ix+nx*iy];
         double rwol = h18171->extension->ddata[ix+nx*iy];
         double rwol2= h18172->extension->ddata[ix+nx*iy];
         double dwol = h18181->extension->ddata[ix+nx*iy];
         double dwol2= h18182->extension->ddata[ix+nx*iy];
         double rm = 0., rs = 0., dm= 0., ds = 0.;

         if ( nwol >= 1.999 && swol > 0. )
         {
            rm = rwol / swol;
            if ( rwol2/swol > rm*rm )
               rs = sqrt((rwol2/swol-rm*rm) * nwol/(nwol-1.));
            dm = dwol / swol;
            if ( dwol2/swol > dm*dm )
               ds = sqrt((dwol2/swol-dm*dm) * nwol/(nwol-1.));
         }

         fill_histogram(h18173,x,y,rm);
         fill_histogram(h18174,x,y,rs);
         fill_histogram(h18183,x,y,dm);
         fill_histogram(h18184,x,y,ds);
      }
   }
}

void fill_ebias_correction(void);

void fill_ebias_correction()
{
   HISTOGRAM *h1, *h2;
   int nbins = 0, nbins2 = 0, i, j;
   double xlow, xhigh, ylow, yhigh;
   if ( (h1 = get_histogram_by_ident(19114)) != NULL )
   {
      if ( (h2 = get_histogram_by_ident(19119)) != NULL )
         free_histogram(h2);
      xlow = h1->specific.real.lower_limit;
      xhigh = h1->specific.real.upper_limit;
      ylow = h1->specific_2d.real.lower_limit;
      yhigh = h1->specific_2d.real.upper_limit;
      nbins = h1->nbins;
      nbins2 = h1->nbins_2d;
      h2 = book_histogram(19119,"Energy bias as lg(E_rec0/E_true)", "D", 1,
            &xlow, &xhigh, &nbins);
      if ( h2 != NULL )
      {
         for ( i=0; i<nbins; i++ )
         {
            double s = 0., all = 0., d, c = 0.;
            for ( j=0; j<nbins2; j++ )
               all += h1->extension->ddata[i+j*nbins];
            for ( j=0; j<nbins2; j++ )
            {
               d = h1->extension->ddata[i+j*nbins];
               if ( s+d > 0.5*all ) /* Median */
               {
                  double f = (0.5*all-s)/d;
                  c = (j+f)/(double)nbins2 * (yhigh - ylow) + ylow;
                  break;
               }
               s += d;
            }
            fill_histogram(h2,xlow+(i+0.5)/(double)nbins*(xhigh-xlow),0.,c);
         }
      }
   }
}

void syntax (char * prgm);

void syntax (char * prgm)
{
   fprintf(stderr,"Syntax: %s [ input_file_name [ output_file_name ]\n", prgm);
   fprintf(stderr,"Default input file name: user_histograms.hdata.gz\n");
   fprintf(stderr,"Default output file name: lookups.hdata.gz\n");
   exit(1);
}

int main(int argc, char **argv)
{
   double xylow[2], xyhigh[2];
   int nbins[2];
   HISTOGRAM *histo, *next_histo, *h;
   const char *histo_file = "user_histograms.hdata.gz";
   const char *lookup_file = "lookups.hdata.gz";
   int tel_type = 0, res_type = 0;
   char *prgm = argv[0];
   int uvers = 0, ntypes = 1;

   if ( argc > 1 )
   {
      if ( strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0 ||
           argv[1][0] == '-'  )
         syntax(prgm);
      histo_file = argv[1];
      if ( argc > 2 )
         lookup_file = argv[2];
   }

   if ( read_histogram_file(histo_file,0) == -1 )
      exit(1);

   sort_histograms();
   if ( get_first_histogram() == 0 )
   {
      fprintf(stderr,"No histograms available.\n");
      exit(1);
   }
   
   if ( (h = get_histogram_by_ident(12099)) == NULL )
   {
      fprintf(stderr, "Cannot identify user parameters.\n");
      exit(1);
   }
   if ( h->extension != NULL && h->extension->ddata != NULL )
   {
      uvers = (int)(h->extension->ddata[0]+0.01);
      if ( uvers < 200 )
         ntypes = 5;
      else
         ntypes = (int)(h->extension->ddata[1]+0.01);
   }
   if ( ntypes < 1 || ntypes > 20 )
   {
      fprintf(stderr,"Unreasonable number of telescope types.\n");
      exit(1);
   }

   for ( tel_type=0; tel_type < ntypes; tel_type++ )
   {
      int ht = 100000*tel_type;
      if ( get_histogram_by_ident(18000+100000*tel_type) == 0 )
         continue;

      h18000 = get_histogram_by_ident(ht+18000);
      h18001 = get_histogram_by_ident(ht+18001);
      h18011 = get_histogram_by_ident(ht+18011);
      h18012 = get_histogram_by_ident(ht+18012);
      h18021 = get_histogram_by_ident(ht+18021);
      h18022 = get_histogram_by_ident(ht+18022);
      h18051 = get_histogram_by_ident(ht+18051);
      h18052 = get_histogram_by_ident(ht+18052);
      
      if ( h18000 == NULL || h18001 == NULL || 
           h18011 == NULL || h18012 == NULL ||
           h18021 == NULL || h18022 == NULL ||
	   h18051 == NULL || h18052 == NULL )
      {
         fprintf(stderr,"Not the right histograms available.\n");
         exit(1);
      }

      xylow[0]  = h18000->specific.real.lower_limit;
      xylow[1]  = h18000->specific_2d.real.lower_limit;
      xyhigh[0] = h18000->specific.real.upper_limit;
      xyhigh[1] = h18000->specific_2d.real.upper_limit;
      nbins[0]  = h18000->nbins;
      nbins[1]  = h18000->nbins_2d;
      if ( (h = get_histogram_by_ident(ht+18100)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18101)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18111)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18112)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18113)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18114)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18121)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18122)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18123)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18124)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18140)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18141)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18151)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18152)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18153)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18154)) != NULL )
         free_histogram(h);
      h18100 = book_histogram(ht+18100, h18000->title, "D", 2, xylow, xyhigh, nbins);
      h18101 = book_histogram(ht+18101, h18001->title, "D", 2, xylow, xyhigh, nbins);
      h18111 = book_histogram(ht+18111, h18011->title, "D", 2, xylow, xyhigh, nbins);
      h18112 = book_histogram(ht+18112, h18012->title, "D", 2, xylow, xyhigh, nbins);
      h18113 = book_histogram(ht+18113, 
            "Log10(image ampl.) versus core distance: mean image width", 
            "D", 2, xylow, xyhigh, nbins);
      h18114 = book_histogram(ht+18114, 
            "Log10(image ampl.) versus core distance: image width r.m.s.", 
            "D", 2, xylow, xyhigh, nbins);
      h18121 = book_histogram(ht+18121, h18021->title, "D", 2, xylow, xyhigh, nbins);
      h18122 = book_histogram(ht+18122, h18022->title, "D", 2, xylow, xyhigh, nbins);
      h18123 = book_histogram(ht+18123,  
            "Log10(image ampl.) versus core distance: mean image length", 
            "D", 2, xylow, xyhigh, nbins);
      h18124 = book_histogram(ht+18124,  
            "Log10(image ampl.) versus core distance: image length r.m.s.", 
            "D", 2, xylow, xyhigh, nbins);
      h18140 = book_histogram(ht+18140, h18000->title, "D", 2, xylow, xyhigh, nbins);
      h18141 = book_histogram(ht+18141, h18001->title, "D", 2, xylow, xyhigh, nbins);
      h18151 = book_histogram(ht+18151, h18051->title, "D", 2, xylow, xyhigh, nbins);
      h18152 = book_histogram(ht+18152, h18052->title, "D", 2, xylow, xyhigh, nbins);
      h18153 = book_histogram(ht+18153,  
            "Log10(image ampl.) versus core distance: mean amp/E", 
            "D", 2, xylow, xyhigh, nbins);
      h18154 = book_histogram(ht+18154,  
            "Log10(image ampl.) versus core distance: amp/E r.m.s.", 
            "D", 2, xylow, xyhigh, nbins);

      h18005 = get_histogram_by_ident(ht+18005);
      h18006 = get_histogram_by_ident(ht+18006);
      h18071 = get_histogram_by_ident(ht+18071);
      h18072 = get_histogram_by_ident(ht+18072);
      h18081 = get_histogram_by_ident(ht+18081);
      h18082 = get_histogram_by_ident(ht+18082);
      if ( (h = get_histogram_by_ident(ht+18105)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18106)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18171)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18172)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18173)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18174)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18181)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18182)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18183)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18184)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18331)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18332)) != NULL )
         free_histogram(h);
      if ( (h = get_histogram_by_ident(ht+18333)) != NULL )
         free_histogram(h);

      if ( h18005 != NULL && h18006 != NULL &&
           h18071 != NULL && h18071 != NULL &&
           h18081 != NULL && h18081 != NULL )
      {
         xylow[0]  = h18005->specific.real.lower_limit;
         xylow[1]  = h18005->specific_2d.real.lower_limit;
         xyhigh[0] = h18005->specific.real.upper_limit;
         xyhigh[1] = h18005->specific_2d.real.upper_limit;
         nbins[0]  = h18005->nbins;
         nbins[1]  = h18005->nbins_2d;

         h18105 = book_histogram(ht+18105, h18005->title, "D", 2, xylow, xyhigh, nbins);
         h18106 = book_histogram(ht+18106, h18006->title, "D", 2, xylow, xyhigh, nbins);
         h18171 = book_histogram(ht+18171, h18071->title, "D", 2, xylow, xyhigh, nbins);
         h18172 = book_histogram(ht+18172, h18072->title, "D", 2, xylow, xyhigh, nbins);
         h18173 = book_histogram(ht+18173,
               "Log10(image ampl.) versus width/length: mean core distance", 
               "D", 2, xylow, xyhigh, nbins);
         h18174 = book_histogram(ht+18174,
               "Log10(image ampl.) versus width/length: core distance r.m.s.", 
               "D", 2, xylow, xyhigh, nbins);
         h18181 = book_histogram(ht+18181, h18081->title, "D", 2, xylow, xyhigh, nbins);
         h18182 = book_histogram(ht+18182, h18082->title, "D", 2, xylow, xyhigh, nbins);
         h18183 = book_histogram(ht+18183,
               "Log10(image ampl.) versus width/length: mean image distance", 
               "D", 2, xylow, xyhigh, nbins);
         h18184 = book_histogram(ht+18184,
               "Log10(image ampl.) versus width/length: image distance r.m.s.", 
               "D", 2, xylow, xyhigh, nbins);
      }

      h18301 = get_histogram_by_ident(ht+18301);
      h18311 = get_histogram_by_ident(ht+18311);
      h18321 = get_histogram_by_ident(ht+18321);
      h18322 = get_histogram_by_ident(ht+18322);

      if ( h18301 != NULL && h18311 != NULL )
      {
         HISTOGRAM *h18331, *h18332, *h18333;
         int i, n;
         xylow[0]  = h18301->specific.real.lower_limit;
         xylow[1]  = h18301->specific_2d.real.lower_limit;
         xyhigh[0] = h18301->specific.real.upper_limit;
         xyhigh[1] = h18301->specific_2d.real.upper_limit;
         nbins[0]  = h18301->nbins;
         nbins[1]  = h18301->nbins_2d;
         n = nbins[0]*nbins[1];

         h18331 = book_histogram(ht+18331, 
            "Log10(image ampl.) versus core distance: fraction of triggered telescopes", 
            "D", 2, xylow, xyhigh, nbins);
         if ( h18331 != NULL )
         {
            for (i=0; i<n; i++ )
               h18331->extension->ddata[i] = ((h18301->extension->ddata[i]==0.) ?
                  -1. : h18311->extension->ddata[i] / h18301->extension->ddata[i]);
         }
         if ( h18321 != NULL )
         {
            h18332 = book_histogram(ht+18332, 
               "Log10(image ampl.) versus core distance: fraction of telescopes with basically usable images", 
               "D", 2, xylow, xyhigh, nbins);
            if ( h18332 != NULL )
            {
               for (i=0; i<n; i++ )
                  h18332->extension->ddata[i] = ((h18301->extension->ddata[i]==0.) ?
                     -1. : h18321->extension->ddata[i] / h18301->extension->ddata[i]);
            }
         }
         if ( h18322 != NULL )
         {
            h18333 = book_histogram(ht+18333, 
               "Log10(image ampl.) versus core distance: fraction of telescopes with fully usable images", 
               "D", 2, xylow, xyhigh, nbins);
            if ( h18333 != NULL )
            {
               for (i=0; i<n; i++ )
                  h18333->extension->ddata[i] = ((h18301->extension->ddata[i]==0.) ?
                     -1. : h18322->extension->ddata[i] / h18301->extension->ddata[i]);
            }
         }
      }

      fill_gaps();

      gen_image_lookups();
   }

   /* Angular resolution versus telescope multiplicity for different cuts. */
   for ( res_type=1; res_type<=7; res_type++ ) 
      if ( (h = get_histogram_by_ident(19001+100*res_type)) != NULL )
      {
         HISTOGRAM *hr68, *hr80;
         int i, j;
         xylow[0]  = h->specific.real.lower_limit - 0.5;
         xyhigh[0] = h->specific.real.upper_limit - 0.5;
         xyhigh[1] = h->specific_2d.real.upper_limit;
         nbins[0]  = h->nbins;
         nbins[1]  = h->nbins_2d;
         if ( (hr68 = get_histogram_by_ident(18900+res_type)) != 0 )
            free_histogram(hr68);
         hr68 = book_histogram(18900+res_type,
            res_type == 1 ?"Angular resolution 68% containment vs. Ntel (shape cuts)" :
            res_type == 2 ?"Angular resolution 68% containment vs. Ntel (shape+dE cuts)" :
            res_type == 3 ?"Angular resolution 68% containment vs. Ntel (shape+dE+dE2+hmax cuts)" :
            res_type == 4 ?"Angular resolution 68% containment vs. Ntel (shape+hmax cuts)" :
            res_type == 5 ?"Angular resolution 68% containment vs. Ntel (shape+dE+hmax cuts)" :
            res_type == 6 ?"Angular resolution 68% containment vs. Ntel (shape+dE2+hmax cuts)" :
            res_type == 7 ?"Angular resolution 68% containment vs. Ntel (shape+dE+dE2 cuts)" :
            "Angular resolution 68% containment vs. Ntel", 
            "D", 1, xylow, xyhigh, nbins);
         if ( (hr80 = get_histogram_by_ident(18910+res_type)) != 0 )
            free_histogram(hr80);
         hr80 = book_histogram(18910+res_type,
            res_type == 1 ?"Angular resolution 80% containment vs. Ntel (shape cuts)" :
            res_type == 2 ?"Angular resolution 80% containment vs. Ntel (shape+dE cuts)" :
            res_type == 3 ?"Angular resolution 80% containment vs. Ntel (shape+dE+dE2+hmax cuts)" :
            res_type == 4 ?"Angular resolution 80% containment vs. Ntel (shape+hmax cuts)" :
            res_type == 5 ?"Angular resolution 80% containment vs. Ntel (shape+dE+hmax cuts)" :
            res_type == 6 ?"Angular resolution 80% containment vs. Ntel (shape+dE2+hmax cuts)" :
            res_type == 7 ?"Angular resolution 80% containment vs. Ntel (shape+dE+dE2 cuts)" :
            "Angular resolution 80% containment vs. Ntel", 
            "D", 1, xylow, xyhigh, nbins);
         if ( hr68 == NULL || hr80 == NULL )
            break;
         for ( i=1; i<nbins[0]; i++ )
         {
            double s100 = 0., s80, s68, s, e;
            for ( j=0; j<h->nbins_2d; j++ )
               s100 += h->extension->ddata[i+j*h->nbins];
            s68 = 0.68*s100;
            s80 = 0.80*s100;
            for ( j=0, s=0.; j<h->nbins_2d; j++ )
            {
               s += (e = h->extension->ddata[i+j*h->nbins]);
               if ( s >= s68 )
               {
                  double bf = (e<=0. ? 0. : (s68-(s-e)) / e);
                  fill_histogram(hr68, (i+0.1)*xyhigh[0]/nbins[0], 0.,
                     (j+bf)*xyhigh[1]/nbins[1]);
                  break;
               }
            }
            for ( j=0, s=0.; j<h->nbins_2d; j++ )
            {
               s += (e = h->extension->ddata[i+j*h->nbins]);
               if ( s >= s80 )
               {
                  double bf = (e<=0. ? 0. : (s80-(s-e)) / e);
                  fill_histogram(hr80, (i+0.1)*xyhigh[0]/nbins[0], 0.,
                     (j+bf)*xyhigh[1]/nbins[1]);
                  break;
               }
            }
         }
      }

   /* Angular resolution versus log energy for different cuts. */
   for ( res_type=1; res_type<=7; res_type++ ) 
      if ( (h = get_histogram_by_ident(19002+100*res_type)) != NULL )
      {
         HISTOGRAM *hr68, *hr80;
         int i, j;
         xylow[0]  = h->specific.real.lower_limit;
         xyhigh[0] = h->specific.real.upper_limit;
         xyhigh[1] = h->specific_2d.real.upper_limit;
         nbins[0]  = h->nbins;
         nbins[1]  = h->nbins_2d;
         if ( (hr68 = get_histogram_by_ident(18920+res_type)) != 0 )
            free_histogram(hr68);
         hr68 = book_histogram(18920+res_type,
            res_type == 1 ?"Angular resolution 68% containment vs. lg E (shape cuts)" :
            res_type == 2 ?"Angular resolution 68% containment vs. lg E (shape+dE cuts)" :
            res_type == 3 ?"Angular resolution 68% containment vs. lg E (shape+dE+dE2+hmax cuts)" :
            res_type == 4 ?"Angular resolution 68% containment vs. lg E (shape+hmax cuts)" :
            res_type == 5 ?"Angular resolution 68% containment vs. lg E (shape+dE+hmax cuts)" :
            res_type == 6 ?"Angular resolution 68% containment vs. lg E (shape+dE2+hmax cuts)" :
            res_type == 7 ?"Angular resolution 68% containment vs. lg E (shape+dE+dE2 cuts)" :
            "Angular resolution 68% containment vs. lg E", 
            "D", 1, xylow, xyhigh, nbins);
         if ( (hr80 = get_histogram_by_ident(18930+res_type)) != 0 )
            free_histogram(hr80);
         hr80 = book_histogram(18930+res_type,
            res_type == 1 ?"Angular resolution 80% containment vs. lg E (shape cuts)" :
            res_type == 2 ?"Angular resolution 80% containment vs. lg E (shape+dE cuts)" :
            res_type == 3 ?"Angular resolution 80% containment vs. lg E (shape+dE+dE2+hmax cuts)" :
            res_type == 4 ?"Angular resolution 80% containment vs. lg E (shape+hmax cuts)" :
            res_type == 5 ?"Angular resolution 80% containment vs. lg E (shape+dE+hmax cuts)" :
            res_type == 6 ?"Angular resolution 80% containment vs. lg E (shape+dE2+hmax cuts)" :
            res_type == 7 ?"Angular resolution 80% containment vs. lg E (shape+dE+dE2 cuts)" :
            "Angular resolution 80% containment vs. lg E", 
            "D", 1, xylow, xyhigh, nbins);
         if ( hr68 == NULL || hr80 == NULL )
            break;
         for ( i=1; i<nbins[0]; i++ )
         {
            double s100 = 0., s80, s68, s, e;
            for ( j=0; j<h->nbins_2d; j++ )
               s100 += h->extension->ddata[i+j*h->nbins];
            s68 = 0.68*s100;
            s80 = 0.80*s100;
            for ( j=0, s=0.; j<h->nbins_2d; j++ )
            {
               s += (e = h->extension->ddata[i+j*h->nbins]);
               if ( s >= s68 )
               {
                  double bf = (e<=0. ? 0. : (s68-(s-e)) / e);
                  fill_histogram(hr68, xylow[0]+(i+0.1)*(xyhigh[0]-xylow[0])/nbins[0], 0.,
                     (j+bf)*xyhigh[1]/nbins[1]);
                  break;
               }
            }
            for ( j=0, s=0.; j<h->nbins_2d; j++ )
            {
               s += (e = h->extension->ddata[i+j*h->nbins]);
               if ( s >= s80 )
               {
                  double bf = (e<=0. ? 0. : (s80-(s-e)) / e);
                  fill_histogram(hr80, xylow[0]+(i+0.1)*(xyhigh[0]-xylow[0])/nbins[0], 0.,
                     (j+bf)*xyhigh[1]/nbins[1]);
                  break;
               }
            }
         }
      }

   fill_ebias_correction();

   for (histo=get_first_histogram(); histo!=NULL; histo=next_histo)
   {
      next_histo = histo->next;
      if ( !(histo->ident%100000 >= 18100 && histo->ident%100000 <= 18184) &&
           !(histo->ident >= 18901 && histo->ident <= 18939) &&
           !(histo->ident%100000 >= 18330 && histo->ident%100000 <= 18339) &&
           histo->ident != 12099 && histo->ident != 19119 )
         free_histogram(histo);
   }

   write_all_histograms(lookup_file);

   fprintf(stderr, 
      "Lookups for %d telescope type%s generated and written to %s\n", 
      tel_type, (tel_type==1)?"":"s", lookup_file);

   return 0;
}
