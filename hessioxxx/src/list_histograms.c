/* ============================================================================

Copyright (C) 2001, 2005, 2009, 2013  Konrad Bernloehr

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

/** @file list_histograms.c
 *  @short Utility program for listing histograms.
@verbatim
Syntax:  list_histograms [ input_file ... ]
@endverbatim
 *  The default input file name is 'testpattern.hdata'.
 *  The histograms may be within multiple I/O blocks of the
 *  input file.
 *
 *  @author Konrad Bernloehr
 *  @date    @verbatim CVS $Date: 2018/09/07 17:43:51 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.4 $ @endverbatim
 */

/** @defgroup list_histograms_c The list_histogram program */
/** @{ */

#include "initial.h"
#include "histogram.h"
#include "io_basic.h"
#include "io_histogram.h"
#include "fileopen.h"


/* --------------------- print_ratio -------------------- */
/**
 *  @short Print ratio of two histograms: fact * histo1 / histo2
 *
 *  @param  histo1 Pointer to histogram 1
 *  @param  histo2 Pointer to histogram 2
 *  @param  fact   Scaling factor
 *
 *  @return  (none)
 */

void print_ratio (HISTOGRAM *histo1, HISTOGRAM *histo2, double fact);

void print_ratio (HISTOGRAM *histo1, HISTOGRAM *histo2, double fact)
{
   double xmin, xmax, xmin2, xmax2;
   double ymin, ymax, ymin2, ymax2;
   double dx, dy;
   int nx=0, ny=1, iy, ix;

   if ( histo1 == (HISTOGRAM *) NULL || histo2 == NULL )
      return; /* Histogram has not been allocated */
   if ( histo1->nbins <= 0 || histo2->nbins != histo1->nbins )
   {
      fprintf(stderr,"Histograms not matching.\n");
      return;
   }
   if ( histo1->nbins_2d != histo2->nbins_2d )
   {
      fprintf(stderr,"Histograms not matching.\n");
      return;
   }
   if ( histo2->type != histo1->type )
   {
      fprintf(stderr,"Histograms not matching.\n");
      return;
   }
   if ( histo1->type == 'F' || histo1->type == 'D' )
   {
      if (histo1->extension == (struct Histogram_Extension *) NULL ||
          histo2->extension == (struct Histogram_Extension *) NULL )
         return;
      if ( (histo1->type=='F' && (histo1->extension)->fdata==(float*)NULL) ||
           (histo1->type=='D' && (histo1->extension)->ddata==(double*)NULL) ||
           (histo2->type=='F' && (histo2->extension)->fdata==(float*)NULL) ||
           (histo2->type=='D' && (histo2->extension)->ddata==(double*)NULL) )
         return;
   }
   else if ( histo1->counts == (unsigned long *) NULL || histo2->counts == (unsigned long *) NULL )
      return; /* This can't happen for dynamically allocated histograms */
   if ( histo1->tentries == 0 && histo2->tentries == 0 )
   {
      fprintf(stderr,"Neither histogram has any entries.\n");
      return;  /* No entries in histogram range */
   }
   if ( histo1->type != 'I' && histo1->type != 'R' &&
        histo1->type != 'i' && histo1->type != 'r' &&
        histo1->type != 'F' && histo1->type != 'D' )
   {
      fprintf(stderr,"Unknown histogram type.\n");
      return; /* Histogram has neither integer nor real type */
   }

   if ( histo1->type == 'I' || histo1->type == 'i' )
   {
      xmin = (double) histo1->specific.integer.lower_limit;
      xmax = (double) histo1->specific.integer.upper_limit;
      xmin2 = (double) histo2->specific.integer.lower_limit;
      xmax2 = (double) histo2->specific.integer.upper_limit;
   }
   else
   {
      xmin = (double) histo1->specific.real.lower_limit;
      xmax = (double) histo1->specific.real.upper_limit;
      xmin2 = (double) histo2->specific.real.lower_limit;
      xmax2 = (double) histo2->specific.real.upper_limit;
   }
   if ( xmin2 != xmin || xmax2 != xmax )
   {
      fprintf(stderr,"Histograms not matching\n");
      return;
   }
   nx = histo1->nbins;
   dx = 0.;
   if ( nx > 0 )
      dx = (xmax-xmin)/nx;

// _HLOCK_
   printf("\n\n");
   printf("# Ratio of two histograms\n");
   printf("# ID1=%ld: ",histo1->ident);
   if ( histo1->title != (char *) NULL )
   {
      printf("%s",histo1->title);
   }
   printf("\n");
   printf("# ID2=%ld: ",histo2->ident);
   if ( histo2->title != (char *) NULL )
   {
      printf("%s",histo2->title);
   }
   printf("\n");
   if ( fact != 1. )
      printf("# Scale factor = %g\n", fact);

   if ( (ny = histo1->nbins_2d) > 0 )
   {
      if ( histo1->type == 'I' || histo1->type == 'i' )
      {
         ymin = (double) histo1->specific_2d.integer.lower_limit;
         ymax = (double) histo1->specific_2d.integer.upper_limit;
         ymin2 = (double) histo2->specific_2d.integer.lower_limit;
         ymax2 = (double) histo2->specific_2d.integer.upper_limit;
      }
      else
      {
         ymin = (double) histo1->specific_2d.real.lower_limit;
         ymax = (double) histo1->specific_2d.real.upper_limit;
         ymin2 = (double) histo2->specific_2d.real.lower_limit;
         ymax2 = (double) histo2->specific_2d.real.upper_limit;
      }
      if ( ymin2 != ymin || ymax2 != ymax )
         return;
      dy = 0.;
      if ( ny > 0 )
         dy = (ymax-ymin)/ny;
      printf("# X = %g to %g, Y = %g to %g, E=%lu/%lu, U=%lu/%lu, O=%lu/%lu\n",
         xmin,xmax,ymin,ymax,
         histo1->entries, histo2->entries,
         histo1->underflow, histo2->underflow,
         histo1->overflow, histo2->overflow);
      if ( histo1->type == 'F' || histo1->type == 'D' )
      {
         printf("# Content: %12.7g/%12.7g total (%12.7g/%12.7g inside range)\n",
            (histo1->extension)->content_all, (histo2->extension)->content_all, 
            (histo1->extension)->content_inside, (histo2->extension)->content_inside);
      }
      printf("\n");
      for ( iy = 0; iy < ny; iy++ )
      {
         printf("\n");
         for ( ix = 0; ix < nx; ix++ )
         {
            double r = 0.;
            if ( histo1->type != 'F' && histo1->type != 'D' )
            {
               if ( histo2->counts[iy*nx+ix] == 0 )
                   printf("   %d %d * %ld %ld\n", ix, iy, 
                     histo1->counts[iy*nx+ix], histo2->counts[iy*nx+ix]);
               else
               {
                  r = fact * (double) histo1->counts[iy*nx+ix] / 
                             (double) histo2->counts[iy*nx+ix];
                  printf("   %d %d %g %ld %ld\n", ix, iy, r, 
                     histo1->counts[iy*nx+ix], histo2->counts[iy*nx+ix]);
               }
            }
            else
            {
               double z1 = 0., z2 = 0.;
               if ( histo1->type == 'F' )
               {
                  z1 = (histo1->extension)->fdata[iy*nx+ix];
                  z2 = (histo2->extension)->fdata[iy*nx+ix];
               }
               else
               {
                  z1 = (histo1->extension)->ddata[iy*nx+ix];
                  z2 = (histo2->extension)->ddata[iy*nx+ix];
               }
               if ( z2 == 0. )
                  printf("   %g\t%g\t*\t%g\t%g\n", 
                     xmin+dx*(ix+0.5), ymin+dy*(iy+0.5), z1, z2);
               else
                  printf("   %g\t%g\t%g\t%g\t%g\n", 
                     xmin+dx*(ix+0.5), ymin+dy*(iy+0.5), fact*z1/z2, z1, z2);
            }
         }
      }
   }
   else
   {
      printf("# X = %g to %g, E=%lu/%lu, U=%lu/%lu, O=%lu/%lu\n",
         xmin,xmax,
         histo1->entries, histo2->entries,
         histo1->underflow, histo2->underflow,
         histo1->overflow, histo2->overflow);
      if ( histo1->type == 'F' || histo1->type == 'D' )
      {
         printf("# Content: %12.7g/%12.7g total (%12.7g/%12.7g inside range)\n",
            (histo1->extension)->content_all,  (histo2->extension)->content_all,
            (histo1->extension)->content_inside, (histo2->extension)->content_inside);
      }
      printf("\n");
      for ( ix = 0; ix < nx; ix++ )
      {
         double r = 0.;
         if ( histo1->type != 'F' && histo1->type != 'D' )
         {
            if ( histo2->counts[ix] == 0 )
                printf("   %d * %ld %ld\n", ix, 
                  histo1->counts[ix], histo2->counts[ix]);
            else
            {
               r = fact * (double) histo1->counts[ix] / 
                          (double) histo2->counts[ix];
               printf("   %d %g %ld %ld\n", ix, r, 
                  histo1->counts[ix], histo2->counts[ix]);
            }
         }
         else
         {
            double z1 = 0., z2 = 0.;
            if ( histo1->type == 'F' )
            {
               z1 = (histo1->extension)->fdata[ix];
               z2 = (histo2->extension)->fdata[ix];
            }
            else
            {
               z1 = (histo1->extension)->ddata[ix];
               z2 = (histo2->extension)->ddata[ix];
            }
            if ( z2 == 0. )
               printf("   %g\t*\t%g\t%g\n", xmin+dx*(ix+0.5), z1, z2);
            else
               printf("   %g\t%g\t%g\t%g\n", xmin+dx*(ix+0.5), fact*z1/z2, z1, z2);
         }
      }
   }

// _HUNLOCK_
}

/** Main program. */

int main (int argc, char **argv)
{
   char *prgm = argv[0];
   int verbose = 0;
   int list_flag = 1;
   int show_id = 0, show_id2 = 0;
   int operation = 0;
   double scale_fact = 1.0;
   const char *ctxt = "";
   
   while ( argc >= 2 )
   {
      if ( strcmp(argv[1],"-V") == 0 || strcmp(argv[1],"--verbose") == 0 )
      {
         argc--;
         argv++;
         verbose++;
      }
      else if ( strcmp(argv[1],"-VV") == 0 )
      {
         argc--;
         argv++;
         verbose+=2;
      }
      else if ( strcmp(argv[1],"-q") == 0 )
      {
         argc--;
         argv++;
         list_flag = 0;
      }
      else if ( strcmp(argv[1],"-h") == 0 && argc > 2 )
      {
         show_id = atoi(argv[2]);
         argc-=2;
         argv+=2;
         list_flag = 0;
         verbose += 2;
      }
      else if ( strcmp(argv[1],"-H") == 0 && argc > 2 )
      {
         show_id = atoi(argv[2]);
         argc-=2;
         argv+=2;
         list_flag = 0;
         verbose += 3;
         ctxt = "# ";
      }
      else if ( strcmp(argv[1],"-R") == 0 && argc > 3 )
      {
         show_id = atoi(argv[2]);
         show_id2 = atoi(argv[3]);
         operation = 1; 
         argc-=3;
         argv+=3;
         list_flag = 0;
         verbose += 3;
         ctxt = "# ";
      }
      else if ( strcmp(argv[1],"--scale-factor") == 0 && argc > 2 )
      {
         scale_fact = atof(argv[2]);
         argc-=2;
         argv+=2;
      }
      else
         break;
   }

   if ( argc >= 2 )
   {
      if ( strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0 )
      {
         fprintf(stderr,"Syntax: %s [-V | -VV ] [ -{h|H} ID ] [ -R ID1 ID2 ] [ input_file ... ]\n",
                  prgm);
         fprintf(stderr,"The '-V'/'-VV' results in more verbose screen output.\n");
         fprintf(stderr,"The '-h' option shows a single histogram in detail.\n");
         fprintf(stderr,"The '-H' option prints contents of a single histogram.\n");
         fprintf(stderr,"The '-R' option prints the ratio of two histogram.\n");
         exit(0);
      }
   }

   {
      int iarg;
      for (iarg=1; iarg<argc; iarg++)
      {
         printf("\n%sFile %s:\n", ctxt, argv[iarg]);
         read_histogram_file(argv[iarg],1+16*list_flag);
      }
   }
   
   if ( operation == 1 && show_id > 0 && show_id2 > 0 )
   {
      sort_histograms();
      HISTOGRAM *h1 = get_histogram_by_ident(show_id);
      HISTOGRAM *h2 = get_histogram_by_ident(show_id2);
      if ( h1 != NULL && h2 != NULL )
         print_ratio(h1,h2,scale_fact);
   }
   else if ( verbose )
   {
      sort_histograms();
      HISTOGRAM *h = get_first_histogram();
      while ( h != NULL )
      {
         if ( show_id == 0 || show_id == h->ident )
         {
            printf("%sHistogram of type %c, ID=%ld, title=\"%s\" is %cD: ", ctxt,
               h->type, h->ident, (h->title!=NULL)?h->title:"(none)",
               (h->nbins_2d > 0) ? '2' : '1');
            if ( h->nbins_2d > 0)
               printf("%d * %d bins",h->nbins,h->nbins_2d);
            else
               printf("%d bins",h->nbins);
            printf(", %ld entries.\n",h->entries);
         }
         h = h->next;
      }

      if ( verbose >= 2 )
      {
         HISTOGRAM *thisto;
         if ( show_id == 0 )
            display_all_histograms();
         else
         {
            if ( (thisto = get_histogram_by_ident(show_id)) != NULL )
            {
               if ( verbose >= 3 )
                  print_histogram(thisto);
               else
                  display_histogram(thisto);
            }
         }
      }
   }

   exit(0);
   return(0);
}

/** @} */
