/* ============================================================================

Copyright (C) 2013  Konrad Bernloehr

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

/** @file add_histograms.c
 *  @short Utility program for adding up matching histograms.
 *
@verbatim
Syntax:  add_histograms [ -x id1,...] input_files ... -o output_file
@endverbatim
 *  The histograms may be within multiple I/O blocks of the
 *  input file. Matching histograms will be added up, unless
 *  set to be excluded with the '-x' option.
 *  Only non-empty histograms are written to output.
 *
 *  @author Konrad Bernloehr
 *  @date    @verbatim CVS $Date: 2014/06/24 14:29:40 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.2 $ @endverbatim
 */

/** @defgroup add_histograms_c The add_histograms program */
/** @{ */

#include "initial.h"
#include "histogram.h"
#include "io_basic.h"
#include "io_histogram.h"
#include "fileopen.h"
#include "straux.h"

void syntax(const char *prgm);

void syntax(const char *prgm)
{
   fprintf(stderr,"Syntax: %s [-V | -VV ] [ -x id1,... ] input_files ... -o output_file\n",
            prgm);
   fprintf(stderr,"The '-V'/'-VV' results in more verbose screen output.\n");
   fprintf(stderr,"The '-x' option excludes histograms from being added up.\n");
   exit(1);
}

/** Main program. */

int main (int argc, char **argv)
{
   char *ofname = NULL, *prgm = argv[0];
   int verbose = 0;
   int iarg;
   long *xcld_ids = NULL;
   int nxcld = 0;

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
      else if ( strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0 )
      {
         syntax(prgm);
      }
      else if ( strcmp(argv[1],"-x") == 0 && argc >= 3 )
      {
         int ipos = 0;
         int xi = 0;
         char word[100];
         while ( getword(argv[2],&ipos,word,99,',','%') > 0 )
         {
            if ( ( xi = atol(word)) > 0 )
            {
               size_t m = (size_t) (nxcld+1) * sizeof(int);
               if ( xcld_ids == NULL )
                  xcld_ids = (long *) malloc(m);
               else
                  xcld_ids = (long *) realloc(xcld_ids,m);
               if ( xcld_ids == NULL )
               {
                  fprintf(stderr,"(Re-)Allocation failed.\n");
                  exit(1);
               }
               xcld_ids[nxcld] = xi;
               nxcld++;
            }
         }
         if ( nxcld > 0 )
         {
            int ix=0;
            fprintf(stderr,"Excluding histogram ID(s): %ld", xcld_ids[ix]);
            for (ix=1; ix<nxcld; ix++)
               fprintf(stderr,", %ld", xcld_ids[ix]);
            fprintf(stderr," (to be replaced instead of added-up)\n");
         }
         argc -= 2;
         argv += 2;
      }
      else if ( argv[1][0] == '-' )
      {
         fprintf(stderr,"Unknown option %s\n", argv[1]);
         syntax(prgm);
      }
      else
         break;
   }

   if ( argc < 4 )
      syntax(prgm);
   if ( strcmp(argv[argc-2],"-o") != 0 )
      syntax(prgm);
   ofname = argv[argc-1];
   argc -= 2;

   for ( iarg=1; iarg<argc; iarg++ )
   {
      if ( argv[iarg][0] == '-' )
         syntax(prgm);
      if ( read_histogram_file_x(argv[iarg],1,xcld_ids,nxcld) < 0 )
         break;
   }

   sort_histograms();
   if ( get_first_histogram() == 0 )
   {
      fprintf(stderr,"No histograms available for conversion.\n");
      exit(1);
   }

   if ( verbose )
   {
      HISTOGRAM *h = get_first_histogram();
      while ( h != NULL )
      {
         printf("Histogram of type %c, ID=%ld, title=\"%s\" is %cD: ",
            h->type, h->ident, (h->title!=NULL)?h->title:"(none)",
            (h->nbins_2d > 0) ? '2' : '1');
         if ( h->nbins_2d > 0)
            printf("%d * %d bins",h->nbins,h->nbins_2d);
         else
            printf("%d bins",h->nbins);
         printf(", %ld entries.\n",h->entries);
         h = h->next;
      }

      if ( verbose >= 2 )
         display_all_histograms(); 
   }
   
   write_all_histograms(ofname);
   return(0);
}

/** @} */
