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
 *  @date    @verbatim CVS $Date: 2013/10/21 12:53:31 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.2 $ @endverbatim
 */

/** @defgroup list_histograms_c The list_histogram program */
/** @{ */

#include "initial.h"
#include "histogram.h"
#include "io_basic.h"
#include "io_histogram.h"
#include "fileopen.h"

/** Main program. */

int main (int argc, char **argv)
{
   char *prgm = argv[0];
   int verbose = 0;
   int list_flag = 1;
   int show_id = 0;
   
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
      }
      else
         break;
   }

   if ( argc >= 2 )
   {
      if ( strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0 )
      {
         fprintf(stderr,"Syntax: %s [-V | -VV ] [ -h ID ] [ input_file ... ]\n",
                  prgm);
         fprintf(stderr,"The '-V'/'-VV' results in more verbose screen output.\n");
         fprintf(stderr,"The '-h' option shows a single histogram in detail.\n");
         fprintf(stderr,"The '-H' option prints contents of a single histogram.\n");
         exit(0);
      }
   }

   {
      int iarg;
      for (iarg=1; iarg<argc; iarg++)
      {
         printf("\nFile %s:\n", argv[iarg]);
         read_histogram_file(argv[iarg],1+16*list_flag);
      }
   }

   if ( verbose )
   {
      sort_histograms();
      HISTOGRAM *h = get_first_histogram();
      while ( h != NULL )
      {
         if ( show_id == 0 || show_id == h->ident )
         {
            printf("Histogram of type %c, ID=%ld, title=\"%s\" is %cD: ",
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
