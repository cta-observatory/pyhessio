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

/** @file check_trgmask.c
 *  @short Check consistency of 'trgmask' files produced with gen_trgmask
 *         for the CTA prod-2 data sets produced in 2013.
 *
@verbatim
Syntax: bin/check_trgmask trgmask-file
@endverbatim
 *
 *  @author Konrad Bernloehr
 */

/** @defgroup check_trgmask_c The check_trgmask program */
/** @{ */


#include "initial.h"  
#include "io_basic.h" 
#include "fileopen.h"
#include "io_trgmask.h"

int main(int argc, char **argv)
{
   struct trgmask_set *tms = (struct trgmask_set *) calloc(1,sizeof(struct trgmask_set));
   struct trgmask_entry *tme = NULL;
   char *in_fname = NULL;
   int n_lst = 0, n_mst = 0, n_dcsst = 0, n_scsst = 0;
   long event = 0;
   size_t i;

   if ( argc < 2 )
      exit(1);
   in_fname = argv[1];

   IO_BUFFER *iobuf = allocate_io_buffer(1000000L);
   if ( iobuf == NULL )
   {
      Error("Cannot allocate I/O buffer");
      exit(1);
   }
   iobuf->max_length = 200000000L;

   IO_ITEM_HEADER item_header;
   iobuf->input_file = fileopen(in_fname,"r");
   if ( iobuf->input_file != NULL )
   {
      if ( find_io_block(iobuf,&item_header) == 0 )
      {
         printf("Found I/O block of type %ld\n",item_header.type);
	 if ( read_io_block(iobuf,&item_header) == 0 )
         {
            // print_trgmask(iobuf); /* Note limitation to 20 or getenv(MAX_PRINT_ARRAY) entries */
            read_trgmask(iobuf,tms);
            printf("Triggermask set for run %ld has %zu entries\n", tms->run, tms->num_entries);
            for (i=0; i<tms->num_entries; i++)
            {
               tme = &tms->mask[i];
               if ( tme->event != event && event > 0 )
               {
                  printf("%d %d %d %d\n", n_lst, n_mst, n_dcsst, n_scsst);
                  n_lst = n_mst = n_dcsst = n_scsst = 0;
               }
               event = tme->event;
               if ( (tme->tel_id >= 1 && tme->tel_id <= 5) || 
                    (tme->tel_id >= 170 && tme->tel_id <= 173) )
                  n_lst++;
               else if ( (tme->tel_id >= 6 && tme->tel_id <= 59) || 
                    tme->tel_id == 169 ||
                    (tme->tel_id >= 174 && tme->tel_id <= 197) )
                  n_mst++;
               else if ( tme->tel_id >= 60 && tme->tel_id <= 96 )
                  n_dcsst++;
               else if ( tme->tel_id >= 97 && tme->tel_id <= 168 )
                  n_scsst++;
            }
            if ( event > 0 )
               printf("%d %d %d %d\n", n_lst, n_mst, n_dcsst, n_scsst);
         }
      }
      fileclose(iobuf->input_file);
      exit(0);
   }

   return 0;
}

/** @} */
