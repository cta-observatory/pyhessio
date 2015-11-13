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

/** @file gen_trgmask.c
 *  @short A utility program for fixing problems with simulation data which
 *         does not have the correct bit pattern of telescope triggers but
 *         the correct pattern can be extracted from the log files.
 *
@verbatim
Syntax: bin/gen_trgmask log-file [ trgmask-file ]
    or: bin/gen_trgmask -l  trgmask-file
The first variant will create a file with a single data block
for the trigger mask patterns recovered from the log file.
The default file name is derived with extension .trgmask.gz
Note that only data for one run per file is supported.
The second variant will list the contents of such a file.
@endverbatim
 *
 *  @author Konrad Bernloehr
 */

/** @defgroup gen_trgmask_c The gen_trgmask program */
/** @{ */

#include "initial.h"      /* This file includes others as required. */
#include "io_basic.h"     /* This file includes others as required. */
#include "fileopen.h"
#include "io_trgmask.h"

void syntax(char *prgname);

void syntax(char *prgname)
{
   fprintf(stderr,"Syntax: %s log-file [ trgmask-file ]\n", prgname);
   fprintf(stderr,"    or: %s -l  trgmask-file\n", prgname);
   fprintf(stderr,"The first variant will create a file with a single data block\n"
                  "for the trigger mask patterns recovered from the log file.\n"
                  "The default file name is derived with extension .trgmask.gz\n"
                  "Note that only data for one run per file is supported.\n"
                  "The second variant will list the contents of such a file.\n");
   exit(1);
}

int main(int argc, char **argv)
{
   struct trgmask_set *tms = (struct trgmask_set *) calloc(1,sizeof(struct trgmask_set));
   // struct trgmask_hash_set *ths = calloc(1,sizeof(struct trgmask_hash_set));
   char *in_fname, *out_fname = NULL;
   char *s;
   char *prgname = argv[0];

   IO_BUFFER *iobuf = allocate_io_buffer(1000000L);
   if ( iobuf == NULL )
   {
      Error("Cannot allocate I/O buffer");
      exit(1);
   }
   iobuf->max_length = 200000000L;

   if ( argc == 2 )
   {
      if ( strcmp(argv[1],"--help") == 0 )
         syntax(prgname);
      in_fname = argv[1];
      out_fname = (char *) malloc(strlen(argv[1])+12);
      strcpy(out_fname, in_fname);
      if ( (s = strstr(out_fname,".log")) != NULL )
         *s = '\0';
      strcat(out_fname,".trgmask.gz");
      if ( strcmp(in_fname,"-") == 0 )
         strcpy(out_fname,"-");
   }
   else if ( argc == 3 )
   {
      in_fname = argv[1];
      out_fname = argv[2];
      if ( strcmp(argv[1],"-l" ) == 0 )
      {
         IO_ITEM_HEADER item_header;
         iobuf->input_file = fileopen(out_fname,"r");
         if ( iobuf->input_file != NULL )
         {
	    if ( find_io_block(iobuf,&item_header) == 0 )
            {
               printf("Found I/O block of type %ld\n",item_header.type);
	       if ( read_io_block(iobuf,&item_header) == 0 )
               {
                  print_trgmask(iobuf); /* Note limitation to 20 or getenv(MAX_PRINT_ARRAY) entries */
                  /* For reading the block instead and setting it up for searches, you would need:
                     read_trgmask(iobuf,tms);
                     trgmask_fill_hashed(tms,ths);
                  */
                  /* For actually searching something in ths, you need:
                     struct trgmask_entry *tme = find_trgmask(ths,event,tel_id);
                  */
               }
            }
            fileclose(iobuf->input_file);
            exit(0);
         }
         else
            perror(out_fname);
         exit(1);
      }
   }
   else
      syntax(prgname);

   if ( trgmask_scan_log(tms,argv[1]) < 0 )
      exit(1);

   // printf("Setting up hashed list ...\n");
   // trgmask_fill_hashed(tms,ths);
   // print_hashed_trgmasks(ths);

   printf("Trigger mask data to be written to file %s\n", out_fname);
   iobuf->output_file = fileopen(out_fname,"w");
   if ( iobuf->output_file != (FILE *) NULL )
   {
      write_trgmask(iobuf,tms);
      fileclose(iobuf->output_file);
   }

   return 0;
}

/** @} */
