/* ============================================================================

Copyright (C) 2003, 2009  Konrad Bernloehr

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

/** @file extract_hess.c
 *  @short Extract part of the H.E.S.S. data from sim_hessarray.
 *
 *  @author Konrad Bernloehr
 *  @date    @verbatim CVS $Date: 2014/10/28 14:23:47 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.7 $ @endverbatim
 */

/** @defgroup extract_hess_c The extract_hess program */
/** @{ */

#include "initial.h"      /* This file includes others as required. */
#include "io_basic.h"     /* This file includes others as required. */
#include "mc_tel.h"
#include "history.h"
#include "io_hess.h"
#include "histogram.h"
#include "io_histogram.h"
#include "fileopen.h"

#include <signal.h>

static int interrupted;

/** Show program syntax */

static void syntax (char *program);

static void syntax (char *program)
{
   printf("Syntax: %s [ options ] [ - | input_fname ... ]\n",program);
   printf("Options:\n");
   printf("   -o fname     (Set output file name.)\n");
   printf("   -t type      (Extract calibration data of given type.)\n");
   printf("      Types:  0 Pedestals with closed lid (default).\n");
   printf("              1 Pedestals with open lid.\n");
   printf("            >=2 LED or Laser events of some amplitude level.\n");
   printf("   -b nb        (Change maximum size of I/O buffers.)\n");
   exit(1);
}

/* -------------------- main program ---------------------- */
/** 
 *  @short Main program 
 *
 *  Main program function of extract_hess.c program.
 */

int main (int argc, char **argv)
{
   IO_BUFFER *iobuf, *iobuf2;
   IO_ITEM_HEADER item_header, item_header2;
   const char *input_fname = NULL, *output_fname = "iact.simhess.extract";
   // int verbose = 0;
   char *program = argv[0];
   int iarg;
   int rc;
   int type_selected = 0;

   if ( argc < 2 )
      input_fname = "iact.simhess";

   if ( (iobuf = allocate_io_buffer(1000000L)) == NULL )
   {
      Error("Cannot allocate I/O buffer");
      exit(1);
   }
   iobuf->max_length = 200000000;
   if ( (iobuf2 = allocate_io_buffer(1000000L)) == NULL )
   {
      Error("Cannot allocate I/O buffer 2");
      exit(1);
   }
   iobuf2->max_length = 200000000;

   for (iarg=1; iarg<argc && input_fname == NULL; iarg++)
   {
       if ( strcmp(argv[iarg],"-o") == 0 && iarg+1 < argc )
	 output_fname = argv[++iarg];
       else if ( strncmp(argv[iarg],"-o",2) == 0 && strlen(argv[iarg]) > 2 )
	 output_fname = argv[iarg]+2;
       // else if ( strcmp(argv[iarg],"-v") == 0 )
       //	  verbose = 1;
       else if ( strcmp(argv[iarg],"-t") == 0 && iarg+1 < argc )
       	  type_selected = atoi(argv[++iarg]);
       else if ( strcmp(argv[iarg],"-b") == 0 && iarg+1 < argc )
       	  iobuf->max_length = iobuf2->max_length = atol(argv[++iarg]);
       else if ( argv[iarg][0] == '-' && argv[iarg][1] != '\0' )
         syntax(program);
       else
	 input_fname = argv[iarg];
    }
    
    fprintf(stderr,"Extracting calibration data of type %d to %s\n",type_selected,output_fname);

    if ( input_fname == 0 )
    {
       Error("No input file.\n");
       syntax(program);
    }
    if ( strcmp(input_fname ,"-") == 0 )
      iobuf->input_file = stdin;
    else if ( (iobuf->input_file = fileopen(input_fname,READ_BINARY)) == NULL )
    {
      perror(input_fname);
      Error("Cannot open input file.");
      exit(1);
    }
    printf("\nInput file '%s' has been opened.\n",input_fname);
    input_fname = NULL;

    if ( (iobuf2->output_file = fileopen(output_fname,WRITE_BINARY)) == NULL )
    {
      perror(output_fname);
      Error("Cannot open output file.");
      exit(1);
    }
    printf("\nOutput file '%s' has been opened.\n",output_fname);

    for (;;) /* Loop over all data in the input file */
    {
      if ( interrupted )
         break;

      /* Find and read the next block of data. */
      /* In case of problems with the data, just give up. */
      if ( find_io_block(iobuf,&item_header) != 0 )
         break;
      if ( read_io_block(iobuf,&item_header) != 0 )
         break;

      /* What did we actually get? */
      switch ( (int) item_header.type )
      {
         /* =================================================== */
         case IO_TYPE_HESS_RUNHEADER:
         case IO_TYPE_HESS_MCRUNHEADER:
	 case IO_TYPE_MC_INPUTCFG:
         case IO_TYPE_HESS_CAMSETTINGS:
         case IO_TYPE_HESS_CAMORGAN:
         case IO_TYPE_HESS_PIXELDISABLE:
         case IO_TYPE_HESS_CAMSOFTSET:
         case IO_TYPE_HESS_POINTINGCOR:
         case IO_TYPE_HESS_TRACKSET:
         case IO_TYPE_MC_TELARRAY:
         case IO_TYPE_HESS_MC_PE_SUM:
         case IO_TYPE_HESS_TEL_MONI:
         case IO_TYPE_HESS_LASCAL:
         case IO_TYPE_HESS_RUNSTAT:
         case IO_TYPE_HESS_MC_RUNSTAT:
            /* Copy it to output. */
            reset_io_block(iobuf2);
            get_item_begin(iobuf,&item_header);
            copy_item_to_io_block(iobuf2,iobuf,&item_header);
            get_item_end(iobuf,&item_header);
            write_io_block(iobuf2);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_CALIBEVENT:
            /* Unwrap it or get rid of it. */
            reset_io_block(iobuf2);
            if ( (rc = get_item_begin(iobuf,&item_header)) < 0 )
               continue;
            if ( item_header.ident == type_selected )
            {
               item_header2.type = IO_TYPE_HESS_EVENT;
               if ( (rc = get_item_begin(iobuf,&item_header2)) < 0 )
                  continue;
               copy_item_to_io_block(iobuf2,iobuf,&item_header2);
               get_item_end(iobuf,&item_header2);
               write_io_block(iobuf2);
            }
            get_item_end(iobuf,&item_header);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_EVENT:
         case IO_TYPE_HESS_MC_EVENT:
         case 70:
         case 100:
            /* Get rid of it */
            break;
      }
   }

   fileclose(iobuf2->output_file);

   return 0;
}

/** @} */
