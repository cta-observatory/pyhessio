/* ============================================================================

   Copyright (C) 1994, 2001, 2007, 2009  Konrad Bernloehr

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

/** @file listio.c
    @short Main function for listing data consisting of eventio blocks.
  
    @author  Konrad Bernloehr
    @date    @verbatim CVS $Date: 2014/06/01 11:33:05 $ @endverbatim
    @version @verbatim CVS $Revision: 1.14 $ @endverbatim

    The item type, version, length and ident are displayed.
    With command line option '-s' all sub-items are shown
    as well.
    Input is from standard input by default, output to standard output.

@verbatim
    Syntax: listio [-s[n]] [-p] [filename]
    List structure of eventio data files.
       -s : also list contained (sub-) items
       -sn: list sub-items up to depth n (n=0,1,...)
       -p : show positions of items in the file
    If no file name given, standard input is used.
@endverbatim

*/

/** @defgroup listio_c The listio program */
/** @{ */

#include "initial.h"
#include "io_basic.h"
#include "fileopen.h"

/** 
 * @short Main function 
 *
 * The main function of the listio program.
 */

int main (int argc, char **argv)
{
   int sub, show_pos;
   IO_BUFFER *iobuf;
   IO_ITEM_HEADER item_header;
#if defined(__USE_LARGEFILE64)
   off64_t pos = 0;
#else
   long pos = 0;
#endif
   FILE *input;
   int verbosity = 0;

#ifdef ALWAYS_WITH_REGISTRY
   /* Use default registry of known types with any compiler. */
   /* Requires linking against the library (not available with */
   /* bare-bone eventio, as for CORSIKA IACT/ATMO extension). */
   /* Not necessary when compiled with GCC. */
   set_ev_reg_std();
#endif

   sub = -1;
   show_pos = 0;
   while ( argc >= 2 )
   {
      if ( argv[1][0] != '-' )
         break;
      if ( strncmp(argv[1],"-s",2) == 0 )
      {
         if ( strlen(argv[1]) > 2 )
            sub = atoi(argv[1]+2);
         else
            sub = 20;
      }
      else if ( strcmp(argv[1],"-p") == 0 )
      {
         show_pos = 1;
      }
      else if ( strcmp(argv[1],"-n") == 0 )
      {
         verbosity = 1;
      }
      else if ( strcmp(argv[1],"-d") == 0 )
      {
         verbosity = 2;
      }
      else
      {
         fprintf(stderr,"Syntax: listio [-s[n]] [-p] [filename]\n");
         fprintf(stderr,"List structure of eventio data files.\n");
         fprintf(stderr,"   -s : also list contained (sub-) items\n");
         fprintf(stderr,"   -sn: list sub-items up to depth n (n=0,1,...)\n");
         fprintf(stderr,"   -p : show positions of items in the file\n");
         fprintf(stderr,"   -n : show type names where known\n");
         fprintf(stderr,"   -d : show type names and descriptions where known\n");
         fprintf(stderr,"If no file name given, standard input is used.\n");
         exit(1);
      }
      argc--;
      argv++;
   }
   if ( (iobuf = allocate_io_buffer(1000)) == (IO_BUFFER *) NULL )
      exit(1);
   iobuf->max_length = 128000000;
   if ( argc > 1 )
   {
      if ( (input = fileopen(argv[1],READ_BINARY)) == (FILE *) NULL )
      {
         perror(argv[1]);
         exit(1);
      }
      iobuf->input_file = input;
   }
   else
      iobuf->input_file = stdin;
   
   if ( show_pos && sub < 0 )
      sub = 0;

   if ( sub < 0 )
      (void) list_io_blocks(iobuf, verbosity);
   else
      while ( find_io_block(iobuf,&item_header) >= 0 )
      {
         if ( show_pos )
	 {
#if defined(__USE_LARGEFILE64)
            pos = ftello64(iobuf->input_file) - 16;
#else
            pos = ftell(iobuf->input_file) - 16;
#endif
         }
         if ( read_io_block(iobuf,&item_header) < 0 )
            break;
         list_sub_items(iobuf,&item_header,sub, verbosity);
         if ( show_pos )
         {
            fflush(NULL);
#if defined(__USE_LARGEFILE64)
            if (sizeof(pos) > sizeof(long))
               printf("(I/O block started at byte offset %lld)\n",(long long)pos);
            else
               printf("(I/O block started at byte offset %ld)\n",(long)pos);
#else
            printf("(I/O block started at byte offset %ld)\n",pos);
#endif
         }
      }

   exit(0);
   return 0;
}

#ifdef OS_OS9
int perror(text)
    char *text;
{
    fprintf(stderr,"%s: Error\n",text);
    return 0;
}
#endif

/** @} */
