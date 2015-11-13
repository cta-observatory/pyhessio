/* ============================================================================

Copyright (C) 2010  Konrad Bernloehr

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

/**
 *  @file fcat.c
 *  @short Trivial test and utility program for the fileopen/fileclose functions. 
 */

/** @defgroup fact_c The fcat program */
/** @{ */

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "fileopen.h"

#define BSIZE 8192

int main (int argc, char **argv)
{
   FILE *f = NULL;
   int iarg;
   size_t n;
   char buffer[BSIZE];

   for (iarg=1; iarg<argc; iarg++)
   {
      if ( (f = fileopen(argv[iarg],"r")) == NULL )
      {
         if ( errno != 0 )
            perror(argv[iarg]);
         fflush(0);
         fprintf(stderr,"(An error opening the file)\n");
         exit(1);
      }
      while ( (n=fread(buffer,1,BSIZE,f)) == BSIZE )
         fwrite(buffer,1,n,stdout);
      if ( n > 0 )
         fwrite(buffer,1,n,stdout);
      if ( ferror(f) )
      {
         if ( errno != 0 )
            perror(argv[iarg]);
         fflush(0);
         fprintf(stderr,"(An error reading from the file)\n");
         exit(1);
      }
      if ( fileclose(f) != 0 )
      {
         if ( errno != 0 )
            perror(argv[iarg]);
         fflush(0);
         fprintf(stderr,"(An error closing the file)\n");
         exit(1);
      }
   }

   return 0;
}

/** @} */
