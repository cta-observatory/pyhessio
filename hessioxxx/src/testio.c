/* ============================================================================

Copyright (C) 1994, 1997, 2000, 2001, 2007, 2009, 2012  Konrad Bernloehr

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

/** @file testio.c
    @short Test program for eventio data format.
    
    @author  Konrad Bernloehr
    @date    1994 to 2010
    @date    @verbatim CVS $Date: 2014/08/04 13:05:28 $ @endverbatim
    @version @verbatim CVS $Revision: 1.22 $ @endverbatim
*/

/** @defgroup testio_c The testio program */
/** @{ */

#include "initial.h"
#include "warning.h"
#include "io_basic.h"
#include "fileopen.h"

struct test_struct
{
   long lvar[2];
   int ilvar[2];
   int isvar[2];
   short svar[3];
   double fvar[2];
   double dvar[2];
   double hvar[2];
   size_t nbvar;
   uint8_t bvar[2];
   int8_t i8var[2];
   uint8_t u8var[2];
   int16_t i16var[2];
   uint16_t u16var[2];
   int32_t i32var[2];
   uint32_t u32var[2];
#ifdef HAVE_64BIT_INT
   int64_t i64var[2];
   uint64_t u64var[2];
#endif
   uint16_t cnt16var[4];
   int16_t scnt16var[10];
   uint32_t cnt32var[6];
   int32_t scnt32var[12];
   size_t cntzvar[6];
   ssize_t scntzvar[12];
#ifdef WITH_UINTMAX_T
   uintmax_t cntvar[8];
#else
   size_t cntvar[8];
#endif
#ifdef WITH_INTMAX_T
   intmax_t scntvar[14];
#else
   ssize_t scntvar[14];
#endif
   char str16var[10];
   char str32var[10];
   char strvvar[10];
};

typedef struct test_struct TEST_DATA;

static int care_long, care_int, care_short;

int datacmp (TEST_DATA *data1, TEST_DATA *data2);
int write_test1 (TEST_DATA *data, IO_BUFFER *iobuf);
int read_test1 (TEST_DATA *data, IO_BUFFER *iobuf);
int write_test2 (TEST_DATA *data, IO_BUFFER *iobuf);
int read_test2 (TEST_DATA *data, IO_BUFFER *iiobuf);
int write_test3 (TEST_DATA *data, IO_BUFFER *iobuf);
int read_test3 (TEST_DATA *data, IO_BUFFER *iiobuf);

/* ------------------------ datacmp ---------------------- */
/**
 *  @short Compare elements of test data structures.
 *  Compare elements of test data structures with the accuracy
 *  relevant to the I/O package.
 *
 *  @param data1 first data structure
 *  @param data2 second data structure
 *
 *  @return 0 (something did not match), 1 (O.K.)
 *
 */

int datacmp(TEST_DATA *data1, TEST_DATA *data2)
{
   size_t i, ok;
   
   ok = 1;
   
   if ( data1->nbvar != data2->nbvar )
   {
      fprintf(stderr,"Bool variables do not match.\n");
      ok = 0;
   }
   else
      for (i=0; i<data1->nbvar; i++)
      {
         size_t m = i/8, n = i%8;
         if ( (data1->bvar[m]&(1<<n)) != (data2->bvar[m]&(1<<n)) )
         {
            fprintf(stderr,"Bool variable %d does not match: %d <--> %d\n",
               (unsigned int)i+1, data1->bvar[m]&(1<<n), data2->bvar[m]&(1<<n));
            ok = 0;
         }
      }

   for (i=0; i<2; i++)
      if ( (int32_t)data1->lvar[i] != (int32_t)data2->lvar[i] )
      {
         fprintf(stderr,"Long variable %d does not match: %08lx <--> %08lx\n",
            (unsigned int)i+1,data1->lvar[i],data2->lvar[i]);
         ok = 0;
      }
      else if ( data1->lvar[i] != data2->lvar[i] )
         care_long = 1;
      
   for (i=0; i<2; i++)
      if ( (int32_t)data1->ilvar[i] != (int32_t)data2->ilvar[i] )
      {
         fprintf(stderr,"Int 'l' variable %d does not match: %08x <--> %08x\n",
            (unsigned int)i+1,data1->ilvar[i],data2->ilvar[i]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( (int16_t)data1->isvar[i] != (int16_t)data2->isvar[i] )
      {
         fprintf(stderr,"Int 's' variable %d does not match: %04x <--> %04x\n",
            (unsigned int)i+1,data1->isvar[i],data2->isvar[i]);
         ok = 0;
      }
      else if ( data1->isvar[i] != data2->isvar[i] )
         care_short = 1;
      
   for (i=0; i<3; i++)
      if ( data1->svar[i] != data2->svar[i] )
      {
         fprintf(stderr,"Short variable %d does not match: %04x <--> %04x\n",
            (unsigned int)i+1,data1->svar[i],data2->svar[i]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
   {
      union f_i { float f; uint32_t i; } f1, f2;
      f1.f = (float) data1->fvar[i];
      f2.f = (float) data2->fvar[i];
      if ( f1.f != f2.f )
      {
         fprintf(stderr,"Float variable %u does not match: %08lx <--> %08lx\n",
            (unsigned int)i+1, (unsigned long)f1.i, (unsigned long)f2.i);
         ok = 0;
      }
   }
   
   for (i=0; i<2; i++)
      if ( data1->dvar[i] != data2->dvar[i] )
      {
         union d_i { double d; uint32_t i[2]; } f1, f2;
	 f1.d = data1->dvar[i];
	 f2.d = data2->dvar[i];
         fprintf(stderr,
	    "Double variable %u does not match: %08lx%08lx <--> %08lx%08lx\n",
            (unsigned int)i+1, (unsigned long)f1.i[0], (unsigned long)f1.i[1], 
            (unsigned long)f2.i[0], (unsigned long)f2.i[1]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1->hvar[i] != data2->hvar[i] && fabs(data1->hvar[i]-data2->hvar[i]) > 0.002*fabs(data1->hvar[i]) )
      {
         union d_i { double d; uint32_t i[2]; } f1, f2;
	 f1.d = data1->hvar[i];
	 f2.d = data2->hvar[i];
         fprintf(stderr,
	    "Sfloat variable %u does not match: %08lx%08lx <--> %08lx%08lx\n",
            (unsigned int)i+1, (unsigned long)f1.i[0], (unsigned long)f1.i[1], 
            (unsigned long)f2.i[0], (unsigned long)f2.i[1]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1->i8var[i] != data2->i8var[i] )
      {
         fprintf(stderr,"Int8 variable %d does not match: %08x <--> %08x\n",
            (unsigned int)i+1,data1->i8var[i],data2->i8var[i]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1->u8var[i] != data2->u8var[i] )
      {
         fprintf(stderr,"UInt8 variable %d does not match: %08x <--> %08x\n",
            (unsigned int)i+1,data1->u8var[i],data2->u8var[i]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1->i16var[i] != data2->i16var[i] )
      {
         fprintf(stderr,"Int16 variable %d does not match: %08x <--> %08x\n",
            (unsigned int)i+1,data1->i16var[i],data2->i16var[i]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1->u16var[i] != data2->u16var[i] )
      {
         fprintf(stderr,"UInt16 variable %d does not match: %08x <--> %08x\n",
            (unsigned int)i+1,data1->u16var[i],data2->u16var[i]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1->i32var[i] != data2->i32var[i] )
      {
         fprintf(stderr,"Int32 variable %d does not match: %08x <--> %08x\n",
            (unsigned int)i+1,data1->i32var[i],data2->i32var[i]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1->u32var[i] != data2->u32var[i] )
      {
         fprintf(stderr,"UInt32 variable %d does not match: %08x <--> %08x\n",
            (unsigned int)i+1,data1->u32var[i],data2->u32var[i]);
         ok = 0;
      }

#ifdef HAVE_64BIT_INT
   for (i=0; i<2; i++)
      if ( data1->i64var[i] != data2->i64var[i] )
      {
#ifdef SIXTY_FOUR_BITS
         fprintf(stderr,"Int64 variable %d does not match: %016lx <--> %016lx\n",
            (unsigned int)i+1,data1->i64var[i],data2->i64var[i]);
#else
         fprintf(stderr,"Int64 variable %d does not match: %016llx <--> %016llx\n",
            (unsigned int)i+1,data1->i64var[i],data2->i64var[i]);
#endif
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1->u64var[i] != data2->u64var[i] )
      {
#ifdef SIXTY_FOUR_BITS
         fprintf(stderr,"UInt64 variable %d does not match: %016lx <--> %016lx\n",
            (unsigned int)i+1,data1->u64var[i],data2->u64var[i]);
#else
         fprintf(stderr,"UInt64 variable %d does not match: %016llx <--> %016llx\n",
            (unsigned int)i+1,data1->u64var[i],data2->u64var[i]);
#endif
         ok = 0;
      }
#endif

   for (i=0; i<4; i++)
      if ( data1->cnt16var[i] != data2->cnt16var[i] )
      {
         fprintf(stderr,
            "Count16 variable %u does not match: %08x <--> %08x\n",
            (unsigned int)i+1, data1->cnt16var[i], data2->cnt16var[i]);
         ok = 0;
      }

   for (i=0; i<6; i++)
      if ( data1->cnt32var[i] != data2->cnt32var[i] )
      {
         fprintf(stderr,
            "Count32 variable %u does not match: %08x <--> %08x\n",
            (unsigned int)i+1, data1->cnt32var[i], data2->cnt32var[i]);
         ok = 0;
      }

   for (i=0; i<6; i++)
      if ( data1->cntzvar[i] != data2->cntzvar[i] )
      {
         fprintf(stderr,
            "Count variable %u does not match: %08zx <--> %08zx\n",
            (unsigned int)i+1, data1->cntzvar[i], data2->cntzvar[i]);
         ok = 0;
      }

   for (i=0; i<8; i++)
      if ( data1->cntvar[i] != data2->cntvar[i] )
      {
         fprintf(stderr,
#if defined(WITH_UINTMAX_T)
            "Long count variable %d does not match: %08jx <--> %08jx\n",
#else
            "Long count variable %d does not match: %08zx <--> %08zx\n",
#endif
            (unsigned int)i+1, data1->cntvar[i], data2->cntvar[i]);
         ok = 0;
      }

   for (i=0; i<10; i++)
      if ( data1->scnt16var[i] != data2->scnt16var[i] )
      {
         fprintf(stderr,
            "Signed count16 variable %u does not match: %08x <--> %08x\n",
            (unsigned int)i+1, data1->scnt16var[i], data2->scnt16var[i]);
         ok = 0;
      }

   for (i=0; i<12; i++)
      if ( data1->scnt32var[i] != data2->scnt32var[i] )
      {
         fprintf(stderr,
            "Signed count32 variable %u does not match: %08x <--> %08x\n",
            (unsigned int)i+1, data1->scnt32var[i], data2->scnt32var[i]);
         ok = 0;
      }

   for (i=0; i<12; i++)
      if ( data1->scntzvar[i] != data2->scntzvar[i] )
      {
         fprintf(stderr,
            "Signed count variable %u does not match: %08zx <--> %08zx\n",
            (unsigned int)i+1, data1->scntzvar[i], data2->scntzvar[i]);
         ok = 0;
      }

   for (i=0; i<14; i++)
      if ( data1->scntvar[i] != data2->scntvar[i] )
      {
         fprintf(stderr,
#if defined(WITH_UINTMAX_T)
            "Signed long count variable %u does not match: %08jx <--> %08jx\n",
#else
            "Signed long count variable %u does not match: %08zx <--> %08zx\n",
#endif
            (unsigned int)i+1, data1->scntvar[i], data2->scntvar[i]);
         ok = 0;
      }

   if ( strcmp(data1->str16var,data2->str16var) != 0 )
   {
      fprintf(stderr,"Strings with 16-bit count do not match: '%s' <--> '%s'\n",
         data1->str16var, data2->str16var);
         ok = 0;
   }
   if ( strcmp(data1->str32var,data2->str32var) != 0 )
   {
      fprintf(stderr,"Strings with 32-bit count do not match: '%s' <--> '%s'\n",
         data1->str32var, data2->str32var);
         ok = 0;
   }
   if ( strcmp(data1->strvvar,data2->strvvar) != 0 )
   {
      fprintf(stderr,"Strings with variable length count do not match: '%s' <--> '%s'\n",
         data1->strvvar, data2->strvvar);
         ok = 0;
   }

   return ok;
}

/* --------------------- write_test1 ---------------------- */
/**
 *  @short Write test data with single-element functions
 *
 *  @param    data   Pointer to test data structure
 *  @param    iobuf  Pointer to I/O buffer
 *
 *  @return  0 (O.K.), <0 (error as for put_item_end())
 *
 */

int write_test1(TEST_DATA *data, IO_BUFFER *iobuf)
{
   IO_ITEM_HEADER item_header;
   int i;
   
   item_header.type = 99;             /* test data */
   item_header.version = 0;          /* Version 0 (test) */
   item_header.ident = 123;
   
   put_item_begin(iobuf,&item_header);

   put_long(data->lvar[0],iobuf);
   put_long(data->lvar[1],iobuf);
   put_long((long)data->ilvar[0],iobuf);
   put_long((long)data->ilvar[1],iobuf);
   put_short(data->isvar[0],iobuf);
   put_short(data->isvar[1],iobuf);
   put_short(data->svar[0],iobuf);
   put_short(data->svar[1],iobuf);
   put_short(data->svar[2],iobuf);
   put_real(data->fvar[0],iobuf);
   put_real(data->fvar[1],iobuf);
   put_double(data->dvar[0],iobuf);
   put_double(data->dvar[1],iobuf);
   put_sfloat(data->hvar[0],iobuf);
   put_sfloat(data->hvar[1],iobuf);
   put_byte(data->i8var[0],iobuf);
   put_byte(data->i8var[1],iobuf);
   put_byte(data->u8var[0],iobuf);
   put_byte(data->u8var[1],iobuf);
   put_short(data->i16var[0],iobuf);
   put_short(data->i16var[1],iobuf);
   put_short(data->u16var[0],iobuf);
   put_short(data->u16var[1],iobuf);
   put_int32(data->i32var[0],iobuf);
   put_int32(data->i32var[1],iobuf);
   put_uint32(data->u32var[0],iobuf);
   put_uint32(data->u32var[1],iobuf);
#ifdef HAVE_64BIT_INT
   put_vector_of_int64(&data->i64var[0],1,iobuf);
   put_vector_of_int64(&data->i64var[1],1,iobuf);
   put_vector_of_uint64(&data->u64var[0],1,iobuf);
   put_vector_of_uint64(&data->u64var[1],1,iobuf);
#endif   
   put_count(data->nbvar,iobuf);
   put_byte(data->bvar[0],iobuf);
   put_byte(data->bvar[1],iobuf);
   for (i=0; i<4; i++)
      put_count16(data->cnt16var[i],iobuf);
   for (i=0; i<6; i++)
      put_count32(data->cnt32var[i],iobuf);
   for (i=0; i<6; i++)
      put_count(data->cntzvar[i],iobuf);
   for (i=0; i<8; i++)
      put_count(data->cntvar[i],iobuf);
   for (i=0; i<10; i++)
      put_scount16(data->scnt16var[i],iobuf);
   for (i=0; i<12; i++)
      put_scount32(data->scnt32var[i],iobuf);
   for (i=0; i<12; i++)
      put_scount(data->scntzvar[i],iobuf);
   for (i=0; i<14; i++)
      put_scount(data->scntvar[i],iobuf);
   put_string(data->str16var,iobuf);
   put_long_string(data->str32var,iobuf);
   put_var_string(data->strvvar,iobuf);
   return(put_item_end(iobuf,&item_header));
}

/* ---------------------- read_test1 ---------------------- */
/**
 *  @short Read test data with single-element functions
 *
 *  @param   data   Pointer to test data structure
 *  @param   iobuf  Pointer to I/O buffer
 *
 *  @return  0 (ok), <0 (error as for get_item_end())
 *
 */

int read_test1(TEST_DATA *data, IO_BUFFER *iobuf)
{
   IO_ITEM_HEADER item_header;
   int i;
   
   item_header.type = 99;             /* test data */
   if ( get_item_begin(iobuf,&item_header) < 0 )
   {
      Warning("Missing or invalid test data block.");
      return -4;
   }

   data->lvar[0]  = get_long(iobuf);
   data->lvar[1]  = get_long(iobuf);
   data->ilvar[0] = get_long(iobuf);
   data->ilvar[1] = get_long(iobuf);
   data->isvar[0] = get_short(iobuf);
   data->isvar[1] = get_short(iobuf);
   data->svar[0]  = get_short(iobuf);
   data->svar[1]  = get_short(iobuf);
   data->svar[2]  = get_short(iobuf);
   data->fvar[0]  = get_real(iobuf);
   data->fvar[1]  = get_real(iobuf);
   data->dvar[0]  = get_double(iobuf);
   data->dvar[1]  = get_double(iobuf);
   data->hvar[0]  = get_sfloat(iobuf);
   data->hvar[1]  = get_sfloat(iobuf);
   data->i8var[0] = get_byte(iobuf);
   data->i8var[1] = get_byte(iobuf);
   data->u8var[0] = get_byte(iobuf);
   data->u8var[1] = get_byte(iobuf);
   data->i16var[0] = get_short(iobuf);
   data->i16var[1] = get_short(iobuf);
   data->u16var[0] = get_short(iobuf);
   data->u16var[1] = get_short(iobuf);
   data->i32var[0] = get_int32(iobuf);
   data->i32var[1] = get_int32(iobuf);
   data->u32var[0] = get_uint32(iobuf);
   data->u32var[1] = get_uint32(iobuf);
#ifdef HAVE_64BIT_INT
   get_vector_of_int64(&data->i64var[0],1,iobuf);
   get_vector_of_int64(&data->i64var[1],1,iobuf);
   get_vector_of_uint64(&data->u64var[0],1,iobuf);
   get_vector_of_uint64(&data->u64var[1],1,iobuf);
#endif   
   data->nbvar = get_count(iobuf);
   data->bvar[0] = get_byte(iobuf);
   data->bvar[1] = get_byte(iobuf);
   for (i=0; i<4; i++)
      data->cnt16var[i] = get_count16(iobuf);
   for (i=0; i<6; i++)
      data->cnt32var[i] = get_count32(iobuf);
   for (i=0; i<6; i++)
      data->cntzvar[i] = get_count(iobuf);
   for (i=0; i<8; i++)
      data->cntvar[i] = get_count(iobuf);
   for (i=0; i<10; i++)
      data->scnt16var[i] = get_scount16(iobuf);
   for (i=0; i<12; i++)
      data->scnt32var[i] = get_scount32(iobuf);
   for (i=0; i<12; i++)
      data->scntzvar[i] = get_scount(iobuf);
   for (i=0; i<14; i++)
      data->scntvar[i] = get_scount(iobuf);
   get_string(data->str16var,sizeof(data->str16var),iobuf);
   get_long_string(data->str32var,sizeof(data->str32var),iobuf);
   get_var_string(data->strvvar,sizeof(data->strvvar),iobuf);

   return(get_item_end(iobuf,&item_header));
}

/* --------------------- write_test2 ---------------------- */
/**
 *  @short Write test data with vector functions as far as possible
 *
 *  @param    data    Pointer to test data structure
 *  @param    iobuf   Pointer to I/O buffer
 *
 *  @return  0 (ok), <0 (error as for put_item_end())
 *
 */

int write_test2(TEST_DATA *data, IO_BUFFER *iobuf)
{
   IO_ITEM_HEADER item_header;
   int i;
   
   item_header.type = 99;             /* test data */
   item_header.version = 0;          /* Version 0 (test) */
   item_header.ident = 123;

   put_item_begin(iobuf,&item_header);

   put_vector_of_long(data->lvar,2,iobuf);
   put_long((long)data->ilvar[0],iobuf);
   put_long((long)data->ilvar[1],iobuf);
   put_vector_of_int(data->isvar,2,iobuf);
   put_vector_of_short(data->svar,3,iobuf);
   put_vector_of_real(data->fvar,2,iobuf);
   put_vector_of_double(data->dvar,2,iobuf);
   put_sfloat(data->hvar[0],iobuf);
   put_sfloat(data->hvar[1],iobuf);
   put_vector_of_byte((uint8_t *)data->i8var,2,iobuf);
   put_vector_of_byte(data->u8var,2,iobuf);
   put_vector_of_short(data->i16var,2,iobuf);
   put_vector_of_short((int16_t *)data->u16var,2,iobuf);
   put_vector_of_int32(data->i32var,2,iobuf);
   put_vector_of_uint32(data->u32var,2,iobuf);
#ifdef HAVE_64BIT_INT
   put_vector_of_int64(data->i64var,2,iobuf);
   put_vector_of_uint64(data->u64var,2,iobuf);
#endif
   put_count(data->nbvar,iobuf);
   put_vector_of_byte(data->bvar,2,iobuf);
   for (i=0; i<4; i++)
      put_count16(data->cnt16var[i],iobuf);
   for (i=0; i<6; i++)
      put_count32(data->cnt32var[i],iobuf);
   for (i=0; i<6; i++)
      put_count(data->cntzvar[i],iobuf);
   for (i=0; i<8; i++)
      put_count(data->cntvar[i],iobuf);
   for (i=0; i<10; i++)
      put_scount16(data->scnt16var[i],iobuf);
   for (i=0; i<12; i++)
      put_scount32(data->scnt32var[i],iobuf);
   for (i=0; i<12; i++)
      put_scount(data->scntzvar[i],iobuf);
   for (i=0; i<14; i++)
      put_scount(data->scntvar[i],iobuf);
   put_string(data->str16var,iobuf);
   put_long_string(data->str32var,iobuf);
   put_var_string(data->strvvar,iobuf);
   
   return(put_item_end(iobuf,&item_header));
}

/* ---------------------- read_test2 ---------------------- */
/**
 *  @short Read test data with vector functions as far as possible
 *
 *  @param   data   Pointer to test data structure
 *  @param   iobuf  Pointer to I/O buffer
 *
 *  @return  0 (ok), <0 (error as for get_item_end())
 *
 */

int read_test2(TEST_DATA *data, IO_BUFFER *iobuf)
{
   IO_ITEM_HEADER item_header;
   int i;
   
   item_header.type = 99;             /* test data */
   if ( get_item_begin(iobuf,&item_header) < 0 )
   {
      Warning("Missing or invalid test data block.");
      return -4;
   }

   get_vector_of_long(data->lvar,2,iobuf);
   data->ilvar[0] = get_long(iobuf);
   data->ilvar[1] = get_long(iobuf);
   get_vector_of_int(data->isvar,2,iobuf);
   get_vector_of_short(data->svar,3,iobuf);
   get_vector_of_real(data->fvar,2,iobuf);
   get_vector_of_double(data->dvar,2,iobuf);
   data->hvar[0] = get_sfloat(iobuf);
   data->hvar[1] = get_sfloat(iobuf);
   get_vector_of_byte((uint8_t *)data->i8var,2,iobuf);
   get_vector_of_byte(data->u8var,2,iobuf);
   get_vector_of_short(data->i16var,2,iobuf);
   get_vector_of_short((int16_t *)data->u16var,2,iobuf);
   get_vector_of_int32(data->i32var,2,iobuf);
   get_vector_of_uint32(data->u32var,2,iobuf);
#ifdef HAVE_64BIT_INT
   get_vector_of_int64(data->i64var,2,iobuf);
   get_vector_of_uint64(data->u64var,2,iobuf);
#endif
   data->nbvar = get_count(iobuf);
   get_vector_of_byte(data->bvar,2,iobuf);
   for (i=0; i<4; i++)
      data->cnt16var[i] = get_count(iobuf);
   for (i=0; i<6; i++)
      data->cnt32var[i] = get_count(iobuf);
   for (i=0; i<6; i++)
      data->cntzvar[i] = get_count(iobuf);
   for (i=0; i<8; i++)
      data->cntvar[i] = get_count(iobuf);
   for (i=0; i<10; i++)
      data->scnt16var[i] = get_scount16(iobuf);
   for (i=0; i<12; i++)
      data->scnt32var[i] = get_scount32(iobuf);
   for (i=0; i<12; i++)
      data->scntzvar[i] = get_scount(iobuf);
   for (i=0; i<14; i++)
      data->scntvar[i] = get_scount(iobuf);
   get_string(data->str16var,sizeof(data->str16var),iobuf);
   get_long_string(data->str32var,sizeof(data->str32var),iobuf);
   get_var_string(data->strvvar,sizeof(data->strvvar),iobuf);

   return(get_item_end(iobuf,&item_header));
}

/* --------------------- write_test3 ---------------------- */
/**
 *  @short Write test data in nested items
 *
 *  @param    data    Pointer to test data structure
 *  @param    iobuf   Pointer to I/O buffer
 *
 *  @return  0 (ok), <0 (error as for put_item_end())
 *
 */

int write_test3(TEST_DATA *data, IO_BUFFER *iobuf)
{
   IO_ITEM_HEADER item_header1, item_header2;
   int i;
   
   item_header1.type = 990;           /* test data */
   item_header1.version = 11;         /* Version 11 */
   item_header1.ident = 1;

   put_item_begin(iobuf,&item_header1);

   item_header2.type = 991;           /* test data */
   item_header2.version = 12;         /* Version 12 */
   item_header2.ident = 1;

   put_item_begin(iobuf,&item_header2);

   put_vector_of_long(data->lvar,2,iobuf);
   put_long((long)data->ilvar[0],iobuf);
   put_long((long)data->ilvar[1],iobuf);
   put_vector_of_int(data->isvar,2,iobuf);
   put_vector_of_short(data->svar,3,iobuf);
   
   put_item_end(iobuf,&item_header2);
   
   item_header2.type = 992;           /* test data */
   item_header2.version = 12;         /* Version 12 */
   item_header2.ident = 2;

   put_item_begin(iobuf,&item_header2);

   put_vector_of_real(data->fvar,2,iobuf);
   put_vector_of_double(data->dvar,2,iobuf);
   put_sfloat(data->hvar[0],iobuf);
   put_sfloat(data->hvar[1],iobuf);
   put_vector_of_byte((uint8_t *)data->i8var,2,iobuf);
   put_vector_of_byte(data->u8var,2,iobuf);
   put_vector_of_short(data->i16var,2,iobuf);
   put_vector_of_short((int16_t *)data->u16var,2,iobuf);
   put_vector_of_int32(data->i32var,2,iobuf);
   put_vector_of_uint32(data->u32var,2,iobuf);
#ifdef HAVE_64BIT_INT
   put_vector_of_int64(data->i64var,2,iobuf);
   put_vector_of_uint64(data->u64var,2,iobuf);
#endif

   put_item_end(iobuf,&item_header2);

   item_header2.type = 993;           /* test data */
   item_header2.version = 12;         /* Version 12 */
   item_header2.ident = 3;

   put_item_begin(iobuf,&item_header2);

   put_count(data->nbvar,iobuf);
   put_vector_of_byte(data->bvar,2,iobuf);
   for (i=0; i<4; i++)
      put_count16(data->cnt16var[i],iobuf);
   for (i=0; i<6; i++)
      put_count32(data->cnt32var[i],iobuf);
   for (i=0; i<6; i++)
      put_count(data->cntzvar[i],iobuf);
   for (i=0; i<8; i++)
      put_count(data->cntvar[i],iobuf);
   for (i=0; i<10; i++)
      put_scount16(data->scnt16var[i],iobuf);
   for (i=0; i<12; i++)
      put_scount32(data->scnt32var[i],iobuf);
   for (i=0; i<12; i++)
      put_scount(data->scntzvar[i],iobuf);
   for (i=0; i<14; i++)
      put_scount(data->scntvar[i],iobuf);

   put_item_end(iobuf,&item_header2);

   item_header2.type = 994;           /* test data */
   item_header2.version = 12;         /* Version 12 */
   item_header2.ident = 4;

   put_item_begin(iobuf,&item_header2);

   put_string(data->str16var,iobuf);
   put_long_string(data->str32var,iobuf);
   put_var_string(data->strvvar,iobuf);

   put_item_end(iobuf,&item_header2);

   return(put_item_end(iobuf,&item_header1));
}

/* ---------------------- read_test3 ---------------------- */
/**
 *  @short Read test data as a nested tree
 *
 *  @param   data   Pointer to test data structure
 *  @param   iobuf  Pointer to I/O buffer
 *
 *  @return  0 (ok), <0 (error as for get_item_end())
 *
 */

int read_test3(TEST_DATA *data, IO_BUFFER *iobuf)
{
   IO_ITEM_HEADER item_header1, item_header2;
   int i;
   
   item_header1.type = 990;             /* test data */
   if ( get_item_begin(iobuf,&item_header1) < 0 )
   {
      Warning("Missing or invalid test data block.");
      return -4;
   }

   item_header2.type = 991;             /* test data */
   if ( get_item_begin(iobuf,&item_header2) < 0 )
   {
      Warning("Missing or invalid test data sub-block 1.");
      return -4;
   }

   get_vector_of_long(data->lvar,2,iobuf);
   data->ilvar[0] = get_long(iobuf);
   data->ilvar[1] = get_long(iobuf);
   get_vector_of_int(data->isvar,2,iobuf);
   get_vector_of_short(data->svar,3,iobuf);
   
   get_item_end(iobuf,&item_header2);

   if ( next_subitem_type(iobuf) != 992 )
   {
      Warning("Failed to look ahead for next sub-item type.");
      return -4;
   }

   item_header2.type = 992;             /* test data */
   if ( get_item_begin(iobuf,&item_header2) < 0 )
   {
      Warning("Missing or invalid test data sub-block 2.");
      return -4;
   }

   get_vector_of_real(data->fvar,2,iobuf);
   get_vector_of_double(data->dvar,2,iobuf);
   data->hvar[0] = get_sfloat(iobuf);
   data->hvar[1] = get_sfloat(iobuf);
   get_vector_of_byte((uint8_t *)data->i8var,2,iobuf);
   get_vector_of_byte(data->u8var,2,iobuf);
   get_vector_of_short(data->i16var,2,iobuf);
   get_vector_of_short((int16_t *)data->u16var,2,iobuf);
   get_vector_of_int32(data->i32var,2,iobuf);
   get_vector_of_uint32(data->u32var,2,iobuf);
#ifdef HAVE_64BIT_INT
   get_vector_of_int64(data->i64var,2,iobuf);
   get_vector_of_uint64(data->u64var,2,iobuf);
#endif

   get_item_end(iobuf,&item_header2);
   
   /* Check that non-sequential access to sub-items also works */

   rewind_item(iobuf,&item_header1);

   item_header2.type = 994;             /* test data */
   if ( search_sub_item(iobuf,&item_header1,&item_header2) < 0 )
   {
      Warning("Cannot find test data sub-block 4.");
      return -4;
   }
   if ( get_item_begin(iobuf,&item_header2) < 0 )
   {
      Warning("Missing or invalid test data sub-block 4.");
      return -4;
   }

   get_string(data->str16var,sizeof(data->str16var),iobuf);
   get_long_string(data->str32var,sizeof(data->str32var),iobuf);
   get_var_string(data->strvvar,sizeof(data->strvvar),iobuf);

   get_item_end(iobuf,&item_header2);

   rewind_item(iobuf,&item_header1);

   item_header2.type = 993;             /* test data */
   if ( search_sub_item(iobuf,&item_header1,&item_header2) < 0 )
   {
      Warning("Cannot find test data sub-block 3.");
      return -4;
   }
   if ( get_item_begin(iobuf,&item_header2) < 0 )
   {
      Warning("Missing or invalid test data sub-block 3.");
      return -4;
   }

   data->nbvar = get_count(iobuf);
   get_vector_of_byte(data->bvar,2,iobuf);
   for (i=0; i<4; i++)
      data->cnt16var[i] = get_count(iobuf);
   for (i=0; i<6; i++)
      data->cnt32var[i] = get_count(iobuf);
   for (i=0; i<6; i++)
      data->cntzvar[i] = get_count(iobuf);
   for (i=0; i<8; i++)
      data->cntvar[i] = get_count(iobuf);
   for (i=0; i<10; i++)
      data->scnt16var[i] = get_scount16(iobuf);
   for (i=0; i<12; i++)
      data->scnt32var[i] = get_scount32(iobuf);
   for (i=0; i<12; i++)
      data->scntzvar[i] = get_scount(iobuf);
   for (i=0; i<14; i++)
      data->scntvar[i] = get_scount(iobuf);

   get_item_end(iobuf,&item_header2);

   return(get_item_end(iobuf,&item_header1));
}

/* ---------------------- perror ------------------------- */
/**
 *  @short Replacement for function missing on OS-9
*/

#ifdef OS_OS9
int perror(text)
    char *text;
{
    fprintf(stderr,"%s: Error\n",text);
    return 0;
}
#endif

void syntax(const char *prg);

void syntax(const char *prg)
{
   fprintf(stderr,"Test basic EventIO write and read functions.\n");
   fprintf(stderr,"Syntax: %s [ -e ] filename\n", prg); 
   fprintf(stderr,"Options:\n");
   fprintf(stderr,"  -e  Use the extension field for all I/O block headers.\n");
   exit(1);
}

/* ------------------------ main ------------------------- */
/**
 *  @short Main function for I/O test program.
 *
 *  First writes a test data structure with the vector
 *  functions, then the same data structure with the
 *  single-element functions. The output file is then
 *  closed and reopened for reading. The first structure
 *  is then read with the single-element functions and
 *  the second with the vector functions (i.e. the other
 *  way as done for writing).
 *  The data from the file is compared with the original
 *  data, taking the relevant accuracy into account.
 *  Note that if an 'int' variable is written via 'put_short()'
 *  and then read again via 'get_short()' not only the
 *  upper two bytes (on a 32-bit machine) are lost but
 *  also the sign bit is propagated from bit 15 to the
 *  upper 16 bits. Similarly, if a 'long' variable is written
 *  via 'put_long()' and later read via 'get_long()' on a
 *  64-bit-machine, not only the upper 4 bytes are lost but
 *  also the sign in bit 31 is propagated to the upper 32 bits.
 */

int main (int argc, char **argv)
{
   IO_BUFFER *iobuf;
   IO_ITEM_HEADER item_header;
   FILE *output = NULL, *input;
   static TEST_DATA tdata, cdata1, cdata2;
   int ok;
   const char *program = argv[0];
   
   tdata.nbvar = 13;
   tdata.bvar[0] = (1<<2) | (1<<7);
   tdata.bvar[1] = (1<<0) | (1<<4);
   tdata.lvar[0]   = 0x01020304L;  tdata.lvar[1]  = 0xf1f2f3f4L;
   tdata.ilvar[0]  = 0x01020304;   tdata.ilvar[1] = 0xf1f2f3f4;
   tdata.isvar[0]  = 0x0102;       tdata.isvar[1] = 0xf1f2;
   tdata.svar[0]   = 0x0102; 
   tdata.svar[1]   = (short) 0xf1f2; 
   tdata.svar[2]   = 0x0a0b;
   tdata.fvar[0]   = 2.38793926059e-38;
   tdata.fvar[1]   = -2.40608939547e+30;
   tdata.dvar[0]   = 2.38793926059674673672e-140;
   tdata.dvar[1]   = -2.40608939547354636548e+180;
   tdata.hvar[0]   = -38189.;
   tdata.hvar[1]   = 0.0782;
   tdata.i8var[0]  = 0x1e;         tdata.i8var[1]  = (int8_t) 0xe1;
   tdata.u8var[0]  = 0x1e;         tdata.u8var[1]  = 0xe1;
   tdata.i16var[0] = 0x1e2e;       tdata.i16var[1] = (int16_t) 0xe2e1;
   tdata.u16var[0] = 0x1e2e;       tdata.u16var[1] = 0xe2e1;
   tdata.i32var[0] = 0x1e2e3e4e;   tdata.i32var[1] = 0xe4e3e2e1;
   tdata.u32var[0] = 0x1e2e3e4e;   tdata.u32var[1] = 0xe4e3e2e1;
#ifdef HAVE_64BIT_INT
# ifdef SIXTY_FOUR_BITS
   tdata.i64var[0] = 0x1a2a3a4a5a6a7a8aL;
   tdata.i64var[1] = 0xa8a7a6a5a4a3a2a1L;
   tdata.u64var[0] = 0x1b2b3b4b5b6b7b8bUL;
   tdata.u64var[1] = 0xb8b7b6b5b4b3b2b1UL;
# else
   tdata.i64var[0] = 0x1a2a3a4a5a6a7a8aLL;
   tdata.i64var[1] = 0xa8a7a6a5a4a3a2a1LL;
   tdata.u64var[0] = 0x1b2b3b4b5b6b7b8bULL;
   tdata.u64var[1] = 0xb8b7b6b5b4b3b2b1ULL;
# endif
#endif
   tdata.cnt16var[0] = 127;
   tdata.cnt16var[1] = 128;
   tdata.cnt16var[2] = 32768;
   tdata.cnt16var[3] = 65535;

   tdata.cnt32var[0] = 127;
   tdata.cnt32var[1] = 129;
   tdata.cnt32var[2] = 32768;
   tdata.cnt32var[3] = 65535;
   tdata.cnt32var[4] = 65536;
   tdata.cnt32var[5] = 123456788;

   tdata.cntzvar[0] = 127;
   tdata.cntzvar[1] = 128;
   tdata.cntzvar[2] = 32768;
   tdata.cntzvar[3] = 65535;
   tdata.cntzvar[4] = 65536;
   tdata.cntzvar[5] = 123456789;

   tdata.cntvar[0] = 127;
   tdata.cntvar[1] = 128;
   tdata.cntvar[2] = 32768;
   tdata.cntvar[3] = 65535;
   tdata.cntvar[4] = 65536;
   tdata.cntvar[5] = 123456789;
   tdata.cntvar[6] = (1UL << 30) + 99UL;
#ifdef HAVE_64BIT_INT
# ifdef SIXTY_FOUR_BITS
   tdata.cntvar[7] = (1UL << 62) + 99UL;
# else
   tdata.cntvar[7] = (1ULL << 62) + 99ULL;
# endif
#else
   tdata.cntvar[7] = (1UL << 31) - 1;
#endif

   tdata.scnt16var[0] = 63;
   tdata.scnt16var[1] = -63;
   tdata.scnt16var[2] = 64;
   tdata.scnt16var[3] = -64;
   tdata.scnt16var[4] = 8191;
   tdata.scnt16var[5] = -8191;
   tdata.scnt16var[6] = 8192;
   tdata.scnt16var[7] = -8192;
   tdata.scnt16var[8] = 32767;
   tdata.scnt16var[9] = -32768;

   tdata.scnt32var[0] = 63;
   tdata.scnt32var[1] = -63;
   tdata.scnt32var[2] = 64;
   tdata.scnt32var[3] = -64;
   tdata.scnt32var[4] = 8191;
   tdata.scnt32var[5] = -8191;
   tdata.scnt32var[6] = 8193;
   tdata.scnt32var[7] = -8195;
   tdata.scnt32var[8] = 32768;
   tdata.scnt32var[9] = 65536;
   tdata.scnt32var[10] = -65536;
   tdata.scnt32var[11] = 123456789;

   tdata.scntzvar[0] = 63;
   tdata.scntzvar[1] = -63;
   tdata.scntzvar[2] = 64;
   tdata.scntzvar[3] = -64;
   tdata.scntzvar[4] = 8191;
   tdata.scntzvar[5] = -8191;
   tdata.scntzvar[6] = 8192;
   tdata.scntzvar[7] = -8192;
   tdata.scntzvar[8] = 32768;
   tdata.scntzvar[9] = 65536;
   tdata.scntzvar[10] = -65536;
   tdata.scntzvar[11] = 123456789;

   tdata.scntvar[0] = 63;
   tdata.scntvar[1] = -63;
   tdata.scntvar[2] = 64;
   tdata.scntvar[3] = -64;
   tdata.scntvar[4] = 8191;
   tdata.scntvar[5] = -8191;
   tdata.scntvar[6] = 8192;
   tdata.scntvar[7] = -8192;
   tdata.scntvar[8] = 32768;
   tdata.scntvar[9] = 65536;
   tdata.scntvar[10] = -65536;
   tdata.scntvar[11] = 123456789;
   tdata.scntvar[12] = -(1L << 30) - 99L;
#ifdef HAVE_64BIT_INT
# ifdef SIXTY_FOUR_BITS
   tdata.scntvar[13] = (1L << 61) + 99L;
# else
   tdata.scntvar[13] = (1LL << 61) + 99LL;
# endif
#else
   tdata.scntvar[13] = (1L << 31) - 1;
#endif
   strcpy(tdata.str16var,"some text");
   strcpy(tdata.str32var,"Some text");
   strcpy(tdata.strvvar,"Some Text");

   if ( (iobuf = allocate_io_buffer((size_t)1000)) == (IO_BUFFER *) NULL )
      exit(1);
   if ( argc > 1 && strcmp(argv[1],"--help") == 0 )
      syntax(program);
   if ( argc > 1 && strcmp(argv[1],"-e") == 0 )
   {
      iobuf->extended = 1;
      Information("Using the extension field for all I/O block headers.");
      argc--;
      argv++;
   }
   if ( argc > 1 )
   {
      if ( (output = fileopen(argv[1],WRITE_BINARY)) == (FILE *) NULL )
      {
         perror(argv[1]);
         exit(1);
      }
      iobuf->output_file = output;
   }
   else
      syntax(program);

   fprintf(stderr,"\nWrite test data to file '%s'.\n",argv[1]);
   fprintf(stderr,"Default byte order, using mainly vector functions.\n");
   write_test2(&tdata,iobuf);
   fprintf(stderr,"Default byte order, using single-element functions.\n");
   write_test1(&tdata,iobuf);
   iobuf->byte_order = 1;
   fprintf(stderr,"Reversed byte order, using single-element functions.\n");
   write_test1(&tdata,iobuf);
   iobuf->byte_order = 0;
   fprintf(stderr,"Normal byte order, using single-element functions.\n");
   write_test1(&tdata,iobuf);
   fprintf(stderr,"Writing as nested item structure.\n");
   write_test3(&tdata,iobuf);
   fprintf(stderr,"Write tests done.\n\n");
   
   fileclose(output);
   iobuf->output_file = output = NULL;
   if ( argc > 2 )
   {
      fprintf(stderr,"\nComparing against test data from external file %s:\n",
         argv[2]);
      if ( (input = fileopen(argv[2],READ_BINARY)) == (FILE *) NULL )
      {
         perror(argv[2]);
         exit(1);
      }
   }
   else if ( (input = fileopen(argv[1],READ_BINARY)) == (FILE *) NULL )
   {
      perror(argv[1]);
      exit(1);
   }
   else
      fprintf(stderr,"Read test data from file '%s'.\n",argv[1]);   

   iobuf->input_file = input;
   ok = 1;

   fprintf(stderr,"Default byte order, using single-element functions.\n");
   memset(&cdata1,0,sizeof(cdata1));
   if ( find_io_block(iobuf,&item_header) < 0 )
   {
      Error("*** Finding I/O block 1 failed");
      exit(1);
   }
   if ( read_io_block(iobuf,&item_header) < 0 )
   {
      Error("*** Reading I/O block 1 failed");
      exit(1);
   }
   if ( read_test1(&cdata1,iobuf) < 0 )
   {
      Error("*** Read test 1 failed");
      exit(1);
   }
   if ( datacmp(&tdata,&cdata1) != 1 )
   {
      Error("*** Data from read test 1 does not match.");
      ok = 0;
   }
   
   fprintf(stderr,"Default byte order, using mainly vector functions.\n");
   memset(&cdata2,0,sizeof(cdata2));
   if ( find_io_block(iobuf,&item_header) < 0 )
   {
      Error("*** Finding I/O block 1 failed");
      exit(1);
   }
   if ( read_io_block(iobuf,&item_header) < 0 )
   {
      Error("*** Reading I/O block 1 failed");
      exit(1);
   }
   if ( read_test2(&cdata2,iobuf) < 0 )
   {
      Error("*** Read test 2 failed");
      exit(1);
   }
   if ( datacmp(&tdata,&cdata2) != 1 )
   {
      Error("*** Data from read test 2 does not match");
      ok = 0;
   }

   fprintf(stderr,"Reversed byte order, using single-element functions.\n");
   memset(&cdata1,0,sizeof(cdata1));
   if ( find_io_block(iobuf,&item_header) < 0 )
   {
      Error("*** Finding I/O block 3 failed");
      exit(1);
   }
   if ( read_io_block(iobuf,&item_header) < 0 )
   {
      Error("*** Reading I/O block 3 failed");
      exit(1);
   }
   if ( read_test1(&cdata1,iobuf) < 0 )
   {
      Error("*** Read test 3 failed");
      exit(1);
   }
   if ( datacmp(&tdata,&cdata1) != 1 )
   {
      Error("*** Data from read test 3 does not match.");
      ok = 0;
   }
   
   fprintf(stderr,"Normal byte order, using single-element functions.\n");
   memset(&cdata2,0,sizeof(cdata2));
   if ( find_io_block(iobuf,&item_header) < 0 )
   {
      Error("*** Finding I/O block 4 failed");
      exit(1);
   }
   if ( read_io_block(iobuf,&item_header) < 0 )
   {
      Error("*** Reading I/O block 4 failed");
      exit(1);
   }
   if ( read_test1(&cdata2,iobuf) < 0 )
   {
      Error("*** Read test 4 failed");
      exit(1);
   }
   if ( datacmp(&tdata,&cdata2) != 1 )
   {
      Error("*** Data from read test 4 does not match");
      ok = 0;
   }
   
   fprintf(stderr,"Reading nested data structure.\n");
   memset(&cdata2,0,sizeof(cdata2));
   if ( find_io_block(iobuf,&item_header) < 0 )
   {
      Error("*** Finding I/O block 5 failed");
      exit(1);
   }
   if ( read_io_block(iobuf,&item_header) < 0 )
   {
      Error("*** Reading I/O block 5 failed");
      exit(1);
   }
   if ( read_test3(&cdata2,iobuf) < 0 )
   {
      Error("*** Read test 5 failed");
      exit(1);
   }
   if ( datacmp(&tdata,&cdata2) != 1 )
   {
      Error("*** Data from read test 5 does not match");
      ok = 0;
   }
   
   Information("Read tests done\n");
   
   if ( ok )
      Information("Everything is ok. Congratulations!\n");
      
   if ( care_long )
   {
      Information("Note: on this machine you should care about the sign propagation");
      Information("of 'LONG' (INT32) data elements to long integer variables.\n");
   }
   
   if ( care_int )
   {
      Information("Note: on this machine you should care about the sign propagation");
      Information("of 'LONG' (INT32) data elements to 32 bit integer variables.\n");
   }
   
   if ( care_short )
   {
      Information("Note: on this machine you should care about the sign propagation");
      Information("of 'SHORT' data elements to integer or long integer variables.\n");
   }
   
#ifdef HAVE_64BIT_INT
   Information("On this machine you can read and write 64-bit integers but you should");
   Information("be aware that this feature is not available on all systems otherwise");
   Information("supported by eventio.");
#ifdef SIXTY_FOUR_BITS
   Information("These 64-bit integers are native types.\n");
#else
   Information("These 64-bit integers are implemented through the C compiler.\n");
#endif
#else
   Information("On this system no 64-bit integers are supported.\n");
#endif

#ifdef WITH_INTMAX_T
   {
      char stmp[200];
      sprintf(stmp,"The longest integers available are %zu bits long.",
         8*sizeof(intmax_t));
      Information(stmp);
   }
#endif
#ifdef WITH_UINTMAX_T
   {
      char stmp[200];
      sprintf(stmp,"The longest unsigned integers available are %zu bits long.\n",
         8*sizeof(uintmax_t));
      Information(stmp);
   }
#endif
#if defined(HAVE_EVENTIO_USER_FLAG) || \
    defined(HAVE_EVENTIO_EXTENDED_LENGTH) || \
    defined(HAVE_EVENTIO_HEADER_LENGTH)
   {
      char stmp[1000];
      sprintf(stmp,"This version of eventio has additional features:\n");
#ifdef HAVE_EVENTIO_USER_FLAG
      strcat(stmp," - User flag bit in the header.\n");
#endif
#ifdef HAVE_EVENTIO_HEADER_LENGTH
      strcat(stmp," - The length of the data field can be seen in the header.\n");
#endif
#ifdef HAVE_EVENTIO_EXTENDED_LENGTH
      if ( sizeof(size_t) == 4 )
         strcat(stmp,
         " - Extended length I/O blocks (only up to 2 GByte on 32-bit-systems).\n");
      else if ( sizeof(size_t) >= 8 )
      {
         strcat(stmp,
         " - Extended length I/O blocks (up to 4 TByte).\n");
         strcat(stmp,
         "   Note: on 32-bit systems only blocks up to 2 GByte are readable.\n");
      }
      strcat(stmp,
         "   Note: Extended length I/O blocks are not readable with\n"
         "   very old eventio versions.\n");
#endif
      Information(stmp);
   }
#endif

   return 0;
}

/** @} */
