/* ============================================================================

Copyright (C) 2003, 2007, 2012 Konrad Bernloehr (Konrad Bernl&ouml;hr)

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

/** @file TestIO.cc
    @short Test program for eventio data format (based on testio.c)

    This file is a re-implementation in C++ that should 
    produce the same output (both on screen and in data file)
    as the original implementation in testio.c.
    Therefore, this is not intended as a masterpiece in object-oriented
    programming but a rather straight-forward C-to-C++ translation.

    The data files produced both in standard and extended [-e] mode
    should be identical to the corresponding data files produced with
    the testio tool. Comparison of the files may serve as an additional test.

    @author  Konrad Bernloehr
    $Date: 2014/11/17 15:44:00 $
    $Revision: 1.25 $
*/

/** @defgroup TestIO_cc The TestIO program */
/** @{ */
#define __STDC_LIMIT_MACROS 1
#include <vector>
#include <valarray>
#include <string>
using std::vector;
#include <iostream>
#include "fileopen.h"
#include "EventIO.hh"

using namespace eventio;

struct test_struct
{
   long lvar[2];
   int ilvar[2];
   int isvar[2];
   short svar[3];
   double fvar[2];
   double dvar[2];
   double hvar[2];
   vector<bool> bvar;
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
   
   void Print(void) const;
};

typedef struct test_struct TEST_DATA;

#define SA(a) std::cout << #a << ": " << a[0]; \
 for (size_t i=1; i<sizeof(a)/sizeof(a[0]); i++) std::cout << ", " << a[i]; \
 std::cout << "\n"
#define SB(a) std::cout << #a << ": " << (int) a[0]; \
 for (size_t i=1; i<sizeof(a)/sizeof(a[0]); i++) std::cout << ", " << (int) a[i]; \
 std::cout << "\n"
#define SS(a) std::cout << #a << ": " << a << "\n";
#define SV(a) std::cout << #a << ": " << a[0]; \
 for (size_t i=1; i<a.size(); i++) std::cout << ", " << a[i]; \
 std::cout << "\n"
 
void test_struct::Print(void) const
{
   SA(lvar);
   SA(ilvar);
   SA(isvar);
   SA(svar);
   SA(fvar);
   SA(dvar);
   SA(hvar);
   SV(bvar);
   SB(i8var);
   SB(u8var);
   SA(i16var);
   SA(u16var);
   SA(i32var);
   SA(u32var);
#ifdef HAVE_64BIT_INT
   SA(i64var);
   SA(u64var);
#endif
   SA(cnt16var);
   SA(scnt16var);
   SA(cnt32var);
   SA(scnt32var);
   SA(cntzvar);
   SA(scntzvar);
   SA(cntvar);
   SA(scntvar);
   SS(str16var);
   SS(str32var);
   SS(strvvar);
}

static int care_long, care_int, care_short;
static bool show_data = false;

int datacmp (TEST_DATA& data1, TEST_DATA& data2);
int write_test1 (TEST_DATA& data, EventIO& iobuf);
int read_test1 (TEST_DATA& data, EventIO& iobuf);
int write_test2 (TEST_DATA& data, EventIO& iobuf);
int read_test2 (TEST_DATA& data, EventIO& iiobuf);

#ifndef WITH_INLINE_WARNING
void Information(const char *text)
{
   std::cerr << text << std::endl;
}
void Warning(const char *text)
{
   std::cerr << text << std::endl;
}
void Error(const char *text)
{
   std::cerr << text << std::endl;
}
#endif

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

int datacmp(TEST_DATA& data1, TEST_DATA& data2)
{
   if ( show_data )
   {
      std::cout << "------------------- data1 -----------------\n\n";
      data1.Print();
      std::cout << "\n";
      std::cout << "------------------- data2 -----------------\n\n";
      data2.Print();
      std::cout << "\n\n";
   }

   size_t i, ok;
   
   ok = 1;
   
   if ( data1.bvar.size() != data2.bvar.size() )
   {
      fprintf(stderr,"Bool variables do not match.\n");
      ok = 0;
   }
   else
      for (i=0; i<data1.bvar.size(); i++)
      {
         if ( data1.bvar[i] != data2.bvar[i] )
         {
            fprintf(stderr,"Bool variable %zu does not match: %d <--> %d\n",
               i+1,int(data1.bvar[i]),int(data2.bvar[i]));
            ok = 0;
         }
      }

   for (i=0; i<2; i++)
      if ( (int32_t)data1.lvar[i] != (int32_t)data2.lvar[i] )
      {
         fprintf(stderr,"Long variable %zu does not match: %08lx <--> %08lx\n",
            i+1,data1.lvar[i],data2.lvar[i]);
         ok = 0;
      }
      else if ( data1.lvar[i] != data2.lvar[i] )
         care_long = 1;
      
   for (i=0; i<2; i++)
      if ( (int32_t)data1.ilvar[i] != (int32_t)data2.ilvar[i] )
      {
         fprintf(stderr,"Int 'l' variable %zu does not match: %08x <--> %08x\n",
            i+1,data1.ilvar[i],data2.ilvar[i]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( (int16_t)data1.isvar[i] != (int16_t)data2.isvar[i] )
      {
         fprintf(stderr,"Int 's' variable %zu does not match: %04x <--> %04x\n",
            i+1,data1.isvar[i],data2.isvar[i]);
         ok = 0;
      }
      else if ( data1.isvar[i] != data2.isvar[i] )
         care_short = 1;
      
   for (i=0; i<3; i++)
      if ( data1.svar[i] != data2.svar[i] )
      {
         fprintf(stderr,"Short variable %zu does not match: %04x <--> %04x\n",
            i+1,data1.svar[i],data2.svar[i]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
   {
      union f_i { float f; uint32_t i; } f1, f2;
      f1.f = (float) data1.fvar[i];
      f2.f = (float) data2.fvar[i];
      if ( f1.f != f2.f )
      {
         fprintf(stderr,"Float variable %zu does not match: %08lx <--> %08lx\n",
            i+1, 
            static_cast<unsigned long>(f1.i), 
            static_cast<unsigned long>(f2.i));
         ok = 0;
      }
   }
   
   for (i=0; i<2; i++)
      if ( data1.dvar[i] != data2.dvar[i] )
      {
         union d_i { double d; uint32_t i[2]; } f1, f2;
	 f1.d = data1.dvar[i];
	 f2.d = data2.dvar[i];
         fprintf(stderr,
	    "Double variable %zu does not match: %08lx%08lx <--> %08lx%08lx\n",
            i+1, 
            static_cast<unsigned long>(f1.i[0]), 
            static_cast<unsigned long>(f1.i[1]),
            static_cast<unsigned long>(f2.i[0]), 
            static_cast<unsigned long>(f2.i[1]));
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1.hvar[i] != data2.hvar[i] && fabs(data1.hvar[i]-data2.hvar[i]) > 0.002*fabs(data1.hvar[i]) )
      {
         union d_i { double d; uint32_t i[2]; } f1, f2;
	 f1.d = data1.dvar[i];
	 f2.d = data2.dvar[i];
         fprintf(stderr,
	    "Sfloat variable %zu does not match: %08lx%08lx <--> %08lx%08lx\n",
            i+1, 
            static_cast<unsigned long>(f1.i[0]), 
            static_cast<unsigned long>(f1.i[1]),
            static_cast<unsigned long>(f2.i[0]), 
            static_cast<unsigned long>(f2.i[1]));
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1.i8var[i] != data2.i8var[i] )
      {
         fprintf(stderr,"Int8 variable %zu does not match: %08x <--> %08x\n",
            i+1,data1.i8var[i],data2.i8var[i]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1.u8var[i] != data2.u8var[i] )
      {
         fprintf(stderr,"UInt8 variable %zu does not match: %08x <--> %08x\n",
            i+1,data1.u8var[i],data2.u8var[i]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1.i16var[i] != data2.i16var[i] )
      {
         fprintf(stderr,"Int16 variable %zu does not match: %08x <--> %08x\n",
            i+1,data1.i16var[i],data2.i16var[i]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1.u16var[i] != data2.u16var[i] )
      {
         fprintf(stderr,"UInt16 variable %zu does not match: %08x <--> %08x\n",
            i+1,data1.u16var[i],data2.u16var[i]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1.i32var[i] != data2.i32var[i] )
      {
         fprintf(stderr,"Int32 variable %zu does not match: %08x <--> %08x\n",
            i+1,data1.i32var[i],data2.i32var[i]);
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1.u32var[i] != data2.u32var[i] )
      {
         fprintf(stderr,"UInt32 variable %zu does not match: %08x <--> %08x\n",
            i+1,data1.u32var[i],data2.u32var[i]);
         ok = 0;
      }

#ifdef HAVE_64BIT_INT
   for (i=0; i<2; i++)
      if ( data1.i64var[i] != data2.i64var[i] )
      {
#ifdef SIXTY_FOUR_BITS
         fprintf(stderr,"Int64 variable %zu does not match: %016lx <--> %016lx\n",
            i+1,data1.i64var[i],data2.i64var[i]);
#else
         fprintf(stderr,"Int64 variable %zu does not match: %016llx <--> %016llx\n",
            i+1,data1.i64var[i],data2.i64var[i]);
#endif
         ok = 0;
      }
      
   for (i=0; i<2; i++)
      if ( data1.u64var[i] != data2.u64var[i] )
      {
#ifdef SIXTY_FOUR_BITS
         fprintf(stderr,"UInt64 variable %zu does not match: %016lx <--> %016lx\n",
            i+1,data1.u64var[i],data2.u64var[i]);
#else
         fprintf(stderr,"UInt64 variable %zu does not match: %016llx <--> %016llx\n",
            i+1,data1.u64var[i],data2.u64var[i]);
#endif
         ok = 0;
      }
#endif
      
   for (i=0; i<4; i++)
      if ( data1.cnt16var[i] != data2.cnt16var[i] )
      {
         fprintf(stderr,
            "Count16 variable %u does not match: %08x <--> %08x\n",
            (unsigned int)i+1, data1.cnt16var[i], data2.cnt16var[i]);
         ok = 0;
      }

   for (i=0; i<6; i++)
      if ( data1.cnt32var[i] != data2.cnt32var[i] )
      {
         fprintf(stderr,
            "Count32 variable %u does not match: %08x <--> %08x\n",
            (unsigned int)i+1, data1.cnt32var[i], data2.cnt32var[i]);
         ok = 0;
      }

   for (i=0; i<6; i++)
      if ( data1.cntzvar[i] != data2.cntzvar[i] )
      {
         fprintf(stderr,
#ifdef WITH_UINTMAX_T
            "Count variable %u does not match: %08jx <--> %08jx\n",
#else
            "Count variable %u does not match: %08zx <--> %08zx\n",
#endif
            (unsigned int)i+1, data1.cntvar[i], data2.cntvar[i]);
         ok = 0;
      }

   for (i=0; i<8; i++)
      if ( data1.cntvar[i] != data2.cntvar[i] )
      {
         fprintf(stderr,
#if defined(WITH_UINTMAX_T)
            "Long count variable %zu does not match: %08jx <--> %08jx\n",
#else
            "Long count variable %zu does not match: %08zx <--> %08zx\n",
#endif
            i+1, data1.cntvar[i], data2.cntvar[i]);
         ok = 0;
      }

   for (i=0; i<10; i++)
      if ( data1.scnt16var[i] != data2.scnt16var[i] )
      {
         fprintf(stderr,
            "Signed count16 variable %u does not match: %08x <--> %08x\n",
            (unsigned int)i+1, data1.scnt16var[i], data2.scnt16var[i]);
         ok = 0;
      }

   for (i=0; i<12; i++)
      if ( data1.scnt32var[i] != data2.scnt32var[i] )
      {
         fprintf(stderr,
            "Signed count32 variable %u does not match: %08x <--> %08x\n",
            (unsigned int)i+1, data1.scnt32var[i], data2.scnt32var[i]);
         ok = 0;
      }

   for (i=0; i<12; i++)
      if ( data1.scntzvar[i] != data2.scntzvar[i] )
      {
         fprintf(stderr,
            "Signed count variable %u does not match: %08zx <--> %08zx\n",
            (unsigned int)i+1, data1.scntzvar[i], data2.scntzvar[i]);
         ok = 0;
      }

   for (i=0; i<14; i++)
      if ( data1.scntvar[i] != data2.scntvar[i] )
      {
         fprintf(stderr,
#if defined(WITH_UINTMAX_T)
            "Signed long count variable %u does not match: %08jx <--> %08jx\n",
#else
            "Signed long count variable %u does not match: %08zx <--> %08zx\n",
#endif
            (unsigned int)i+1, data1.cntvar[i], data2.cntvar[i]);
         ok = 0;
      }

   if ( strcmp(data1.str16var,data2.str16var) != 0 )
   {
      fprintf(stderr,"Strings with 16-bit count do not match: '%s' <--> '%s'\n",
         data1.str16var, data2.str16var);
         ok = 0;
   }
   if ( strcmp(data1.str32var,data2.str32var) != 0 )
   {
      fprintf(stderr,"Strings with 16-bit count do not match: '%s' <--> '%s'\n",
         data1.str32var, data2.str32var);
         ok = 0;
   }
   if ( strcmp(data1.strvvar,data2.strvvar) != 0 )
   {
      fprintf(stderr,"Strings do not match: '%s' <--> '%s'\n",
         data1.strvvar, data2.strvvar);
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

int write_test1(TEST_DATA& data, EventIO& iobuf)
{
   EventIO::Item item(iobuf,"put",99,0,123);
   if ( item.Status() != 0 )
      return item.Status();
   
   item.PutInt32(data.lvar[0]);
   item.PutInt32(data.lvar[1]);
   item.PutInt32((long)data.ilvar[0]);
   item.PutInt32((long)data.ilvar[1]);
   item.PutInt16(data.isvar[0]);
   item.PutInt16(data.isvar[1]);
   item.PutInt16(data.svar[0]);
   item.PutInt16(data.svar[1]);
   item.PutInt16(data.svar[2]);
   item.PutReal(data.fvar[0]);
   item.PutReal(data.fvar[1]);
   item.PutDouble(data.dvar[0]);
   item.PutDouble(data.dvar[1]);
   item.PutSfloat(data.hvar[0]);
   item.PutSfloat(data.hvar[1]);
   item.PutUint8(data.i8var[0]);
   item.PutUint8(data.i8var[1]);
   item.PutUint8(data.u8var[0]);
   item.PutUint8(data.u8var[1]);
   item.PutInt16(data.i16var[0]);
   item.PutInt16(data.i16var[1]);
   item.PutInt16(data.u16var[0]);
   item.PutInt16(data.u16var[1]);
   item.PutInt32(data.i32var[0]);
   item.PutInt32(data.i32var[1]);
   item.PutUint32(data.u32var[0]);
   item.PutUint32(data.u32var[1]);
#ifdef HAVE_64BIT_INT
   item.PutInt64(&data.i64var[0],1);
   item.PutInt64(&data.i64var[1],1);
   item.PutUint64(&data.u64var[0],1);
   item.PutUint64(&data.u64var[1],1);
#endif
   item.PutBool(data.bvar);
   size_t i;
   for (i=0; i<4; i++)
      item.PutCount16(data.cnt16var[i]);
   for (i=0; i<6; i++)
      item.PutCount32(data.cnt32var[i]);
   for (i=0; i<6; i++)
      item.PutCount(data.cntzvar[i]);
   for (i=0; i<8; i++)
      item.PutCount(data.cntvar[i]);
   for (i=0; i<10; i++)
      item.PutSCount16(data.scnt16var[i]);
   for (i=0; i<12; i++)
      item.PutSCount32(data.scnt32var[i]);
   for (i=0; i<12; i++)
      item.PutSCount(data.scntzvar[i]);
   for (i=0; i<14; i++)
      item.PutSCount(data.scntvar[i]);
   item.PutUint16(strlen(data.str16var));
   for (i=0; i<strlen(data.str16var); i++)
      item.PutUint8(data.str16var[i]);
   item.PutUint32((uint32_t)strlen(data.str32var));
   for (i=0; i<strlen(data.str32var); i++)
      item.PutUint8(data.str32var[i]);
   item.PutString(data.strvvar);

   item.Done(); 
   return (item.Status()==1?0:item.Status());
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

int read_test1(TEST_DATA& data, EventIO& iobuf)
{
   EventIO::Item item(iobuf,"get",99);
   if ( item.Status() < 0 )
   {
      Warning("Missing or invalid test data block.");
      return -4;
   }

   data.lvar[0]   = item.GetInt32();
   data.lvar[1]   = item.GetInt32();
   data.ilvar[0]  = item.GetInt32();
   data.ilvar[1]  = item.GetInt32();
   data.isvar[0]  = item.GetInt16();
   data.isvar[1]  = item.GetInt16();
   data.svar[0]   = item.GetInt16();
   data.svar[1]   = item.GetInt16();
   data.svar[2]   = item.GetInt16();
   data.fvar[0]   = item.GetReal();
   data.fvar[1]   = item.GetReal();
   data.dvar[0]   = item.GetDouble();
   data.dvar[1]   = item.GetDouble();
   data.hvar[0]   = item.GetSfloat();
   data.hvar[1]   = item.GetSfloat();
   data.i8var[0]  = item.GetUint8();
   data.i8var[1]  = item.GetUint8();
   data.u8var[0]  = item.GetUint8();
   data.u8var[1]  = item.GetUint8();
   data.i16var[0] = item.GetInt16();
   data.i16var[1] = item.GetInt16();
   data.u16var[0] = item.GetInt16();
   data.u16var[1] = item.GetInt16();
   data.i32var[0] = item.GetInt32();
   data.i32var[1] = item.GetInt32();
   data.u32var[0] = item.GetUint32();
   data.u32var[1] = item.GetUint32();
#ifdef HAVE_64BIT_INT
   item.GetInt64(&data.i64var[0],1);
   item.GetInt64(&data.i64var[1],1);
   item.GetUint64(&data.u64var[0],1);
   item.GetUint64(&data.u64var[1],1);
#endif
   item.GetBool(data.bvar);
   size_t i; 
   for (i=0; i<4; i++)
      data.cnt16var[i] = item.GetCount16();
   for (i=0; i<6; i++)
      data.cnt32var[i] = item.GetCount32();
   for (i=0; i<6; i++)
      data.cntzvar[i] = item.GetCount();
   for (i=0; i<8; i++)
      data.cntvar[i] = item.GetCount();
   for (i=0; i<10; i++)
      data.scnt16var[i] = item.GetSCount16();
   for (i=0; i<12; i++)
      data.scnt32var[i] = item.GetSCount32();
   for (i=0; i<12; i++)
      data.scntzvar[i] = item.GetSCount();
   for (i=0; i<14; i++)
      data.scntvar[i] = item.GetSCount();
   size_t n16 = item.GetUint16();
   for (i=0; i<n16 && i+1<sizeof(data.str16var); i++)
      data.str16var[i] = item.GetUint8();
   data.str16var[i] = '\0';
   for ( ; i<n16; i++)
      (void) item.GetUint8();
   size_t n32 = item.GetUint32();
   for (i=0; i<n32 && i+1<sizeof(data.str32var); i++)
      data.str32var[i] = item.GetUint8();
   data.str32var[i] = '\0';
   for ( ; i<n32; i++)
      (void) item.GetUint8();
   item.GetString(data.strvvar,sizeof(data.strvvar));

   item.Done();
   return (item.Status()==1?0:item.Status());
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

int write_test2(TEST_DATA& data, EventIO& iobuf)
{
   EventIO::Item item(iobuf,"put",99,0,123);
   if ( item.Status() != 0 )
      return item.Status();
   
   item.PutInt32(data.lvar,2);
   item.PutInt32((long)data.ilvar[0]);
   item.PutInt32((long)data.ilvar[1]);
   item.PutInt16(data.isvar,2);
   item.PutInt16(data.svar,3);
   item.PutReal(data.fvar,2);
   item.PutDouble(data.dvar,2);
   item.PutSfloat(data.hvar,2);
   item.PutUint8((uint8_t *)data.i8var,2);
   item.PutUint8(data.u8var,2);
   item.PutInt16(data.i16var,2);
   item.PutInt16((int16_t *)data.u16var,2);
   item.PutInt32(data.i32var,2);
   item.PutUint32(data.u32var,2);
#ifdef HAVE_64BIT_INT
   item.PutInt64(data.i64var,2);
   item.PutUint64(data.u64var,2);
#endif   
   item.PutBool(data.bvar);
   item.PutCount(data.cnt16var,4);
   item.PutCount(data.cnt32var,6);
   item.PutCount(data.cntzvar,6);
   item.PutCount(data.cntvar,8);
   item.PutSCount(data.scnt16var,10);
   item.PutSCount(data.scnt32var,12);
   item.PutSCount(data.scntzvar,12);
   item.PutSCount(data.scntvar,14);
   item.PutUint16(strlen(data.str16var));
   item.PutUint8((uint8_t *)data.str16var,strlen(data.str16var));
   item.PutUint32((uint32_t)strlen(data.str32var));
   item.PutUint8((uint8_t *)data.str32var,strlen(data.str32var));
   item.PutString(data.strvvar);
   
   item.Done();
   return (item.Status()==1?0:item.Status());
}

/* ---------------------- read_test2 ---------------------- */
/**
 *  @short Read test data with vector functions as far as possible
 *
 *  @param   data   Test data structure
 *  @param   iobuf  I/O buffer
 *
 *  @return  0 (ok), <0 (error as for get_item_end())
 *
 */

int read_test2(TEST_DATA& data, EventIO& iobuf)
{
   EventIO::Item item(iobuf,"get",99);
   if ( item.Status() < 0 )
   {
      Warning("Missing or invalid test data block.");
      return -4;
   }

   item.GetInt32(data.lvar,2);
   data.ilvar[0] = item.GetInt32();
   data.ilvar[1] = item.GetInt32();
   item.GetInt16(data.isvar,2);
   item.GetInt16(data.svar,3);
   item.GetReal(data.fvar,2);
   item.GetDouble(data.dvar,2);
   item.GetSfloat(data.hvar,2);
   item.GetUint8((uint8_t *)data.i8var,2);
   item.GetUint8(data.u8var,2);
   item.GetInt16(data.i16var,2);
   item.GetInt16((int16_t *)data.u16var,2);
   item.GetInt32(data.i32var,2);
   item.GetUint32(data.u32var,2);
#ifdef HAVE_64BIT_INT
   item.GetInt64(data.i64var,2);
   item.GetUint64(data.u64var,2);
#endif   
   item.GetBool(data.bvar);
   item.GetCount(data.cnt16var,4);
   item.GetCount(data.cnt32var,6);
   item.GetCount(data.cntzvar,6);
   item.GetCount(data.cntvar,8);
   item.GetSCount(data.scnt16var,10);
   item.GetSCount(data.scnt32var,12);
   item.GetSCount(data.scntzvar,12);
   item.GetSCount(data.scntvar,14);
   size_t i;
   size_t n16 = item.GetUint16();
   for (i=0; i<n16 && i+1<sizeof(data.str16var); i++)
      data.str16var[i] = item.GetUint8();
   data.str16var[i] = '\0';
   for ( ; i<n16; i++)
      (void) item.GetUint8();
   size_t n32 = item.GetUint32();
   for (i=0; i<n32 && i+1<sizeof(data.str32var); i++)
      data.str32var[i] = item.GetUint8();
   data.str32var[i] = '\0';
   for ( ; i<n32; i++)
      (void) item.GetUint8();
   item.GetString(data.strvvar,sizeof(data.strvvar));

   item.Done();
   return (item.Status()==1?0:item.Status());
}

/* --------------------- write_test3 ---------------------- */
/**
 *  @short Write test data in nested items
 *
 *  @param    data    Test data structure
 *  @param    iobuf   I/O buffer
 *
 *  @return  0 (ok), <0 (error as for put_item_end())
 *
 */

int write_test3(TEST_DATA& data, EventIO& iobuf)
{
   EventIO::Item item1(iobuf,"put",990,11,1);
   if ( item1.Status() != 0 )
      return item1.Status();

   EventIO::Item item2(item1,"put",991,12,1);
   if ( item2.Status() != 0 )
      return item2.Status();

   item2.PutInt32(data.lvar,2);
   item2.PutInt32((long)data.ilvar[0]);
   item2.PutInt32((long)data.ilvar[1]);
   item2.PutInt16(data.isvar,2);
   item2.PutInt16(data.svar,3);
   
   item2.Done();

   EventIO::Item item3(item1,"put",992,12,2);
   if ( item3.Status() != 0 )
      return item3.Status();

   item3.PutReal(data.fvar,2);
   item3.PutDouble(data.dvar,2);
   item3.PutSfloat(data.hvar,2);
   item3.PutUint8((uint8_t *)data.i8var,2);
   item3.PutUint8(data.u8var,2);
   item3.PutInt16(data.i16var,2);
   item3.PutInt16((int16_t *)data.u16var,2);
   item3.PutInt32(data.i32var,2);
   item3.PutUint32(data.u32var,2);
#ifdef HAVE_64BIT_INT
   item3.PutInt64(data.i64var,2);
   item3.PutUint64(data.u64var,2);
#endif   

   item3.Done();

   EventIO::Item item4(item1,"put",993,12,3);
   if ( item4.Status() != 0 )
      return item4.Status();

   item4.PutBool(data.bvar);
   item4.PutCount(data.cnt16var,4);
   item4.PutCount(data.cnt32var,6);
   item4.PutCount(data.cntzvar,6);
   item4.PutCount(data.cntvar,8);
   item4.PutSCount(data.scnt16var,10);
   item4.PutSCount(data.scnt32var,12);
   item4.PutSCount(data.scntzvar,12);
   item4.PutSCount(data.scntvar,14);

   item4.Done();

   EventIO::Item item5(item1,"put",994,12,4);
   if ( item5.Status() != 0 )
      return item5.Status();

   item5.PutUint16(strlen(data.str16var));
   item5.PutUint8((uint8_t *)data.str16var,strlen(data.str16var));
   item5.PutUint32((uint32_t)strlen(data.str32var));
   item5.PutUint8((uint8_t *)data.str32var,strlen(data.str32var));
   item5.PutString(data.strvvar);

   item5.Done();
   
   item1.Done();
   return (item1.Status()==1?0:item1.Status());
}

/* ---------------------- read_test3 ---------------------- */
/**
 *  @short Read test data as a nested tree
 *
 *  @param   data   Test data structure
 *  @param   iobuf  I/O buffer
 *
 *  @return  0 (ok), <0 (error as for get_item_end())
 *
 */

int read_test3(TEST_DATA& data, EventIO& iobuf)
{
   EventIO::Item item1(iobuf,"get",990);
   if ( item1.Status() < 0 )
   {
      Warning("Missing or invalid test data block.");
      return -4;
   }

   EventIO::Item item2(item1,"get",991);
   if ( item2.Status() < 0 )
   {
      Warning("Missing or invalid test data sub-block 1.");
      return -4;
   }

   item2.GetInt32(data.lvar,2);
   data.ilvar[0] = item2.GetInt32();
   data.ilvar[1] = item2.GetInt32();
   item2.GetInt16(data.isvar,2);
   item2.GetInt16(data.svar,3);
   
   item2.Done();

   if ( item1.NextSubItemType() != 992 )
   {
      Warning("Failed to look ahead for next sub-item type.");
      return -4;
   }

   EventIO::Item item3(item1,"get",992);
   if ( item3.Status() < 0 )
   {
      Warning("Missing or invalid test data sub-block 2.");
      return -4;
   }

   item3.GetReal(data.fvar,2);
   item3.GetDouble(data.dvar,2);
   item3.GetSfloat(data.hvar,2);
   item3.GetUint8((uint8_t *)data.i8var,2);
   item3.GetUint8(data.u8var,2);
   item3.GetInt16(data.i16var,2);
   item3.GetInt16((int16_t *)data.u16var,2);
   item3.GetInt32(data.i32var,2);
   item3.GetUint32(data.u32var,2);
#ifdef HAVE_64BIT_INT
   item3.GetInt64(data.i64var,2);
   item3.GetUint64(data.u64var,2);
#endif   

   item3.Done();
   
   /* Check that non-sequential access to sub-items also works */

   item1.Rewind();

   if ( item1.Search(994) < 0 )
   {
      Warning("Cannot find test data sub-block 4.");
      return -4;
   }

   EventIO::Item item5(item1,"get",994);
   if ( item5.Status() < 0 )
   {
      Warning("Missing or invalid test data sub-block 4.");
      return -4;
   }

   size_t i;
   size_t n16 = item5.GetUint16();
   for (i=0; i<n16 && i+1<sizeof(data.str16var); i++)
      data.str16var[i] = item5.GetUint8();
   data.str16var[i] = '\0';
   for ( ; i<n16; i++)
      (void) item5.GetUint8();
   size_t n32 = item5.GetUint32();
   for (i=0; i<n32 && i+1<sizeof(data.str32var); i++)
      data.str32var[i] = item5.GetUint8();
   data.str32var[i] = '\0';
   for ( ; i<n32; i++)
      (void) item5.GetUint8();
   item5.GetString(data.strvvar,sizeof(data.strvvar));

   item5.Done();

   item1.Rewind();

   if ( item1.Search(993) < 0 )
   {
      Warning("Cannot find test data sub-block 3.");
      return -4;
   }

   EventIO::Item item4(item1,"get",993);
   if ( item4.Status() < 0 )
   {
      Warning("Missing or invalid test data sub-block 3.");
      return -4;
   }

   item4.GetBool(data.bvar);
   item4.GetCount(data.cnt16var,4);
   item4.GetCount(data.cnt32var,6);
   item4.GetCount(data.cntzvar,6);
   item4.GetCount(data.cntvar,8);
   item4.GetSCount(data.scnt16var,10);
   item4.GetSCount(data.scnt32var,12);
   item4.GetSCount(data.scntzvar,12);
   item4.GetSCount(data.scntvar,14);

   item4.Done();

   item1.Done();
   return (item1.Status()==1?0:item1.Status());
}

int write_test_ex(EventIO& iobuf)
{
   std::cerr << "Expecting an exception for PutUint16(99999) ...\n";
   
   EventIO::Item item(iobuf,"put",99,0,123);
   if ( item.Status() != 0 )
      return item.Status();
   
   size_t n = 99999;
   bool caught = false;
   try 
   {
      item.PutUint16(n);
   }
   catch (const std::exception& error)
   {
      std::cerr << "Exception caught: " << error.what() << std::endl;
      std::cerr << "That looks fine.\n";
      caught = true;
   }
   if ( ! caught )
      std::cerr << "No exception caught. Hmmm ...\n";
   
   item.Done();
   return (item.Status()==1?0:item.Status());
}

void syntax(const char *prg)
{
   fprintf(stderr,"Test basic EventIO write and read functions (C++ version).\n");
   fprintf(stderr,"Syntax: %s [ -s ] [ -e ] filename [ filename2 ]\n", prg); 
   fprintf(stderr,"Options:\n");
   fprintf(stderr,"  -s  Show all data values as they are compared.\n");
   fprintf(stderr,"  -e  Use the extension field for all I/O block headers.\n");
   fprintf(stderr,"The optional second file is only opened for reading, and can\n");
   fprintf(stderr,"can be used to check for consistency between C and C++ versions.\n");
   fprintf(stderr,"The files produced by both versions should be indentical.\n");
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
   FILE *input;
   static TEST_DATA tdata, cdata1, cdata2;
   int ok;
   char *prg = argv[0];
   
   if ( getenv("SHOWDATA") != 0 )
      show_data = true;

   tdata.bvar.resize(13); 
   for (size_t i=0; i<tdata.bvar.size(); ++i) tdata.bvar[i] = false;
   tdata.bvar[2] = tdata.bvar[7] = tdata.bvar[8] = tdata.bvar[12] = true;
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
   tdata.i8var[0]  = 0x1e;         tdata.i8var[1]  = 0xe1;
   tdata.u8var[0]  = 0x1e;         tdata.u8var[1]  = 0xe1;
   tdata.i16var[0] = 0x1e2e;       tdata.i16var[1] = 0xe2e1;
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

   EventIO iobuf;
   if ( argc > 1 && strcmp(argv[1],"--help") == 0 )
      syntax(prg);
   if ( argc > 1 && strcmp(argv[1],"-s") == 0 )
   {
      show_data = true;
      argc--;
      argv++;
   }
   if ( argc > 1 && strcmp(argv[1],"-e") == 0 )
   {
      iobuf.SetExtended();
      Information("Using the extension field for all I/O block headers.");
      argc--;
      argv++;
   }

   if ( argc > 1 )
   {
      if ( strcmp(argv[1],"--help") == 0 )
      	 syntax(prg);
      iobuf.OpenOutput(argv[1]);
   }
   else
      syntax(prg);

   write_test_ex(iobuf);
   iobuf.OpenOutput(argv[1]);
   fprintf(stderr,"\nWrite test data to file '%s'.\n",argv[1]);
   fprintf(stderr,"Default byte order, using mainly vector functions.\n");
   write_test2(tdata,iobuf);
   fprintf(stderr,"Default byte order, using single-element functions.\n");
   write_test1(tdata,iobuf);
   iobuf.Buffer()->byte_order = 1;
   fprintf(stderr,"Reversed byte order, using single-element functions.\n");
   write_test1(tdata,iobuf);
   iobuf.Buffer()->byte_order = 0;
   fprintf(stderr,"Normal byte order, using single-element functions.\n");
   write_test1(tdata,iobuf);
   fprintf(stderr,"Writing as nested item structure.\n");
   write_test3(tdata,iobuf);
   fprintf(stderr,"Write tests done.\n\n");
   
   iobuf.CloseOutput();

   if ( argc > 2 )
   {
      fprintf(stderr,"\nComparing against data from external file %s:\n",
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

   iobuf.OpenInput(input);
   ok = 1;

   fprintf(stderr,"Default byte order, using single-element functions.\n");
   memset(&cdata1,0,sizeof(cdata1));
   if ( iobuf.Find() < 0 )
   {
      Error("*** Finding I/O block 1 failed");
      exit(1);
   }
   if ( iobuf.Read() < 0 )
   {
      Error("*** Reading I/O block 1 failed");
      exit(1);
   }
   if ( read_test1(cdata1,iobuf) < 0 )
   {
      Error("*** Read test 1 failed");
      exit(1);
   }
   if ( datacmp(tdata,cdata1) != 1 )
   {
      Error("*** Data from read test 1 does not match.");
      ok = 0;
   }
   
   fprintf(stderr,"Default byte order, using mainly vector functions.\n");
   memset(&cdata2,0,sizeof(cdata2));
   if ( iobuf.Find() < 0 )
   {
      Error("*** Finding I/O block 2 failed");
      exit(1);
   }
   if ( iobuf.Read() < 0 )
   {
      Error("*** Reading I/O block 2 failed");
      exit(1);
   }
   if ( read_test2(cdata2,iobuf) < 0 )
   {
      Error("*** Read test 2 failed");
      exit(1);
   }
   if ( datacmp(tdata,cdata2) != 1 )
   {
      Error("*** Data from read test 2 does not match");
      ok = 0;
   }

   fprintf(stderr,"Reversed byte order, using single-element functions.\n");
   memset(&cdata1,0,sizeof(cdata1));
   if ( iobuf.Find() < 0 )
   {
      Error("*** Finding I/O block 3 failed");
      exit(1);
   }
   if ( iobuf.Read() < 0 )
   {
      Error("*** Reading I/O block 3 failed");
      exit(1);
   }
   if ( read_test1(cdata1,iobuf) < 0 )
   {
      Error("*** Read test 3 failed");
      exit(1);
   }
   if ( datacmp(tdata,cdata1) != 1 )
   {
      Error("*** Data from read test 3 does not match.");
      ok = 0;
   }
   
   fprintf(stderr,"Normal byte order, using single-element functions.\n");
   memset(&cdata2,0,sizeof(cdata2));
   if ( iobuf.Find() < 0 )
   {
      Error("*** Finding I/O block 4 failed");
      exit(1);
   }
   if ( iobuf.Read() < 0 )
   {
      Error("*** Reading I/O block 4 failed");
      exit(1);
   }
   if ( read_test1(cdata2,iobuf) < 0 )
   {
      Error("*** Read test 4 failed");
      exit(1);
   }
   if ( datacmp(tdata,cdata2) != 1 )
   {
      Error("*** Data from read test 4 does not match");
      ok = 0;
   }
   
   fprintf(stderr,"Reading nested data structure.\n");
   memset(&cdata2,0,sizeof(cdata2));
   if ( iobuf.Find() < 0 )
   {
      Error("*** Finding I/O block 5 failed");
      exit(1);
   }
   if ( iobuf.Read() < 0 )
   {
      Error("*** Reading I/O block 5 failed");
      exit(1);
   }
   if ( read_test3(cdata2,iobuf) < 0 )
   {
      Error("*** Read test 5 failed");
      exit(1);
   }
   if ( datacmp(tdata,cdata2) != 1 )
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
# ifdef SIXTY_FOUR_BITS
   Information("These 64-bit integers are native types.\n");
# else
   Information("These 64-bit integers are implemented through the C compiler.\n");
# endif
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
      char stmp[1200];
      sprintf(stmp,"This version of eventio has additional features:\n");
#ifdef HAVE_EVENTIO_USER_FLAG
      strcat(stmp," - User flag bit in the header.\n");
#endif
#ifdef HAVE_EVENTIO_HEADER_LENGTH
      strcat(stmp," - The length of the data field can also be seen in the C language header.\n");
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
