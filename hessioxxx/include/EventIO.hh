/* ============================================================================

   Copyright (C) 2003, 2005, 2009, 2010, 2016 Konrad Bernloehr (Konrad Bernl&ouml;hr)

   This file is part of the eventio/hessio library.

   The eventio/hessio library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this library. If not, see <http://www.gnu.org/licenses/>.

============================================================================ */

/** @file EventIO.hh
    @short C++ interface to the eventio data format.
    
    @author  Konrad Bernloehr
    @date Initial release: April 2003
    $Date: 2016/11/18 12:30:51 $
    $Revision: 1.44 $
    
    Header file for C++ interface to the eventio data format.
    In constrast to the C interface, the C++ interface can
    take care of the completion of "put" and "get" operations
    at the end of the lifetime of an object.
    Nested data objects are much more straight-forward
    than in C.

    Overloading of methods is used to have the same name
    for all "put" and "get" operations of the same
    basic data type. The operations are optionally
    available also for std::vector and std::valarray
    arguments of (selected) suitable types - but only
    if the includes were done before this file is included.
    For standards conforming compilers the macros
    _CPP_VECTOR, _CPP_VALARRAY, and _CPP_STRING are
    used as indicators. For gcc 2.96, the macros
    __SGI_STL_VECTOR, __STD_VALARRAY__, and __BASTRING__
    seem to be corresponding indicators. For GCC version 3.4
    this was changed to _GLIBCXX_VECTOR, _GLIBCXX_VALARRAY, 
    and _GLIBCXX_STRING.
*/

#ifndef EVENTIO_HH_INCLUDED

#define EVENTIO_HH_INCLUDED 1

#define HAVE_64BIT_INT 1
#include "io_basic.h"
#include "eventio_registry.h"
#include <stdexcept>
#include <limits.h>

#include <vector>
#include <valarray>
#include <string>

#if defined(_CPP_VECTOR) || defined(__SGI_STL_VECTOR) || defined(_GLIBCXX_VECTOR) || defined(_LIBCPP_VECTOR)
#define HAVE_STD_VECTOR 1
#endif

#if defined(_CPP_VALARRAY) || defined(__STD_VALARRAY__) || defined(_GLIBCXX_VALARRAY) || defined(_LIBCPP_VALARRAY)
#define HAVE_STD_VALARRAY 1
#endif

#if defined(_CPP_STRING) || defined(__BASTRING__) || defined(_GLIBCXX_STRING) || defined(_LIBCPP_STRING)
#define HAVE_STD_STRING
#endif

#if defined(UINTMAX_MIN) && defined(HAVE_64BIT_INT) && !defined(SIXTY_FOUR_BITS)
#define WITH_INTMAX_T 1
#define WITH_UINTMAX_T 1
#endif

// Optional support for half-procession floating point header, if not included yet.
#ifndef HALF_HALF_HPP
# ifdef HAVE_HALF
// Explicitly requested; assume it is available.
#  include "half.hpp"
using half_float::half;
// #pragma message "Including half.hpp(1)"
# else
// When not explicitly requested, the preprocessor can perhaps find out if the half.hpp header file is available.
#  if defined(__has_include)
#   if __has_include("half.hpp")
#    include "half.hpp"
using half_float::half;
// #pragma message "Including half.hpp(2)"
#    define HAVE_HALF 1
#   endif
#  endif
# endif
#else
# ifndef HAVE_HALF
#  define HAVE_HALF 1
# endif
using half_float::half;
// #pragma message "Already included half.hpp"
#endif

/** The classes of this interface belong to the namespace "eventio". */

namespace eventio 
{
   class EventIO;  // Forward declaration.
   // class EventIO::Item;   // Forward declaration.
   
   /// This class provides the embedded buffer, 
   /// the file I/O interface and the toplevel item access.

   class EventIO
   {
      public:

      /// This (sub-) class provides all the interfaces for putting
      /// and getting basic data to and from the embedded buffer.

      class Item
      {
         friend class eventio::EventIO;

         private:
            IO_BUFFER *iobuf;
            IO_ITEM_HEADER item_header;
            int rc;
            bool put_flag;
            bool throw_on_error;
            bool done;
            bool active; ///< From {get|put}_item_begin to {get|put}_item_end
            Item *child;
            Item *orig_parent;
            EventIO *toplevel_parent;
            
            inline void check_throw(const char *op="unknown") const
            { if ( rc != 0 && throw_on_error )
                 throw std::logic_error(std::string("I/O buffer in error state at ")
                    + op + " operation"); }

         public:
            Item(EventIO &ev, const char *method, 
               size_t type=0, size_t version=0, long ident=0, 
               bool user_flag=false, bool extended_flag=false);
            Item(Item &parent, const char *method, 
               size_t type=0, size_t version=0, long ident=0, 
               bool user_flag=false, bool extended_flag=false);
            Item(const Item& item);
            ~Item(void);
            Item& operator= (const Item& item);
            int Done(void);
            /// Returns <0 in case of error, 0 for active item, 1 for finished item. 
            int Status(void) const { return rc; }
            /// For more tricky things you can still access the underlying I/O buffer.
            IO_BUFFER *Buffer(void) { return iobuf; }
            /// Return the item type
            size_t Type(void) const { return size_t(item_header.type); }
            /// Return the version number
            size_t Version(void) const { return size_t(item_header.version); }
            /// Return the item ID
            long Ident(void) const { return item_header.ident; }
            /// Return the nesting level depth
            int Depth(void) const { return item_header.level; }
            /// Return the data size (excluding header) 
            size_t size(void) const { int l = item_header.level;
               if ( l<0 || l>iobuf->item_level || l>MAX_IO_ITEM_LEVEL ) return 0;
               return iobuf->item_length[l]; }
            /// Return the length of the data area in the current item
            size_t Length(void) const { return size(); }
            /// Return whether an item can be searched for sub-items
            bool IsSearchable(void) const { return item_header.can_search; }
            /// Search for a specific sub-item
            int Search(size_t sub);
            /// Rewind to beginning of data area
            int Rewind(void);
            /// Completely undo getting this item, i.e. rewind to beginning of its header.
            /// The item should not be used any more afterwards.
            int Unget(void);
            /// Completely undo putting this item, as if never started.
            /// The item should not be used any more afterwards.
            int Unput(void);
            /// Skip a sub-item starting at the current  position
            int Skip(void);
            /// Is the user flag set?
            bool UserFlag(void) const { return item_header.user_flag; }
            bool IsExtended(void) const { return item_header.use_extension; }
            /// List item type, length, ID, ... + sub-items
            int List(int maxlevel=0, int verbosity=0);
            /// Registered name for this type of I/O block
            const char *TypeName(void);
            /// Registered description for this type of I/O block
            const char *Description(void);

            /// Access to information of next sub-item at current reading position
            int NextSubItemType(void) const;
            size_t NextSubItemLength(void) const;
            long NextSubItemIdent(void) const;

            // ------------- Get operations ---------------

            /// Get operations for bytes (only unsigned flavour implemented).
            /// Even for std::vector and std::valarray the number of
            /// elements to get from the input buffer has to be specified
            /// explicitly.

            uint8_t GetUint8(void) { if ( rc == 0 ) return get_byte(iobuf); 
               else if ( throw_on_error ) 
                  throw std::logic_error("I/O buffer in error state at GetUint8 operation");
               else return 0; }
            void GetUint8(uint8_t& v) { v=GetUint8(); }
            void GetUint8(uint8_t *vec, size_t num) { if ( rc == 0 )
               get_vector_of_uint8(vec,num,iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetUint8 operation"); }
#ifdef HAVE_STD_VECTOR
            void GetUint8(std::vector<uint8_t>& vec, size_t num);
            void GetUint8(std::vector<uint8_t>& vec);
#endif
#ifdef HAVE_STD_VALARRAY
            void GetUint8(std::valarray<uint8_t>& vec, size_t num);
            void GetUint8(std::valarray<uint8_t>& vec);
#endif

            /// Individual 'bool' elements are stored as unsigned chars holding 
            /// up to 8 bools.
            /// Note that retrieving a vector of n bools is not the same as
            /// n times retrieving a single bool.

            bool GetBool(void) { if ( rc == 0 ) return GetUint8(); 
               else if ( throw_on_error ) 
                  throw std::logic_error("I/O buffer in error state at GetBool operation");
               else return false; }
            void GetBool(bool& v) { v = GetBool(); }
#ifdef HAVE_STD_VECTOR
            void GetBool(std::vector<bool>& vec, size_t num);
            void GetBool(std::vector<bool>& vec);
#endif
#ifdef HAVE_STD_VALARRAY
            void GetBool(std::valarray<bool>& vec, size_t num);
            void GetBool(std::valarray<bool>& vec);
#endif

            /// Get 'count' type data (unsigned integer stored with variable lenghth).

            uint16_t GetCount16(void) { if ( rc == 0 ) return get_count16(iobuf);
               else if ( throw_on_error ) 
                  throw std::logic_error("I/O buffer in error state at GetCount16 operation");
               else return 0; }
            uint32_t GetCount32(void) { if ( rc == 0 ) return get_count32(iobuf);
               else if ( throw_on_error ) 
                  throw std::logic_error("I/O buffer in error state at GetCount32 operation");
               else return 0; }
            uintmax_t GetCount(void) { if ( rc == 0 ) return get_count(iobuf);
               else if ( throw_on_error ) 
                  throw std::logic_error("I/O buffer in error state at GetCount operation");
               else return 0; }
            void GetCount(uint16_t *vec, size_t num);
#ifdef SIXTY_FOUR_BITS
            void GetCount(uint32_t *vec, size_t num);
#endif
            void GetCount(size_t *vec, size_t num);
            void GetCount(uint8_t& v) { v = int8_t(GetCount16()); }
            void GetCount(uint16_t& v) { v = int16_t(GetCount16()); }
#ifdef SIXTY_FOUR_BITS
            void GetCount(uint32_t& v) { v = uint32_t(GetCount32()); }
#endif
            void GetCount(size_t& v) { v = size_t(GetCount()); }
#if defined(WITH_UINTMAX_T) && !defined(SIXTY_FOUR_BITS)
            void GetCount(uintmax_t *vec, size_t num);
            void GetCount(uintmax_t& v) { v = GetCount(); }
#endif
#ifdef HAVE_STD_VECTOR
            void GetCount(std::vector<uint16_t>& vec, size_t num);
            void GetCount(std::vector<uint16_t>& vec);
#ifdef SIXTY_FOUR_BITS
            void GetCount(std::vector<uint32_t>& vec, size_t num);
            void GetCount(std::vector<uint32_t>& vec);
#endif
            void GetCount(std::vector<size_t>& vec, size_t num);
            void GetCount(std::vector<size_t>& vec);
# if defined(WITH_UINTMAX_T) && !defined(SIXTY_FOUR_BITS)
            void GetCount(std::vector<uintmax_t>& vec, size_t num);
            void GetCount(std::vector<uintmax_t>& vec);
# endif
#endif
#ifdef HAVE_STD_VALARRAY
            void GetCount(std::valarray<uint16_t>& vec, size_t num);
            void GetCount(std::valarray<uint16_t>& vec);
#ifdef SIXTY_FOUR_BITS
            void GetCount(std::valarray<uint32_t>& vec, size_t num);
            void GetCount(std::valarray<uint32_t>& vec);
#endif
            void GetCount(std::valarray<size_t>& vec, size_t num);
            void GetCount(std::valarray<size_t>& vec);
# if defined(WITH_UINTMAX_T) && !defined(SIXTY_FOUR_BITS)
            void GetCount(std::valarray<uintmax_t>& vec, size_t num);
            void GetCount(std::valarray<uintmax_t>& vec);
# endif
#endif

            /// Get 'scount' type data (signed integer stored with variable lenghth).

            int16_t GetSCount16(void) { if ( rc == 0 ) return get_scount16(iobuf);
               else if ( throw_on_error ) 
                  throw std::logic_error("I/O buffer in error state at GetSCount16 operation");
               else return 0; }
            int32_t GetSCount32(void) { if ( rc == 0 ) return get_scount32(iobuf);
               else if ( throw_on_error ) 
                  throw std::logic_error("I/O buffer in error state at GetSCount32 operation");
               else return 0; }
            intmax_t GetSCount(void) { if ( rc == 0 ) return get_scount(iobuf);
               else if ( throw_on_error ) 
                  throw std::logic_error("I/O buffer in error state at GetSCount operation");
               else return 0; }
            void GetSCount(int16_t *vec, size_t num);
            void GetDiffSCount(int16_t *vec, size_t num);
#ifdef SIXTY_FOUR_BITS
            void GetSCount(int32_t *vec, size_t num);
            void GetDiffSCount(int32_t *vec, size_t num);
#endif
            void GetSCount(ssize_t *vec, size_t num);
            void GetDiffSCount(ssize_t *vec, size_t num);
            void GetSCount(int8_t& v) { v = int8_t(GetSCount16()); }
            void GetDiffSCount(int8_t& v) { v = int8_t(GetSCount16()); }
            void GetSCount(int16_t& v) { v = int16_t(GetSCount16()); }
            void GetDiffSCount(int16_t& v) { v = int16_t(GetSCount16()); }
#ifdef SIXTY_FOUR_BITS
            void GetSCount(int32_t& v) { v = int32_t(GetSCount32()); }
            void GetDiffSCount(int32_t& v) { v = int32_t(GetSCount32()); }
#endif
            void GetSCount(ssize_t& v) { v = ssize_t(GetSCount()); }
            void GetDiffSCount(ssize_t& v) { v = ssize_t(GetSCount()); }
#if defined(WITH_INTMAX_T) && !defined(SIXTY_FOUR_BITS)
            void GetSCount(intmax_t *vec, size_t num);
            void GetDiffSCount(intmax_t *vec, size_t num);
            void GetSCount(intmax_t& v) { v = GetSCount(); }
            void GetDiffSCount(intmax_t& v) { v = GetSCount(); }
#endif
#ifdef HAVE_STD_VECTOR
            void GetSCount(std::vector<int16_t>& vec, size_t num);
            void GetDiffSCount(std::vector<int16_t>& vec, size_t num);
            void GetSCount(std::vector<int16_t>& vec);
            void GetDiffSCount(std::vector<int16_t>& vec);
#ifdef SIXTY_FOUR_BITS
            void GetSCount(std::vector<int32_t>& vec, size_t num);
            void GetDiffSCount(std::vector<int32_t>& vec, size_t num);
            void GetSCount(std::vector<int32_t>& vec);
            void GetDiffSCount(std::vector<int32_t>& vec);
#endif
            void GetSCount(std::vector<ssize_t>& vec, size_t num);
            void GetDiffSCount(std::vector<ssize_t>& vec, size_t num);
            void GetSCount(std::vector<ssize_t>& vec);
            void GetDiffSCount(std::vector<ssize_t>& vec);
# if defined(WITH_INTMAX_T) && !defined(SIXTY_FOUR_BITS)
            void GetSCount(std::vector<intmax_t>& vec, size_t num);
            void GetDiffSCount(std::vector<intmax_t>& vec, size_t num);
            void GetSCount(std::vector<intmax_t>& vec);
            void GetDiffSCount(std::vector<intmax_t>& vec);
# endif
#endif
#ifdef HAVE_STD_VALARRAY
            void GetSCount(std::valarray<int16_t>& vec, size_t num);
            void GetDiffSCount(std::valarray<int16_t>& vec, size_t num);
            void GetSCount(std::valarray<int16_t>& vec);
            void GetDiffSCount(std::valarray<int16_t>& vec);
#ifdef SIXTY_FOUR_BITS
            void GetSCount(std::valarray<int32_t>& vec, size_t num);
            void GetDiffSCount(std::valarray<int32_t>& vec, size_t num);
            void GetSCount(std::valarray<int32_t>& vec);
            void GetDiffSCount(std::valarray<int32_t>& vec);
#endif
            void GetSCount(std::valarray<ssize_t>& vec, size_t num);
            void GetDiffSCount(std::valarray<ssize_t>& vec, size_t num);
            void GetSCount(std::valarray<ssize_t>& vec);
            void GetDiffSCount(std::valarray<ssize_t>& vec);
# if defined(WITH_INTMAX_T) && !defined(SIXTY_FOUR_BITS)
            void GetSCount(std::valarray<intmax_t>& vec, size_t num);
            void GetDiffSCount(std::valarray<intmax_t>& vec, size_t num);
            void GetSCount(std::valarray<intmax_t>& vec);
            void GetDiffSCount(std::valarray<intmax_t>& vec);
# endif
#endif

            /// Get operations for 'Int16' item type (signed 16-bit integers, 'SHORT')

            int16_t GetInt16(void) { if ( rc == 0 ) return get_short(iobuf);
               else if ( throw_on_error ) 
                  throw std::logic_error("I/O buffer in error state at GetInt16 operation");
               else return 0; }
            void GetInt16(int16_t& v) { v=GetInt16(); }
            void GetInt16(int& v) { v=GetInt16(); }
            void GetInt16(int16_t *vec, size_t num) { if ( rc == 0 )
               get_vector_of_int16(vec,num,iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetInt16 operation"); }
            void GetInt16(int *vec, size_t num) { if ( rc == 0 )
               get_vector_of_int(vec,num,iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetInt16 operation"); }
#ifdef HAVE_STD_VECTOR
            void GetInt16(std::vector<int16_t>& vec, size_t num);
            void GetInt16(std::vector<int>& vec, size_t num);
            void GetInt16(std::vector<int16_t>& vec);
            void GetInt16(std::vector<int>& vec);
#endif
#ifdef HAVE_STD_VALARRAY
            void GetInt16(std::valarray<int16_t>& vec, size_t num);
            void GetInt16(std::valarray<int>& vec, size_t num);
            void GetInt16(std::valarray<int16_t>& vec);
            void GetInt16(std::valarray<int>& vec);
#endif

            /// Get operations for 'Uint16' item type (unsigned 16-bit integers)

            uint16_t GetUint16(void) { if ( rc == 0 )
                  return static_cast<uint16_t>(get_short(iobuf)); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetUint16 operation");
               else return 0; }
            void GetUint16(uint16_t &v) { v=GetUint16(); }
            void GetUint16(size_t &v) { v=GetUint16(); }
            void GetUint16(uint16_t *vec, size_t num) { if ( rc == 0 )
               get_vector_of_uint16(vec,num,iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetUint16 operation"); }
#ifdef HAVE_STD_VECTOR
            void GetUint16(std::vector<uint16_t>& vec, size_t num);
            void GetUint16(std::vector<unsigned int>& vec, size_t num);
            void GetUint16(std::vector<uint16_t>& vec);
            void GetUint16(std::vector<unsigned int>& vec);
#endif
#ifdef HAVE_STD_VALARRAY
            void GetUint16(std::valarray<uint16_t>& vec, size_t num);
            void GetUint16(std::valarray<unsigned int>& vec, size_t num);
            void GetUint16(std::valarray<uint16_t>& vec);
            void GetUint16(std::valarray<unsigned int>& vec);
#endif

            /// Get operations for 'Int32' item type (signed 32-bit integers, 'LONG')

            int32_t GetInt32(void) { if ( rc == 0 ) return get_int32(iobuf);
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetInt32 operation");
               else return 0; }
            void GetInt32(int32_t &v) { v = GetInt32(); }
            void GetInt32(size_t &v) { v = GetInt32(); }
            void GetInt32(long &v) { v = GetInt32(); }
            void GetInt32(int32_t *vec, size_t num) { if ( rc == 0 )
               get_vector_of_int32(vec,num,iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetInt32 operation"); }
            void GetInt32(long *vec, size_t num) { if ( rc == 0 )
               get_vector_of_long(vec,num,iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetInt32 operation"); }
#ifdef HAVE_STD_VECTOR
            void GetInt32(std::vector<int32_t>& vec, size_t num);
            void GetInt32(std::vector<long>& vec, size_t num);
            void GetInt32(std::vector<int32_t>& vec);
            void GetInt32(std::vector<long>& vec);
#endif
#ifdef HAVE_STD_VALARRAY
            void GetInt32(std::valarray<int32_t>& vec, size_t num);
            void GetInt32(std::valarray<long>& vec, size_t num);
            void GetInt32(std::valarray<int32_t>& vec);
            void GetInt32(std::valarray<long>& vec);
#endif

            /// Get operations for 'Uint32' item type (unsigned 32-bit integers)

            uint32_t GetUint32(void) { if ( rc == 0 ) return get_uint32(iobuf);
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetUint32 operation");
               else return 0; }
            void GetUint32(uint32_t& v) { v = GetUint32(); }
            void GetUint32(unsigned long& v) { v = GetUint32(); }
            void GetUint32(uint32_t *vec, size_t num) { if ( rc == 0 )
               get_vector_of_uint32(vec,num,iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetUnt32 operation"); }
#ifdef HAVE_STD_VECTOR
            void GetUint32(std::vector<uint32_t>& vec, size_t num);
            void GetUint32(std::vector<unsigned long>& vec, size_t num);
            void GetUint32(std::vector<uint32_t>& vec);
            void GetUint32(std::vector<unsigned long>& vec);
#endif
#ifdef HAVE_STD_VALARRAY
            void GetUint32(std::valarray<uint32_t>& vec, size_t num);
            void GetUint32(std::valarray<unsigned long>& vec, size_t num);
            void GetUint32(std::valarray<uint32_t>& vec);
            void GetUint32(std::valarray<unsigned long>& vec);
#endif
#ifdef HAVE_64BIT_INT

            /// Get operations for 'Int64' item type (signed 64-bit integers)

            int64_t GetInt64(void) { if ( rc == 0 )
               { int64_t i; get_vector_of_int64(&i,1,iobuf);
               return i; }
               else if ( throw_on_error )
                 throw std::logic_error("I/O buffer in error state at GetInt64 operation");
               else return 0; }
            void GetInt64(int64_t& v) { v = GetInt64(); }
            void GetInt64(int64_t *vec, size_t num) { if ( rc == 0 )
               get_vector_of_int64(vec,num,iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetInt64 operation"); }
#ifdef HAVE_STD_VECTOR
            void GetInt64(std::vector<int64_t>& vec, size_t num);
            void GetInt64(std::vector<int64_t>& vec);
#endif
#ifdef HAVE_STD_VALARRAY
            void GetInt64(std::valarray<int64_t>& vec, size_t num);
            void GetInt64(std::valarray<int64_t>& vec);
#endif

            /// Get operations for 'Uint64' item type (unsigned 64-bit integers)

            uint64_t GetUint64(void) { if ( rc == 0 )
               { uint64_t i; get_vector_of_uint64(&i,1,iobuf);
               return i; }
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetUint64 operation");
               else return 0; }
            void GetUint64(uint64_t& v) { v = GetUint64(); }
            void GetUint64(uint64_t *vec, size_t num) { if ( rc == 0 )
               get_vector_of_uint64(vec,num,iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetUint64 operation"); }
#ifdef HAVE_STD_VECTOR
            void GetUint64(std::vector<uint64_t>& vec, size_t num);
            void GetUint64(std::vector<uint64_t>& vec);
#endif
#ifdef HAVE_STD_VALARRAY
            void GetUint64(std::valarray<uint64_t>& vec, size_t num);
            void GetUint64(std::valarray<uint64_t>& vec);
#endif
#endif

            /// Get operations for 'Real' item type (32-bit IEEE float)

            double GetReal(void) { if ( rc == 0 ) return get_real(iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetReal operation");
               else return 0.; }
            void GetReal(float& v) { v = (float)GetReal(); }
            void GetReal(double& v) { v = GetReal(); }
            void GetReal(float *vec, size_t num) { if ( rc == 0 )
               get_vector_of_float(vec,num,iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetReal operation"); }
            void GetReal(double *vec, size_t num) { if ( rc == 0 )
               get_vector_of_real(vec,num,iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetReal operation"); }
#ifdef HAVE_STD_VECTOR
            void GetReal(std::vector<float>& vec, size_t num);
            void GetReal(std::vector<double>& vec, size_t num);
            void GetReal(std::vector<float>& vec);
            void GetReal(std::vector<double>& vec);
#endif
#ifdef HAVE_STD_VALARRAY
            void GetReal(std::valarray<float>& vec, size_t num);
            void GetReal(std::valarray<double>& vec, size_t num);
            void GetReal(std::valarray<float>& vec);
            void GetReal(std::valarray<double>& vec);
#endif

            /// Get operations for 'Double' item type (64-bit IEEE double)

            double GetDouble(void) { if ( rc == 0 ) return get_double(iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetDouble operation");
               else return 0.; }
            void GetDouble(double &v) { v = GetDouble(); }
            void GetDouble(double *vec, size_t num) { if ( rc == 0 )
               get_vector_of_double(vec,num,iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetDouble operation"); }
#ifdef HAVE_STD_VECTOR
            void GetDouble(std::vector<double>& vec, size_t num);
            void GetDouble(std::vector<double>& vec);
#endif
#ifdef HAVE_STD_VALARRAY
            void GetDouble(std::valarray<double>& vec, size_t num);
            void GetDouble(std::valarray<double>& vec);
#endif

            /// Get operations for 'Sfloat' item type (16-bit OpenGL/IEEE 754-2008 float)

            double GetSfloat(void) { if ( rc == 0 ) return get_sfloat(iobuf); 
               else if ( throw_on_error ) 
                 throw std::logic_error("I/O buffer in error state at GetSfloat operation");
               else return 0.; }
            void GetSfloat(float &v) { v = (float)GetSfloat(); }
            void GetSfloat(double &v) { v = GetSfloat(); }
            void GetSfloat(double *vec, size_t num) { if ( rc == 0 )
               {
                  for ( size_t i=0; i<num; i++ )
                     vec[i] = get_sfloat(iobuf);
               }
               else if ( throw_on_error )
                 throw std::logic_error("I/O buffer in error state at GetSfloat operation"); }
            void GetSfloat(float *vec, size_t num) { if ( rc == 0 )
               {
                  for ( size_t i=0; i<num; i++ )
                     vec[i] = (float)get_sfloat(iobuf);
               }
               else if ( throw_on_error )
                 throw std::logic_error("I/O buffer in error state at GetSfloat operation"); }
#ifdef HAVE_STD_VECTOR
            void GetSfloat(std::vector<float>& vec, size_t num);
            void GetSfloat(std::vector<double>& vec, size_t num);
            void GetSfloat(std::vector<float>& vec);
            void GetSfloat(std::vector<double>& vec);
#endif
#ifdef HAVE_STD_VALARRAY
            void GetSfloat(std::valarray<float>& vec, size_t num);
            void GetSfloat(std::valarray<double>& vec, size_t num);
            void GetSfloat(std::valarray<float>& vec);
            void GetSfloat(std::valarray<double>& vec);
#endif

#ifdef HAVE_HALF
            /// Alternate access to 16-bit floats ('half' has uint16_t storage representation, same as 'sfloat')
            half GetHalf(void) { uint16_t u = GetUint16(); return *((half *) (&u)); }
            void GetHalf(half &v) { v = GetHalf(); }
            void GetHalf(float &v) { v = (float) GetHalf(); }
            void GetHalf(double &v) { v = (double) GetHalf(); }
            void GetHalf(half *vec, size_t num) { if ( rc == 0 )
               {
                  uint16_t *p = (uint16_t *) (&(vec[0]));
                  for ( size_t i=0; i<num; i++ )
                     p[i] = get_short(iobuf);
               }
               else if ( throw_on_error )
                 throw std::logic_error("I/O buffer in error state at GetHalf operation"); }
#ifdef HAVE_STD_VECTOR
            void GetHalf(std::vector<half>& vec, size_t num);
            void GetHalf(std::vector<half>& vec);
#endif
#ifdef HAVE_STD_VALARRAY
            void GetHalf(std::valarray<half>& vec, size_t num);
            void GetHalf(std::valarray<half>& vec);
#endif
#else
/// Without HAVE_HALF you can still use GetUint16() and use dbl_from_sfloat() or such to
/// convert to other types.
#endif


            /// Get operations for strings with length encoded as Int16.

            int GetString(char *s, size_t nmax);
#ifdef HAVE_STD_STRING
            int GetString16(std::string& s);
            std::string GetString16() { std::string s; GetString16(s); return s; }
            int GetString32(std::string& s);
            std::string GetString32() { std::string s; GetString32(s); return s; }
            int GetString(std::string& s);
            std::string GetString() { std::string s; GetString(s); return s; }
#ifdef HAVE_STD_VECTOR
            void GetString(std::vector<std::string>& vec, size_t num);
            void GetString(std::vector<std::string>& vec);
#endif
#ifdef HAVE_STD_VALARRAY
            void GetString(std::valarray<std::string>& vec, size_t num);
            void GetString(std::valarray<std::string>& vec);
#endif
#endif

            // ------------------------ Put operations ----------------------

            /// Put operations for bytes (only unsigned flavour implemented).
            /// For std::vector and std::valarray sloppy variants are available
            /// to put all elements.

            void PutUint8(uint8_t b) { if ( rc == 0 ) put_byte(b,iobuf);
               else check_throw("Put"); }
            void PutUint8(const uint8_t *vec, size_t num) { if ( rc == 0 )
               put_vector_of_uint8(const_cast<uint8_t *>(vec),num,iobuf);
               else check_throw("Put"); }
#ifdef HAVE_STD_VECTOR
            void PutUint8(const std::vector<uint8_t>& vec, size_t num);
            void PutUint8(const std::vector<uint8_t>& vec)
               { PutCount(vec.size()); PutUint8(vec,vec.size()); }
#endif
#ifdef HAVE_STD_VALARRAY
            void PutUint8(const std::valarray<uint8_t>& vec, size_t num);
            void PutUint8(const std::valarray<uint8_t>& vec)
               { PutCount(vec.size()); PutUint8(vec,vec.size()); }
#endif
            void PutBool(bool flag) { PutUint8(flag?1:0); }
#ifdef HAVE_STD_VECTOR
            void PutBool(const std::vector<bool>& vec, size_t num);
            void PutBool(const std::vector<bool>& vec)
               { PutCount(vec.size()); PutBool(vec,vec.size()); }
#endif
#ifdef HAVE_STD_VALARRAY
            void PutBool(const std::valarray<bool>& vec, size_t num);
            void PutBool(const std::valarray<bool>& vec)
               { PutCount(vec.size()); PutBool(vec,vec.size()); }
#endif

            /// Put operation for 'Count' item type (unsigned integer of variable length);

            void PutCount16(size_t n) { if ( rc == 0 ) put_count16(n,iobuf);
               else check_throw("Put"); }
            void PutCount(const uint16_t *vec, size_t num) { if ( rc == 0 )
               { for (size_t i=0; i<num; i++) put_count16(vec[i],iobuf); }
               else check_throw("Put"); }
            void PutCount32(uint32_t n) { if ( rc == 0 ) put_count32(n,iobuf);
               else check_throw("Put"); }
#ifdef SIXTY_FOUR_BITS
//            void PutCount32(size_t n) { if ( rc == 0 ) put_count32(n,iobuf);
//               else check_throw("Put"); }
            void PutCount(const uint32_t *vec, size_t num) { if ( rc == 0 )
               { for (size_t i=0; i<num; i++) put_count32(vec[i],iobuf); }
               else check_throw("Put"); }
#endif
            void PutCount(const size_t *vec, size_t num) { if ( rc == 0 )
               { for (size_t i=0; i<num; i++) put_count(vec[i],iobuf); }
               else check_throw("Put"); }
#if defined(WITH_UINTMAX_T) && !defined(SIXTY_FOUR_BITS)
            void PutCount(uintmax_t n) { if ( rc == 0 ) put_count(n,iobuf);
               else check_throw("Put"); }
            void PutCount(const uintmax_t *vec, size_t num) { if ( rc == 0 )
               { for (size_t i=0; i<num; i++) put_count(vec[i],iobuf); }
               else check_throw("Put"); }
#else
            void PutCount(size_t n) { if ( rc == 0 ) put_count(n,iobuf);
               else check_throw("Put"); }
#endif
#ifdef HAVE_STD_VECTOR
            void PutCount(const std::vector<uint16_t>& vec, size_t num);
            void PutCount(const std::vector<uint16_t>& vec)
               { PutCount(vec.size()); PutCount(vec,vec.size()); }
#ifdef SIXTY_FOUR_BITS
            void PutCount(const std::vector<uint32_t>& vec, size_t num);
            void PutCount(const std::vector<uint32_t>& vec)
               { PutCount(vec.size()); PutCount(vec,vec.size()); }
#endif
            void PutCount(const std::vector<size_t>& vec, size_t num);
            void PutCount(const std::vector<size_t>& vec)
               { PutCount(vec.size()); PutCount(vec,vec.size()); }
#if defined(WITH_UINTMAX_T) && !defined(SIXTY_FOUR_BITS)
            void PutCount(const std::vector<uintmax_t>& vec, size_t num);
            void PutCount(const std::vector<uintmax_t>& vec)
               { PutCount(vec.size()); PutCount(vec,vec.size()); }
#endif
#endif
#ifdef HAVE_STD_VALARRAY
            void PutCount(const std::valarray<uint16_t>& vec, size_t num);
            void PutCount(const std::valarray<uint16_t>& vec)
               { PutCount(vec.size()); PutCount(vec,vec.size()); }
#ifdef SIXTY_FOUR_BITS
            void PutCount(const std::valarray<uint32_t>& vec, size_t num);
            void PutCount(const std::valarray<uint32_t>& vec)
               { PutCount(vec.size()); PutCount(vec,vec.size()); }
#endif
            void PutCount(const std::valarray<size_t>& vec, size_t num);
            void PutCount(const std::valarray<size_t>& vec)
               { PutCount(vec.size()); PutCount(vec,vec.size()); }
#if defined(WITH_UINTMAX_T) && !defined(SIXTY_FOUR_BITS)
            void PutCount(const std::valarray<uintmax_t>& vec, size_t num);
            void PutCount(const std::valarray<uintmax_t>& vec)
               { PutCount(vec.size()); PutCount(vec,vec.size()); }
#endif
#endif

            /// Put operation for 'SCount' item type (signed integer of variable length);

            void PutSCount16(int16_t n) { if ( rc == 0 ) put_scount16(n,iobuf);
               else check_throw("Put"); }
            void PutSCount(const int16_t *vec, size_t num) { if ( rc == 0 )
               { for (size_t i=0; i<num; i++) put_scount16(vec[i],iobuf); }
               else check_throw("Put"); }
            void PutSCount32(int32_t n) { if ( rc == 0 ) put_scount32(n,iobuf);
               else check_throw("Put"); }
#ifdef SIXTY_FOUR_BITS
            void PutSCount32(ssize_t n) { if ( rc == 0 ) put_scount32(n,iobuf);
               else check_throw("Put"); }
            void PutSCount(const int32_t *vec, size_t num) { if ( rc == 0 )
               { for (size_t i=0; i<num; i++) put_scount32(vec[i],iobuf); }
               else check_throw("Put"); }
#endif
            void PutSCount(const ssize_t *vec, size_t num) { if ( rc == 0 )
               { for (size_t i=0; i<num; i++) put_scount(vec[i],iobuf); }
               else check_throw("Put"); }
#if defined(WITH_INTMAX_T) && !defined(SIXTY_FOUR_BITS)
            void PutSCount(intmax_t n) { if ( rc == 0 ) put_scount(n,iobuf);
               else check_throw("Put"); }
            void PutSCount(const intmax_t *vec, size_t num) { if ( rc == 0 )
               for (size_t i=0; i<num; i++) put_scount(vec[i],iobuf); }
#else
            void PutSCount(ssize_t n) { if ( rc == 0 ) put_scount(n,iobuf);
               else check_throw("Put"); }
#endif
#ifdef HAVE_STD_VECTOR
            void PutSCount(const std::vector<int16_t>& vec, size_t num);
            void PutDiffSCount(const std::vector<int16_t>& vec, size_t num);
            void PutSCount(const std::vector<int16_t>& vec)
               { PutCount(vec.size()); PutSCount(vec,vec.size()); }
            void PutDiffSCount(const std::vector<int16_t>& vec)
               { PutCount(vec.size()); PutDiffSCount(vec,vec.size()); }
#ifdef SIXTY_FOUR_BITS
            void PutSCount(const std::vector<int32_t>& vec, size_t num);
            void PutDiffSCount(const std::vector<int32_t>& vec, size_t num);
            void PutSCount(const std::vector<int32_t>& vec)
               { PutCount(vec.size()); PutSCount(vec,vec.size()); }
            void PutDiffSCount(const std::vector<int32_t>& vec)
               { PutCount(vec.size()); PutDiffSCount(vec,vec.size()); }
#endif
            void PutSCount(const std::vector<ssize_t>& vec, size_t num);
            void PutDiffSCount(const std::vector<ssize_t>& vec, size_t num);
            void PutSCount(const std::vector<ssize_t>& vec)
               { PutCount(vec.size()); PutSCount(vec,vec.size()); }
            void PutDiffSCount(const std::vector<ssize_t>& vec)
               { PutCount(vec.size()); PutDiffSCount(vec,vec.size()); }
#if defined(WITH_INTMAX_T) && !defined(SIXTY_FOUR_BITS)
            void PutSCount(const std::vector<intmax_t>& vec, size_t num);
            void PutDiffSCount(const std::vector<intmax_t>& vec, size_t num);
            void PutSCount(const std::vector<intmax_t>& vec)
               { PutCount(vec.size()); PutSCount(vec,vec.size()); }
            void PutDiffSCount(const std::vector<intmax_t>& vec)
               { PutCount(vec.size()); PutDiffSCount(vec,vec.size()); }
#endif
#endif
#ifdef HAVE_STD_VALARRAY
            void PutSCount(const std::valarray<int16_t>& vec, size_t num);
            void PutDiffSCount(const std::valarray<int16_t>& vec, size_t num);
            void PutSCount(const std::valarray<int16_t>& vec)
               { PutCount(vec.size()); PutSCount(vec,vec.size()); }
            void PutDiffSCount(const std::valarray<int16_t>& vec)
               { PutCount(vec.size()); PutDiffSCount(vec,vec.size()); }
#ifdef SIXTY_FOUR_BITS
            void PutSCount(const std::valarray<int32_t>& vec, size_t num);
            void PutDiffSCount(const std::valarray<int32_t>& vec, size_t num);
            void PutSCount(const std::valarray<int32_t>& vec)
               { PutCount(vec.size()); PutSCount(vec,vec.size()); }
            void PutDiffSCount(const std::valarray<int32_t>& vec)
               { PutCount(vec.size()); PutDiffSCount(vec,vec.size()); }
#endif
            void PutSCount(const std::valarray<ssize_t>& vec, size_t num);
            void PutDiffSCount(const std::valarray<ssize_t>& vec, size_t num);
            void PutSCount(const std::valarray<ssize_t>& vec)
               { PutCount(vec.size()); PutSCount(vec,vec.size()); }
            void PutDiffSCount(const std::valarray<ssize_t>& vec)
               { PutCount(vec.size()); PutDiffSCount(vec,vec.size()); }
#if defined(WITH_INTMAX_T) && !defined(SIXTY_FOUR_BITS)
            void PutSCount(const std::valarray<intmax_t>& vec, size_t num);
            void PutDiffSCount(const std::valarray<intmax_t>& vec, size_t num);
            void PutSCount(const std::valarray<intmax_t>& vec)
               { PutCount(vec.size()); PutSCount(vec,vec.size()); }
            void PutDiffSCount(const std::valarray<intmax_t>& vec)
               { PutCount(vec.size()); PutDiffSCount(vec,vec.size()); }
#endif
#endif

            /// Put operations for 'Int16' item type (signed 16-bit integers, 'SHORT')

            void PutInt16(int16_t s) { if ( rc == 0 )
               put_vector_of_int16(&s,1,iobuf);
               else check_throw("Put"); }
            void PutInt16(int s) { if ( rc == 0 ) put_short(s,iobuf); }
            void PutInt16(const int16_t *vec, size_t num) { if ( rc == 0 )
               put_vector_of_int16(const_cast<int16_t *>(vec),num,iobuf);
               else check_throw("Put"); }
            void PutInt16(const int *vec, size_t num) { if ( rc == 0 )
               put_vector_of_int(const_cast<int *>(vec),num,iobuf);
               else check_throw("Put"); }
#ifdef HAVE_STD_VECTOR
            void PutInt16(const std::vector<int16_t>& vec, size_t num);
            void PutInt16(const std::vector<int16_t>& vec)
               { PutCount(vec.size()); PutInt16(vec,vec.size()); }
            void PutInt16(const std::vector<int>& vec, size_t num);
            void PutInt16(const std::vector<int>& vec)
               { PutCount(vec.size()); PutInt16(vec,vec.size()); }
#endif
#ifdef HAVE_STD_VALARRAY
            void PutInt16(const std::valarray<int16_t>& vec, size_t num);
            void PutInt16(const std::valarray<int16_t>& vec)
               { PutCount(vec.size()); PutInt16(vec,vec.size()); }
            void PutInt16(const std::valarray<int>& vec, size_t num);
            void PutInt16(const std::valarray<int>& vec)
               { PutCount(vec.size()); PutInt16(vec,vec.size()); }
#endif

            /// Put operations for 'Uint16' item type (unsigned 16-bit integers)

            void PutUint16(uint16_t us) { if ( rc == 0 )
               { put_vector_of_uint16(&us,1,iobuf); }
               else { check_throw("Put"); } }
            /// The 'size_t' to uint16_t conversion gets extra checking since it
            /// is used in many vector length values and silently chopping off here
            /// would easily lead to data corruption.
            inline void PutUint16(size_t us) { if (us > 65535) 
               { throw std::range_error("PutUint16(size_t): value too large for uint16_t data"); }
               uint16_t us16 = us; put_vector_of_uint16(&us16,1,iobuf); }
            void PutUint16(const uint16_t *vec, size_t num) { if ( rc == 0 )
               { put_vector_of_uint16(const_cast<uint16_t *>(vec),num,iobuf); }
               else { check_throw("Put"); } }
#ifdef HAVE_STD_VECTOR
            void PutUint16(const std::vector<uint16_t>& vec, size_t num);
            void PutUint16(const std::vector<uint16_t>& vec)
               { PutCount(vec.size()); PutUint16(vec,vec.size()); }
            void PutUint16(const std::vector<unsigned int>& vec, size_t num);
            void PutUint16(const std::vector<unsigned int>& vec)
               { PutCount(vec.size()); PutUint16(vec,vec.size()); }
#endif
#ifdef HAVE_STD_VALARRAY
            void PutUint16(const std::valarray<uint16_t>& vec, size_t num);
            void PutUint16(const std::valarray<uint16_t>& vec)
               { PutCount(vec.size()); PutUint16(vec,vec.size()); }
            void PutUint16(const std::valarray<unsigned int>& vec, size_t num);
            void PutUint16(const std::valarray<unsigned int>& vec)
               { PutCount(vec.size()); PutUint16(vec,vec.size()); }
#endif

            /// Put operations for 'Int32' item type (signed 32-bit integers, 'LONG')

            void PutInt32(int32_t i) { if ( rc == 0 )
               { put_int32(i,iobuf); }
               else { check_throw("Put"); } }
            void PutInt32(const int32_t *vec, size_t num) { if ( rc == 0 )
               { put_vector_of_int32(const_cast<int32_t *>(vec),num,iobuf); }
               else { check_throw("Put"); } }
            void PutInt32(const long *vec, size_t num) { if ( rc == 0 )
               { put_vector_of_long(const_cast<long *>(vec),num,iobuf); }
               else { check_throw("Put"); } }
#ifdef HAVE_STD_VECTOR
            void PutInt32(const std::vector<int32_t>& vec, size_t num);
            void PutInt32(const std::vector<int32_t>& vec)
               { PutCount(vec.size()); PutInt32(vec,vec.size()); }
            void PutInt32(const std::vector<long>& vec, size_t num);
            void PutInt32(const std::vector<long>& vec)
               { PutCount(vec.size()); PutInt32(vec,vec.size()); }
#endif
#ifdef HAVE_STD_VALARRAY
            void PutInt32(const std::valarray<int32_t>& vec, size_t num);
            void PutInt32(const std::valarray<int32_t>& vec)
               { PutCount(vec.size()); PutInt32(vec,vec.size()); }
            void PutInt32(const std::valarray<long>& vec, size_t num);
            void PutInt32(const std::valarray<long>& vec)
               { PutCount(vec.size()); PutInt32(vec,vec.size()); }
#endif

            /// Put operations for 'Uint32' item type (unsigned 32-bit integers)

            void PutUint32(uint32_t ui) { if ( rc == 0 )
               { put_uint32(ui,iobuf); }
               else { check_throw("Put"); } }
            void PutUint32(const uint32_t *vec, size_t num) { if ( rc == 0 )
               { put_vector_of_uint32(const_cast<uint32_t *>(vec),num,iobuf); }
               else { check_throw("Put"); } }
#ifdef HAVE_64BIT_INT
            /// The 'uint_t' to uint32_t conversion gets extra checking 
            /// on machines with 64 bit integers since it is used in many
            /// vector length values (as a size_t) and silently 
            /// chopping off here would easily lead to data corruption
            /// (well actually only with really huge data blocks).
            inline void PutUint32(uint64_t us) { if (us > UINT_MAX) 
               { throw std::range_error("PutUint32(unit64_t): value too large for uint32_t data"); }
               uint32_t us32 = us; put_vector_of_uint32(&us32,1,iobuf); }
#endif
#ifdef HAVE_STD_VECTOR
            void PutUint32(const std::vector<uint32_t>& vec, size_t num);
            void PutUint32(const std::vector<uint32_t>& vec)
               { PutCount(vec.size()); PutUint32(vec,vec.size()); }
            void PutUint32(const std::vector<unsigned long>& vec, size_t num);
            void PutUint32(const std::vector<unsigned long>& vec)
               { PutCount(vec.size()); PutUint32(vec,vec.size()); }
#endif
#ifdef HAVE_STD_VALARRAY
            void PutUint32(const std::valarray<uint32_t>& vec, size_t num);
            void PutUint32(const std::valarray<uint32_t>& vec)
               { PutCount(vec.size()); PutUint32(vec,vec.size()); }
            void PutUint32(const std::valarray<unsigned long>& vec, size_t num);
            void PutUint32(const std::valarray<unsigned long>& vec)
               { PutCount(vec.size()); PutUint32(vec,vec.size()); }
#endif

#ifdef HAVE_64BIT_INT

            /// Put operations for 'Int64' item type (signed 64-bit integers)

            void PutInt64(int64_t ll) { if ( rc == 0 )
               { put_vector_of_int64(&ll,1,iobuf); }
               else { check_throw("Put"); } }
            void PutInt64(const int64_t *vec, size_t num) { if ( rc == 0 )
               { put_vector_of_int64(const_cast<int64_t *>(vec),num,iobuf); }
               else { check_throw("Put"); } }
#ifdef HAVE_STD_VECTOR
            void PutInt64(const std::vector<int64_t>& vec, size_t num);
            void PutInt64(const std::vector<int64_t>& vec)
               { PutCount(vec.size()); PutInt64(vec,vec.size()); }
#endif
#ifdef HAVE_STD_VALARRAY
            void PutInt64(const std::valarray<int64_t>& vec, size_t num);
            void PutInt64(const std::valarray<int64_t>& vec)
               { PutCount(vec.size()); PutInt64(vec,vec.size()); }
#endif

            /// Put operations for 'Uint64' item type (unsigned 64-bit integers)

            void PutUint64(uint64_t ull) { if ( rc == 0 )
               { put_vector_of_uint64(&ull,1,iobuf); }
               else { check_throw("Put"); } }
            void PutUint64(const uint64_t *vec, size_t num) { if ( rc == 0 )
               { put_vector_of_uint64(const_cast<uint64_t *>(vec),num,iobuf); }
               else { check_throw("Put"); } }
#ifdef HAVE_STD_VECTOR
            void PutUint64(const std::vector<uint64_t>& vec, size_t num);
            void PutUint64(const std::vector<uint64_t>& vec)
               { PutCount(vec.size()); PutUint64(vec,vec.size()); }
#endif
#ifdef HAVE_STD_VALARRAY
            void PutUint64(const std::valarray<uint64_t>& vec, size_t num);
            void PutUint64(const std::valarray<uint64_t>& vec)
               { PutCount(vec.size()); PutUint64(vec,vec.size()); }
#endif
#endif

            /// Put operations for 'Real' item type (32-bit IEEE float)

            void PutReal(float f) { if ( rc == 0 )
               { put_real(f,iobuf); }
               else { check_throw("Put"); } }
            void PutReal(const float *vec, size_t num) { if ( rc == 0 )
               { put_vector_of_float(const_cast<float *>(vec),num,iobuf); }
               else { check_throw("Put"); } }
#ifdef HAVE_STD_VECTOR
            void PutReal(const std::vector<float>& vec, size_t num);
            void PutReal(const std::vector<float>& vec)
               { PutCount(vec.size()); PutReal(vec,vec.size()); }
#endif
#ifdef HAVE_STD_VALARRAY
            void PutReal(const std::valarray<float>& vec, size_t num);
            void PutReal(const std::valarray<float>& vec)
               { PutCount(vec.size()); PutReal(vec,vec.size()); }
#endif
            void PutReal(double d) { if ( rc == 0 )
               { put_real(d,iobuf); }
               else { check_throw("Put"); } }
            void PutReal(const double *vec, size_t num) { if ( rc == 0 )
               { put_vector_of_real(const_cast<double *>(vec),num,iobuf); }
               else { check_throw("Put"); } }
#ifdef HAVE_STD_VECTOR
            void PutReal(const std::vector<double>& vec, size_t num);
            void PutReal(const std::vector<double>& vec)
               { PutCount(vec.size()); PutReal(vec,vec.size()); }
#endif
#ifdef HAVE_STD_VALARRAY
            void PutReal(const std::valarray<double>& vec, size_t num);
            void PutReal(const std::valarray<double>& vec)
               { PutCount(vec.size()); PutReal(vec,vec.size()); }
#endif

            /// Put operations for 'Double' item type (64-bit IEEE double)

            void PutDouble(double d) { if ( rc == 0 )
               { put_double(d,iobuf); }
               else { check_throw("Put"); } }
            void PutDouble(const double *vec, size_t num) { if ( rc == 0 )
               { put_vector_of_double(const_cast<double *>(vec),num,iobuf); }
               else { check_throw("Put"); } }
#ifdef HAVE_STD_VECTOR
            void PutDouble(const std::vector<double>& vec, size_t num);
            void PutDouble(const std::vector<double>& vec)
               { PutCount(vec.size()); PutDouble(vec,vec.size()); }
#endif
#ifdef HAVE_STD_VALARRAY
            void PutDouble(const std::valarray<double>& vec, size_t num);
            void PutDouble(const std::valarray<double>& vec)
               { PutCount(vec.size()); PutDouble(vec,vec.size()); }
#endif

            /// Put operations for 'Sfloat' item type (16-bit OpenGL/IEEE 754-2008 float)

            void PutSfloat(double d) { if ( rc == 0 )
               { put_sfloat(d,iobuf); }
               else { check_throw("Put");}  }
            void PutSfloat(const float *vec, size_t num) { if ( rc == 0 )
               {
                  for ( size_t i=0; i<num; i++ )
                     put_sfloat(vec[i],iobuf);
               }
               else { check_throw("Put"); } }
            void PutSfloat(const double *vec, size_t num) { if ( rc == 0 )
               {
                  for ( size_t i=0; i<num; i++ )
                     put_sfloat(vec[i],iobuf);
               }
               else { check_throw("Put"); } }
#ifdef HAVE_STD_VECTOR
            void PutSfloat(const std::vector<float>& vec, size_t num);
            void PutSfloat(const std::vector<double>& vec, size_t num);
            void PutSfloat(const std::vector<float>& vec)
               { PutCount(vec.size()); PutSfloat(vec,vec.size()); }
            void PutSfloat(const std::vector<double>& vec)
               { PutCount(vec.size()); PutSfloat(vec,vec.size()); }
#endif
#ifdef HAVE_STD_VALARRAY
            void PutSfloat(const std::valarray<float>& vec, size_t num);
            void PutSfloat(const std::valarray<double>& vec, size_t num);
            void PutSfloat(const std::valarray<float>& vec)
               { PutCount(vec.size()); PutSfloat(vec,vec.size()); }
            void PutSfloat(const std::valarray<double>& vec)
               { PutCount(vec.size()); PutSfloat(vec,vec.size()); }
#endif


            /// Put operations for 'Half' type (same as 'Sfloat', 16-bit OpenGL/IEEE 754-2008 float)

#ifdef HAVE_HALF
            void PutHalf(half h) { if ( rc == 0 )
               { uint16_t *p = (uint16_t *) &h; put_vector_of_uint16(p, 1, iobuf); }
               else { check_throw("Put");}  }
            void PutHalf(const half *vec, size_t num) { if ( rc == 0 )
               {
                  const uint16_t *p = (const uint16_t *) vec;
                  put_vector_of_uint16(p, num, iobuf);
               }
               else { check_throw("Put"); } }
#ifdef HAVE_STD_VECTOR
            void PutHalf(const std::vector<half>& vec, size_t num);
            void PutHalf(const std::vector<half>& vec)
               { PutCount(vec.size()); PutHalf(vec,vec.size()); }
#endif
#ifdef HAVE_STD_VALARRAY
            void PutHalf(const std::valarray<half>& vec, size_t num);
            void PutHalf(const std::valarray<half>& vec)
               { PutCount(vec.size()); PutHalf(vec,vec.size()); }
#endif
#endif


            /// Put operations for strings of any length.
            /// (Note: encoding different from old C method put_string().)

            void PutString(const char *s);
            void PutString16(const char *s);
            void PutString32(const char *s);
#ifdef HAVE_STD_STRING
            void PutString(const std::string& s);
            void PutString16(const std::string& s);
            void PutString32(const std::string& s);
#ifdef HAVE_STD_VECTOR
            void PutString(const std::vector<std::string>& vec, size_t num);
            void PutString(const std::vector<std::string>& vec)
               { PutCount(vec.size()); PutString(vec,vec.size()); }
#endif
#ifdef HAVE_STD_VALARRAY
            void PutString(const std::valarray<std::string>& vec, size_t num);
            void PutString(const std::valarray<std::string>& vec)
               { PutCount(vec.size()); PutString(vec,vec.size()); }
#endif
#endif

      };  // End of Item (sub-) class

      friend class eventio::EventIO::Item;

      private:
         IO_BUFFER *iobuf;
         eventio::EventIO::Item *toplevel;
         IO_ITEM_HEADER search_header;
         std::string input_fname;
         std::string output_fname;
         bool local_input;
         bool local_output;
         bool throw_on_error;
         bool external_buffer;

      public:
         EventIO(size_t initial_size=65536, size_t max_size=100000000);
         EventIO(IO_BUFFER *bf);
         EventIO(const EventIO& eventio);
         ~EventIO(void);
         EventIO& operator= (const EventIO& eventio);

         /// For more tricky things you can still access the underlying I/O buffer.
         IO_BUFFER *Buffer(void) { return iobuf; }
         const IO_BUFFER *Buffer(void) const { return iobuf; }
         /// Return data size (including 16 or 20 byte header)
         size_t size(void) const { if ( iobuf->item_level < 0 ) return 0;
             return iobuf->item_length[0] + 16 + (iobuf->item_extension[0]?4:0); }
         size_t Length(void) const { return size(); }
         /// Open input file
         int OpenInput(const char *fname);
         int OpenInput(const std::string& fname) { return OpenInput(fname.c_str()); }
         int OpenInput(FILE *f);
         /// Open Output file
         int OpenOutput(const char *fname, const char *mode);
         int OpenOutput(const std::string& fname, const std::string& mode)
         { return OpenOutput(fname.c_str(), mode.c_str()); }
         int OpenOutput(const char *fname, bool append=false)
         { return OpenOutput(fname, append?"a":"w"); }
         int OpenOutput(const std::string& fname, bool append=false)
         { return OpenOutput(fname.c_str(), append); }
         int OpenOutput(FILE *f);
         int CloseInput(void);
         int CloseOutput(void);
         /// Input via user function rather than file.
         int OpenFunction(IO_USER_FUNCTION ufunc);
         int CloseFunction(void);
         bool HaveOutput(void) const;
         bool HaveInput(void) const;

         /// Append one full I/O block at the current writing position into another one.
         int Append(const EventIO& ev2);
         /// Copy a sub-item to another I/O buffer as top-level item.
         int Copy(const EventIO::Item& item);

         /// The item type found by the last Find() call.
         unsigned long ItemType(void) const { return search_header.type; }
         /// The item ident found by the last Find() call.
         unsigned int ItemVersion(void) const { return search_header.version; }
         /// The item ident found by the last Find() call.
         long ItemIdent(void) const { return search_header.ident; }
         /// The length of the data in the item found by the last Find() call.
         size_t ItemLength(void) const { if (iobuf != 0) { return iobuf->item_length[0]; }
               else { return 0; } }

         int Find(void);
         int Read(void);
         int Skip(void);
         int List(int verbosity=0);
         int Write(void);

         bool SetThrow(bool on=true) { bool t=throw_on_error; 
            throw_on_error=on; return t; }
         bool SetExtended(bool on=true) { bool t=iobuf->extended;
            iobuf->extended=on; return t; }
   };

   typedef class EventIO::Item EventIO_Item;

}  // End of namespace eventio

#endif

