/* ============================================================================

   Copyright (C) 2003, 2009, 2010 Konrad Bernloehr (Konrad Bernl&ouml;hr)

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

/** @file EventIO.cc
    @short Implementation of methods for the C++ interface to the eventio data format.
    
    @author  Konrad Bernloehr
    @date Initial release: April 2003
    $Date: 2014/11/17 15:44:00 $
    $Revision: 1.33 $
*/

#include <string>
#include <valarray>
#include <vector>
#include <cstdio>
#include <cstring>
#include <iostream>
#include "fileopen.h"
#include "EventIO.hh"

namespace eventio
{

// ===================== EventIO methods ==========================

/// Constructor with initial and maximum sizes for underlying buffer.
///
/// @param initial_size The initial size of the I/O buffer in bytes.
/// @param max_size The maximum size to which the I/O buffer can be expanded.

EventIO::EventIO (size_t initial_size, size_t max_size) :
   iobuf(0), toplevel(0), input_fname(""), output_fname(""),
   local_input(false), local_output(false), throw_on_error(false),
   external_buffer(false)
{
   if ( max_size < initial_size )
      max_size = initial_size;
   iobuf = allocate_io_buffer(initial_size);
   iobuf->max_length = max_size;
   search_header.type = 0;
   search_header.version = 0;
   search_header.can_search = 0;
   search_header.level = 0;
   search_header.ident = 0;
   search_header.use_extension = 0;
}

/// Constructor that wraps around a pre-existing IO_BUFFER.
/// 
/// It is certainly a bad idea to do this with IO_BUFFERs that are
/// being processed by C (or C-style) code since the chain of parent
/// and top-level items would be incomplete. Thus this constructor
/// should better be applied to IO_BUFFERs at the top-level and
/// the resulting EventIO either destroyed or returned to top-level
/// (e.g. via Done()) before C(-style) code takes over again.

EventIO::EventIO (IO_BUFFER *bf) : iobuf(bf), toplevel(0), 
   input_fname(""), output_fname(""),
   local_input(false), local_output(false), throw_on_error(false),
   external_buffer(true)
{
}

/// The copy constructor is something that might work to some extend, but
/// sooner or later may result in corrupted data, in particular when
/// multi-threaded programs try to write to the same file.

EventIO::EventIO (const EventIO& eventio) :
   iobuf(0), toplevel(0), input_fname(""), output_fname(""),
   local_input(false), local_output(false), throw_on_error(false),
   external_buffer(false)
{ 
   throw std::logic_error(std::string("I/O buffers should not be copy-constructed.\n"));
}

/// The same restrictions as for the copy constructor also applies to
/// the assignment operator.

EventIO& EventIO::operator= (const EventIO& eventio)
{
   throw std::logic_error(std::string("I/O buffers should not be copied.\n")); 
   return *this;
}

/// Destructor takes care of finishing any active item, 
/// closing input and output files and releasing the
/// underlying I/O buffer (unless the I/O buffer already existed 
/// beforehand).

EventIO::~EventIO (void) 
{
   if ( toplevel != 0 )
   {
      // Tell the toplevel Item and all of its descendants that it is done.
      toplevel->Done();
      // Even though the Item objects may persist, they have no buffer any more.
      for (EventIO::Item *next=toplevel; next!=0; next=next->child)
         next->iobuf = 0;
   }
   toplevel = 0;
   // Close any remaining input and output files.
   CloseInput();
   CloseOutput();
   // Free the underlying I/O buffer.
   if ( iobuf != 0 && !external_buffer )
      free_io_buffer(iobuf);
}

/// Open input file locally (this object takes responsibility)
/// for closing it afterwards.
///
/// @param fname The name of the input file or "-" for standard input.
/// @return 0 (OK), -1 (opening the input file failed)

int EventIO::OpenInput (const char *fname) 
{
   CloseInput();
   local_input = true;
   if ( fname == 0 )
      return 0; // ???
   input_fname = fname;
   if ( strcmp(fname,"-") == 0 )
   {
      iobuf->input_file = stdin;
      local_input = false;
      return 0;
   }
   iobuf->input_file = fileopen(fname,"r"); 
   if ( iobuf->input_file == 0 )
   {
      perror(fname);
      std::cerr << "EventIO::OpenInput(\"" << fname << "\",\"" << "r" << "\") failed.\n";
      return -1;
   }
   return 0; 
}

/// Use externally opened input file or pipe (to be closed externally afterwards).
///
/// @param f A FILE pointer for the input file.
/// @return 0 (OK), -1 (NULL file pointer passed)

int EventIO::OpenInput (FILE *f)
{
   CloseInput();
   local_input = false;
   if ( f == 0 )
      return -1;
   iobuf->input_file = f;
   return 0;
}

/// Close any still open input file for which this object is responsible.
///
/// @return 0 (OK), -1 (error on closing file or pipe)

int EventIO::CloseInput (void)
{
   int rc = 0;
   if ( iobuf->user_function != 0 )
      CloseFunction();
   if ( local_input && iobuf->input_file != 0 )
   {
      rc = fileclose(iobuf->input_file);
      iobuf->input_file = 0;
      local_input = false;
   }
   input_fname = "";
   return rc;
}

/// Open output file locally (this object takes responsibility)
/// for closing it afterwards.
///
/// @param fname The name of the file to be used for output or "-" for standard output.
/// @param mode The output mode like "w" or "a" (see 'man fopen').
/// @return 0 (OK), -1 (opening the output file failed)

int EventIO::OpenOutput (const char *fname, const char *mode) 
{
   CloseOutput();
   local_output = true;
   if ( fname == 0 )
      return 0; // ???
   output_fname = fname;
   if ( strcmp(fname,"-") == 0 )
   {
      iobuf->output_file = stdout;
      local_output = false;
      return 0;
   }
   iobuf->output_file = fileopen(fname,mode);
   if ( iobuf->output_file == 0 )
   {
      perror(fname);
      std::cerr << "EventIO::OpenOutput(\"" << fname << "\",\"" << mode << "\") failed.\n";
      return -1;
   }
   return 0; 
}

/// Use externally opened output file or pipe (to be closed externally afterwards).
///
/// @param f A FILE pointer for the output file.
/// @return 0 (OK), -1 (NULL file pointer passed)

int EventIO::OpenOutput (FILE *f)
{
   CloseOutput();
   local_output = false;
   if ( f == 0 )
      return -1;
   iobuf->output_file = f;
   return 0;
}

/// Close any still open output file or pipe for which this object is responsible.
///
/// @return 0 (OK), -1 (error on closing file or pipe)

int EventIO::CloseOutput (void)
{
   int rc = 0;
   if ( iobuf->user_function != 0 )
      CloseFunction();
   if ( local_output && iobuf->output_file != 0 )
   {
      rc = fileclose(iobuf->output_file);
      iobuf->output_file = 0;
      local_output = false;
   }
   output_fname = "";
   return rc;
}

/// Use the external function method for input and/or output.
///
/// @param ufunc A user function providing the necessary functionality.
/// @return 0 (always)

int EventIO::OpenFunction (IO_USER_FUNCTION ufunc)
{
   CloseInput();
   CloseOutput();
   iobuf->user_function = ufunc;
   return 0;
}

/// Finish with the extern function method for input and/or output.
///
///  @return 0 (OK), -1 (error)

int EventIO::CloseFunction (void)
{
   if ( iobuf->user_function != 0 )
   {
      iobuf->user_function(0,0L,1); // Flush
      iobuf->user_function = 0;
   }
   return 0;
}

/// Check if we have any kind of output assigned.
///
/// @return true if any file or pipe or function is assigned for output.

bool EventIO::HaveOutput(void) const
{
   if ( iobuf->output_file != 0 ||
        iobuf->output_fileno != -1 ||
        iobuf->user_function != 0 )
      return true;
   else
      return false;
}

/// Check if we have any kind of input assigned.
///
/// @return true if any file or pipe or function is assigned for input.

bool EventIO::HaveInput(void) const
{
   if ( iobuf->input_file != 0 ||
        iobuf->input_fileno != -1 ||
        iobuf->user_function != 0 )
      return true;
   else
      return false;
}

/// Find the next toplevel item in the input stream (file etc.)

int EventIO::Find (void)
{
   if ( toplevel != 0 )
      toplevel->Done();
   return find_io_block(iobuf,&search_header);
}

/// Read the next toplevel item from the input stream to the buffer.

int EventIO::Read (void)
{
   return read_io_block(iobuf,&search_header);
}

/// Skip the next toplevel item on the input stream.

int EventIO::Skip (void)
{
   return skip_io_block(iobuf,&search_header);
}

/// List the next toplevel item on the input stream.

int EventIO::List (int verbosity)
{
   return list_io_blocks(iobuf, verbosity);
}

/// Explicitly write the current toplevel item to the output stream, if any.

int EventIO::Write (void)
{
   return write_io_block(iobuf);
}

/// Append one full I/O block at the current writing position into another one.

int EventIO::Append(const EventIO& ev2)
{
   IO_ITEM_HEADER dummy;
   return append_io_block_as_item(iobuf,&dummy,ev2.Buffer()->buffer,ev2.size());
}

/// Copy a sub-item to another I/O buffer as top-level item.

int EventIO::Copy(const EventIO::Item& item)
{
   if ( iobuf == 0 || item.iobuf == 0 )
      return -1;
   return copy_item_to_io_block (iobuf,item.iobuf,&item.item_header);
}

// ===================== EventIO::Item methods ===========================

/// Item constuctor for toplevel item takes the EventIO buffer as first argument.
///
/// There are three methods available here:
/// a) "get" for reading from top-level of an I/O buffer,
/// b) "put" for writing to an I/O level at top level (after finishing
///     anything that may still be active in it currently),
/// c) "append" for appending to an I/O buffer at its current position.
///    If the I/O buffer has nothing active, this is identical to "put".

EventIO::Item::Item (EventIO &ev, 
   const char *method, size_t type, size_t version, long ident, 
      bool user_flag, bool extended) :
   iobuf(ev.Buffer()), rc(-1), 
   put_flag(false), throw_on_error(ev.throw_on_error), 
   done(false), active(false),
   child(0), orig_parent(0), toplevel_parent(0)
{
   item_header.type = 0;
   if ( method == 0 )
      return;
   if ( ev.toplevel != 0 )
   {
      if ( strcmp(method,"append") != 0 )
      {
         ev.toplevel->Done();
         ev.toplevel = 0;
      }
      else
      {
         EventIO::Item *parent = ev.toplevel;
         while ( parent->child != 0 && !parent->child->done )
            parent = parent->child;
         parent->child = this;
      }
   }
   else if ( strcmp(method,"append") == 0 )
      method = "put";
   toplevel_parent = &ev;
   if ( iobuf == 0 || 
         ( strcmp(method,"get") != 0 && strcmp(method,"put") != 0 &&
           strcmp(method,"append") != 0 ) )
      return;

   item_header.type = type;
   item_header.version = version;
   item_header.ident = ident;
   item_header.level = 0;
   item_header.can_search = 0;
   item_header.user_flag = 0;
   item_header.use_extension = 0;

   if ( strcmp(method,"get") == 0 )
      rc = get_item_begin(iobuf,&item_header);
   else
   {
      rc = put_item_begin_with_flags(iobuf,&item_header,user_flag,extended);
      if ( rc == 0 ) // Only set item into "put" mode if put_item_begin was OK.
         put_flag = true;
   }
   done = false;
   // Start enforcing a more systematic version number check than
   // currently done in the eventio low-level engine.
   // If asking for version '0' that would mean 'any version'.
   // If asking for a version above zero, all encountered version numbers
   // below that would be assumed OK (item handlers taking care of it).
   // But with higher version numbers than desired, your software
   // is probably out of date with respect to the data being read
   // and should be updated.
   if ( rc == 0 && version > 0 && item_header.version > version )
   {
      std::cerr << "Item type " << type << " version "
          << item_header.version << " is beyond max. version " << version
          << " (the last version for which this software was built).\n";
   }
   else if ( rc == 0 )
      active = true;
   else
      std::cerr << "Error " << rc << " on item " << method << ".\n";
   if ( strcmp(method,"append") != 0 )
      ev.toplevel = this;

   check_throw("item constructor");
}

/// Item constructor for sub-items takes the parent item as first argument.
/// There are two methods available here:
/// a) "get" for reading from an I/O buffer at the level of the parent,
///    after any active children are finished.
/// b) "put" for writing to an I/O level at the level of the parent,
///    after any active children are finished.

EventIO::Item::Item (EventIO::Item &parent, 
   const char *method, size_t type, size_t version, long ident, 
   bool user_flag, bool extended) :
   iobuf(parent.Buffer()), rc(-1), 
   put_flag(false), throw_on_error(parent.throw_on_error), 
   done(false), active(false),
   child(0), orig_parent(&parent), toplevel_parent(0)
{
   item_header.type = 0;
   if ( parent.child != 0 )
   {
      parent.child->Done();
      parent.child = 0;
   }
   if ( iobuf == 0 || (strcmp(method,"get") != 0 && strcmp(method,"put") != 0) )
      return;
   item_header.type = type;
   item_header.version = version;
   item_header.ident = ident;
   item_header.can_search = 0;
   item_header.user_flag = 0;
   item_header.use_extension = 0;
   
   if ( strcmp(method,"get") == 0 )
      rc = get_item_begin(iobuf,&item_header);
   else
   {
      rc = put_item_begin_with_flags(iobuf,&item_header,user_flag,extended);
      if ( rc == 0 ) // Only set item into "put" mode if put_item_begin was OK.
         put_flag = true;
   }
   done = false;
   if ( rc == 0 )
      active = true;
   else
      std::cerr << "Error " << rc << " on item " << method << ".\n";
   parent.child = this;

   check_throw("item constructor");
}

/// The item copy constructor is something that might work to some extend, but
/// sooner or later may result in corrupted data, in particular when
/// multi-threaded programs try to write to the same file.

EventIO::Item::Item (const EventIO::Item& eventio) :
   iobuf(0), rc(-1), 
   put_flag(false), throw_on_error(false), 
   done(false), active(false),
   child(0), orig_parent(0), toplevel_parent(0)
{
   throw std::logic_error(std::string("I/O items should not be copy-constructed.\n"));
}

/// Restrictions to the copy constructor also apply to the assignment operator

EventIO::Item& EventIO::Item::operator= (const EventIO::Item& item /* unused */ )
{ throw std::logic_error(std::string("I/O items should not be copied.\n")); return *this; }

/// Item destructor takes care that all sub-items are finished first.
/// In case of a toplevel item, it unregisters from the EventIO buffer.

EventIO::Item::~Item (void)
{
   // fprintf(stderr,"Delete header for item type %lu\n",item_header.type);
   Done();
   if ( orig_parent && orig_parent->child == this )
      orig_parent->child = 0;
   if ( toplevel_parent != 0 )
      toplevel_parent->toplevel = 0;
}

/// All put/get operations for this item and all sub-items are done.

int EventIO::Item::Done (void)
{
   if ( done || !active )
      return 1;

   if ( child != 0 )
   {
      child->Done();
      child = 0;
   }
   
   done = true;

   if ( iobuf == 0 )
   {
      std::cerr << "No I/O buffer available in Done operation. Hmm...\n";
      return rc = -4;
   }

   if ( rc == 0 && iobuf != 0 )
   {
      // fprintf(stderr,"Done with header for item type %lu\n",item_header.type);
      rc = 1;
      if ( put_flag )
         rc = put_item_end(iobuf,&item_header);
      else
         rc = get_item_end(iobuf,&item_header);
      active = false;
      return rc;
   }
   return 0;
}

/// Search the current item, starting at the current reading position,
/// for the next sub-item of a specific type.
///
/// If the item is marked as searchable, it will check each
/// sub-item (at the next level only) if it is of the desired type.
/// That sub-item can then be immediately opened as a new Eventio::Item.
///
/// @param sub The type of sub-item that we are searching for.
/// @return  0 (O.k., sub-item was found),
///         -1 (error),
///         -2 (no such sub-item),
///         -3 (cannot skip sub-items).

int EventIO::Item::Search (size_t sub)
{
   IO_ITEM_HEADER sub_item_header;
   sub_item_header.type = sub;
   if ( iobuf == 0 )
      return rc = -1;
   return rc = search_sub_item(iobuf,&item_header,&sub_item_header);
}

/// Rewind to beginning of the data area of the current item.
///
/// You can restart searching for specific sub-items again.
///
///  @return  0 (ok), -1 (error)

int EventIO::Item::Rewind ()
{
   if ( iobuf == 0 )
      return rc = -1;
   return rc = rewind_item(iobuf,&item_header);
}

/// Skip a sub-item starting at the current position.
///
/// If the sub-item at the current position is of no further interest,
/// you can skip it without having to decode its contents and go
/// to the next item, if there is any.
///
///  @return  0 (ok), -1 (error)

int EventIO::Item::Skip ()
{
   if ( iobuf == 0 )
      return rc = -1;
   return rc = skip_subitem(iobuf);
}

/// Completely undo getting this item, i.e. rewind to beginning of its header.
/// The item should not be used any more afterwards.
///
///  @return  0 (ok), -1 (error)

int EventIO::Item::Unget ()
{
   child = 0; // Children get irrelevant.
   done = true; // A subsequent Done() call or destructor has nothing to do.
   if ( iobuf == 0 )
      return rc = -1;
   return rc = unget_item(iobuf,&item_header);
}

/// Completely undo putting this item, as if never started.
/// The item should not be used any more afterwards.
///
///  @return  0 (ok), -1 (error)

int EventIO::Item::Unput ()
{
   child = 0; // Children get irrelevant.
   done = true; // A subsequent Done() call or destructor has nothing to do.
   if ( iobuf == 0 )
      return rc = -1;
   return rc = unput_item(iobuf,&item_header);
}

/// Access to type information of next sub-item at current reading position.
///
/// @return  >= 0 (O.k., sub-item type),  -1 (error), -2 (end-of-buffer), -3 (no buffer)

int EventIO::Item::NextSubItemType(void) const
{
   if ( iobuf != 0 )
      return next_subitem_type(iobuf);
   else
      return -3;
}

/// Access to length information of next sub-item at current reading position.
///
/// @return  >= 0 (O.k., sub-item length),  -1 (error), -2 (end-of-buffer), -3 (no buffer)

size_t EventIO::Item::NextSubItemLength(void) const
{
   if ( iobuf != 0 )
      return next_subitem_length(iobuf);
   else
      return 0;
}

/// Access to ID information of next sub-item at current reading position.
///
/// @return  >= 0 (O.k., sub-item ident),  -1 (error), -2 (end-of-buffer), -3 (no buffer)

long EventIO::Item::NextSubItemIdent(void) const
{
   if ( iobuf != 0 )
      return next_subitem_ident(iobuf);
   else
      return -1;
}

/// Show I/O block information to standard output

int EventIO::Item::List(int maxlevel, int verbosity)
{
   return list_sub_items(iobuf, &item_header, maxlevel, verbosity);
}

/// Access to the registered name for this type of I/O block, if available

const char *EventIO::Item::TypeName()
{
   return eventio_registered_typename(item_header.type);
}

/// Access to the registered description for this type of I/O block, if available

const char *EventIO::Item::Description()
{
   return eventio_registered_description(item_header.type);
}

// ------------------------ Get operations ----------------------

/// Get operations for bytes (only unsigned flavour implemented).

void EventIO::Item::GetUint8(std::vector<uint8_t>& vec, size_t num)
{
   check_throw("GetUint8");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint8();
}

/// Get operations for bytes (only unsigned flavour implemented).

void EventIO::Item::GetUint8(std::vector<uint8_t>& vec)
{
   check_throw("GetUint8");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint8();
}

/// Get operations for bytes (only unsigned flavour implemented).

void EventIO::Item::GetUint8(std::valarray<uint8_t>& vec, size_t num)
{
   check_throw("GetUint8");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint8();
}

/// Get operations for bytes (only unsigned flavour implemented).

void EventIO::Item::GetUint8(std::valarray<uint8_t>& vec)
{
   check_throw("GetUint8");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint8();
}

/// Get operations for bools.

void EventIO::Item::GetBool(std::vector<bool>& vec, size_t num)
{
   check_throw("GetBool");

   if ( vec.size() < num )
      vec.resize(num);

   uint8_t v = 0;

   for (size_t i=0, j=0; i<num; ++i, ++j)
   {
      if ( (i&0x07) == 0 )
      {
         v = GetUint8();
         j = 0;
      }
      if ( (v&(1<<j)) != 0 )
         vec[i] = true;
      else
         vec[i] = false;
   }
}

/// Get operations for bools.

void EventIO::Item::GetBool(std::vector<bool>& vec)
{
   check_throw("GetBool");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   uint8_t v = 0;

   for (size_t i=0, j=0; i<num; ++i, ++j)
   {
      if ( (i&0x07) == 0 )
      {
         v = GetUint8();
         j = 0;
      }
      if ( (v&(1<<j)) != 0 )
         vec[i] = true;
      else
         vec[i] = false;
   }
}

/// Get operations for bools.

void EventIO::Item::GetBool(std::valarray<bool>& vec, size_t num)
{
   check_throw("GetBool");

   if ( vec.size() < num )
      vec.resize(num);

   uint8_t v = 0;

   for (size_t i=0, j=0; i<num; ++i, ++j)
   {
      if ( (i&0x07) == 0 )
      {
         v = GetUint8();
         j = 0;
      }
      if ( (v&(1<<j)) != 0 )
         vec[i] = true;
      else
         vec[i] = false;
   }
}

/// Get operations for bools.

void EventIO::Item::GetBool(std::valarray<bool>& vec)
{
   check_throw("GetBool");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   uint8_t v = 0;

   for (size_t i=0, j=0; i<num; ++i, ++j)
   {
      if ( (i&0x07) == 0 )
      {
         v = GetUint8();
         j = 0;
      }
      if ( (v&(1<<j)) != 0 )
         vec[i] = true;
      else
         vec[i] = false;
   }
}

/// Get operations for counts (16 bit variant)

void EventIO::Item::GetCount(uint16_t *vec, size_t num)
{
   check_throw("GetCount");

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount16();
}

#ifdef SIXTY_FOUR_BITS

/// Get operations for counts (32 bit variant)

void EventIO::Item::GetCount(uint32_t *vec, size_t num)
{
   check_throw("GetCount");

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount32();
}

#endif

/// Get operations for counts

void EventIO::Item::GetCount(size_t *vec, size_t num)
{
   check_throw("GetCount");

   for (size_t i=0; i<num; ++i)
      vec[i] = size_t(GetCount());
}

/// Get operations for counts

#if defined(WITH_UINTMAX_T) && !defined(SIXTY_FOUR_BITS)
void EventIO::Item::GetCount(uintmax_t *vec, size_t num)
{
   check_throw("GetCount");

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount();
}
#endif

/// Get operations for counts (16 bit variant).

void EventIO::Item::GetCount(std::vector<uint16_t>& vec, size_t num)
{
   check_throw("GetCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount16();
}

#ifdef SIXTY_FOUR_BITS

/// Get operations for counts (32 bit variant).

void EventIO::Item::GetCount(std::vector<uint32_t>& vec, size_t num)
{
   check_throw("GetCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount32();
}
#endif

/// Get operations for counts.

void EventIO::Item::GetCount(std::vector<size_t>& vec, size_t num)
{
   check_throw("GetCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = size_t(GetCount());
}

/// Get operations for counts.

#if defined(WITH_UINTMAX_T) && !defined(SIXTY_FOUR_BITS)
void EventIO::Item::GetCount(std::vector<uintmax_t>& vec, size_t num)
{
   check_throw("GetCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount();
}

#endif

/// Get operations for counts (16 bit variant).

void EventIO::Item::GetCount(std::vector<uint16_t>& vec)
{
   check_throw("GetCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount16();
}

#ifdef SIXTY_FOUR_BITS

/// Get operations for counts (32 bit variant).

void EventIO::Item::GetCount(std::vector<uint32_t>& vec)
{
   check_throw("GetCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount32();
}

#endif

/// Get operations for counts.

void EventIO::Item::GetCount(std::vector<size_t>& vec)
{
   check_throw("GetCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = size_t(GetCount());
}

/// Get operations for counts.

#if defined(WITH_UINTMAX_T) && !defined(SIXTY_FOUR_BITS)
void EventIO::Item::GetCount(std::vector<uintmax_t>& vec)
{
   check_throw("GetCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount();
}
#endif

/// Get operations for counts (16 bit variant).

void EventIO::Item::GetCount(std::valarray<uint16_t>& vec, size_t num)
{
   check_throw("GetCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount16();
}

#ifdef SIXTY_FOUR_BITS

/// Get operations for counts (32 bit variant).

void EventIO::Item::GetCount(std::valarray<uint32_t>& vec, size_t num)
{
   check_throw("GetCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount32();
}

#endif

/// Get operations for counts.

void EventIO::Item::GetCount(std::valarray<size_t>& vec, size_t num)
{
   check_throw("GetCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = size_t(GetCount());
}

/// Get operations for counts.

#if defined(WITH_UINTMAX_T) && !defined(SIXTY_FOUR_BITS)
void EventIO::Item::GetCount(std::valarray<uintmax_t>& vec, size_t num)
{
   check_throw("GetCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount();
}
#endif

/// Get operations for counts (16 bit variant).

void EventIO::Item::GetCount(std::valarray<uint16_t>& vec)
{
   check_throw("GetCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount16();
}

#ifdef SIXTY_FOUR_BITS

/// Get operations for counts (32 bit variant).

void EventIO::Item::GetCount(std::valarray<uint32_t>& vec)
{
   check_throw("GetCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount32();
}

#endif

/// Get operations for counts.

void EventIO::Item::GetCount(std::valarray<size_t>& vec)
{
   check_throw("GetCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount();
}

/// Get operations for counts.

#if defined(WITH_UINTMAX_T) && !defined(SIXTY_FOUR_BITS)
void EventIO::Item::GetCount(std::valarray<uintmax_t>& vec)
{
   check_throw("GetCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetCount();
}
#endif

/// Get operations for signed counts (16 bit variant)

void EventIO::Item::GetSCount(int16_t *vec, size_t num)
{
   check_throw("GetSCount");

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount16();
}

void EventIO::Item::GetDiffSCount(int16_t *vec, size_t num)
{
   check_throw("GetDiffSCount");

   if ( num > 0 )
      vec[0] = GetSCount16();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount16();
}

#ifdef SIXTY_FOUR_BITS

/// Get operations for signed counts (32 bit variant)

void EventIO::Item::GetSCount(int32_t *vec, size_t num)
{
   check_throw("GetSCount");

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount32();
}

void EventIO::Item::GetDiffSCount(int32_t *vec, size_t num)
{
   check_throw("GetDiffSCount");

   if ( num > 0 )
      vec[0] = GetSCount32();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount32();
}

#endif

/// Get operations for signed counts

void EventIO::Item::GetSCount(ssize_t *vec, size_t num)
{
   check_throw("GetSCount");

   for (size_t i=0; i<num; ++i)
      vec[i] = ssize_t(GetSCount());
}

void EventIO::Item::GetDiffSCount(ssize_t *vec, size_t num)
{
   check_throw("GetDiffSCount");

   if ( num > 0 )
      vec[0] = ssize_t(GetSCount());
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + ssize_t(GetSCount());
}

/// Get operations for signed counts

#if defined(WITH_INTMAX_T) && !defined(SIXTY_FOUR_BITS)
void EventIO::Item::GetSCount(intmax_t *vec, size_t num)
{
   check_throw("GetSCount");

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount();
}

void EventIO::Item::GetDiffSCount(intmax_t *vec, size_t num)
{
   check_throw("GetDiffSCount");

   if ( num > 0 )
      vec[0] = GetSCount();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount();
}
#endif

/// Get operations for signed counts (16 bit variant).

void EventIO::Item::GetSCount(std::vector<int16_t>& vec, size_t num)
{
   check_throw("GetSCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount16();
}

void EventIO::Item::GetDiffSCount(std::vector<int16_t>& vec, size_t num)
{
   check_throw("GetDiffSCount");

   if ( vec.size() < num )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = GetSCount16();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount16();
}

#ifdef SIXTY_FOUR_BITS

/// Get operations for signed counts (32 bit variant).

void EventIO::Item::GetSCount(std::vector<int32_t>& vec, size_t num)
{
   check_throw("GetSCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount32();
}

void EventIO::Item::GetDiffSCount(std::vector<int32_t>& vec, size_t num)
{
   check_throw("GetDiffSCount");

   if ( vec.size() < num )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = GetSCount32();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount32();
}

#endif

/// Get operations for signed counts.

void EventIO::Item::GetSCount(std::vector<ssize_t>& vec, size_t num)
{
   check_throw("GetSCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = ssize_t(GetSCount());
}

void EventIO::Item::GetDiffSCount(std::vector<ssize_t>& vec, size_t num)
{
   check_throw("GetDiffSCount");

   if ( vec.size() < num )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = ssize_t(GetSCount());
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + ssize_t(GetSCount());
}

/// Get operations for signed counts.

#if defined(WITH_INTMAX_T) && !defined(SIXTY_FOUR_BITS)
void EventIO::Item::GetSCount(std::vector<intmax_t>& vec, size_t num)
{
   check_throw("GetSCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount();
}

void EventIO::Item::GetDiffSCount(std::vector<intmax_t>& vec, size_t num)
{
   check_throw("GetDiffSCount");

   if ( vec.size() < num )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = GetSCount();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount();
}
#endif

/// Get operations for signed counts (16 bit variant).

void EventIO::Item::GetSCount(std::vector<int16_t>& vec)
{
   check_throw("GetSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount16();
}

void EventIO::Item::GetDiffSCount(std::vector<int16_t>& vec)
{
   check_throw("GetDiffSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = GetSCount16();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount16();
}

#ifdef SIXTY_FOUR_BITS

/// Get operations for signed counts (32 bit variant).

void EventIO::Item::GetSCount(std::vector<int32_t>& vec)
{
   check_throw("GetSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount32();
}

void EventIO::Item::GetDiffSCount(std::vector<int32_t>& vec)
{
   check_throw("GetDiffSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = GetSCount32();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount32();
}

#endif

/// Get operations for signed counts.

void EventIO::Item::GetSCount(std::vector<ssize_t>& vec)
{
   check_throw("GetSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = ssize_t(GetSCount());
}

void EventIO::Item::GetDiffSCount(std::vector<ssize_t>& vec)
{
   check_throw("GetDiffSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = ssize_t(GetSCount());
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + ssize_t(GetSCount());
}

/// Get operations for signed counts.

#if defined(WITH_INTMAX_T) && !defined(SIXTY_FOUR_BITS)
void EventIO::Item::GetSCount(std::vector<intmax_t>& vec)
{
   check_throw("GetSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount();
}

void EventIO::Item::GetDiffSCount(std::vector<intmax_t>& vec)
{
   check_throw("GetDiffSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = GetSCount();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount();
}
#endif

/// Get operations for signed counts (16 bit variant).

void EventIO::Item::GetSCount(std::valarray<int16_t>& vec, size_t num)
{
   check_throw("GetSCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount16();
}

void EventIO::Item::GetDiffSCount(std::valarray<int16_t>& vec, size_t num)
{
   check_throw("GetDiffSCount");

   if ( vec.size() < num )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = GetSCount16();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount16();
}

#ifdef SIXTY_FOUR_BITS

/// Get operations for signed counts (32 bit variant).

void EventIO::Item::GetSCount(std::valarray<int32_t>& vec, size_t num)
{
   check_throw("GetSCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount32();
}

void EventIO::Item::GetDiffSCount(std::valarray<int32_t>& vec, size_t num)
{
   check_throw("GetDiffSCount");

   if ( vec.size() < num )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = GetSCount32();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount32();
}

#endif

/// Get operations for signed counts.

void EventIO::Item::GetSCount(std::valarray<ssize_t>& vec, size_t num)
{
   check_throw("GetSCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = ssize_t(GetSCount());
}

void EventIO::Item::GetDiffSCount(std::valarray<ssize_t>& vec, size_t num)
{
   check_throw("GetDiffSCount");

   if ( vec.size() < num )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = ssize_t(GetSCount());
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + ssize_t(GetSCount());
}

/// Get operations for signed counts.

#if defined(WITH_INTMAX_T) && !defined(SIXTY_FOUR_BITS)
void EventIO::Item::GetSCount(std::valarray<intmax_t>& vec, size_t num)
{
   check_throw("GetSCount");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount();
}

void EventIO::Item::GetDiffSCount(std::valarray<intmax_t>& vec, size_t num)
{
   check_throw("GetDiffSCount");

   if ( vec.size() < num )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = GetSCount();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount();
}
#endif

/// Get operations for signed counts (16 bit variant).

void EventIO::Item::GetSCount(std::valarray<int16_t>& vec)
{
   check_throw("GetSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount16();
}

void EventIO::Item::GetDiffSCount(std::valarray<int16_t>& vec)
{
   check_throw("GetDiffSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = GetSCount16();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount16();
}

#ifdef SIXTY_FOUR_BITS

/// Get operations for signed counts (32 bit variant).

void EventIO::Item::GetSCount(std::valarray<int32_t>& vec)
{
   check_throw("GetSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount32();
}

void EventIO::Item::GetDiffSCount(std::valarray<int32_t>& vec)
{
   check_throw("GetDiffSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = GetSCount32();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount32();
}

#endif

/// Get operations for signed counts.

void EventIO::Item::GetSCount(std::valarray<ssize_t>& vec)
{
   check_throw("GetSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount();
}

void EventIO::Item::GetDiffSCount(std::valarray<ssize_t>& vec)
{
   check_throw("GetSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = GetSCount();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount();
}

/// Get operations for signed counts.

#if defined(WITH_INTMAX_T) && !defined(SIXTY_FOUR_BITS)
void EventIO::Item::GetSCount(std::valarray<intmax_t>& vec)
{
   check_throw("GetSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSCount();
}

void EventIO::Item::GetDiffSCount(std::valarray<intmax_t>& vec)
{
   check_throw("GetDiffSCount");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   if ( num > 0 )
      vec[0] = GetSCount();
   // Except for the first elements, all elements are stored as differences
   for (size_t i=1; i<num; ++i)
      vec[i] = vec[i-1] + GetSCount();
}
#endif

/// Get operations for 'Int16' item type (signed 16-bit integers, 'SHORT')

void EventIO::Item::GetInt16(std::vector<int16_t>& vec, size_t num)
{
   check_throw("GetInt16");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt16();
}

/// Get operations for 'Int16' item type (signed 16-bit integers, 'SHORT')

void EventIO::Item::GetInt16(std::vector<int16_t>& vec)
{
   check_throw("GetInt16");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt16();
}

/// Get operations for 'Int16' item type (signed 16-bit integers, 'SHORT')

void EventIO::Item::GetInt16(std::vector<int>& vec, size_t num)
{
   check_throw("GetInt16");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt16();
}

/// Get operations for 'Int16' item type (signed 16-bit integers, 'SHORT')

void EventIO::Item::GetInt16(std::vector<int>& vec)
{
   check_throw("GetInt16");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt16();
}

/// Get operations for 'Int16' item type (signed 16-bit integers, 'SHORT')

void EventIO::Item::GetInt16(std::valarray<int16_t>& vec, size_t num)
{
   check_throw("GetInt16");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt16();
}

/// Get operations for 'Int16' item type (signed 16-bit integers, 'SHORT')

void EventIO::Item::GetInt16(std::valarray<int16_t>& vec)
{
   check_throw("GetInt16");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt16();
}

/// Get operations for 'Int16' item type (signed 16-bit integers, 'SHORT')

void EventIO::Item::GetInt16(std::valarray<int>& vec, size_t num)
{
   check_throw("GetInt16");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt16();
}

/// Get operations for 'Int16' item type (signed 16-bit integers, 'SHORT')

void EventIO::Item::GetInt16(std::valarray<int>& vec)
{
   check_throw("GetInt16");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt16();
}

/// Get operations for 'Uint16' item type (unsigned 16-bit integers)

void EventIO::Item::GetUint16(std::vector<uint16_t>& vec, size_t num)
{
   check_throw("GetUint16");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint16();
}

/// Get operations for 'Uint16' item type (unsigned 16-bit integers)

void EventIO::Item::GetUint16(std::vector<uint16_t>& vec)
{
   check_throw("GetUint16");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint16();
}

/// Get operations for 'Uint16' item type (unsigned 16-bit integers)

void EventIO::Item::GetUint16(std::valarray<uint16_t>& vec, size_t num)
{
   check_throw("GetUint16");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint16();
}

/// Get operations for 'Uint16' item type (unsigned 16-bit integers)

void EventIO::Item::GetUint16(std::valarray<uint16_t>& vec)
{
   check_throw("GetUint16");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint16();
}

/// Get operations for 'Int32' item type (signed 32-bit integers, 'LONG')

void EventIO::Item::GetInt32(std::vector<int32_t>& vec, size_t num)
{
   check_throw("GetInt32");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt32();
}

/// Get operations for 'Int32' item type (signed 32-bit integers, 'LONG')

void EventIO::Item::GetInt32(std::vector<int32_t>& vec)
{
   check_throw("GetInt32");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt32();
}

/// Get operations for 'Int32' item type (signed 32-bit integers, 'LONG')

void EventIO::Item::GetInt32(std::vector<long>& vec, size_t num)
{
   check_throw("GetInt32");

   if ( vec.size() < num )
      vec.resize(num);
   
   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt32();
}

/// Get operations for 'Int32' item type (signed 32-bit integers, 'LONG')

void EventIO::Item::GetInt32(std::vector<long>& vec)
{
   check_throw("GetInt32");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);
   
   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt32();
}

/// Get operations for 'Int32' item type (signed 32-bit integers, 'LONG')

void EventIO::Item::GetInt32(std::valarray<int32_t>& vec, size_t num)
{
   check_throw("GetInt32");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt32();
}

/// Get operations for 'Int32' item type (signed 32-bit integers, 'LONG')

void EventIO::Item::GetInt32(std::valarray<int32_t>& vec)
{
   check_throw("GetInt32");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt32();
}

/// Get operations for 'Int32' item type (signed 32-bit integers, 'LONG')

void EventIO::Item::GetInt32(std::valarray<long>& vec, size_t num)
{
   check_throw("GetInt32");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt32();
}

/// Get operations for 'Int32' item type (signed 32-bit integers, 'LONG')

void EventIO::Item::GetInt32(std::valarray<long>& vec)
{
   check_throw("GetInt32");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt32();
}

/// Get operations for 'Uint32' item type (unsigned 32-bit integers)

void EventIO::Item::GetUint32(std::vector<uint32_t>& vec, size_t num)
{
   check_throw("GetUint32");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint32();
}

/// Get operations for 'Uint32' item type (unsigned 32-bit integers)

void EventIO::Item::GetUint32(std::vector<uint32_t>& vec)
{
   check_throw("GetUint32");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint32();
}

/// Get operations for 'Uint32' item type (unsigned 32-bit integers)

void EventIO::Item::GetUint32(std::valarray<uint32_t>& vec, size_t num)
{
   check_throw("GetUint32");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint32();
}

/// Get operations for 'Uint32' item type (unsigned 32-bit integers)

void EventIO::Item::GetUint32(std::valarray<uint32_t>& vec)
{
   check_throw("GetUint32");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint32();
}

#ifdef HAVE_64BIT_INT

/// Get operations for 'Int64' item type (signed 64-bit integers)

void EventIO::Item::GetInt64(std::vector<int64_t>& vec, size_t num)
{
   check_throw("GetInt64");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt64();
}

/// Get operations for 'Int64' item type (signed 64-bit integers)

void EventIO::Item::GetInt64(std::vector<int64_t>& vec)
{
   check_throw("GetInt64");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt64();
}

/// Get operations for 'Int64' item type (signed 64-bit integers)

void EventIO::Item::GetInt64(std::valarray<int64_t>& vec, size_t num)
{
   check_throw("GetInt64");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt64();
}

/// Get operations for 'Int64' item type (signed 64-bit integers)

void EventIO::Item::GetInt64(std::valarray<int64_t>& vec)
{
   check_throw("GetInt64");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetInt64();
}

/// Get operations for 'Uint64' item type (unsigned 64-bit integers)

void EventIO::Item::GetUint64(std::vector<uint64_t>& vec, size_t num)
{
   check_throw("GetUint64");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint64();
}

/// Get operations for 'Uint64' item type (unsigned 64-bit integers)

void EventIO::Item::GetUint64(std::vector<uint64_t>& vec)
{
   check_throw("GetUint64");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint64();
}

/// Get operations for 'Uint64' item type (unsigned 64-bit integers)

void EventIO::Item::GetUint64(std::valarray<uint64_t>& vec, size_t num)
{
   check_throw("GetUint64");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint64();
}

/// Get operations for 'Uint64' item type (unsigned 64-bit integers)

void EventIO::Item::GetUint64(std::valarray<uint64_t>& vec)
{
   check_throw("GetUint64");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetUint64();
}
#endif

/// Get operations for 'Real' item type (32-bit IEEE float)

void EventIO::Item::GetReal(std::vector<float>& vec, size_t num)
{
   check_throw("GetReal");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetReal();
}

/// Get operations for 'Real' item type (32-bit IEEE float)

void EventIO::Item::GetReal(std::vector<float>& vec)
{
   check_throw("GetReal");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetReal();
}

/// Get operations for 'Real' item type (32-bit IEEE float)

void EventIO::Item::GetReal(std::valarray<float>& vec, size_t num)
{
   check_throw("GetReal");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetReal();
}

/// Get operations for 'Real' item type (32-bit IEEE float)

void EventIO::Item::GetReal(std::valarray<float>& vec)
{
   check_throw("GetReal");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetReal();
}

/// Get operations for 'Real' item type (32-bit IEEE float)

void EventIO::Item::GetReal(std::vector<double>& vec, size_t num)
{
   check_throw("GetReal");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetReal();
}

/// Get operations for 'Real' item type (32-bit IEEE float)

void EventIO::Item::GetReal(std::vector<double>& vec)
{
   check_throw("GetReal");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetReal();
}

/// Get operations for 'Real' item type (32-bit IEEE float)

void EventIO::Item::GetReal(std::valarray<double>& vec, size_t num)
{
   check_throw("GetReal");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetReal();
}

/// Get operations for 'Real' item type (32-bit IEEE float)

void EventIO::Item::GetReal(std::valarray<double>& vec)
{
   check_throw("GetReal");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetReal();
}

/// Get operations for 'Double' item type (64-bit IEEE double)

void EventIO::Item::GetDouble(std::vector<double>& vec, size_t num)
{
   check_throw("GetDouble");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetDouble();
}

/// Get operations for 'Double' item type (64-bit IEEE double)

void EventIO::Item::GetDouble(std::vector<double>& vec)
{
   check_throw("GetDouble");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetDouble();
}

/// Get operations for 'Double' item type (64-bit IEEE double)

void EventIO::Item::GetDouble(std::valarray<double>& vec, size_t num)
{
   check_throw("GetDouble");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetDouble();
}

/// Get operations for 'Double' item type (64-bit IEEE double)

void EventIO::Item::GetDouble(std::valarray<double>& vec)
{
   check_throw("GetDouble");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetDouble();
}


/// Get operations for 'Sfloat' item type (16-bit OpenGL float)

void EventIO::Item::GetSfloat(std::vector<float>& vec, size_t num)
{
   check_throw("GetSfloat");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = (float)GetSfloat();
}

void EventIO::Item::GetSfloat(std::vector<double>& vec, size_t num)
{
   check_throw("GetSfloat");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSfloat();
}

/// Get operations for 'Sfloat' item type (16-bit OpenGL float)

void EventIO::Item::GetSfloat(std::vector<float>& vec)
{
   check_throw("GetSfloat");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = (float)GetSfloat();
}

void EventIO::Item::GetSfloat(std::vector<double>& vec)
{
   check_throw("GetSfloat");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSfloat();
}

/// Get operations for 'Sfloat' item type (16-bit OpenGL float)

void EventIO::Item::GetSfloat(std::valarray<float>& vec, size_t num)
{
   check_throw("GetSfloat");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = (float)GetSfloat();
}

void EventIO::Item::GetSfloat(std::valarray<double>& vec, size_t num)
{
   check_throw("GetSfloat");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSfloat();
}

/// Get operations for 'Sfloat' item type (16-bit OpenGL float)

void EventIO::Item::GetSfloat(std::valarray<float>& vec)
{
   check_throw("GetSfloat");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = (float)GetSfloat();
}

void EventIO::Item::GetSfloat(std::valarray<double>& vec)
{
   check_throw("GetSfloat");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      vec[i] = GetSfloat();
}


/// Get operation for old-style strings with 16-bit length prefix.

int EventIO::Item::GetString16(std::string& s)
{
   check_throw("GetString");

   s = "";
   size_t n = GetUint16();
   s.resize(n);
   for (size_t i=0; i<n; ++i)
      s[i] = GetUint8();
   return int(n);
}

/// Get operation for old-style strings with 32-bit length prefix.

int EventIO::Item::GetString32(std::string& s)
{
   check_throw("GetString");

   s = "";
   size_t n = GetUint32();
   s.resize(n);
   for (size_t i=0; i<n; ++i)
      s[i] = GetUint8();
   return int(n);
}

/// Get operation for strings of any length.

int EventIO::Item::GetString(std::string& s)
{
   check_throw("GetString");

   s = "";
   size_t n = GetCount();
   s.resize(n);
   for (size_t i=0; i<n; ++i)
      s[i] = GetUint8();
   return int(n);
}

/// @short Get operation for C strings of any length.
///
/// Note: This is not the same as the old C method get_string().
/// To emulate that (e.g. for reading data written from old C code), use
/// @code
///    size_t n = item.GetUint16();
///    item.GetUint8((const char *) s, n);
/// @endcode
/// In the same way, the old C get_long_string() can be emulated by
/// using a GetUint32() in the above code.
/// The C language equivalent of the method here is called get_var_string().

int EventIO::Item::GetString(char *s, size_t nmax)
{
   check_throw("GetString");

   return get_var_string(s,nmax,iobuf);
}

/// Get operation for vectors of strings of any length.

void EventIO::Item::GetString(std::vector<std::string>& vec)
{
   check_throw("GetString");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      GetString(vec[i]);
}

/// Get operation for vectors of strings of any length.

void EventIO::Item::GetString(std::vector<std::string>& vec, size_t num)
{
   check_throw("GetString");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      GetString(vec[i]);
}

/// Get operation for valarrays of strings of any length.

void EventIO::Item::GetString(std::valarray<std::string>& vec)
{
   check_throw("GetString");

   size_t num = GetCount();
   if ( num != vec.size() )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      GetString(vec[i]);
}

/// Get operation for valarrays of strings of any length.

void EventIO::Item::GetString(std::valarray<std::string>& vec, size_t num)
{
   check_throw("GetString");

   if ( vec.size() < num )
      vec.resize(num);

   for (size_t i=0; i<num; ++i)
      GetString(vec[i]);
}

// ------------------------ Put operations ----------------------

/// Put operations for bytes (only unsigned flavour implemented).

void EventIO::Item::PutUint8(const std::vector<uint8_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutUint8(vec[i]);
   for ( ; i<num; ++i)
      PutUint8(0);
}

/// Put operations for bytes (only unsigned flavour implemented).

void EventIO::Item::PutUint8(const std::valarray<uint8_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutUint8(vec[i]);
   for ( ; i<num; ++i)
      PutUint8(0);
}

/// Put operations for bools.

void EventIO::Item::PutBool(const std::vector<bool>& vec, size_t num)
{
   size_t i, j, n=vec.size();
   
   uint8_t v = 0;

   for (i=0, j=0; i<n && i<num; ++i, ++j)
   {
      if ( (i&0x07) == 0 )
      {
         if ( i>0 )
            PutUint8(v);
         v = 0;
         j = 0;
      }
      if ( vec[i] )
         v |= (1<<j);
   }
   for ( ; i<num; ++i, ++j)
   {
      if ( (i&0x07) == 0 )
      {
         if ( i>0 )
            PutUint8(v);
         v = 0;
         j = 0;
      }
   }
   if ( i > 0 )
      PutUint8(v);
}

/// Put operations for bools.

void EventIO::Item::PutBool(const std::valarray<bool>& vec, size_t num)
{
   size_t i, j, n=vec.size();
   
   uint8_t v = 0;

   for (i=0, j=0; i<n && i<num; ++i, ++j)
   {
      if ( (i&0x07) == 0 )
      {
         if ( i>0 )
            PutUint8(v);
         v = 0;
         j = 0;
      }
      if ( vec[i] )
         v |= (1<<j);
   }
   for ( ; i<num; ++i, ++j)
   {
      if ( (i&0x07) == 0 )
      {
         if ( i>0 )
            PutUint8(v);
         v = 0;
         j = 0;
      }
   }
   if ( i > 0 )
      PutUint8(v);
}

/// Put operations for counts.

void EventIO::Item::PutCount(const std::vector<uint16_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutCount16(vec[i]);
   for ( ; i<num; ++i)
      PutCount16(0);
}

#ifdef SIXTY_FOUR_BITS

void EventIO::Item::PutCount(const std::vector<uint32_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutCount32(vec[i]);
   for ( ; i<num; ++i)
      PutCount32(0);
}

#endif

void EventIO::Item::PutCount(const std::vector<size_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutCount(vec[i]);
   for ( ; i<num; ++i)
      PutCount(0);
}

/// Put operations for counts.

#if defined(WITH_UINTMAX_T) && !defined(SIXTY_FOUR_BITS)
void EventIO::Item::PutCount(const std::vector<uintmax_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutCount(vec[i]);
   for ( ; i<num; ++i)
      PutCount(0);
}
#endif

/// Put operations for counts.

void EventIO::Item::PutCount(const std::valarray<uint16_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutCount16(vec[i]);
   for ( ; i<num; ++i)
      PutCount16(0);
}

#ifdef SIXTY_FOUR_BITS

void EventIO::Item::PutCount(const std::valarray<uint32_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutCount32(vec[i]);
   for ( ; i<num; ++i)
      PutCount32(0);
}

#endif

void EventIO::Item::PutCount(const std::valarray<size_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutCount(vec[i]);
   for ( ; i<num; ++i)
      PutCount(0);
}

/// Put operations for counts.

#if defined(WITH_UINTMAX_T) && !defined(SIXTY_FOUR_BITS)
void EventIO::Item::PutCount(const std::valarray<uintmax_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutCount(vec[i]);
   for ( ; i<num; ++i)
      PutCount(0);
}
#endif

/// Put operations for scounts.

void EventIO::Item::PutSCount(const std::vector<int16_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutSCount16(vec[i]);
   for ( ; i<num; ++i)
      PutSCount16(0);
}

void EventIO::Item::PutDiffSCount(const std::vector<int16_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   // The first element is written as-is, all others as differences.
   if ( n>0 && num>0 )
      PutSCount16(vec[0]);
   for (i=1; i<n && i<num; ++i)
      PutSCount16(vec[i]-vec[i-1]);
   // If more elements requested than available, fill with effectively zeroes.
   if ( i < num )
   {
      if ( i>0 )
      {
         PutSCount16(-vec[i-1]);
         i++;
      }
      for ( ; i<num; ++i)
         PutSCount16(0);
   }
}

#ifdef SIXTY_FOUR_BITS

void EventIO::Item::PutSCount(const std::vector<int32_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutSCount32(vec[i]);
   for ( ; i<num; ++i)
      PutSCount32(0);
}

void EventIO::Item::PutDiffSCount(const std::vector<int32_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   // The first element is written as-is, all others as differences.
   if ( n>0 && num>0 )
      PutSCount32(vec[0]);
   for (i=1; i<n && i<num; ++i)
      PutSCount32(vec[i]-vec[i-1]);
   // If more elements requested than available, fill with effectively zeroes.
   if ( i < num )
   {
      if ( i>0 )
      {
         PutSCount32(-vec[i-1]);
         i++;
      }
      for ( ; i<num; ++i)
         PutSCount32(0);
   }
}

#endif

void EventIO::Item::PutSCount(const std::vector<ssize_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutSCount(vec[i]);
   for ( ; i<num; ++i)
      PutSCount(0);
}

void EventIO::Item::PutDiffSCount(const std::vector<ssize_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   // The first element is written as-is, all others as differences.
   if ( n>0 && num>0 )
      PutSCount(vec[0]);
   for (i=1; i<n && i<num; ++i)
      PutSCount(vec[i]-vec[i-1]);
   // If more elements requested than available, fill with effectively zeroes.
   if ( i < num )
   {
      if ( i>0 )
      {
         PutSCount(-vec[i-1]);
         i++;
      }
      for ( ; i<num; ++i)
         PutSCount(0);
   }
}

/// Put operations for scounts.

#if defined(WITH_INTMAX_T) && !defined(SIXTY_FOUR_BITS)
void EventIO::Item::PutSCount(const std::vector<intmax_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutSCount(vec[i]);
   for ( ; i<num; ++i)
      PutSCount(0);
}

void EventIO::Item::PutDiffSCount(const std::vector<intmax_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   // The first element is written as-is, all others as differences.
   if ( n>0 && num>0 )
      PutSCount(vec[0]);
   for (i=1; i<n && i<num; ++i)
      PutSCount(vec[i]-vec[i-1]);
   // If more elements requested than available, fill with effectively zeroes.
   if ( i < num )
   {
      if ( i>0 )
      {
         PutSCount(-vec[i-1]);
         i++;
      }
      for ( ; i<num; ++i)
         PutSCount(0);
   }
}
#endif

/// Put operations for scounts.

void EventIO::Item::PutSCount(const std::valarray<int16_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutSCount16(vec[i]);
   for ( ; i<num; ++i)
      PutSCount16(0);
}

void EventIO::Item::PutDiffSCount(const std::valarray<int16_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   // The first element is written as-is, all others as differences.
   if ( n>0 && num>0 )
      PutSCount(vec[0]);
   for (i=1; i<n && i<num; ++i)
      PutSCount16(vec[i]-vec[i-1]);
   // If more elements requested than available, fill with effectively zeroes.
   if ( i < num )
   {
      if ( i>0 )
      {
         PutSCount16(-vec[i-1]);
         i++;
      }
      for ( ; i<num; ++i)
         PutSCount16(0);
   }
}

#ifdef SIXTY_FOUR_BITS

void EventIO::Item::PutSCount(const std::valarray<int32_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutSCount32(vec[i]);
   for ( ; i<num; ++i)
      PutSCount32(0);
}

void EventIO::Item::PutDiffSCount(const std::valarray<int32_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   // The first element is written as-is, all others as differences.
   if ( n>0 && num>0 )
      PutSCount(vec[0]);
   for (i=1; i<n && i<num; ++i)
      PutSCount32(vec[i]-vec[i-1]);
   // If more elements requested than available, fill with effectively zeroes.
   if ( i < num )
   {
      if ( i>0 )
      {
         PutSCount32(-vec[i-1]);
         i++;
      }
      for ( ; i<num; ++i)
         PutSCount32(0);
   }
}

#endif

void EventIO::Item::PutSCount(const std::valarray<ssize_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutSCount(vec[i]);
   for ( ; i<num; ++i)
      PutSCount(0);
}

void EventIO::Item::PutDiffSCount(const std::valarray<ssize_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   // The first element is written as-is, all others as differences.
   if ( n>0 && num>0 )
      PutSCount(vec[0]);
   for (i=1; i<n && i<num; ++i)
      PutSCount(vec[i]-vec[i-1]);
   // If more elements requested than available, fill with effectively zeroes.
   if ( i < num )
   {
      if ( i>0 )
      {
         PutSCount(-vec[i-1]);
         i++;
      }
      for ( ; i<num; ++i)
         PutSCount(0);
   }
}

/// Put operations for counts.

#if defined(WITH_INTMAX_T) && !defined(SIXTY_FOUR_BITS)
void EventIO::Item::PutSCount(const std::valarray<intmax_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutSCount(vec[i]);
   for ( ; i<num; ++i)
      PutSCount(0);
}

void EventIO::Item::PutDiffSCount(const std::valarray<intmax_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   // The first element is written as-is, all others as differences.
   if ( n>0 && num>0 )
      PutSCount(vec[0]);
   for (i=1; i<n && i<num; ++i)
      PutSCount(vec[i]-vec[i-1]);
   // If more elements requested than available, fill with effectively zeroes.
   if ( i < num )
   {
      if ( i>0 )
      {
         PutSCount(-vec[i-1]);
         i++;
      }
      for ( ; i<num; ++i)
         PutSCount(0);
   }
}
#endif

/// Put operations for 'Int16' item type (signed 16-bit integers, 'SHORT')

void EventIO::Item::PutInt16(const std::vector<int16_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutInt16(vec[i]);
   for ( ; i<num; ++i)
      PutInt16(0);
}

/// Put operations for 'Int16' item type (signed 16-bit integers, 'SHORT')

void EventIO::Item::PutInt16(const std::vector<int>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutInt16(vec[i]);
   for ( ; i<num; ++i)
      PutInt16(0);
}

/// Put operations for 'Int16' item type (signed 16-bit integers, 'SHORT')

void EventIO::Item::PutInt16(const std::valarray<int16_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutInt16(vec[i]);
   for ( ; i<num; ++i)
      PutInt16(0);
}

/// Put operations for 'Int16' item type (signed 16-bit integers, 'SHORT')

void EventIO::Item::PutInt16(const std::valarray<int>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutInt16(vec[i]);
   for ( ; i<num; ++i)
      PutInt16(0);
}

/// Put operations for 'Uint16' item type (unsigned 16-bit integers)

void EventIO::Item::PutUint16(const std::vector<uint16_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutUint16(vec[i]);
   for ( ; i<num; ++i)
      PutUint16(uint16_t(0));
}

/// Put operations for 'Uint16' item type (unsigned 16-bit integers)

void EventIO::Item::PutUint16(const std::valarray<uint16_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutUint16(vec[i]);
   for ( ; i<num; ++i)
      PutUint16(uint16_t(0));
}

/// Put operations for 'Int32' item type (signed 32-bit integers, 'LONG')

void EventIO::Item::PutInt32(const std::vector<int32_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutInt32(vec[i]);
   for ( ; i<num; ++i)
      PutInt32(0);
}

/// Put operations for 'Int32' item type (signed 32-bit integers, 'LONG')

void EventIO::Item::PutInt32(const std::vector<long>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutInt32(vec[i]);
   for ( ; i<num; ++i)
      PutInt32(0);
}

/// Put operations for 'Int32' item type (signed 32-bit integers, 'LONG')

void EventIO::Item::PutInt32(const std::valarray<int32_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutInt32(vec[i]);
   for ( ; i<num; ++i)
      PutInt32(0);
}

/// Put operations for 'Int32' item type (signed 32-bit integers, 'LONG')

void EventIO::Item::PutInt32(const std::valarray<long>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutInt32(vec[i]);
   for ( ; i<num; ++i)
      PutInt32(0);
}

/// Put operations for 'Uint32' item type (unsigned 32-bit integers)

void EventIO::Item::PutUint32(const std::vector<uint32_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutUint32(vec[i]);
   for ( ; i<num; ++i)
      PutUint32(0U);
}

/// Put operations for 'Uint32' item type (unsigned 32-bit integers)

void EventIO::Item::PutUint32(const std::valarray<uint32_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutUint32(vec[i]);
   for ( ; i<num; ++i)
      PutUint32(0U);
}

#ifdef HAVE_64BIT_INT

/// Put operations for 'Int64' item type (signed 64-bit integers)

void EventIO::Item::PutInt64(const std::vector<int64_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutInt64(vec[i]);
   for ( ; i<num; ++i)
      PutInt64(0);
}

/// Put operations for 'Int64' item type (signed 64-bit integers)

void EventIO::Item::PutInt64(const std::valarray<int64_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutInt64(vec[i]);
   for ( ; i<num; ++i)
      PutInt64(0);
}

/// Put operations for 'Uint64' item type (unsigned 64-bit integers)

void EventIO::Item::PutUint64(const std::vector<uint64_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutUint64(vec[i]);
   for ( ; i<num; ++i)
      PutUint64(0);
}

/// Put operations for 'Uint64' item type (unsigned 64-bit integers)

void EventIO::Item::PutUint64(const std::valarray<uint64_t>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutUint64(vec[i]);
   for ( ; i<num; ++i)
      PutUint64(0);
}
#endif

/// Put operations for 'Real' item type (32-bit IEEE float)

void EventIO::Item::PutReal(const std::vector<float>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutReal(vec[i]);
   for ( ; i<num; ++i)
      PutReal(0.);
}

/// Put operations for 'Real' item type (32-bit IEEE float)

void EventIO::Item::PutReal(const std::valarray<float>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutReal(vec[i]);
   for ( ; i<num; ++i)
      PutReal(0.);
}

/// Put operations for 'Real' item type (32-bit IEEE float)

void EventIO::Item::PutReal(const std::vector<double>& vec, size_t num)
{
   size_t i, n=vec.size();
   
   for (i=0; i<n && i<num; ++i)
      PutReal(vec[i]);
   for ( ; i<num; ++i)
      PutReal(0.);
}

/// Put operations for 'Real' item type (32-bit IEEE float)

void EventIO::Item::PutReal(const std::valarray<double>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutReal(vec[i]);
   for ( ; i<num; ++i)
      PutReal(0.);
}

/// Put operations for 'Double' item type (64-bit IEEE double)

void EventIO::Item::PutDouble(const std::vector<double>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutDouble(vec[i]);
   for ( ; i<num; ++i)
      PutDouble(0.);
}

/// Put operations for 'Double' item type (64-bit IEEE double)

void EventIO::Item::PutDouble(const std::valarray<double>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutDouble(vec[i]);
   for ( ; i<num; ++i)
      PutDouble(0.);
}


/// Put operations for 'Sfloat' item type (16-bit OpenGL float)

void EventIO::Item::PutSfloat(const std::vector<float>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutSfloat(vec[i]);
   for ( ; i<num; ++i)
      PutSfloat(0.);
}

void EventIO::Item::PutSfloat(const std::vector<double>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutSfloat(vec[i]);
   for ( ; i<num; ++i)
      PutSfloat(0.);
}

/// Put operations for 'Sfloat' item type (16-bit OpenGL float)

void EventIO::Item::PutSfloat(const std::valarray<float>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutSfloat(vec[i]);
   for ( ; i<num; ++i)
      PutSfloat(0.);
}

void EventIO::Item::PutSfloat(const std::valarray<double>& vec, size_t num)
{
   size_t i, n=vec.size();

   for (i=0; i<n && i<num; ++i)
      PutSfloat(vec[i]);
   for ( ; i<num; ++i)
      PutSfloat(0.);
}

/// Put operations for strings of any length

void EventIO::Item::PutString(const std::string& s)
{
   size_t l = s.size();
   PutCount(l);
   for (size_t i=0; i<l; ++i)
      PutUint8(s[i]);
}

void EventIO::Item::PutString(const char *s)
{
   size_t l = strlen(s);
   PutCount(l);
   for (size_t i=0; i<l; ++i)
      PutUint8(s[i]);
}

/// Put operations for strings with old-style 16 bit length prefix.

void EventIO::Item::PutString16(const std::string& s)
{
   size_t l = s.size();
   if ( l>32767UL )
      l = 32767UL; // Need to truncate
   uint16_t l16 = (uint16_t) l;

   PutUint16(l16);
   for (size_t i=0; i<l; ++i)
      PutUint8(s[i]);
}

void EventIO::Item::PutString16(const char *s)
{
   size_t l = strlen(s);
   if ( l>32767UL )
      l = 32767UL; // Need to truncate
   uint16_t l16 = (uint16_t) l;

   PutUint16(l16);
   for (size_t i=0; i<l; ++i)
      PutUint8(s[i]);
}

/// Put operations for strings with old-style 32 bit length prefix.

void EventIO::Item::PutString32(const std::string& s)
{
   size_t l = s.size();
   if ( l>2147483647UL )
      l = 2147483647UL; // Need to truncate
   uint32_t l32 = (uint32_t) l;

   PutUint32(l32);
   for (size_t i=0; i<l; ++i)
      PutUint8(s[i]);
}

void EventIO::Item::PutString32(const char *s)
{
   size_t l = strlen(s);
   if ( l>2147483647UL )
      l = 2147483647UL; // Need to truncate
   uint32_t l32 = (uint32_t) l;

   PutUint32(l32);
   for (size_t i=0; i<l; ++i)
      PutUint8(s[i]);
}

/// Put operations for vectors of strings of any length.

void EventIO::Item::PutString(const std::vector<std::string>& vec, size_t num)
{
   size_t i, n=vec.size();
   
   for (i=0; i<n && i<num; ++i)
      PutString(vec[i]);
   for ( ; i<num; ++i)
      PutString("");
}

/// Put operations for valarrays of strings of any length.

void EventIO::Item::PutString(const std::valarray<std::string>& vec, size_t num)
{
   size_t i, n=vec.size();
   
   for (i=0; i<n && i<num; ++i)
      PutString(vec[i]);
   for ( ; i<num; ++i)
      PutString("");
}

} // End of namespace eventio

