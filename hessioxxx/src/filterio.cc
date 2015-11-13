/* ============================================================================

Copyright (C) 2014 Konrad Bernloehr (Konrad Bernl&ouml;hr)

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

/** @file filterio.cc
 *  @short A program for filtering eventio data blocks.
 *
@verbatim
Filter (or filter out) some EventIO blocks in given files.
Syntax: bin/filterio [ options ] [ block types ... ] filename [ ... ]
Options:
  -s  Show statistics of input and output.
  -v  Verbose output.
  -q  Quiet mode, not even showing total blocks count.
  -n  Not accepting any blocks of the listed types.
      The default is to accept only the listed types.
  -o fname Output goes to given file instead of standard output.
@endverbatim
 *
 *  @author Konrad Bernloehr
 *
 */

/** @defgroup filterio_cc The filterio program */
/** @{ */

#define __STDC_LIMIT_MACROS 1
#include <map>
using std::map;
#include <vector>
using std::vector;
#include <iostream>
#include <string>
#include "fileopen.h"
#include "hconfig.h"
#include "EventIO.hh"

using namespace eventio;
using std::string;
using std::cout;
using std::cerr;

struct iostats
{
   size_t count;
   uint64_t bytes;
   unsigned version_low;
   unsigned version_high;
   
   iostats() {count=0; bytes=0; version_low=0; version_high=0; }
   iostats(size_t c, uint64_t b, unsigned v) { count=c; bytes=b;
      version_low = version_high = v; }
   ~iostats() {}
};

void syntax(const char *prg)
{
   fprintf(stderr,"Filter (or filter out) some EventIO blocks in given files.\n");
   fprintf(stderr,"Syntax: %s [ options ] [ block types ... ] filename [ ... ]\n", prg); 
   fprintf(stderr,"Options:\n");
   fprintf(stderr,"  -s  Show statistics of input and output.\n");
   fprintf(stderr,"  -v  Verbose output.\n");
   fprintf(stderr,"  -q  Quiet mode, not even showing total blocks count.\n");
   fprintf(stderr,"  -n  Not accepting any blocks of the listed types.\n");
   fprintf(stderr,"      The default is to accept only the listed types.\n");
   fprintf(stderr,"  -o fname Output goes to given file instead of standard output.\n");
   exit(1);
}

int main (int argc, char **argv)
{
   EventIO iobuf(10000000,1000000000);
   iobuf.OpenOutput(stdout);
   bool verbose = false; ///< Verbose output (for stats)
   bool negate = false;  ///< Negate the list of block types.
   bool show_stats = false; ///< Show statistics of input and output.
   bool quiet = false; ///< Absolute no output other than data.
   char *prg = argv[0];
   size_t nf = 0;
   size_t ni = 0, no = 0;

   vector<unsigned long> type_list;
   size_t nl = 0;

   map<unsigned long,iostats> sm_in, sm_out;

#ifdef ALWAYS_WITH_REGISTRY
   /* Use default registry of known types with any compiler. */
   /* Requires linking against the library (not available with */
   /* bare-bone eventio, as for CORSIKA IACT/ATMO extension). */
   /* Not necessary when compiled with GCC. */
   set_ev_reg_std();
#endif

   for (int iarg=1; iarg<argc; iarg++)
   {
      if ( string(argv[iarg]) == "-v" || string(argv[iarg]) == "--verbose" )
      {
         verbose = true;
         continue;
      } 
      if ( string(argv[iarg]) == "-n" ||  string(argv[iarg]) == "--not" )
      {
         negate = true;
         continue;
      } 
      if ( string(argv[iarg]) == "-q" ||  string(argv[iarg]) == "--quiet" )
      {
         quiet = true;
         continue;
      } 
      if ( string(argv[iarg]) == "-s" ||  string(argv[iarg]) == "--statistics" )
      {
         show_stats = true;
         continue;
      } 
      if ( string(argv[iarg]) == "-o"  && iarg+1<argc )
      {
         iobuf.OpenOutput(argv[++iarg]);
         continue;
      }

      if ( string(argv[iarg]) == "--help" )
      {
         syntax(prg);
      }
      
      // Unsigned numeric arguments denote block types.
      if ( is_unsigned_number(argv[iarg]) )
      {
         unsigned long tp = 0;
         if ( sscanf(argv[iarg],"%lu",&tp) == 1 )
         {
            type_list.push_back(tp);
            nl = type_list.size();
         }
         continue;
      }
   
      // Anything else should be a file name.
#ifdef READ_BINARY
      FILE *f = fileopen(argv[iarg],READ_BINARY);
#else
      FILE *f = fileopen(argv[iarg],"r");
#endif
      if ( f == 0 )
      {
         perror(argv[iarg]);
         continue;
      }
      nf++;
      iobuf.OpenInput(f);
      for (;;)
      {
         if ( iobuf.Find() < 0 )
            break;
         unsigned long it = iobuf.ItemType();
         unsigned iv = iobuf.ItemVersion();
         size_t ln = iobuf.size();
         if ( iobuf.Read() < 0 )
            break;
         ni++;
         if ( show_stats )
         {
            map<unsigned long,iostats>::iterator ti = sm_in.find(it);
            if ( ti != sm_in.end() )
            {
               ti->second.count++;
               ti->second.bytes += ln;
               if ( iv < ti->second.version_low )
                   ti->second.version_low = iv;
               if ( iv > ti->second.version_high )
                   ti->second.version_high = iv;
            }
            else
            {
               sm_in[it] = iostats(1,ln,iv);
            }
         }
         bool in_list = false;
         for ( size_t i=0; i<nl; i++ )
         {
            if ( type_list[i] == it )
            {
               in_list = true;
               break;
            }
         }
         if ( negate )
            in_list = !in_list;
         if ( !in_list )
            continue;
         if ( iobuf.Write() != 0 )
            break;
         no++;
         if ( show_stats )
         {
            map<unsigned long,iostats>::iterator ti = sm_out.find(it);
            if ( ti != sm_out.end() )
            {
               ti->second.count++;
               ti->second.bytes += ln;
               if ( iv < ti->second.version_low )
                   ti->second.version_low = iv;
               if ( iv > ti->second.version_high )
                   ti->second.version_high = iv;
            }
            else
            {
               sm_out[it] = iostats(1,ln,iv);
            }
         }
      }
      iobuf.CloseInput();
   }

   if ( nf == 0 )
   {
      syntax(prg);
   }

   if ( show_stats )
   {
      cerr << "\n================== Input =========================\n\n";

      size_t total_c = 0;
      uint64_t total_b = 0;
      if ( sm_in.size() > 0 && ! verbose )
         cerr << "Type\tBlocks\tBytes\t    Version(s)\tName\n";

      for (map<unsigned long,iostats>::iterator ti = sm_in.begin(); 
            ti != sm_in.end(); ++ti)
      {
         const char *name = eventio_registered_typename(ti->first);
         if ( verbose )
         {
            const char *desc = eventio_registered_description(ti->first);
            cerr << "Type " << ti->first << ": " 
                 << ti->second.count << " blocks with " 
                 << ti->second.bytes << " bytes";
            if ( ti->second.version_low == ti->second.version_high )
               cerr << " (version " << ti->second.version_low << ")";
            else
               cerr << " (versions " << ti->second.version_low 
                    << " to " << ti->second.version_high << ")";
            if ( name != NULL && *name != '\0' )
               cerr << "\t[" << name << "] " << desc << "\n";
            else
               cerr << "\n";
         }
         else
         {
            cerr << ti->first << "\t" << ti->second.count << "\t" 
                 << ti->second.bytes << "\t"
                 << (ti->second.bytes < 10000000 ? "\t" : "")
                 << ti->second.version_low;
            if ( ti->second.version_low != ti->second.version_high )
               cerr << "-" << ti->second.version_high;
            if ( name != NULL && *name != '\0' )
               cerr << "\t[" << name << "]\n";
            else
               cerr << "\n";
         }
         total_c += ti->second.count;
         total_b += ti->second.bytes;
      }

      if ( verbose )
      {
         cerr << "Total:" << "\t" << total_c << " blocks containing " << total_b << " bytes";
         if ( total_b >= (1UL << 30) )
            cerr << " (" << total_b/(double)(1UL << 30) << " GiB)";
         cerr << "\n";
      }
      else
         cerr << "Total:" << "\t" << total_c << "\t" << total_b << "\n";

      cerr << "\n================== Output ========================\n\n";

      total_c = 0;
      total_b = 0;
      if ( sm_out.size() > 0 && ! verbose )
         cerr << "Type\tBlocks\tBytes\t    Version(s)\tName\n";

      for (map<unsigned long,iostats>::iterator ti = sm_out.begin(); 
            ti != sm_out.end(); ++ti)
      {
         const char *name = eventio_registered_typename(ti->first);
         if ( verbose )
         {
            const char *desc = eventio_registered_description(ti->first);
            cerr << "Type " << ti->first << ": " 
                 << ti->second.count << " blocks with " 
                 << ti->second.bytes << " bytes";
            if ( ti->second.version_low == ti->second.version_high )
               cerr << " (version " << ti->second.version_low << ")";
            else
               cerr << " (versions " << ti->second.version_low 
                    << " to " << ti->second.version_high << ")";
            if ( name != NULL && *name != '\0' )
               cerr << "\t[" << name << "] " << desc << "\n";
            else
               cerr << "\n";
         }
         else
         {
            cerr << ti->first << "\t" << ti->second.count << "\t" 
                 << ti->second.bytes << "\t"
                 << (ti->second.bytes < 10000000 ? "\t" : "")
                 << ti->second.version_low;
            if ( ti->second.version_low != ti->second.version_high )
               cerr << "-" << ti->second.version_high;
            if ( name != NULL && *name != '\0' )
               cerr << "\t[" << name << "]\n";
            else
               cerr << "\n";
         }
         total_c += ti->second.count;
         total_b += ti->second.bytes;
      }

      if ( verbose )
      {
         cerr << "Total:" << "\t" << total_c << " blocks containing " << total_b << " bytes";
         if ( total_b >= (1UL << 30) )
            cerr << " (" << total_b/(double)(1UL << 30) << " GiB)";
         cerr << "\n";
      }
      else
         cerr << "Total:" << "\t" << total_c << "\t" << total_b << "\n";

      cerr << "\n==================================================\n\n";
   }
   else if ( ! quiet )
   {
      cerr << ni << (ni==1?" block in, ":" blocks in, ") 
           << no << (no==1?" block out\n":" blocks out\n");
   }

   return 0;
}

/** @} */
