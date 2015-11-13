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

/** @file io_trgmask.h
 *  @short EventIO plus helper functions for trigger type bit patterns extracted
 *         from sim_telarray log files (only relevant for simulations with
 *         multiple trigger types using sim_telarray versions before mid-2013). 
 */

#ifndef IO_TRGMASK_LOADED

#define IO_TRGMASK_LOADED 1

#ifdef __cplusplus
extern "C" {
#endif

/** Extra (or external - not in normal data file) trigger mask data block type. */
#define IO_TYPE_HESS_XTRGMASK 2090

struct trgmask_entry
{
   long event;                 ///< The event number.
   int tel_id;                 ///< The telescope ID number.
   int trg_mask;               ///< The trigger mask bit pattern which got messed up in data files.
   struct trgmask_entry *next; ///< Can be used in arrays but also in linked lists
};

struct trgmask_set
{
   long run;
   size_t num_entries;
   struct trgmask_entry *mask;
};

// #define TRGMASK_PRIME 10067
#define TRGMASK_PRIME 15269

#define TRGMASK_HASH(ev,ti) (((ti)*10000+(ev))%TRGMASK_PRIME)

struct trgmask_hash_set
{
   long run;
   struct trgmask_entry *h_e[TRGMASK_PRIME]; ///< Start of linked list for each possible hash value
};

int trgmask_scan_log (struct trgmask_set *tms, const char *fname);
int write_trgmask (IO_BUFFER *iobuf, struct trgmask_set *tms);
int print_trgmask (IO_BUFFER *iobuf);
int read_trgmask (IO_BUFFER *iobuf, struct trgmask_set *tms);
int trgmask_fill_hashed(struct trgmask_set *tms, struct trgmask_hash_set *ths);
struct trgmask_entry *find_trgmask(struct trgmask_hash_set *ths, long event, int tel_id);
void print_hashed_trgmasks(struct trgmask_hash_set *ths);

#ifdef __cplusplus
}
#endif

#endif

