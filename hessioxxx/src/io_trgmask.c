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

/** @file io_trgmask.c
 *  @short EventIO plus helper functions for trigger type bit patterns extracted
 *         from sim_telarray log files (only relevant for simulations with
 *         multiple trigger types using sim_telarray versions before mid-2013). 
 */

#include "initial.h"      /* This file includes others as required. */
#include "io_basic.h"     /* This file includes others as required. */
#include "fileopen.h"
#include "io_trgmask.h"

#define TMS_ALLOCS 100

/* ---------------------- trgmask_scan_log -------------------------------- */
/**
 *  @short Scan a sim_telarray log file for lines related to trigger type mask bit patterns.
 *
 *  @param tms The trigger mask structure into which results should be filled in.
 *  @param fname The name of the log file to be opened.
 *  @return 0 (OK), -1 (invalid parameters or file not found), -2 (allocation error, partially filled)
 */

int trgmask_scan_log (struct trgmask_set *tms, const char *fname)
{
   size_t last_alloc = 0;
   char line [10000];
   FILE *f = NULL;
   long event = 0;

   if ( tms == NULL || fname == NULL )
      return -1;
   tms->run = 0;
   tms->num_entries = 0;
   if ( tms->mask != NULL )
      free(tms->mask);

   tms->mask = (struct trgmask_entry *) calloc(TMS_ALLOCS,sizeof(struct trgmask_entry));
   if ( tms->mask == NULL )
   {
      fprintf(stderr,"Trigger mask allocation failed\n");
      return -1;
   }
   last_alloc = TMS_ALLOCS;

   if ( strcmp(fname,"-") == 0 )
      f = stdin;
   else
      f = fileopen(fname,"r");
   if ( f == NULL )
   {
      fprintf(stderr,"Cannot open log file %s as source for trigger mask patterns.\n",fname);
      return -1;
   }

   while ( fgets(line,sizeof(line)-1,f) != NULL )
   {
      if ( strncmp(line,"Run ",4) == 0 && strstr(line,"(started") != NULL )
      {
         long run = atol(line+4);
         if ( tms->run != 0 )
         {
            fprintf(stderr,"Warning: run number changed in log file from %ld to %ld.\n", tms->run, run);
         }
         else
         {
            printf("Processing log file for run %ld\n", run);
         }
         tms->run = run;
      }
      else if ( strncmp(line,"Event ",6) == 0 && strstr(line," has triggered") != NULL )
      {
         event = atol(line+6);
      }
      else if ( strncmp(line,"Telescope ",10) == 0 && strstr(line," triggered ") != NULL &&
                strstr(line,"mask") != NULL )
      {
         int tel_id = atoi(line+10);
         char *s = strstr(line,"mask");
         int mask = (s!=NULL ? atoi(s+5) : 0);
         if ( tms->num_entries >= last_alloc )
         {
            struct trgmask_entry *old_mask = tms->mask;
            if ( (tms->mask = (struct trgmask_entry *) realloc(tms->mask,(last_alloc+TMS_ALLOCS)*sizeof(struct trgmask_entry)))
                  == NULL )
            {
               fprintf(stderr,"Reallocation of trigger mask entries failed after %ld entries.\n",
                  (long)last_alloc);
               tms->mask = old_mask;
               return -2;
            }
            last_alloc += TMS_ALLOCS;
         }

         tms->mask[tms->num_entries].event = event;
         tms->mask[tms->num_entries].tel_id = tel_id;
         tms->mask[tms->num_entries].trg_mask = mask;
         tms->mask[tms->num_entries].next = (struct trgmask_entry *) NULL; // Not needed here
         // We could also set up the linked list, although we do not really need it.
         // if ( tms->num_entries > 0 )
         //   tms->mask[tms->num_entries-1] = &tms->mask[tms->num_entries].trg_mask;
         tms->num_entries++;

//         printf("Run %ld, event %ld, telescope %d has mask %d\n",
//            tms->run, event, tel_id, mask);

//if ( tms->num_entries >= 20 ) // For debugging
//   return 0;
      }
   }
   if ( strcmp(fname,"-") != 0 )
      fileclose(f);
   return 0;
}

/* -------------------------- write_trgmask ---------------------------- */
/**
 *  Write the accumulated trigger mask bit patterns as an I/O block.
 */

int write_trgmask (IO_BUFFER *iobuf, struct trgmask_set *tms)
{
   IO_ITEM_HEADER item_header;
   size_t i;
   
   if ( iobuf == (IO_BUFFER *) NULL || tms == NULL )
      return -1;
   if ( tms->mask == (struct trgmask_entry *) NULL )
      return -1;
   
   item_header.type = IO_TYPE_HESS_XTRGMASK;  /* Data type */
   item_header.version = 1;                    /* Version 2 now */
   item_header.ident = tms->run;
   put_item_begin(iobuf,&item_header);
   
   put_count(tms->num_entries,iobuf);
   for (i=0; i<tms->num_entries; i++)
   {
      put_count(tms->mask[i].event,iobuf);
      put_scount32(tms->mask[i].tel_id,iobuf);
      put_scount32(tms->mask[i].trg_mask,iobuf);
   }

   return put_item_end(iobuf,&item_header);
}

/* -------------------------- print_trgmask ---------------------------- */
/**
 *  Print the trigger mask bit patterns contained in an I/O block.
 */

int print_trgmask (IO_BUFFER *iobuf)
{
   IO_ITEM_HEADER item_header;
   size_t num_entries = 0, i;
   long run, event;
   int tel_id, trg_mask;
   int rc;
   int maxprt = 20;
   if ( getenv("MAX_PRINT_ARRAY") != NULL )
      maxprt = atoi(getenv("MAX_PRINT_ARRAY"));
   
   if ( iobuf == (IO_BUFFER *) NULL )
      return -1;

   item_header.type = IO_TYPE_HESS_XTRGMASK;  /* Data type */
   if ( (rc = get_item_begin(iobuf,&item_header)) < 0 )
      return rc;
   if ( item_header.version > 1 )
   {
      fprintf(stderr,"Unsupported trigger mask version: %d.\n",
         item_header.version);
      get_item_end(iobuf,&item_header);
      return -1;
   }
   
   run = item_header.ident;
   num_entries = (size_t) get_count(iobuf);
   printf("Trigger mask block for run %ld with %zu entries:\n", run, num_entries);
   for ( i=0; i<num_entries; i++ )
   {
      event = get_count(iobuf);
      tel_id = get_scount32(iobuf);
      trg_mask = get_scount32(iobuf);
      if ( (int) i < maxprt )
         printf("   Event %ld, telescope %d has mask %d\n", event, tel_id, trg_mask);
      else if ( (int) i == maxprt )
         printf("   ...\n");
   }

   return get_item_end(iobuf,&item_header);
}

/* -------------------------- print_trgmask ---------------------------- */
/**
 *  Read the trigger mask bit patterns contained in an I/O block.
 */

int read_trgmask (IO_BUFFER *iobuf, struct trgmask_set *tms)
{
   IO_ITEM_HEADER item_header;
   size_t i;
   int rc;

   if ( iobuf == (IO_BUFFER *) NULL || tms == (struct trgmask_set *) NULL )
      return -1;
   if ( tms->mask != (struct trgmask_entry *) NULL )
      free(tms->mask);

   item_header.type = IO_TYPE_HESS_XTRGMASK;  /* Data type */
   if ( (rc = get_item_begin(iobuf,&item_header)) < 0 )
      return rc;
   if ( item_header.version > 1 )
   {
      fprintf(stderr,"Unsupported trigger mask version: %d.\n",
         item_header.version);
      get_item_end(iobuf,&item_header);
      return -1;
   }

   tms->run = item_header.ident;
   tms->num_entries = (size_t) get_count(iobuf);
   tms->mask = (struct trgmask_entry *) calloc(tms->num_entries, sizeof(struct trgmask_entry));
   if ( tms->mask == (struct trgmask_entry *) NULL )
   {
      fprintf(stderr,"Trigger mask allocation error for %zu entries\n", tms->num_entries);
      tms->num_entries = 0;
      return -1;
   }
   for ( i=0; i<tms->num_entries; i++ )
   {
      tms->mask[i].event = get_count(iobuf);
      tms->mask[i].tel_id = get_scount32(iobuf);
      tms->mask[i].trg_mask = get_scount32(iobuf);
   }

   return get_item_end(iobuf,&item_header);
}

// #define TRGMASK_HASH(ev,ti) (((ti)*10000+(ev))%TRGMASK_PRIME)

/* ------------------------- trgmask_fill_hashed --------------------------- */
/**
 *  Fill an array of linked lists of trgmask entries, suitable for hashing.
 *  Hash collisions are handled by linear search through the linked list at each hash entry.
 */

int trgmask_fill_hashed (struct trgmask_set *tms, struct trgmask_hash_set *ths)
{
   size_t i, hi;
   if ( tms == (struct trgmask_set *) NULL || ths == (struct trgmask_hash_set *) NULL )
      return -1;
   if ( tms->mask == (struct trgmask_entry *) NULL )
      return -1;

   ths->run = tms->run;
   for ( hi=0; hi<TRGMASK_PRIME; hi++ )
      ths->h_e[hi] = NULL;

   for ( i=0; i<tms->num_entries; i++ )
   {
      hi = TRGMASK_HASH(tms->mask[i].event,tms->mask[i].tel_id);
      if ( hi < 0 || hi >= TRGMASK_PRIME )
      {
         fprintf(stderr,
            "Not filling invalid trigger bit pattern entry for event %ld telescope ID %d.\n",
            tms->mask[i].event, tms->mask[i].tel_id);
            continue;
      }
      struct trgmask_entry *h = ths->h_e[hi];
      if ( h == (struct trgmask_entry *) NULL )
      {
         ths->h_e[hi] = &tms->mask[i];
      }
      else
      {
         while ( h->next != NULL )
            h = h->next;
         h->next = &tms->mask[i];
      }
      tms->mask[i].next = (struct trgmask_entry *) NULL;
   }

   return 0;
}

/* ------------------------- find_trgmask --------------------------- */
/**
 *  Find the trgmask entry for a given event and telescope in the hashed list.
 *  Hash collisions are handled by linear search through the linked list at each hash entry.
 *
 *  @param ths The trgmask hash set.
 *  @param event The event number in the search.
 *  @param tel_id The telescope ID in the search.
 *  @return A pointer to the trgmask entry searched for, or NULL for not found.
 */

struct trgmask_entry *find_trgmask(struct trgmask_hash_set *ths, long event, int tel_id)
{
   long hi = TRGMASK_HASH(event,tel_id);
   struct trgmask_entry *h;
   if ( hi < 0 || hi >= TRGMASK_PRIME )
   {
      fprintf(stderr,
         "Request for invalid trigger pattern for event %ld, telescope ID %d.\n",
         event, tel_id);
      return NULL;
   }
   h = ths->h_e[hi];
   if ( h == (struct trgmask_entry *) NULL )
      return NULL;
   while ( h != NULL )
   {
      if ( event == h->event && tel_id == h->tel_id )
         return h;
      h = h->next;
   }
   /* No matching entry found */
   return NULL;
}

/* -------------------- print_hashed_trgmasks ------------------------- */
/**
 *  Print the collected trgmask entries in the order as hashed.
 *  Also show the maximum number of colliding entries under one hash value.
 */

void print_hashed_trgmasks(struct trgmask_hash_set *ths)
{
   long hi;
   size_t np=0, nx=0;
   printf("Trigger patterns as hashed for run %ld:\n", ths->run);
   for (hi=0; hi<TRGMASK_PRIME; hi++)
   {
      struct trgmask_entry *h = ths->h_e[hi];
      np = 0;
      while ( h != (struct trgmask_entry *) NULL )
      {
         printf("   Event %ld, telescope %d has mask %d\n", h->event, h->tel_id, h->trg_mask);
         h = h->next;
         np++;
      }
      if ( np > nx )
         nx = np;
   }
   printf("Longest list of hashed entries: %zu\n", nx);
}

#ifdef TEST_TRGMASK
int main(int argc, char **argv)
{
   struct trgmask_set *tms = calloc(1,sizeof(struct trgmask_set));
   struct trgmask_hash_set *ths = calloc(1,sizeof(struct trgmask_hash_set));

   if ( argc > 1 )
      if ( trgmask_scan_log(tms,argv[1]) < 0 )
         exit(1);

   printf("Setting up hashed list ...\n");
   trgmask_fill_hashed(tms,ths);
   print_hashed_trgmasks(ths);
   
   /* Avoid cppcheck false positives: */
   free(tms);
   free(ths);

   return 0;
}
#endif

