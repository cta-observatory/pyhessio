/* ============================================================================

   Copyright (C) 2014  Konrad Bernloehr

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

/** @file eventio_registry.c
 *  @short Register and enquire about well-known I/O block types.
 *
 *
 *  @author  Konrad Bernloehr 
 *  @date    2014
 *  @date    CVS $Date: 2018/09/19 12:11:36 $
 *  @version CVS $Revision: 1.5 $
 */

#include "initial.h"
#include "eventio_registry.h"
#include "fileopen.h"

#if 0
/**
 *  An entry in the registry of well-known data block types.
 */

struct ev_reg_entry
{
   unsigned long type; /**< The data block type number */
   char *name;         /**< The data block name (short) */
   char *description;  /**< Optional longer description of the data block */
};
#endif

/** Allocate a new entry for the registry */

struct ev_reg_entry *new_reg_entry(unsigned long t, const char *n, const char *d);

struct ev_reg_entry *new_reg_entry(unsigned long t, const char *n, const char *d)
{
   struct ev_reg_entry *e;
   if ( t<=0 || n==NULL )
      return NULL;
   if ( (e = (struct ev_reg_entry *) malloc(sizeof(struct ev_reg_entry))) == NULL )
      return NULL;
   e->type = t;
   e->name = strdup(n);
   if ( d != NULL )
      e->description = strdup(d);
   else
      e->description = NULL;
   return e;
}

/** Use a double-linked list for the registry */

struct ev_reg_chain
{
   struct ev_reg_entry *entry;  /**< The current entry */
   struct ev_reg_chain *prev;
   struct ev_reg_chain *next;
};

static struct ev_reg_chain *ev_reg_start = NULL;

/** Read the type names and descriptions into the registry. */
/** Note: this will only be done once. */

int read_eventio_registry(const char *fname)
{
   char line[2048];
   FILE *f;
   struct ev_reg_chain *ev_reg_pos = NULL;
   static int init_done = 0;
   if ( init_done ) /* Just try once */
      return 1;
   init_done = 1;
   if ( ev_reg_start != NULL )
      return -1;
   f = fileopen(fname,"r");
   if ( f == NULL )
   {
      perror(fname);
      return -1;
   }
   while ( fgets(line,sizeof(line)-1,f) != NULL )
   {
      size_t l;
      struct ev_reg_entry *e = NULL;
      char *c = strchr(line,'#');
      unsigned long type = 0;
      char *tp = NULL, *np = NULL, *dp = NULL;
      struct ev_reg_chain *ev_reg_next = NULL;
      if ( c != NULL )
         *c = '\0';
      l = strlen(line);
      if ( l>0 )
      for (l=l-1;;l--)
      {
         if ( line[l] == '\n' || line[l] == ' ' || 
              line[l] == '\t' || line[l] == '#' )
            line[l] = '\0';
         else
            break;
         if ( l == 0 )
            break;
      }
      if ( line[0] == '\0' )
         continue;
      for ( tp=line; *tp == ' ' || *tp == '\t'; tp++ )
      { /* continue; */
      }
      if ( *tp == ':' || *tp == '\0' )
         continue; /* No type number */
      if ( sscanf(tp,"%lu",&type) != 1 )
         continue; /* No type number */
      for ( np=tp; *np != ':' && *np != '\0'; np++ )
      { /* continue; */
      }
      if ( *np != ':' )
         continue; /* No name */
      np++;
      for ( dp=np; *dp != ':' && *dp != '\0'; dp++ )
      { /* continue; */
      }
      if ( *dp == ':' )
      {
         *dp = '\0';
         dp++;
      }
      else
         dp = NULL;
      e = new_reg_entry(type,np,dp);
      if ( e == NULL )
         continue;
      if ( (ev_reg_next = (struct ev_reg_chain *) calloc(sizeof(struct ev_reg_chain),1)) == NULL )
         break;
      if ( ev_reg_start == NULL )
      {
         ev_reg_start = ev_reg_pos = ev_reg_next;
      }
      else
      {
         ev_reg_pos->next = ev_reg_next;
         ev_reg_next->prev = ev_reg_pos;
         ev_reg_pos = ev_reg_next;
      }
      ev_reg_pos->entry = e;
   }
   fileclose(f);
   return 0;
}

/** By default the registry contents will be searched in a few places */

static void read_default_registry(void);

static void read_default_registry()
{
   if ( getenv("EVENTIO_REGISTRY") != NULL )
      (void) read_eventio_registry(getenv("EVENTIO_REGISTRY"));
   else if ( access("EventioRegisteredNames.dat", R_OK) == 0 ) 
      (void) read_eventio_registry("EventioRegisteredNames.dat");
   else 
   {
      char home_reg[1024];
      *home_reg = '\0';
      if ( getenv("HOME") != NULL )
      {
         strncpy(home_reg,getenv("HOME"),sizeof(home_reg)-3-strlen(".EventioRegisteredNames"));
         strcat(home_reg,"/");
      }
      strcat(home_reg,".EventioRegisteredNames");
      if ( access(home_reg, R_OK) == 0 )
         (void) read_eventio_registry(home_reg);
      else
      {
         fprintf(stderr,"No eventIO types registry\n");
         ev_reg_start = (struct ev_reg_chain *) calloc(sizeof(struct ev_reg_chain),1);
      }
   }
}

/** Find an entry for a given type number in the registry */
/** This is the standard implementation being used by default where available. */

/* struct ev_reg_entry *find_ev_reg_std(unsigned long t); */

struct ev_reg_entry *find_ev_reg_std(unsigned long t)
{
   struct ev_reg_chain *ev_reg_pos;
   if ( ev_reg_start == NULL )
   {
      read_default_registry();
      if ( ev_reg_start == NULL )
         return NULL;
   }
   /* FIXME: For long lists, this linear search could be accelerated */
   for ( ev_reg_pos = ev_reg_start; ev_reg_pos != NULL; ev_reg_pos = ev_reg_pos->next )
   {
      if ( ev_reg_pos->entry == NULL )
         continue;
      if ( ev_reg_pos->entry->type == t )
         return ev_reg_pos->entry;
   }
   return NULL;
}

/** Set the default registry search function */
/** At least with GCC we can do this without explicitly calling it. */

#ifdef __GNUC__
void __attribute__ ((constructor)) set_ev_reg_std()
#else
void set_ev_reg_std()
#endif
{
   set_eventio_registry_hook(&find_ev_reg_std);
}

#if 0
int main()
{
   struct ev_reg_chain *ev_reg_pos;
   read_default_registry();
   for ( ev_reg_pos = ev_reg_start; ev_reg_pos != NULL; ev_reg_pos = ev_reg_pos->next )
   {
      if ( ev_reg_pos->entry != NULL )
      {
         printf("%lu : '%s' : '%s'\n", ev_reg_pos->entry->type,
            ev_reg_pos->entry->name, ev_reg_pos->entry->description);
      }
   }
   printf("Type 1000 has name '%s'\n", eventio_registered_typename(1000));
   printf("Type 1001 has name '%s'\n", eventio_registered_typename(1001));
   return 0;
}
#endif
