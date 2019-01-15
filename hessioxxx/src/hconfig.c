/* ============================================================================

   Copyright (C) 2001, 2006, 2008, 2010  Konrad Bernloehr

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

/* ================================================================ */
/**
   @file hconfig.c
   @brief Configuration control and procedure call interface.

   @author  Konrad Bernloehr
   @date    $Date: 2018/05/11 11:54:29 $
   @version $Revision: 1.21 $


     This is the module controlling all configuration except
     that a function has to be supplied that collects input line
     for line. Most functions in this file are for internal use only
     and are given a 'static' modifier. The only functions to be
     called by the user are
@verbatim   
	build_config()
	init_config()
	reconfig().
@endverbatim

     In order to set up the configuration, one or several calls to
     build_config() should be done, each with a list of 'configuration
     items' ('CONFIG_ITEM *items') terminated by a NULL_CONFIG_ITEM
     as an end marker. The list must be of 'static' or global/'extern'
     type and none of its entries must be modified by the user in any
     way, once they have been passed to build_config.

     Such a list might look like the following example:
@verbatim   
	static CONFIG_ITEM cfg_list[] =
	{
	   { "ANY_Numbers", "Int", 30, iarray },
	   { "ANY_Function", "function", -1, NULL, some_function },
	   { "REAL_Number", "R", 10, dblarray, NULL, "0-9: 99.9", 
	       "10", "100", CFG_REQUIRE_ALL_DATA | CFG_REJECT_MODIFICATION },
	   { "DYnAllocArray", "i", 100, NULL, NULL },
	   { NULL_CONFIG_ITEM }
	}
@endverbatim

     The components of each item are:
@verbatim   
	1) The name, consisting of letters, digits, and '_'.
	   In external data the items are referenced by their name
	   which may be abbreviated and is case-insensitive. However,
	   the name used for the definition is case-sensitive in the
	   current implementation. The first lowercase letter indicates
	   the minimum length of accepted abbrevations. In the example
	   above "ANY_Numbers" may be abbreviated as "any_n", "any_nu",
	   and so on, "DYnAllocArray" as "dy", "dyn", and so on.
	   It is the user's responsibility the avoid conflicts of the
	   accepted abbreviations of any two items.
	2) The type which may be an abbreviation of one of the following:
	      "Character", "Short", "Integer", "Long" (signed integer types),
	      "UCharacter", "UShort", "UInteger, "ULong" (unsigned types)
	      "FLoat", "Real", "Double" (floating point, "Real" == "Double"),
	      "Text" (simple text, character string),
	      "FUnction" (a function reference, not a data reference).
	3) The number of data element. Must be -1 for "FUnction" type.
	   The terminating '\0' in characters strings should be included.
	4) A data pointer of any type. Must be NULL for "FUnction type.
	   If the data should be dynamically allocated by the configuration
	   software it should be a pointer to the pointer that should
	   be set. Allocated data is initialized with '0's.
	5) A function pointer. Must not be NULL except for "FUnction" type
	   and is optional (may be NULL) for data type entries.
	   For the "Function" type, the data (normally a character string)
	   is passed as the only argument. For data type entries,
	   the associated functions are called with an extended
	   calling syntax.
	6) A pointer to a character string with the default initialization
	   values or NULL.
	7) A pointer to a character string with a lower bound value or NULL.
	8) A pointer to a character string with an upper bound value or NULL.
	9) An integer where any of the following flags may be combined
	   by a bitwise OR '|':
	      CFG_REQUIRE_DATA
	      CFG_REQUIRE_ALL_DATA
	      CFG_REJECT_MODIFICATION
       10) Reserved. In multi_threaded mode, use
	      CFG_MUTEX(&some_pthread_mutex)
	   if the associated function is not fully reentrant or
	   if a set of functions should only be called one at a time.
       11) Reserved. Do not modify. Is 1 if reconfigured.
@endverbatim
     Components not specified are automatically initialized to NULL or 0.

     The reason why build_config may be called several times (with
     different configuration items each time) is that this way
     the configuration items for each more or less independent part of
     a program may be defined separately and there is no need for global
     data sharing. You only need to call a 'configuration definition
     function' for each part which has its items defined and only calls
     build_config().

     Once the whole configuration items from all parts have been passed
     to build_config(), a single call to init_config() is required to
     make the configuration effective. init_config() first sets those
     initial values declared in the items (if any) and then tries to
     get external data line by line from a function passed to init_config(),
     unless a NULL pointer is passed instead of a function pointer.
     This user-defined function (declared 'char *user_function(void);')
     should return the address of the first character of each line read
     from a configuration file, the command line, or anywhere else, until
     the end of input which the function must indicate by returning a
     NULL pointer. Input lines can be of any length up to 10240 bytes
     and may include a linefeed character as read by fgets().
     Note that there used to be a problem with semicolons in comments,
     which should be fixed now - but beware of possible side-effects.

     Later, configuration data can be changed by calling reconfig() with
     a line of input passed as argument. Configuration data marked as
     'not to be modified' will not be changed. If a configuration item
     is of 'function' type that function will be called with the remaining
     line (after extracting the item name and processing special characters)
     passed as argument.
*/
/* ====================================================================== */

#include "initial.h"
#include "io_basic.h"
#include <ctype.h>
#include "warning.h"
#include <errno.h>
#include <strings.h>

#ifdef _REENTRANT
#include <pthread.h>
#ifndef PTHREAD_ONCE_INIT
#define PTHREAD_ONCE_INIT pthread_once_init
#define PTHREAD_MUTEX_INITIALIZER 0
#endif
#endif

#include "hconfig.h"
#include "io_basic.h"
#include "history.h"

/** Configuration is organized in sections. CONFIG_BLOCK used for bookkeeping of that. */

struct ConfigBlockStruct
{
   const char *section;
   struct ConfigItemStruct *items;
   struct ConfigBlockStruct *next;
   int flag;
};
typedef struct ConfigBlockStruct CONFIG_BLOCK;

static int do_config (CONFIG_ITEM *item, CONST char *data);
static void config_syntax_error (const char *name, const char *text);
static void config_info (const char *name, const char *text);
static int set_config_values (CONFIG_ITEM *item, int first,
   int last, char *text);
static void display_config_current (CONFIG_ITEM *item);
static void display_config_item (CONFIG_ITEM *item);
static int do_reset_func (const char *text);
static int signed_config_val (const char *name, const char *text, const char *lbound,
   const char *ubound, int strict, long *ival);
static int unsigned_config_val (const char *name, const char *text, const char *lbound,
   const char *ubound, int strict, unsigned long *uval);
static int hex_config_val (const char *name, const char *text, const char *lbound,
   const char *ubound, int strict, unsigned long *uval);
static int real_config_val (const char *name, const char *text, const char *lbound,
   const char *ubound, int strict, double *rval);

static int f_show_config (const char *name, CONFIG_VALUES *val);
static int f_lock_config (const char *name, CONFIG_VALUES *val);
static int f_unlock_config (const char *name, CONFIG_VALUES *val);
static int f_limit_config (const char *name, CONFIG_VALUES *val);
static int f_status_config (const char *name, CONFIG_VALUES *val);
static int f_list_config (const char *name, CONFIG_VALUES *val);
static int f_get_config (const char *name, CONFIG_VALUES *val);
static int f_echo (const char *name, CONFIG_VALUES *val);
static int f_warning (const char *name, CONFIG_VALUES *val);
static int f_error (const char *name, CONFIG_VALUES *val);

#ifdef WITH_HCONFIG_BINARY
static int do_binary_config (CONFIG_ITEM *item, IO_BUFFER *iobuf);
static int config_dummy_io_function (unsigned char *buf, long n, int code);
#endif
static int save_config_values (CONFIG_ITEM *item, int first, int last);
static int restore_config_values (CONFIG_ITEM *item, int first, int last);


#ifdef _REENTRANT
static pthread_mutex_t mlock_hconfig = PTHREAD_MUTEX_INITIALIZER;
#endif

/**
  * Internal functions of the hconfig package.
 */
 
static CONFIG_ITEM default_config[] =
{
   /* All missing elements are initialized to zero. So there shouldn't be any unitialized variables,
      even if your compiler wants to tell you something else. */
   { "SHOW",   "FUN", -1, NULL, f_show_config,   NULL, NULL, NULL, 0, NULL, CFG_MUTEX(mlock_hconfig) },
   { "LOCK",   "FUN", -1, NULL, f_lock_config,   NULL, NULL, NULL, 0, NULL, CFG_MUTEX(mlock_hconfig) },
   { "UNLOCK", "FUN", -1, NULL, f_unlock_config, NULL, NULL, NULL, 0, NULL, CFG_MUTEX(mlock_hconfig) },
   { "LIMITS", "FUN", -1, NULL, f_limit_config,  NULL, NULL, NULL, 0, NULL, CFG_MUTEX(mlock_hconfig) },
   { "STATUS", "FUN", -1, NULL, f_status_config, NULL, NULL, NULL, 0, NULL, CFG_MUTEX(mlock_hconfig) },
   { "LIST",   "FUN", -1, NULL, f_list_config,   NULL, NULL, NULL, 0, NULL, CFG_MUTEX(mlock_hconfig) },
   { "GET",    "FUN", -1, NULL, f_get_config,    NULL, NULL, NULL, 0, NULL, CFG_MUTEX(mlock_hconfig) },
   { "ECHO",   "FUN", -1, NULL, f_echo,          NULL, NULL, NULL, 0, NULL, CFG_MUTEX(mlock_hconfig) },
   { "WARNING","FUN", -1, NULL, f_warning,       NULL, NULL, NULL, 0, NULL, CFG_MUTEX(mlock_hconfig) },
   { "ERROR",  "FUN", -1, NULL, f_error,         NULL, NULL, NULL, 0, NULL, CFG_MUTEX(mlock_hconfig) },
   { NULL_CONFIG_ITEM }
};
static CONFIG_BLOCK first_config_block =
   { "_internal_", default_config, (CONFIG_BLOCK *) NULL, 0 };
static int internal_unhooked = 0;

static PFITI history_function;
int config_level;

struct config_specific_data
{
   char default_section[65];
};

static struct config_specific_data config_defaults =
{ 
   "_internal_"
};

struct Binary_Interface_Chain
{
   struct Config_Binary_Item_Interface *interface;
   struct Binary_Interface_Chain *next;
};

static struct Binary_Interface_Chain *bin_chain_root;

#ifdef _REENTRANT

static int delete_config_specific(void);
static void config_destructor (void *whatever);
static void config_func_once(void);
static int create_config_specific(void);
static struct config_specific_data *get_config_specific(void);

/** Global key for thread-specific data. */
static pthread_key_t config_tsd_key;
static pthread_once_t config_key_once = PTHREAD_ONCE_INIT;
/* static pthread_mutex_t mlock_config = PTHREAD_MUTEX_INITIALIZER; */
#define mlock_config mlock_hconfig

static int delete_config_specific()
{
   void *specific;
   
   if ( (specific = pthread_getspecific(config_tsd_key)) == NULL )
      return 0;
   free(specific);
   if ( pthread_setspecific(config_tsd_key,NULL) < 0 )
      return -1;
   return 0;
}

static void config_destructor (void *whatever)
{
   fprintf(stderr,"delete_config_specific() called\n");
   delete_config_specific();
}

static void config_func_once()
{
   fprintf(stderr,"Doing one-time thread-specific config data initialization.\n");
#ifdef OS_LYNX
   pthread_keycreate(&config_tsd_key,config_destructor);
#else
   pthread_key_create(&config_tsd_key,config_destructor);
#endif
}

static int create_config_specific()
{
   void *specific;
   /* Make sure that the key has been set up */
   if ( pthread_once(&config_key_once,config_func_once) != 0 )
   {
      fprintf(stderr,"Thread specific one-time initialization failed.\n");
      return -1;
   }
   /* Any prior data ? This would be a memory leak. */
   if ( (specific = pthread_getspecific(config_tsd_key)) != NULL )
   {
      fprintf(stderr,"Prior thread-specific config data being deleted.\n");
      free(specific);
      pthread_setspecific(config_tsd_key,NULL);
   }
   
   /* New data allocated and initialized with defaults */
   if ( (specific = calloc(1,sizeof(struct config_specific_data))) == NULL )
      return -1;
   memcpy(specific,&config_defaults,sizeof(struct config_specific_data));
   if ( pthread_setspecific(config_tsd_key,specific) != 0 )
      return -1;
   return 0;
}

static struct config_specific_data *get_config_specific()
{
   void *specific;
   
   if ( (specific = pthread_getspecific(config_tsd_key)) == NULL )
   {
#if 0
      fprintf(stderr,"Dynamically creating config thread-specific data.\n");
      fprintf(stderr,"Failure to call config_delete_specific() at thread termination may result in memory leaks.\n");
#endif
      create_config_specific();
      specific = pthread_getspecific(config_tsd_key);
   }
   
   return (struct config_specific_data *) specific;
}

#else
#define get_config_specific() (&config_defaults)
#endif

/* ---------------------- build_config ------------------------- */
/**
 *  Build up the configuration by adding another section of
 *  configuration definitions.
 *
 *  @param  items  Vector of configuration items, which is
 *			   terminated by a NULL_CONFIG_ITEM
 *  @param  section  Name of this configuration section.
 *
 *  @return   0 (O.k.),  -1 (memory allocation failed), -2 (other error)
 *
*/

int build_config (CONFIG_ITEM *items, const char *section)
{
   CONFIG_BLOCK *block = 0;

   /* Nothing to add ? */
   if ( items == (CONFIG_ITEM *) NULL )
      return 0;
      
   if ( section != NULL )
   {
      if ( strcasecmp(section,"all") == 0 )
      {
      	 Warning("Reserved name 'all' is invalid as configuration section");
	 return -2;
      }
      if ( strcasecmp(section,"_internal_") == 0 )
      {
      	 Warning("Reserved name '_internal_' is invalid as configuration section");
	 return -2;
      }
   }

   /* Scan the allocated linked list to find the last block. */
   for ( block=(&first_config_block); block->items != (CONFIG_ITEM *) NULL &&
         block->next != (CONFIG_BLOCK *) NULL; block=block->next )
      ;  /* Just find the end of the block list */

   /* Except for the very first time, we have to allocate the new block. */
   if ( block->items != (CONFIG_ITEM *) NULL )
   {
      if ( (block->next = (CONFIG_BLOCK *) calloc(1,sizeof(CONFIG_BLOCK))) ==
           (CONFIG_BLOCK *) NULL )
      {
         Warning("Insufficient memory for building configuration setup.");
         return -1;
      }
      block = block->next;
      block->next = (CONFIG_BLOCK *) NULL;
   }
   block->items = items;
   if ( section != (char *) NULL )
   {
      char *new_section;
      if ( (new_section = (char *) malloc(strlen(section)+1)) ==
           (char *) NULL )
      {
         Warning("Insufficient memory for building configuration setup.");
         return -1;
      }
      strcpy(new_section,section);
      block->section = new_section;
   }
   else
      block->section = (char *) NULL;
   return 0;
}

/* ---------------------- init_config --------------------------- */
/**
 *  @short Initialize the configuration after all build_config() calls.
 *
 *  Initialize the configuration after all sections
 *  have been supplied via build_config(). A function
 *  may be specified for reading external configuration data
 *  after the internal specifications have been processed.
 *  This function may be called only once.
 *
 *  @param   fptr	   Pointer to function that returns a
 *			   string pointer as long as external
 *			   configuration data is available, and
 *			   NULL when no more data is available.
 *			   fptr may be NULL if no such function
 *			   should be called.
 *
 *  @return  0 (O.k.),   -1 (called a second time or
 *	   invalid configuration data)
 *
 */

int init_config (char *(*fptr) (void))
{
   CONFIG_BLOCK *block;
   CONFIG_ITEM *itemlist, *item;
   int item_no, i;
   char c;
   int bad_syntax, bad_config, bad_init;
   char *cfg_data;
   static int done = 0;

   if ( first_config_block.items == (CONFIG_ITEM *) NULL )
      return 0;

#ifdef _REENTRANT
   pthread_mutex_lock(&mlock_config);
#endif

   if ( done )
   {
      Warning("Invalid attempt to re-initialize configuration data");
#ifdef _REENTRANT
      pthread_mutex_unlock(&mlock_config);
#endif
      return -1;
   }

   done = 1;
#ifdef _REENTRANT
   pthread_mutex_unlock(&mlock_config);
#endif
   bad_syntax = 0;
   for ( block=(&first_config_block); block != (CONFIG_BLOCK *) NULL;
         block=block->next )
   {
      if ( (itemlist=block->items) == (CONFIG_ITEM *) NULL )
         break;
      for (item_no=0; itemlist[item_no].name != (char *) NULL; item_no++)
      {
         item = &itemlist[item_no];
         item->flags &= 31;
         for ( i=0; (c=item->name[i]) != '\0'; i++ )
            if ( !(tolower((int)c) >= 'a' && tolower((int)c) <= 'z') &&
                 !isdigit((int)c) && c != '_' )
            {
	       char message[1024];
               snprintf(message,sizeof(message),
    "Invalid character 0x%02x in configuration item '%s' (section '%s').",
                   c,item->name,block->section);
               Warning(message);
               bad_syntax = 1;
               continue;
            }
	 if ( block->section != NULL && 
	      strcmp(block->section,"_internal_") != 0 )
	   if ( strcasecmp(item->name,"all") == 0 ||
		strcasecmp(item->name,"set") == 0 ||
		strcasecmp(item->name,"reset") == 0 ||
		strcasecmp(item->name,"lock") == 0 ||
		strcasecmp(item->name,"unlock") == 0 ||
		strcasecmp(item->name,"list") == 0 ||
		strcasecmp(item->name,"limits") == 0 ||
		strcasecmp(item->name,"status") == 0 ||
		strcasecmp(item->name,"get") == 0 )
	 {
	    char message[1024];
            snprintf(message,sizeof(message),
	       "Name '%s' is invalid (reserved name)\n",item->name);
	    Warning(message);
            bad_syntax = 1;
            continue;
	 }
         if ( item->type == (char *) NULL )
         {
	    char message[1024];
            snprintf(message,sizeof(message),
	        "Type missing in configuration item '%s' of section [%s].",
                item->name, block->section);
            Warning(message);
            bad_syntax = 1;
            continue;
         }
         if ( abbrev(item->type,"FUnction") )
         {
            item->internal.itype = 50;
            item->type = "FUnc";
            if ( item->size != -1 || item->data != (void *) NULL ||
                 item->lbound != (char *) NULL ||
                 item->ubound != (char *) NULL )
            {
	       char message[1024];
               snprintf(message,sizeof(message),
                  "Invalid function type configuration item '%s' of section [%s].",
                  item->name, block->section);
               Warning(message);
               bad_syntax = 1;
               continue;
            }
	    item->internal.values.itype = item->internal.itype;
	    item->internal.values.name = item->name;
	    item->internal.values.section = block->section;
	    continue;
         }
         else if ( abbrev(item->type,"Short") )
         {
            item->internal.itype = 1;
	    item->internal.elem_size = sizeof(short);
            item->type = "Short";
         }
         else if ( abbrev(item->type,"UShort") )
         {
            item->internal.itype = 11;
	    item->internal.elem_size = sizeof(short);
            item->type = "UShort";
         }
         else if ( abbrev(item->type,"XShort") )
         {
            item->internal.itype = 21;
	    item->internal.elem_size = sizeof(short);
            item->type = "XShort";
         }
         else if ( abbrev(item->type,"Integer") )
         {
            item->internal.itype = 2;
	    item->internal.elem_size = sizeof(int);
            item->type = "Int";
         }
         else if ( abbrev(item->type,"UInteger") )
         {
            item->internal.itype = 12;
	    item->internal.elem_size = sizeof(int);
            item->type = "UInt";
         }
         else if ( abbrev(item->type,"XInteger") )
         {
            item->internal.itype = 22;
	    item->internal.elem_size = sizeof(int);
            item->type = "XInt";
         }
         else if ( abbrev(item->type,"Long") )
         {
            item->internal.itype = 3;
	    item->internal.elem_size = sizeof(long);
            item->type = "Long";
         }
         else if ( abbrev(item->type,"ULong") )
         {
            item->internal.itype = 13;
	    item->internal.elem_size = sizeof(long);
            item->type = "ULong";
         }
         else if ( abbrev(item->type,"XLong") )
         {
            item->internal.itype = 23;
	    item->internal.elem_size = sizeof(long);
            item->type = "XLong";
         }
         else if ( abbrev(item->type,"Character") )
         {
            item->internal.itype = 4;
	    item->internal.elem_size = sizeof(char);
            item->type = "Char";
         }
         else if ( abbrev(item->type,"UCharacter") )
         {
            item->internal.itype = 14;
	    item->internal.elem_size = sizeof(char);
            item->type = "UChar";
         }
         else if ( abbrev(item->type,"XCharacter") )
         {
            item->internal.itype = 24;
	    item->internal.elem_size = sizeof(char);
            item->type = "UChar";
         }
         else if ( abbrev(item->type,"FLoat") )
         {
            item->internal.itype = 31;
	    item->internal.elem_size = sizeof(float);
            item->type = "FLoat";
         }
         else if ( abbrev(item->type,"Real") ||
               abbrev(item->type,"Double") )
         {
            item->internal.itype = 32;
	    item->internal.elem_size = sizeof(double);
            item->type = "Double";
         }
         else if ( abbrev(item->type,"Text") )
         {
            item->internal.itype = 30;
	    item->internal.elem_size = sizeof(char);
            item->type = "Text";
         }
	 else if ( atof(item->type) > 0 )
	 {
	    item->internal.itype = -atof(item->type);
	    if ( (item->internal.bin_interface =
	          find_config_binary_interface(-item->internal.itype)) 
		  == NULL )
	    {
	       char message[1024];
               snprintf(message,sizeof(message),
	       "No binary interface for type '%s' of configuration item '%s' of section [%s].",
        	  item->type,item->name, block->section);
               Warning(message);
               bad_syntax = 1;
               continue;
	    }
	    item->internal.elem_size = item->internal.bin_interface->elem_size;
	 }
         else
         {
	    char message[1024];
            snprintf(message,sizeof(message),
	       "Invalid type '%s' of configuration item '%s' of section [%s].",
               item->type, item->name, block->section);
            Warning(message);
            bad_syntax = 1;
            continue;
         }
         if ( item->size <= 0 && item->internal.itype != 50 )
         {
	    char message[1024];
            snprintf(message,sizeof(message),
	       "Invalid data size for configuration item '%s' of section [%s].",
               item->name, block->section);
            Warning(message);
            bad_syntax = 1;
            continue;
         }
	 item->internal.bin_alloc_elements = 0;
         if ( item->data == (void *) NULL && item->internal.itype != 50 )
         {
	    char message[1024];
	    
	    /* Complex item types may need some initialization. */
	    if ( item->internal.itype < 0 && 
	         item->internal.bin_interface != NULL )
	    {
	       if ( item->internal.bin_interface->new_func != NULL )
	          item->data = 
	          (   *item->internal.bin_interface->new_func)(item->size,
		        item->internal.bin_interface->io_item_type);
	       else
	          item->data = calloc(item->size,item->internal.elem_size);
	       item->internal.bin_alloc_elements = item->size;
	    }
	    else /* Standard types are simply set to 0. */
	       item->data = calloc(item->size,item->internal.elem_size);
	    if ( item->data == NULL )
	    {
               snprintf(message,sizeof(message),
                "Dynamic allocation for configuration item '%s' of section [%s] failed.",
                item->name, block->section);
               Warning(message);
               bad_syntax = 1;
               continue;
	    }
         }
      	 if ( item->internal.itype != 30 )
	    item->internal.values.max_mod = item->size;
	 else
	    item->internal.values.max_mod = 1;
	 item->internal.values.list_mod = 
	    (int *) calloc(item->internal.values.max_mod,sizeof(int));
	 item->internal.values.mod_flag =
	    (unsigned char *) calloc(item->internal.values.max_mod,sizeof(char));
	 if ( item->internal.itype < 0 && 
	      item->internal.bin_interface != NULL )
	 {
	    if ( item->internal.bin_interface->new_func != NULL )
	       item->internal.values.data_saved = 
		  (*item->internal.bin_interface->new_func)(item->size,
	             item->internal.bin_interface->io_item_type);
	    else
	       item->internal.values.data_saved =
		  calloc(item->size,item->internal.elem_size);
	 }
      	 else
	    item->internal.values.data_saved =
	       calloc(item->size,item->internal.elem_size);
	 item->internal.values.data_changed = item->data;
	 
	 if ( item->internal.values.list_mod == NULL ||
	      item->internal.values.mod_flag == NULL ||
	      item->internal.values.data_saved == NULL ||
	      item->internal.values.data_changed == NULL )
	 {
	    char message[1024];
            snprintf(message,sizeof(message),
               "Dynamic allocation for configuration item '%s' of section [%s] failed.",
               item->name, block->section);
            Warning(message);
            bad_syntax = 1;
            continue;
	 }
	 item->internal.values.itype = item->internal.itype;
	 item->internal.values.name = item->name;
	 item->internal.values.section = block->section;
	 item->internal.values.elements = item->size;
	 item->internal.values.elem_size = item->internal.elem_size;
#ifdef _REENTRANT
      	 pthread_mutex_init (&item->internal.mutex_setval,NULL);
      	 pthread_mutex_init (&item->internal.mutex_item,NULL);
#endif
      }
      if ( bad_syntax )
         block->flag = -1;
   }
   if ( bad_syntax )
   {
      first_config_block.items = (CONFIG_ITEM *) NULL;
      first_config_block.next = (CONFIG_BLOCK *) NULL;
      Warning("Configuration data ignored due to syntax or other error(s)");
      return -1;
   }

   for ( block=(&first_config_block); block != (CONFIG_BLOCK *) NULL;
         block=block->next )
   {
      if ( (itemlist=block->items) == (CONFIG_ITEM *) NULL )
         break;
      for (item_no=0; itemlist[item_no].name != (char *) NULL; item_no++)
      {
         item = &itemlist[item_no];
         item->flags &= ~CFG_INITIALIZED;
         item->flags &= ~CFG_ALL_INITIALIZED;
         item->flags &= ~CFG_NOT_INITIAL;
         if ( item->initial != (char *) NULL )
            if ( do_config(item,item->initial) != 0 )
         {
	    char message[1024];
            snprintf(message,sizeof(message),
               "Invalid internal configuration data for item '%s' of section [%s].",
               item->name, block->section);
            Warning(message);
	    bad_init = 1;
            block->flag = -2;
         }
      }
   }

   bad_config = 0;
   config_level = 0;
   if ( fptr != NULL )
   {
      while ( (cfg_data = (*fptr)()) != (char *) NULL )
      {
         if ( reconfig(cfg_data) != 0 )
            bad_config = 1;
      }
   }
   config_level = 1;

   bad_init = 0;
   for ( block=(&first_config_block); block != (CONFIG_BLOCK *) NULL;
         block=block->next )
   {
      if ( (itemlist=block->items) == (CONFIG_ITEM *) NULL )
         break;
      for (item_no=0; itemlist[item_no].name != (char *) NULL; item_no++)
      {
         item = &itemlist[item_no];
         item->flags |= CFG_NOT_INITIAL;
         if ( (item->flags & CFG_REQUIRE_DATA) &&
              !(item->flags & CFG_INITIALIZED) )
         {
	    char message[1024];
            snprintf(message,sizeof(message),
	       "No configuration data for item '%s' of section [%s].",
               item->name, block->section);
            Warning(message);
            bad_init = 1;
            block->flag = -3;
            continue;
         }
         if ( (item->flags & CFG_REQUIRE_ALL_DATA) &&
              !(item->flags & CFG_ALL_INITIALIZED) )
         {
	    char message[1024];
            snprintf(message,sizeof(message),
               "No complete configuration data for item '%s' of section [%s].",
	       item->name, block->section);
            Warning(message);
            bad_init = 1;
            block->flag = -3;
            continue;
         }
      }
      if ( block->flag < 0 )
      {
	 char message[1024];
         snprintf(message,sizeof(message),
      	    "Bad items in configuration section [%s].",block->section);
         Warning(message);
         bad_init = 1;
      }
   }
   if ( bad_config || bad_init )
   {
      Warning("Error in configuration data or incomplete configuration");
      return -1;
   }

   return 0;
}

/* ----------------------- unhook_internal ---------------------- */
/**
    Disable access to internal functions via configuration.
*/

void unhook_internal()
{
#ifdef _REENTRANT
   pthread_mutex_lock(&mlock_hconfig);
#endif
   internal_unhooked = 1;
#ifdef _REENTRANT
   pthread_mutex_unlock(&mlock_hconfig);
#endif
}

/* ----------------------- rehook_internal ---------------------- */
/**
    Enable access again to internal functions via configuration.
*/

void rehook_internal()
{
#ifdef _REENTRANT
   pthread_mutex_lock(&mlock_hconfig);
#endif
   internal_unhooked = 0;
#ifdef _REENTRANT
   pthread_mutex_unlock(&mlock_hconfig);
#endif
}

/* ------------------------ do_reset_func ----------------------- */

static int do_reset_func (const char *text)
{
   char keyword[130];
   int ipos=0;
   
   Warning("Reset: Not yet implemented");
   while ( getword(text,&ipos,keyword,129,',','%') > 0 )
   {
      if ( keyword[0] == '[' && keyword[strlen(keyword)-1] == ']' )
      {
         /* Reset a whole section */
      }
      else
      {
      	 /* Reset a single item */
      }
   }
   return 0;
}

/* ------------------------ reload_config ----------------------- */
/**
 *  Reload some configuration using the file name/preprocessor
 *  as set up for init_config() or with different file etc.
 *
 *  @param   fptr	   Pointer to function that returns a
 *			   string pointer as long as external
 *			   configuration data is available, and
 *			   NULL when no more data is available.
 *
 *  @return  0 (O.k.),   -1 (invalid configuration data)
 *
 */

int reload_config (char *(*fptr) (void))
{
   int bad_config = 0;
   char *cfg_data;

   if ( fptr == NULL )
   {
      Warning("Cannot reload configuration without input function.");
      return -1;
   }
 
   while ( (cfg_data = (*fptr)()) != (char *) NULL )
   {
      if ( reconfig(cfg_data) != 0 )
         bad_config = 1;
   }

   if ( bad_config )
   {
      Warning("Error in configuration data");
      return -1;
   }

   return 0;
}

/* ------------------------ find_config_item --------------------- */
/**
 *  Find a configuration item by its name (mainly for internal usage).
 *
 *  @param name Item name or block:name
 *
 *  @return Pointer to (first) configuration item found or NULL.
*/

CONFIG_ITEM *find_config_item (const char *name)
{
   CONFIG_BLOCK *block;
   CONFIG_ITEM *itemlist, *item;
   int item_no, i;
   char section[65];
   const char *cp;

   if ( (cp = strchr(name,':')) != (char *) NULL )
   {
      for (i=0; i<64 && name+i<cp; i++)
         section[i] = name[i];
      section[i] = '\0';
      name = cp+1;
   }
   else
      section[0] = '\0';

   for ( block=(&first_config_block); block != (CONFIG_BLOCK *) NULL;
         block=block->next )
   {
      if ( internal_unhooked && block == &first_config_block )
      	 continue;
      if ( (itemlist=block->items) == (CONFIG_ITEM *) NULL )
         break;
      if ( *section )
      {
         if ( block->section == (char *) NULL )
            continue;
         if ( strcmp(block->section,section) != 0 )
            continue;
      }
      for (item_no=0; itemlist[item_no].name != (char *) NULL; item_no++)
      {
         item = &itemlist[item_no];
         if ( abbrev(name,item->name) )
            return(item);
      }
   }

   return((CONFIG_ITEM *) NULL);
}

/* ---------------------- verify_config_section --------------------- */

int verify_config_section (char *section)
{
   CONFIG_BLOCK *block;

   if ( section == (char *) NULL )
      return -1;

   for ( block=(&first_config_block); block != (CONFIG_BLOCK *) NULL;
         block=block->next )
   {
      if ( block->items == (CONFIG_ITEM *) NULL )
         break;
      if ( *section )
      {
         if ( block->section == (char *) NULL )
            continue;
         if ( strcmp(block->section,section) == 0 )
            return(block->flag);
      }
   }
   return -1;
}

/* ------------------------- set_config_history --------------------- */
/**
 *  Set a function for recording the history of the configuration settings.
 *
 *  @param  fptr  --  Pointer to function of type
 *			  'int fptr(char *text,int flag)' where
 *			  'text' is the configuration line and
 *			  flag is 0 for configuration file processing
 *			  and 1 for latre reconfiguration.
 *
 *  @return  0
 *
*/

int set_config_history (PFITI fptr)
{
   history_function = fptr;
   return 0;
}

/* ------------------------------ reconfig -------------------------- */
/**
 *  Modify the configuration after init_config() has been called.
 *
 *  @param  text    String consisting of configuration keyword
 *			  (separated by a blank or '=' from the rest)
 *			  and the corresponding data.
 *
 *  @return  0 (O.k.),  -1 (invalid or undefined configuration
 *		   keyword or error in the data)
 *
*/

int reconfig (char *text)
{
   char keyword[129], keyword2[198];
   int i, j, ipos, rc, lk;
   CONFIG_BLOCK *block;
   CONFIG_ITEM *item;
   char section[65];
   char *default_section;
   char *s;
   struct config_specific_data *spec = get_config_specific();

   section[0] = '\0';
   if ( spec == NULL )
      default_section = section;
   else
      default_section = spec->default_section;

   /* A leading word 'SET' (case insensitive) is ignored */
   for ( s=text; *s==' ' || *s=='\t'; s++ )
      ;
   if ( strncasecmp(s,"SET ",4) == 0 || strncasecmp(s,"SET\t",4) == 0 )
      text = s+4;
      
   /* Reset must be handled right away because it calls back to reconfig() */
   if ( strncasecmp(s,"RESET ",6) == 0 || strncasecmp(s,"RESET\t",6) == 0 )
      return do_reset_func(s+6);
   else if ( strcasecmp(s,"RESET") == 0 )
      return do_reset_func(s+5);

   /* Extract the first keyword */
   ipos = 0;
   if ( (rc = getword(text,&ipos,keyword,128,'=','%')) == 0 )
      return 0;  /* An empty or comment line */
   if ( rc < 0 )
   {
      Warning("Invalid or missing configuration keyword");
      return -1;
   }
   lk = strlen(keyword);
   /* A line for setting the default section? */
   if ( rc == 2 && (keyword[lk-1] == ':' || 
      (keyword[0] == '[' && keyword[lk-1] == ']')) )
   {
      default_section[0] = '\0';
      j = 0;
      if ( keyword[lk-1] == ']' && keyword[0] == '[' )
      	 j++;
      for ( i=0; i<64 && keyword[i+j]!=':' && keyword[i+j]!=']' && 
      	    keyword[i+j]!='\0'; i++ )
         default_section[i] = keyword[i+j];
      default_section[i] = '\0';
      for ( block=(&first_config_block); block != (CONFIG_BLOCK *) NULL;
            block=block->next )
      {
	 if ( block->items == (CONFIG_ITEM *) NULL )
            break;
         if ( block->section == (char *) NULL )
            continue;
         if ( strcmp(block->section,default_section) != 0 )
            continue;

      	 if ( history_function != NULL )
      	    push_config_history(keyword,config_level);
      	 return 0;
      }
      {
      	 char message[1024];
	 snprintf(message,1024,"No such section: %s",default_section);
	 Warning(message);
	 *default_section = '\0';
	 return -1;
      }
   }

   /* If the keyword has no explicit section name and a default */
   /* section is available, try to locate the item in the default */
   /* section first. */
   if ( strchr(keyword,':') == (char *) NULL && *default_section != '\0' )
   {
      strncpy(keyword2,default_section,sizeof(keyword2)-3);
      strcat(keyword2,":");
      lk = strlen(keyword2);
      strncpy(keyword2+lk,keyword,sizeof(keyword2)-lk-1);
      item = find_config_item(keyword2);
   }
   else
      item = (CONFIG_ITEM *) NULL;

   /* If nothing appropriate in the default section, search in all sections */
   if ( item == (CONFIG_ITEM *) NULL )
      if ( (item = find_config_item(keyword)) == (CONFIG_ITEM *) NULL )
      {
	 char message[1024];
         snprintf(message,sizeof(message),
	    "Unknown configuration item '%s'.",keyword);
         Warning(message);
         return -1;
      }

   if ( item->internal.locked )
   {
      char message[1024];
      snprintf(message,sizeof(message),
	 "Configuration item '%s' is locked. Unlock before using it.",keyword);
      Warning(message);
      return -1;
   }

   if ( history_function != NULL )
   {
      char message[1024];
      snprintf(message,sizeof(message),"%s %s",item->name,text+ipos);
      push_config_history(message,config_level);
   }
   
   rc = 0;
   if ( do_config(item,text+ipos) != 0 )
      rc = -1;

   return rc;
}

/* --------------------------- do_config ---------------------------- */
/**
 *   Internal configuration function.
*/

static int do_config (CONFIG_ITEM *item, CONST char *line)
{
   int ichar, len, nextnum, init, ipos, jpos, kpos, rc;
   int first, last, i, j, ilist;
   char c;
   char list[128], value[65], range[33], lidx[21], uidx[21], copy[10241];
   char *text, *tc, *xcopy;

   if ( (item->flags & CFG_REJECT_MODIFICATION) &&
        (item->flags & CFG_NOT_INITIAL) )
   {
      char message[1024];
      snprintf(message,sizeof(message),
      	 "Modification of configuration item '%s' rejected",
         item->name);
      Warning(message);
      return -1;
   }

   if ( item->internal.itype < 0 )
   {
      char message[1024];
      snprintf(message,sizeof(message),
         "Item '%s' is binary (type %d) and cannot be accessed in text mode.",
	 item->name, -item->internal.itype);
      Warning(message);
      return -1;
   }

   len = strlen(line);
   /* Short lines are copied to the stack. */
   if ( (size_t) len < sizeof(copy) )
   {
      strcpy(copy,line);
      xcopy = NULL;
      text = copy;
   }
   else /* For longer lines we take the extra effort of allocating */
   {
      if ( (xcopy = (char *) malloc(len+1)) == NULL )
      {
      	 Warning("Memory allocation failed in do_config()\n");
	 return -1;
      }
      strcpy(xcopy,line);
      text = xcopy;
   }

#ifdef _REENTRANT
   /* There should be only one return beyond this point. */
   if ( item->mutex_pointer != NULL )
      pthread_mutex_lock(item->mutex_pointer);
   /* We never want to mix up values from two calls or show values */
   /* just while we change values. */
   pthread_mutex_lock(&item->internal.mutex_item);
#endif

   /* Skip leading whitespace */
   while (*text == ' ' || *text == '\t' )
      text++;
   len = strlen(text);
   /* Process text input to the end */
   for ( i=j=0; i<len; i++, j++ )
   {
      c = text[i];
      if ( c == '\\' ) /* Check for backslash-escaped characters */
      {
         if ( text[i+1] == '\0' || text[i+1] == '\n' || text[i+1] == '\r' )
         {
            text[j] = '\0';
            break;
         }
         else if ( text[i+1] == 'n' )
            text[j] = '\n';
         else if ( text[i+1] == 'r' )
            text[j] = '\r';
         else if ( text[i+1] == 't' )
            text[j] = '\t';
         else if ( text[i+1] == '0' )
            text[j] = '\0';
         else
            text[j] = text[i+1];
         for (ichar=1; (text[j+ichar] = text[i+1+ichar]) != '\0'; ichar++)
            ;  /* No-op */
      }
      else if ( c=='\n' || c=='\r' || c=='%' )
      {
         text[j] = '\0';
         break;
      }
   }
   /* Squash all trailing whitespace */
   for (ichar=len=strlen(text)-1; ichar>=0 &&
        (text[ichar] == '\t' || text[ichar] == ' '); ichar--)
      text[ichar] = '\0';

   if ( abbrev(item->type,"Text") || abbrev(item->type,"FUnction") )
   {
      len = strlen(text);
      if ( (text[0] == '\'' && text[len-1] == '\'') ||
           (text[0] == '"' && text[len-1] == '"') )
      {
         text[len-1] = '\0';
         text++;
         len -= 2;
      }
      if ( abbrev(item->type,"FUnction") && item->function != NULL )
      {
#if defined(_REENTRANT) && !defined(STRICT_HCONFIG_LOCKING)
      	 /* Unless a user-defined mutex item->mutex_pointer is present, */
	 /* the mutex is temporarily released since the function is */
	 /* apparently meant to be thread-safe. */
      	 pthread_mutex_unlock(&item->internal.mutex_item);
#endif
         rc = (*item->function)(text,(CONFIG_VALUES *) NULL);
#if defined(_REENTRANT) && !defined(STRICT_HCONFIG_LOCKING)
         /* Regain the mutex just to release it again at the end. */
      	 pthread_mutex_lock(&item->internal.mutex_item);
#endif
      	 if ( rc < 0 )
	    goto error_return;
      }
      else /* Text item */
      {
         if ( len >= item->size )
            text[item->size-1] = '\0';
#ifdef _REENTRANT
      	 pthread_mutex_lock(&item->internal.mutex_setval);
#endif
	 strncpy((char *)item->internal.values.data_saved,(char *)item->data,item->size);
         strcpy((char *)item->data,text);
#ifdef _REENTRANT
      	 pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
	 if ( item->function != NULL )
	 {
	    item->internal.values.nmod = 1;
	    item->internal.values.list_mod[0] = 0;
	    item->internal.values.mod_flag[0] = 1;
	    item->internal.values.elem_size = strlen((char *)item->data);
            /* Note that unlike for pure functions temporary release */
	    /* of the mutex is not thread-safe (even if the function */
	    /* itself is thread-safe) because item->internal.values */
	    /* elements could be modified by one thread while the other */
	    /* one just works with them. */
	    rc = (*item->function)(text,&item->internal.values);
	    item->internal.values.nmod = 0;
	    item->internal.values.mod_flag[0] = 0;
	    item->internal.values.elem_size = 1;
	    if ( rc < 0 )
	    {
#ifdef _REENTRANT
      	       pthread_mutex_lock(&item->internal.mutex_setval);
#endif
	       strcpy((char *) item->data,
	          (char *)item->internal.values.data_saved);
#ifdef _REENTRANT
      	       pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
	       goto error_return;
	    }
	 }
      }
      item->flags |= CFG_INITIALIZED | CFG_ALL_INITIALIZED;
      goto good_return;
   }

   item->internal.values.nmod = 0;

   nextnum = init = 0;
   ipos = 0;
   for (;;)
   {
      *list = '\0';
      /* Compound range lists must be in parenthesis, like '(1-4,7,13-20):' */
      if ( text[ipos] == '(' )
      {
         ipos++;
         for ( ilist=0; ilist<127 && text[ipos]!='\0' && text[ipos]!=')';
               ilist++, ipos++ )
            list[ilist] = text[ipos];
         if ( ilist>=127 || text[ipos]!=')' || text[ipos+1]!=':' )
         {
            config_syntax_error(item->name,"Wrong usage of parentheses");
	    goto error_return;
         }
         list[ilist] = '\0';
         ipos += 2;
         while ( text[ipos] == ' ' || text[ipos] == '\t' )
            ipos++;
      }
      if ( (rc = getword(text,&ipos,value,65,',','%')) < 0 )
      {
         config_syntax_error(item->name,"String invalid or too long");
         goto error_return;
      }
      else if ( rc == 0 )
         break;
      if ( *list )
         ; /* No-op for list in parenthesises */
      /* Simple range lists like '1-4:' don't need parenthesises. */
      else if ( (tc=strchr(value,':')) != (char *) NULL )
      {
         *tc = '\0';
         if ( strcasecmp(value,"all") == 0 )
            sprintf(list,"0-%d",item->size-1);
         else
            strcpy(list,value);
         len = strlen(tc+1);
         for (ichar=0; ichar<=len; ichar++)
            value[ichar] = *(tc+1+ichar);
         if ( !*value )
         {
            if ( (rc = getword(text,&ipos,value,65,',','%')) < 0 )
            {
               config_syntax_error(item->name,"String invalid or too long");
	       goto error_return;
            }
            else if ( rc == 0 )
               break;
         }
      }
      /* Repetition factor instead of range list, like '6*2.46' */
      else if ( (tc=strchr(value,'*')) != (char *) NULL )
      {
         *tc = '\0';
         if ( (ilist = atoi(value)) <= 0 )
         {
            config_syntax_error(item->name,"Invalid repetition factor");
	    goto error_return;
         }
         len = strlen(tc+1);
         for (ichar=0; ichar<=len; ichar++)
            value[ichar] = *(tc+1+ichar);
         if ( !*value )
         {
            if ( (rc = getword(text,&ipos,value,65,',','%')) < 0 )
            {
               config_syntax_error(item->name,"String invalid or too long");
	       goto error_return;
            }
            else if ( rc == 0 )
               break;
         }
         /* Convert to a range list. */
         sprintf(list,"%d-%d",nextnum,nextnum+ilist-1);
      }

#ifdef _REENTRANT
      pthread_mutex_lock(&item->internal.mutex_setval);
#endif

      /* If a list of ranges is present, assign the value to each element */
      if ( *list )
      {
         jpos = 0;
         while ( (rc = getword(list,&jpos,range,32,',','%')) > 0 )
         {
            if ( is_unsigned_number(range) )
            {
               nextnum = atoi(range);
               if ( set_config_values(item,nextnum,nextnum,value) != 0 )
	       {
#ifdef _REENTRANT
      	          pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
	          goto error_return;
	       }
	       if ( !init )
                  init = 1;
            }
            else
            {
               *uidx = *lidx = '\0';
               kpos = 0;
               if ( getword(range,&kpos,lidx,20,'-','%') == 1 )
                  if ( is_unsigned_number(lidx) )
                     if ( getword(range,&kpos,uidx,20,'-','%') == 2 )
                        if ( is_unsigned_number(uidx) )
                        {
                           first = atoi(lidx);
                           last = atoi(uidx);
                           nextnum = last+1;
                           if ( set_config_values(item,first,last,value) != 0 )
			   {
#ifdef _REENTRANT
      	             	      pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
			      goto error_return;
			   }
                           if ( first == 0 && last == item->size-1 )
                              item->flags |= CFG_ALL_INITIALIZED;
	             	   if ( !init )
                              init = 1;
                           continue;
                        }
               config_syntax_error(item->name,"Invalid index range");
#ifdef _REENTRANT
      	       pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
	       goto error_return;
            }
         }
         if ( rc < 0 )
         {
            config_syntax_error(item->name,"Invalid index");
#ifdef _REENTRANT
      	    pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
	    goto error_return;
         }
      }
      /* Otherwise, the value is assigned to just one element. */
      else
      {
         if ( set_config_values(item,nextnum,nextnum,value) != 0 )
	 {
#ifdef _REENTRANT
      	    pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
	    goto error_return;
	 }
         nextnum++;
      }

#ifdef _REENTRANT
   pthread_mutex_unlock(&item->internal.mutex_setval);
#endif

   }
   
   /* There may have been overlapping regions, resulting in more */
   /* modifications indicated than we have elements. */
   if ( item->internal.values.nmod > item->internal.values.max_mod )
   {
      item->internal.values.nmod = 0;
      for ( i=0; i<item->internal.values.max_mod; i++)
      	 if ( item->internal.values.mod_flag[i] )
	    item->internal.values.list_mod[item->internal.values.nmod++] = 1;
   }
   

   if ( item->function != NULL )
   {
      /* Note that unlike for pure functions temporary release */
      /* of the mutex is not thread-safe (even if the function */
      /* itself is thread-safe) because item->internal.values */
      /* elements could be modified by one thread while the other */
      /* one just works with them. */
      rc = (*item->function)(NULL,&item->internal.values);
      if ( rc < 0 )
      {
      	 /* Modifications were declined: restore old values */
	 for ( i=0; i<item->internal.values.nmod; i++)
	 {
	    j = item->internal.values.list_mod[i];
#ifdef _REENTRANT
      	    pthread_mutex_lock(&item->internal.mutex_setval);
#endif
      	    restore_config_values(item,j,j);
#ifdef _REENTRANT
      	    pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
	    item->internal.values.mod_flag[j] = 0;
	 }
	 item->internal.values.nmod = 0;
      }
   }
   
   /* Indicate that some values have been initialized */
   if ( init )
      item->flags |= CFG_INITIALIZED;
   /* Indicate if all values have been initialized (if so required) */
   if ( item->flags & CFG_REQUIRE_ALL_DATA && !(item->flags & CFG_ALL_INITIALIZED) )
   {
      if ( item->internal.values.nmod < item->internal.values.max_mod )
      	 init = 0;
      else
	 for ( i=0; i<item->internal.values.max_mod; i++)
      	    if ( !item->internal.values.mod_flag[i] )
	    {
	       if ( init )
	          init = 0;
	       break;
	    }
      if ( init )
         item->flags |= CFG_ALL_INITIALIZED;
   }
   
   /* Clean up modifications flags for next time */
   if ( item->internal.values.nmod > item->internal.values.max_mod/2 )
   {
      for ( i=0; i<item->internal.values.max_mod; i++ )
      	 item->internal.values.mod_flag[i] = 0;
   }
   else
   {
      for ( i=0; i<item->internal.values.nmod; i++ )
      	 item->internal.values.mod_flag[item->internal.values.list_mod[i]] = 0;
   }
   item->internal.values.nmod = 0;

good_return:
   rc = 0;
   goto common_return;
   
error_return:
   rc = -1;
   
common_return:

#ifdef _REENTRANT
   pthread_mutex_unlock(&item->internal.mutex_item);
   if ( item->mutex_pointer != NULL )
      pthread_mutex_unlock(item->mutex_pointer);
#endif

   if ( xcopy != NULL )
      free(xcopy);

   return rc;
}

#ifdef WITH_HCONFIG_BINARY

/* ------------------- config_dummy_io_function --------------------- */

static int config_dummy_io_function (unsigned char *buf, long n, int code)
{
   if ( code != 2 && code != 3 )
      return -1;
   return 0;
}

/* ------------------------ reconfig_binary ------------------------- */

int reconfig_binary (char *buffer, size_t buflen)
{
   char keyword[130], keyword2[198];
   struct config_specific_data *spec = get_config_specific();
   CONFIG_BLOCK *block;
   CONFIG_ITEM *item;
   char section[65];
   char *default_section;
   IO_BUFFER *iobuf;
   IO_ITEM_HEADER item_header;
   uint32_t marker;
   int type, rc, lk;

   section[0] = '\0';
   if ( spec == NULL )
      default_section = section;
   else
      default_section = spec->default_section;

   if ( buffer == NULL || buflen < 16 )
      return -1;
   marker = *((uint32_t *) buffer);
   if ( marker != (uint32_t) 0xD41F8A37 && 
        marker != (uint32_t) 0x378A1FD4 )
   {
      Warning("Unsupported binary data.");
      return -1;
   }
   if ( (iobuf = allocate_io_buffer(16)) == NULL )
   {
      Warning("Memory allocation failed.");
      return -1;
   }
   free(iobuf->buffer);
   iobuf->buffer = buffer;
   iobuf->is_allocated = 0;
   iobuf->buflen = iobuf->min_length = iobuf->max_length = buflen;
   iobuf->user_function = config_dummy_io_function;
   if ( find_io_block(iobuf,&item_header) < 0 )
   {
      Warning("Error in binary configuration data.");
      goto bad_return;
   }
   if ( item_header.type != IO_TYPE_HCONFIG_ENVELOPE )
   {
      Warning("Wrong type of binary configuration data.");
      goto bad_return;
   }
   if ( read_io_block(iobuf,&item_header) < 0 )
   {
      Warning("Error in binary configuration data.");
      goto bad_return;
   }
   item_header.type = IO_TYPE_HCONFIG_ENVELOPE;
   if ( get_item_begin(iobuf,&item_header) < 0 )
      goto bad_return;
   if ( item_header.version != 0 )
      goto bad_return;
   type = next_subitem_type(iobuf);
   /* If the contents are just text, pass it over to the text mode function. */
   if ( type == IO_TYPE_HCONFIG_TEXT )
   {
      char text[4096];
      if ( config_binary_read_text(iobuf,keyword,sizeof(keyword)) < 0 )
      	 goto bad_return;
      if ( (rc=reconfig(text)) < 0 )
      	 goto all_return;
      goto good_return;
   }
   /* For true binary data, we expect first the item name or section:item. */
   if ( type != IO_TYPE_HCONFIG_NAME )
      goto bad_return;
   if ( config_binary_read_text(iobuf,keyword,sizeof(keyword)) < 0 )
      goto bad_return;
   if ( (lk=strlen(keyword)) <= 0 )
      goto bad_return;
   /* Is it just to set the default section? */
   if ( keyword[0] == '[' && keyword[lk-1] == ']' )
   {
      keyword[lk-1] = '\0';
      strncpy(default_section,keyword+1,sizeof(default_section)-1);
      keyword[lk-1] = ']';
      for ( block=(&first_config_block); block != (CONFIG_BLOCK *) NULL;
            block=block->next )
      {
	 if ( block->items == (CONFIG_ITEM *) NULL )
            break;
         if ( block->section == (char *) NULL )
            continue;
         if ( strcmp(block->section,default_section) != 0 )
            continue;

      	 if ( history_function != NULL )
      	    push_config_history(keyword,config_level);
      	 goto good_return;
      }
      {
      	 char message[1024];
	 snprintf(message,1024,"No such section: %s",default_section);
	 Warning(message);
	 *default_section = '\0';
	 goto bad_return;
      }
   }

   /* If the keyword has no explicit section name and a default */
   /* section is available, try to locate the item in the default */
   /* section first. */
   if ( strchr(keyword,':') == (char *) NULL && *default_section != '\0' )
   {
      strncpy(keyword2,default_section,sizeof(keyword2)-3);
      strcat(keyword2,":");
      lk = strlen(keyword2);
      strncpy(keyword2+lk,keyword,sizeof(keyword2)-lk-1);
      item = find_config_item(keyword2);
   }
   else
      item = (CONFIG_ITEM *) NULL;

   /* If nothing appropriate in the default section, search in all sections */
   if ( item == (CONFIG_ITEM *) NULL )
      if ( (item = find_config_item(keyword)) == (CONFIG_ITEM *) NULL )
      {
	 char message[1024];
         snprintf(message,sizeof(message),
	    "Unknown configuration item '%s'.",keyword);
         Warning(message);
         goto bad_return;
      }

   /* If the item is locked, deny any modifications and function calls. */
   if ( item->internal.locked )
   {
      char message[1024];
      snprintf(message,sizeof(message),
	 "Configuration item '%s' of section [%s] is locked. Unlock before using it.",
	 item->name, item->internal.values.section);
      Warning(message);
      goto bad_return;
   }

   /* For binary data, we make very little attempt of writing the history. */
   if ( history_function != NULL )
   {
      char message[1024];
      snprintf(message,sizeof(message),"%s:%s [binary]",
      	 item->internal.values.section, item->name);
      push_config_history(message,config_level);
   }
   
   /* The real work is done by do_binary_config(). */
   if ( (rc = do_binary_config(item,iobuf)) != 0 )
      goto all_return;

good_return:
   rc = 0;
   goto all_return;

bad_return:
   rc = -1;
   
all_return:
   free_io_buffer(iobuf);
   return rc;
}

/* --------------------- do_binary_config ------------------------ */

static int do_binary_config (CONFIG_ITEM *item, IO_BUFFER *iobuf)
{
#define MAX_INDEX_LIST 100
   int nidx, idx_low[MAX_INDEX_LIST], idx_high[MAX_INDEX_LIST];
   int i, j, type, rc, len;
   
   if ( (item->flags & CFG_REJECT_MODIFICATION) &&
        (item->flags & CFG_NOT_INITIAL) )
   {
      char message[1024];
      snprintf(message,sizeof(message),
      	 "Modification of configuration item '%s' rejected",
         item->name);
      Warning(message);
      return -1;
   }

#ifdef _REENTRANT
   /* There should be only one return beyond this point. */
   if ( item->mutex_pointer != NULL )
      pthread_mutex_lock(item->mutex_pointer);
   /* We never want to mix up values from two calls or show values */
   /* just while we change values. */
   pthread_mutex_lock(&item->internal.mutex_item);
#endif

   item->internal.values.binary_config = 1;

   /* If index missing at all, we start at offset 0. */
   idx_low[0] = idx_high[0] = -1;
   
   /* There may many data blocks (optionally preceeded by index lists). */
   for (;;)
   {
      /* If there is no index list, just advance to the next element. */
      if ( next_subitem_type(iobuf) != IO_TYPE_HCONFIG_INDEX )
      {
      	 nidx = 1;
	 idx_low[0] = idx_high[0] + 1;
	 idx_high[0] = idx_low[0];
      }
      else if ( config_binary_read_index(iobuf,&nidx,
	        idx_low,idx_high,MAX_INDEX_LIST) < 0 )
      {
	 Warning("Error in index list");
	 goto bad_return;
      }
      
      /* Check that all index pairs are in the valid range. */
      for ( i=0; i<nidx; i++ )
      {
      	 /* A range from -1 to -1 is a wildcard like '*:' in text mode. */
      	 if ( idx_low[i] == -1 && idx_high[i] == -1 )
	 {
	    idx_low[i] = 0;
	    idx_high[i] = item->internal.values.max_mod-1;
	 }
      	 if ( idx_low[i] < 0 || idx_high[i] >= item->internal.values.max_mod )
	 {
	    Warning("Invalid element index");
	    goto bad_return;
	 }
      }
      
      type = next_subitem_type(iobuf);
      /* Binary-only items. */
      if ( item->internal.itype < 0 )
      {
      	 if ( item->internal.bin_interface == NULL )
	 {
	    Warning("Binary interface missing");
	    goto bad_return;
	 }
      	 /* Binary-only items must match exactly. */
	 if ( item->internal.bin_interface->io_item_type != type )
	 {
	    Warning("Wrong data type for binary data.");
	    goto bad_return;
	 }
	 if ( item->internal.bin_interface->read_func == NULL )
	 {
	    Warning("No function defined for reading this item type");
	    goto bad_return;
	 }
#ifdef _REENTRANT
      	 pthread_mutex_lock(&item->internal.mutex_setval);
#endif
	 /* Save current values of all elements to be modified. */
      	 for ( i=0; i<nidx; i++ )
	    save_config_values(item,idx_low[i],idx_high[i]);
	 /* Read the data to the first 'value' on our list. */
	 rc = (*item->internal.bin_interface->read_func)(
      	    (char *)item->internal.values.data_changed +
      	    idx_low[0]*item->internal.values.elem_size,
	    iobuf,
	    item->internal.bin_interface->io_item_type);
	 /* Did it fail? */
	 if ( rc < 0 )
	 {
	    /* Restore every 'value' marked as modified. */
	    restore_config_values(item,0,item->internal.values.max_mod);
	    item->internal.values.nmod = 0;
#ifdef _REENTRANT
      	    pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
	    goto bad_return;
	 }
	 /* Then copy this 'value' to all others on the list. */
      	 for ( i=0; i<nidx; i++ )
	    for ( j=idx_low[i]; j<=idx_high[i]; j++ )
	    {
	       /* Never copy onto itself. */
	       if ( j == idx_low[0] )
	          continue;
	       /* This may be using a dedicated copy function ... */
	       if ( item->internal.bin_interface->copy_func != NULL )
		  (*item->internal.bin_interface->copy_func)(
		     (char *)item->internal.values.data_changed +
      		     j*item->internal.values.elem_size,
      		     (char *)item->internal.values.data_changed +
      		     idx_low[0]*item->internal.values.elem_size,
		     item->internal.bin_interface->io_item_type);
	       else
	          /* ... or the standard memory copy function. */
      		  memcpy((char *)item->internal.values.data_changed +
      		     j*item->internal.values.elem_size,
      		     (char *)item->internal.values.data_changed +
      		     idx_low[0]*item->internal.values.elem_size,
      		     item->internal.values.elem_size);	       
	    }
#ifdef _REENTRANT
      	 pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
      }
      /* Ordinary element types. Data conversion is supported. */
      else
      {	 
	 /* Save current values of all elements to be modified. */
      	 for ( i=0; i<nidx; i++ )
	    save_config_values(item,idx_low[i],idx_high[i]);
	 
         /* For consistency with text-mode interface the 'FUNCTION' and */
	 /* 'TEXT' type items are only allowed to get text arguments. */
	 /* While other data types have their functions only called after */
	 /* all values are set, the 'FUNCTION' and 'TEXT' type items */
	 /* have their function called for each data block coming. */
	 /* On the other hand any index-ranges are ignored. */
	 if ( (item->internal.itype == 30 || item->internal.itype == 50) &&
	       type != IO_TYPE_HCONFIG_TEXT )
	 {
	    Warning("Unsupported data type. Text required.");
	    continue;
	 }

	 switch ( item->internal.itype )
	 {
	    case 50: /* Function */
	    {
      	       char *text;
	       len = config_binary_text_length(iobuf);
	       if ( len < 0 )
	          goto bad_return;
	       if ( (text = malloc(len+1)) == NULL )
	          goto bad_return;
      	       if ( config_binary_read_text(
		     iobuf,text,len+1) < 0 )
	       {
	          free(text);
      	          break;
	       }
#if defined(_REENTRANT) && !defined(STRICT_HCONFIG_LOCKING)
      	       /* Unless a user-defined mutex item->mutex_pointer is present, */
	       /* the mutex is temporarily released since the function is */
	       /* apparently meant to be thread-safe. */
      	       pthread_mutex_unlock(&item->internal.mutex_item);
#endif
               rc = (*item->function)(text,(CONFIG_VALUES *) NULL);
#if defined(_REENTRANT) && !defined(STRICT_HCONFIG_LOCKING)
               /* Regain the mutex just to release it again at the end. */
      	       pthread_mutex_lock(&item->internal.mutex_item);
#endif
      	       free(text);
	    }
	       break;
	    case 30: /* Text */
#ifdef _REENTRANT
      	       pthread_mutex_lock(&item->internal.mutex_setval);
#endif
      	       /* Save old values. */
	       strncpy((char *)item->internal.values.data_saved,item->data,
	          item->size);
	       /* Get new values. */
	       config_binary_read_text(iobuf,
	          item->internal.values.data_changed, item->size);
#ifdef _REENTRANT
      	       pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
	       if ( item->function != NULL )
	       {
		  item->internal.values.nmod = 1;
		  item->internal.values.list_mod[0] = 0;
		  item->internal.values.mod_flag[0] = 1;
		  item->internal.values.elem_size = strlen(item->data)+1;
        	  /* Note that unlike for pure functions temporary release */
		  /* of the mutex is not thread-safe (even if the function */
		  /* itself is thread-safe) because item->internal.values */
		  /* elements could be modified by one thread while the other */
		  /* one just works with them. */
		  rc = (*item->function)(item->data,&item->internal.values);
		  item->internal.values.nmod = 0;
		  item->internal.values.mod_flag[0] = 0;
		  item->internal.values.elem_size = 1;
		  if ( rc < 0 )
		  {
#ifdef _REENTRANT
        	     pthread_mutex_lock(&item->internal.mutex_setval);
#endif
      	             /* Restore old values. */
		     strcpy((char *) item->data,
	        	(char *)item->internal.values.data_saved);
#ifdef _REENTRANT
      		     pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
		     goto bad_return;
		  }
	       }
	       break;
	       
	    /* For the other types we have two options: */
	    /* binary values or text values. */
	    case 1:  /* Short */
	    case 11: /* UShort */
	    case 21: /* XShort */
	    case 2:  /* Integer */
	    case 12: /* UInteger */
	    case 22: /* XInteger */
	    case 3:  /* Long */
	    case 13: /* ULong */
	    case 23: /* XLong */
	    case 4:  /* Character */
	    case 14: /* Ucharacter */
	    case 24: /* XCharacter */
	    case 31: /* Float */
	    case 32: /* Double */
	       /* Save current values of all elements to be modified. */
#ifdef _REENTRANT
               pthread_mutex_lock(&item->internal.mutex_setval);
#endif
	       for ( i=0; i<nidx; i++ )
		  save_config_values(item,idx_low[i],idx_high[i]);
#ifdef _REENTRANT
               pthread_mutex_unlock(&item->internal.mutex_setval);
#endif

	       if ( type == IO_TYPE_HCONFIG_TEXT )
	       {
	          Warning("Not yet implemented");
	          /* Not yet implemented. */
		  goto bad_return;
	       }
	       else if ( type == IO_TYPE_HCONFIG_NUMBERS )
	       {
	          union ConfigDataPointer indata, locdata;
		  int ntype, nsize, zero_opt, ltype, lsize, itype;
	          int32_t num, n, k;
		  long zero[2];

      	          itype = item->internal.itype;
		  ltype = lsize = 0;

                  /* Our hconfig data types are categorized in the same */
		  /* way like the received data. */
		  if ( itype < 30 )
	          {
		     if ( itype % 10 == 1 )
		     	lsize = sizeof(short);
		     else if ( itype % 10 == 2 )
		     	lsize = sizeof(int);
		     else if ( itype % 10 == 3 )
		     	lsize = sizeof(long);
		     else if ( itype % 10 == 4 )
		     	lsize = sizeof(char);
#ifdef HAVE_64BIT_INT
		     else if ( itype % 10 == 5 )
		     	lsize = sizeof(int64_t);
#endif
      	             if ( itype < 10 )
		     	ltype = 1;
		     else
		     	ltype = 2;
		  }
		  else if ( itype == 31 )
		  {
		     lsize = sizeof(float);
		     ltype = 3;
		  }
		  else if ( itype == 32 )
		  {
		     lsize = sizeof(double);
		     ltype = 3;
		  }
		  if ( lsize == 0 || ltype == 0 )
		     goto bad_return;

	          if ( config_binary_inquire_numbers(iobuf,
		     	&ntype,&nsize,&num,&zero_opt) < 0 )
		     continue;
		  if ( (indata.anything = calloc(nsize,num)) == NULL )
		     goto bad_return;
		  if ( config_binary_read_numbers(iobuf,
		     indata.anything,(size_t) (nsize*num)) < 0 )
		  {
		     free(indata.anything);
		     goto bad_return;
		  }
		  locdata.anything = item->internal.values.data_changed;
		  n = 0;
#ifdef _REENTRANT
        	  pthread_mutex_lock(&item->internal.mutex_setval);
#endif
		  for (i=0; i<nidx; i++)
		  {
		     int ld = idx_high[i]-idx_low[i];

		     /* We might be lucky and have only to copy the bytes. */
		     if ( ( ltype == ntype || (ntype <= 2 && ltype <= 2) ) && 
		     	  lsize == nsize && n + ld < num )
		     {
		     	memcpy(locdata.cdata+idx_low[i]*lsize,
			       indata.cdata+n*nsize, ld*nsize);
			continue;
		     }

		     for (j=idx_low[i]; j<=idx_high[i]; j++)
		     {
		     	if ( n >= num )
			{
			   if ( zero_opt )
			   {
			      zero[0] = zero[1] = 0;
			      k = 0;
			      locdata.ldata = zero;
			   }
			   else
			     k = num-1;
			}
			else
			   k = n;
			if ( config_binary_convert_data(locdata.cdata+j*lsize, 
			      ltype, lsize,
			      indata.cdata + k*nsize, ntype, nsize) < 0 )
			{
#ifdef _REENTRANT
        	     	   pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
			   goto bad_return;
			}
		     }
		  }
#ifdef _REENTRANT
        	  pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
	       }
	       break;
	    
	    default:
	       Warning("Invalid item type.");
	       goto bad_return;
	 }
      }
   }

   if ( item->function != NULL )
   {
      /* Note that unlike for pure functions temporary release */
      /* of the mutex is not thread-safe (even if the function */
      /* itself is thread-safe) because item->internal.values */
      /* elements could be modified by one thread while the other */
      /* one just works with them. */
      rc = (*item->function)(NULL,&item->internal.values);
      if ( rc < 0 )
      {
      	 /* Modifications were declined: restore old values */
	 for ( i=0; i<item->internal.values.nmod; i++)
	 {
	    j = item->internal.values.list_mod[i];
#ifdef _REENTRANT
      	    pthread_mutex_lock(&item->internal.mutex_setval);
#endif
      	    restore_config_values(item,j,j);
#ifdef _REENTRANT
      	    pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
	    item->internal.values.mod_flag[j] = 0;
	 }
	 item->internal.values.nmod = 0;
      }
   }
   
   /* Clean up modifications flags for next time */
   if ( item->internal.values.nmod > item->internal.values.max_mod/2 )
   {
      for ( i=0; i<item->internal.values.max_mod; i++ )
      	 item->internal.values.mod_flag[i] = 0;
   }
   else
   {
      for ( i=0; i<item->internal.values.nmod; i++ )
      	 item->internal.values.mod_flag[item->internal.values.list_mod[i]] = 0;
   }
   item->internal.values.nmod = 0;

/*
good_return:
*/
   rc = 0;
   goto all_return;

bad_return:
   Warning("Setting values or calling the function failed.");
   rc = -1;
   
all_return:
#ifdef _REENTRANT
   pthread_mutex_unlock(&item->internal.mutex_item);
   if ( item->mutex_pointer != NULL )
      pthread_mutex_unlock(item->mutex_pointer);
#endif
   return rc;
}

#endif /* WITH_HCONFIG_BINARY */

/* -------------------- save_config_values ------------------------ */

static int save_config_values (CONFIG_ITEM *item, int first, int last)
{
   int inum;
   
   /* The function-only type has no saved values. */
   if ( item->internal.itype == 50 )
      return 0;

   /* Save previous values */
   for ( inum=first; inum<=last; inum++ )
   {
      /* Copy only if not yet modified */
      if ( !item->internal.values.mod_flag[inum] )
      {
      	 if ( item->internal.itype < 0 && 
	      item->internal.bin_interface != NULL &&
	      item->internal.bin_interface->copy_func != NULL )
	    (*item->internal.bin_interface->copy_func)(
	       (char *)item->internal.values.data_saved +
      	       inum*item->internal.values.elem_size,
      	       (char *)item->internal.values.data_changed +
      	       inum*item->internal.values.elem_size,
	       item->internal.bin_interface->io_item_type);
	 else
      	    memcpy((char *)item->internal.values.data_saved +
      	       inum*item->internal.values.elem_size,
      	       (char *)item->internal.values.data_changed +
      	       inum*item->internal.values.elem_size,
      	       item->internal.values.elem_size);
      	 item->internal.values.mod_flag[inum] = 1;
      	 if ( item->internal.values.nmod < item->internal.values.max_mod )
      	 {
      	    item->internal.values.list_mod[item->internal.values.nmod++] = inum;
      	 }
      }
   }
   return 0;
}

/* -------------------- restore_config_values ------------------------- */

static int restore_config_values (CONFIG_ITEM *item, int first, int last)
{
   int inum;

   /* The function-only type has no saved values. */
   if ( item->internal.itype == 50 )
      return 0;

   /* Restore previous values */
   for ( inum=first; inum<=last; inum++ )
   {
      if ( inum < 0 || inum >= item->internal.values.max_mod )
      	 continue;
      /* Copy it only if it was modified */
      if ( !item->internal.values.mod_flag[inum] )
      	 continue;
      else
      {
      	 if ( item->internal.itype < 0 && 
	      item->internal.bin_interface != NULL &&
	      item->internal.bin_interface->copy_func != NULL )
	    (*item->internal.bin_interface->copy_func)(
	       (char *)item->internal.values.data_changed +
      	       inum*item->internal.values.elem_size,
      	       (char *)item->internal.values.data_saved +
      	       inum*item->internal.values.elem_size,
	       item->internal.bin_interface->io_item_type);
	 else
      	    memcpy((char *)item->internal.values.data_changed +
      	       inum*item->internal.values.elem_size,
      	       (char *)item->internal.values.data_saved +
      	       inum*item->internal.values.elem_size,
      	       item->internal.values.elem_size);
      	 item->internal.values.mod_flag[inum] = 0;
      }
   }
   return 0;
}

/* ---------------------- set_config_values ------------------------- */
/**
 *  Set configuration values (internal usage only).
*/

static int set_config_values (CONFIG_ITEM *item, int first, int last,
   char *text)
{
   int rc, inum;
   union ConfigDataPointer data;
   long ival=0;
   unsigned long uval=0;
   double rval=0.;
   int strict = 0;

   if ( (item->flags & CFG_HARD_BOUND) )
      strict = 1;
   if ( (item->flags & CFG_STRICT_BOUND) )
      strict += 2;

   if ( first < 0 || last >= item->size )
   {
      config_syntax_error(item->name,"Index out of range");
      return -1;
   }
   data.anything = item->data; /* Same as item->internal.values.data_changed */

   if ( item->internal.itype < 0 )
      rc = -1; /* We cannot set binary types yet */
   if ( item->internal.itype < 10 )
      rc = signed_config_val(item->name,text,item->lbound,
         item->ubound,strict,&ival);
   else if ( item->internal.itype < 20 )
      rc = unsigned_config_val(item->name,text,item->lbound,
         item->ubound,strict,&uval);
   else if ( item->internal.itype < 30 )
      rc = hex_config_val(item->name,text,item->lbound,
         item->ubound,strict,&uval);
   else
      rc = real_config_val(item->name,text,item->lbound,
         item->ubound,strict,&rval);
   if ( rc != 0 )
      return -1;

   /* Save previous values */
   save_config_values(item,first,last);

   switch ( item->internal.itype )
   {
      case 1:
         for (inum=first; inum<=last; inum++)
            data.sdata[inum] = (short) ival;
         break;
      case 2:
         for (inum=first; inum<=last; inum++)
            data.idata[inum] = (int) ival;
         break;
      case 3:
         for (inum=first; inum<=last; inum++)
            data.ldata[inum] = ival;
         break;
      case 4:
         for (inum=first; inum<=last; inum++)
            data.cdata[inum] = (char) ival;
         break;
      case 11:
         for (inum=first; inum<=last; inum++)
            data.usdata[inum] = (unsigned short) uval;
         break;
      case 12:
         for (inum=first; inum<=last; inum++)
            data.uidata[inum] = (unsigned int) uval;
         break;
      case 13:
         for (inum=first; inum<=last; inum++)
            data.uldata[inum] = uval;
         break;
      case 14:
         for (inum=first; inum<=last; inum++)
            data.ucdata[inum] = (unsigned char) uval;
         break;
      case 21:
         for (inum=first; inum<=last; inum++)
            data.usdata[inum] = (unsigned short) uval;
         break;
      case 22:
         for (inum=first; inum<=last; inum++)
            data.uidata[inum] = (unsigned int) uval;
         break;
      case 23:
         for (inum=first; inum<=last; inum++)
            data.uldata[inum] = uval;
         break;
      case 24:
         for (inum=first; inum<=last; inum++)
            data.ucdata[inum] = (unsigned char) uval;
         break;
      case 31:
         for (inum=first; inum<=last; inum++)
            data.fdata[inum] = (float) rval;
         break;
      case 32:
         for (inum=first; inum<=last; inum++)
            data.ddata[inum] = rval;
         break;
      default:
         config_syntax_error(item->name,"Invalid data type");
         return -1;
   }

   return 0;
}

/* -------------------- lock_unlock_status ------------------------ */

static int lock_unlock_status (const char *name, int lock)
{
   CONFIG_ITEM *item, *itemlist;
   CONFIG_BLOCK *block;
   const char *s;
   char section[65];
   const char *list = name;
   int l, found, rc = 0, item_no;

   /* Lock a whole section or list of sections ? */
   if ( (name[0] == '[' && name[strlen(name)-1] == ']') || 
      	 strcmp(name,"all") == 0 )
   {
      s = list;
      if ( name[0] == '[' )
      	 s++;
      for (;;)
      {
	 while ( *s == ' ' || *s == '\t' )
      	    s++;
	 for ( l=0; l<64 && (isalnum(*s) || *s=='_'); l++ )
      	    section[l] = *s++;
	 section[l] = '\0';
	 if ( l == 0 )
	    continue;

	 for ( block = (&first_config_block), found = 0; 
	       block != (CONFIG_BLOCK *) NULL;
               block = block->next )
	 {
	    if ( (itemlist=block->items) == (CONFIG_ITEM *) NULL )
               break;
            if ( block->section == (char *) NULL )
               continue;
            if ( strcmp(block->section,section) == 0 ||
		 strcasecmp("all",section) == 0 )
            {
	       found = 1;
	       for ( item_no = 0; 
	             itemlist[item_no].name != (char *) NULL; 
		     item_no++ )
	       {
        	  item = &itemlist[item_no];
		  if ( lock == -1 )
		  {
	             char message[1024];
		     if ( item->internal.values.section != NULL &&
		     	  item->internal.values.name != NULL )
		     {
	        	snprintf(message, 1024, "%s:%s (%s) is %s.\n",
			   item->internal.values.section,
			   item->internal.values.name,
			   item->name,
			   (item->internal.locked==0)?"unlocked":"locked");
			Output(message);
		     }
		     else
		     {
			snprintf(message, 1024, "%s is %s.\n",
			   item->name,
			   (item->internal.locked==0)?"unlocked":"locked");
			Output(message);
		     }
		  }
		  else if ( strcmp(block->section,"_internal_") == 0 )
		  {
	             if ( strcasecmp("all",section) != 0 )
		     {
	        	Warning("Cannot lock or unlock internal configuration functions");
			rc = -1;
		     }
		     continue;
		  }
		  else
		     item->internal.locked = lock; /* Lock or unlock */
	       }
	    }
	 }
	 if ( !found )
	 {
	    char message[1024];
	    sprintf(message,
	       "Cannot %s configuration section %s: no such section",
	       (lock==0)?"unlock":"lock", section);
	    rc = -1;
	 }

	 while ( *s == ' ' || *s == '\t' )
      	    s++;
	 if ( *s++ != ',' )
	    break;
      }
   }
   if ( (item = find_config_item(name)) == (CONFIG_ITEM *) NULL )
   {
      char message[1024];
      snprintf(message,1024,"%s failed: no item '%s' known",
      	 ((lock==-1)?"Status inquiry":(lock==0)?"Unlocking":"Locking"), name);
      Warning(message);
      
      return -1;
   }
   
   if ( strcmp(item->internal.values.section,"_internal_") == 0 && lock != -1 )
   {
      char message[1024];
      snprintf(message,1024,"Internal configuration functions cannot be %s.",
      	 (lock==0)?"unlocked":"locked");
      Warning(message);
      return -1;
   }

   if ( lock == -1 )
   {
      char message[1024];
      if ( item->internal.values.section != NULL &&
      	   item->internal.values.name != NULL )
      {
	 snprintf(message, 1024, "%s:%s (%s) is %s.\n",
	    item->internal.values.section,
	    item->internal.values.name,
	    item->name,
	    (item->internal.locked==0)?"unlocked":"locked");
	 Output(message);
      }
      else
      {
	 snprintf(message, 1024, "%s is %s.\n",
	    item->name,
	    (item->internal.locked==0)?"unlocked":"locked");
	 Output(message);
      }
   }
   else
      item->internal.locked = lock;

   return rc;
}

/* ----------------------- f_lock_config ------------------------------ */

static int f_lock_config (const char *name, CONFIG_VALUES *dummy)
{
   return lock_unlock_status(name,1);
}

/* ----------------------- f_unlock_config ------------------------------ */

static int f_unlock_config (const char *name, CONFIG_VALUES *dummy)
{
   return lock_unlock_status(name,0);
}

/* ----------------------- f_limit_config ------------------------------ */

static int f_limit_config (const char *name, CONFIG_VALUES *dummy)
{
   return 0;
}

/* ----------------------- f_status_config ------------------------------ */

static int f_status_config (const char *name, CONFIG_VALUES *dummy)
{
   return lock_unlock_status(name,-1);
}

/* ----------------------- f_list_config ------------------------------ */

static int f_list_config (const char *name, CONFIG_VALUES *dummy)
{
   CONFIG_BLOCK *block;
   CONFIG_ITEM *itemlist, *item;
   int item_no;
   
   char *li_style = getenv("HCONFIG_LIST_STYLE");
   char *li_prefix = getenv("HCONFIG_LIST_PREFIX");
   int ltx_style = 0;
   if ( li_style != NULL )
      if ( strcmp(li_style,"latex") == 0 )
         ltx_style = 1;

   for ( block=(&first_config_block); block != (CONFIG_BLOCK *) NULL;
         block=block->next )
   {
      if ( (itemlist=block->items) == (CONFIG_ITEM *) NULL )
         break;
      for (item_no=0; itemlist[item_no].name != (char *) NULL; item_no++)
      {
         int show_it = 0;
         item = &itemlist[item_no];
         if ( item == (CONFIG_ITEM *) NULL )
            continue;
         if ( name[0] == '\0' || strcmp(name,"*") == 0 || strcmp(name,"all") == 0 )
            show_it = 1;
         else if ( strcmp(name,item->name) == 0 )
            show_it = 1;
         else if ( strcmp(name,"unlocked") == 0 && item->internal.locked == 0 )
            show_it = 1;
         else if ( strcmp(name,"locked") == 0 && item->internal.locked != 0 )
            show_it = 1;
         else if ( strcmp(name,"modified") == 0 && item->internal.values.nmod > 0 )
            show_it = 1;
         if ( ! show_it )
            continue;
#ifdef _REENTRANT
         pthread_mutex_lock(&item->internal.mutex_setval);
#endif
         if ( li_prefix != NULL )
         { 
            Output(li_prefix);
            if ( ltx_style )
               Output(" & ");
            else
               Output("\t");
         }
         Output(item->name);
         if ( ltx_style )
            Output(" & ");
         else
            Output(" = ");
         display_config_current(item);
         if ( ltx_style )
            Output(" \\\\\n");
         else
            Output("\n");
#ifdef _REENTRANT
         pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
      }
   }

   return 0;
}

/* ----------------------- f_get_config ------------------------------ */

static int f_get_config (const char *name, CONFIG_VALUES *dummy)
{
   return 0;
}

/* ----------------------- f_echo ------------------------------ */

static int f_echo (const char *name, CONFIG_VALUES *dummy)
{
   fprintf(stdout,"%s\n", name);
   return 0;
}

/* ----------------------- f_warning ------------------------------ */

static int f_warning (const char *name, CONFIG_VALUES *dummy)
{
   fprintf(stderr,"%s\n", name);
   return 0;
}

/* ----------------------- f_error ------------------------------ */

static int f_error (const char *name, CONFIG_VALUES *dummy)
{
   fprintf(stderr,"%s\n", name);
   return -1;
}

/* ----------------------- f_show_config ------------------------------ */
/**
 *  Display the current configuration status (internal usage only).
*/

static int f_show_config (const char *name, CONFIG_VALUES *val)
{
   CONFIG_BLOCK *block;
   CONFIG_ITEM *itemlist, *item;
   int item_no, i;
   char section[65];
   const char *cp;
   const char *org_name;

   if ( name == (char *) NULL )
      return -1;

   org_name = name;
   if ( (cp = strchr(name,':')) != (char *) NULL )
   {
      for (i=0; i<64 && name+i<cp; i++)
         section[i] = name[i];
      section[i] = '\0';
      name = cp+1;
   }
   else
      section[0] = '\0';

   Output(
"\n   Element name      Type    Elements  Lower b. Upper b.  Initial  Current\n");
   Output(
"---------------------------------------------------------------------------\n");
   if ( name[0] == '\0' || strcmp(name,"*") == 0 || strcmp(name,"all") == 0 )
   {
      for ( block=(&first_config_block); block != (CONFIG_BLOCK *) NULL;
            block=block->next )
      {
         if ( *section )
         {
            if ( block->section == (char *) NULL )
               continue;
            if ( strcmp(block->section,section) != 0 )
               continue;
         }
         if ( (itemlist=block->items) == (CONFIG_ITEM *) NULL )
            break;
         if ( block->section != (char *) NULL )
         {
	    char message[1024];
            snprintf(message,sizeof(message)-1,
	       "Section '%s':\n",block->section);
            Output(message);
         }
         else
            Output("Unnamed section:\n");
         for (item_no=0; itemlist[item_no].name != (char *) NULL; item_no++)
         {
            item = &itemlist[item_no];
            display_config_item(item);
         }
      }
   }
   else if ( (item = find_config_item(org_name)) != (CONFIG_ITEM *) NULL )
      display_config_item(item);
   else
   {
      char message[1024];
      snprintf(message,1024,"   %s    \t(Not found)\n",org_name);
      Output(message);
   }

   return 0;
}

/* --------------------- display_config_current ------------------------ */
/**
 *  Display current values of a single configuration item (internal usage only).
*/

static void display_config_current (CONFIG_ITEM *item)
{
   union ConfigDataPointer data;
   char message[10240];

   data.anything = item->data;

   if ( item->internal.itype > 0 && item->internal.itype != 50 )
   {
      if ( item->internal.itype != 30 )
      {
         int all_shown = 0, inum = 0;
         if ( item->size > 1 )
         {
            all_shown = 1;
            for (inum=1; inum<item->size; inum++)
            {
               switch ( item->internal.itype )
               {
               case 1:
                  if ( data.sdata[inum] != data.sdata[0] )
                     all_shown = 0;
                  break;
               case 2:
                  if ( data.idata[inum] != data.idata[0] )
                     all_shown = 0;
                  break;
               case 3:
                  if ( data.ldata[inum] != data.ldata[0] )
                     all_shown = 0;
                  break;
               case 4:
                  if ( data.cdata[inum] != data.cdata[0] )
                     all_shown = 0;
                  break;
               case 11:
               case 21:
                  if ( data.usdata[inum] != data.usdata[0] )
                     all_shown = 0;
                  break;
               case 12:
               case 22:
                  if ( data.uidata[inum] != data.uidata[0] )
                     all_shown = 0;
                  break;
               case 13:
               case 23:
                  if ( data.uldata[inum] != data.uldata[0] )
                     all_shown = 0;
                  break;
               case 14:
               case 24:
                  if ( data.ucdata[inum] != data.ucdata[0] )
                     all_shown = 0;
                  break;
               case 31:
                  if ( data.fdata[inum] != data.fdata[0] )
                     all_shown = 0;
                  break;
               case 32:
                  if ( data.ddata[inum] != data.ddata[0] )
                     all_shown = 0;
                  break;
               }
               if ( all_shown == 0 )
                  break;
            }
         }
         if ( all_shown )
            Output("all: ");
         for (inum=0; inum<item->size; inum++)
         {
            switch ( item->internal.itype )
            {
            case 1:
               sprintf(message,"%d ",data.sdata[inum]);
               break;
            case 2:
               sprintf(message,"%d ",data.idata[inum]);
               break;
            case 3:
               sprintf(message,"%ld ",data.ldata[inum]);
               break;
            case 4:
               sprintf(message,"%d ",data.cdata[inum]);
               break;
            case 11:
               sprintf(message,"%u ",data.usdata[inum]);
               break;
            case 12:
               sprintf(message,"%u ",data.uidata[inum]);
               break;
            case 13:
               sprintf(message,"%lu ",data.uldata[inum]);
               break;
            case 14:
               sprintf(message,"%u ",data.ucdata[inum]);
               break;
            case 21:
               sprintf(message,"%x ",data.usdata[inum]);
               break;
            case 22:
               sprintf(message,"%x ",data.uidata[inum]);
               break;
            case 23:
               sprintf(message,"%lx ",data.uldata[inum]);
               break;
            case 31:
               sprintf(message,"%g ",(double)data.fdata[inum]);
               break;
            case 32:
               sprintf(message,"%g ",data.ddata[inum]);
               break;
            }
            Output(message);
            if ( all_shown )
               break;
         }
      }
      else /* itype==30 */
         Output((char *)item->data);
   }
}

/* --------------------- display_config_item ------------------------ */
/**
 *  Display a single configuration item (internal usage only).
*/

static void display_config_item (CONFIG_ITEM *item)
{
   char message[1024];

   if ( item == (CONFIG_ITEM *) NULL )
      return;
#ifdef _REENTRANT
   pthread_mutex_lock(&item->internal.mutex_setval);
#endif
   if (strlen(item->name) <= 17 )
      sprintf(message,"  %c %-17s %-7s ",
      	 item->internal.locked?'#':' ',item->name,item->type);
   else
      sprintf(message,"  %c %s\n    %-17s %-7s ",
      	 item->internal.locked?'#':' ',item->name," ",item->type);
   Output(message);
   if ( item->internal.itype > 0 && item->internal.itype != 50 )
   {
      snprintf(message,sizeof(message),"%s%-8d %-8s %-8s ",
          /* If the data pointer was allocated, mark it with a '*' */
          (item->function != NULL)?"*":" ",
            item->size,
                 (item->lbound != (char *) NULL)?item->lbound:" ",
                      (item->ubound != (char *) NULL)?item->ubound:" ");
      if ( item->initial != (char *) NULL )
         if ( strlen(item->initial) > 8 )
            snprintf(message+strlen(message),sizeof(message)-strlen(message)-1,
	       "\n    %-17s Initial:  "," ");
      strcat(message,(item->flags & CFG_REJECT_MODIFICATION)?"*":" ");
      Output(message);
      if ( strlen((item->initial != (char *) NULL)?item->initial:" ") <= 8 )
      {
         sprintf(message,"%-8s ",
            (item->initial != (char *) NULL)?item->initial:" ");
         Output(message);
      }
      else
         Output((item->initial != (char *) NULL)?item->initial:" ");
      if ( item->size > 1 )
      {
         snprintf(message,sizeof(message),"\n    %-17s Current:   "," ");
         Output(message);
      }
      display_config_current(item);
   }
   else if ( item->internal.itype < 0 )
   {
      snprintf(message,sizeof(message),"%s %d [of %d bytes each]",
       (item->function != NULL)?"*":" ", item->size, 
       item->internal.values.elem_size);
      Output(message);
   }
   Output("\n");
#ifdef _REENTRANT
   pthread_mutex_unlock(&item->internal.mutex_setval);
#endif
}

/* ---------------------- signed_config_val ------------------------- */

static int signed_config_val (const char *name, const char *text, 
   const char *lbound, const char *ubound, int strict, long *ival)
{
   long value, lower, upper;

   lower = upper = 0;
   if ( text[0] == '0' && text[1] == 'x' )
   {
      if ( !is_hex_number(text+2) )
      {
	 config_syntax_error(name,"Value is not an integer");
	 return -1;
      }
      sscanf(text+2,"%lx",&value);
   }
   else if ( text[0] == '0' && text[1] == 'b' )
   {
      if ( !is_bin_number(text+2) )
      {
	 config_syntax_error(name,"Value is not an integer");
	 return -1;
      }
      value = decode_bin_number(text+2);
   }
   else if ( !is_signed_number(text) )
   {
      config_syntax_error(name,"Value is not an integer");
      return -1;
   }
   else
      value = atol(text);

   if ( lbound != (char *) NULL )
   {
      if ( !is_signed_number(lbound) )
      {
         config_syntax_error(name,"Lower bound is not an integer");
         return -1;
      }
      lower = atol(lbound);
      if ( value < lower )
      {
fprintf(stderr,"%s: integer lower bound violation for %ld >= %ld\n", name, value, lower);
         if ( strict ) 
         {
            config_syntax_error(name,"Violates lower bound");
            return -1;
         }
         value = lower;
         config_info(name,"Set to lower bound");
      }
   }

   if ( ubound != (char *) NULL )
   {
      if ( !is_signed_number(ubound) )
      {
         config_syntax_error(name,"Upper bound is not an integer");
         return -1;
      }
      upper = atol(ubound);
      if ( value > upper )
      {
fprintf(stderr,"%s: integer upper bound violation for %ld <= %ld\n", name, value, upper);
         if ( strict ) 
         {
            config_syntax_error(name,"Violates upper bound");
            return -1;
         }
         value = upper;

         config_info(name,"Set to upper bound");
      }
   }

   if ( lbound != (char *) NULL && ubound != (char *) NULL )
      if ( lower > upper )
      {
         config_syntax_error(name,"Upper bound is below lower bound");
         return -1;
      }

   *ival = value;
   return 0;
}

/* --------------------- unsigned_config_val ------------------------ */

static int unsigned_config_val (const char *name, const char *text, 
   const char *lbound, const char *ubound, int strict, unsigned long *uval)
{
   unsigned long value, lower, upper;

   value = 2;
   lower = 1;
   upper = 0;
   if ( text[0] == '0' && text[1] == 'x' )
   {
      if ( !is_hex_number(text+2) )
      {
	 config_syntax_error(name,"Value is not an unsigned integer");
	 return -1;
      }
      sscanf(text+2,"%lx",&value);
   }
   else if ( text[0] == '0' && text[1] == 'b' )
   {
      if ( !is_bin_number(text+2) )
      {
	 config_syntax_error(name,"Value is not an unsigned integer");
	 return -1;
      }
      value = decode_bin_number(text+2);
   }
   else if ( !is_unsigned_number(text) )
   {
      config_syntax_error(name,"Value is not an unsigned integer");
      return -1;
   }
   else
   {
#if ( defined(ANSI_C) && !defined(OS_LYNX) )
      sscanf(text,"%lu",&value);
#else
      sscanf(text,"%ld",&value);
#endif
   }

   if ( lbound != (char *) NULL )
   {
      if ( !is_unsigned_number(lbound) )
      {
         config_syntax_error(name,"Lower bound is not an unsigned integer");
         return -1;
      }
#if ( defined(ANSI_C) && !defined(OS_LYNX) )
      sscanf(lbound,"%lu",&lower);
#else
      sscanf(lbound,"%ld",&lower);
#endif
      if ( value < lower )
      {
fprintf(stderr,"%s: unsigned lower bound violation for %lu >= %lu\n", name, value, lower);
         if ( strict ) 
         {
            config_syntax_error(name,"Violates lower bound");
            return -1;
         }
         config_info(name,"Set to lower bound");
         value = lower;
         config_info(name,"Set to lower bound");
      }
   }

   if ( ubound != (char *) NULL )
   {
      if ( !is_unsigned_number(ubound) )
      {
         config_syntax_error(name,"Upper bound is not an unsigned integer");
         return -1;
      }
#if ( defined(ANSI_C) && !defined(OS_LYNX) )
      sscanf(ubound,"%lu",&upper);
#else
      sscanf(ubound,"%ld",&upper);
#endif
      if ( value > upper )
      {
fprintf(stderr,"%s: unsigned upper bound violation for %lu <= %lu\n", name, value, upper);
         if ( strict ) 
         {
            config_syntax_error(name,"Violates upper bound");
            return -1;
         }
         value = upper;
         config_info(name,"Set to upper bound");
      }
   }

   if ( lbound != (char *) NULL && ubound != (char *) NULL )
      if ( lower > upper )
      {
         config_syntax_error(name,"Upper bound is below lower bound");
         return -1;
      }

   *uval = value;
   return 0;
}

/* ----------------------- hex_config_val -------------------------- */

static int hex_config_val (const char *name, const char *text, 
   const char *lbound, const char *ubound, int strict, unsigned long *uval)
{
   unsigned long value, lower, upper;

   value = 2;
   lower = 1;
   upper = 0;
   if ( !is_hex_number(text) )
   {
      config_syntax_error(name,"Value is not a hexadecimal integer");
      return -1;
   }
   sscanf(text,"%lx",&value);
   if ( lbound != (char *) NULL )
   {
      if ( !is_hex_number(lbound) )
      {
         config_syntax_error(name,"Lower bound is not a hexadecimal integer");
         return -1;
      }
      sscanf(lbound,"%lx",&lower);
      if ( value < lower )
      {
         if ( strict ) 
         {
            config_syntax_error(name,"Violates lower bound");
            return -1;
         }
         value = lower;
         config_info(name,"Set to lower bound");
      }
   }
   if ( ubound != (char *) NULL )
   {
      if ( !is_hex_number(ubound) )
      {
         config_syntax_error(name,"Upper bound is not a hexadecimal integer");
         return -1;
      }
      sscanf(ubound,"%lx",&upper);
      if ( value > upper )
      {
         if ( strict ) 
         {
            config_syntax_error(name,"Violates upper bound");
            return -1;
         }
         value = upper;
         config_info(name,"Set to upper bound");
      }
   }

   if ( lbound != (char *) NULL && ubound != (char *) NULL )
      if ( lower > upper )
      {
         config_syntax_error(name,"Upper bound is below lower bound");
         return -1;
      }

   *uval = value;
   return 0;
}

/* ----------------------- real_config_val -------------------------- */

static int real_config_val (const char *name, const char *text, 
   const char *lbound, const char *ubound, int strict, double *rval)
{
   double value, lower, upper;

   lower = upper = 0;
   if ( !is_real_number(text) )
   {
      config_syntax_error(name,"Value is not a real number");
      return -1;
   }
   value = atof(text);
   if ( lbound != (char *) NULL )
   {
      if ( !is_real_number(lbound) )
      {
         config_syntax_error(name,"Lower bound is not a real number");
         return -1;
      }
      lower = atof(lbound);
      if ( value < lower )
      {
fprintf(stderr,"%s: lower bound violation for %f >= %f\n", name, value, lower);
         if ( strict ) 
         {
            config_syntax_error(name,"Violates lower bound");
            return -1;
         }
         value = lower;
         config_info(name,"Set to lower bound");
      }
   }
   if ( ubound != (char *) NULL )
   {
      if ( !is_real_number(ubound) )
      {
         config_syntax_error(name,"Upper bound is not a real number");
         return -1;
      }
      upper = atof(ubound);
      if ( value > upper )
      {
fprintf(stderr,"%s: upper bound violation for %f (%s) <= %f (%s)\n", name, value, text, upper, ubound);
         if ( strict ) 
         {
            config_syntax_error(name,"Violates upper bound");
            return -1;
         }
         value = upper;
         config_info(name,"Set to upper bound");
      }
   }

   if ( lbound != (char *) NULL && ubound != (char *) NULL )
      if ( lower > upper )
      {
         config_syntax_error(name,"Upper bound is below lower bound");
         return -1;
      }

   *rval = value;
   return 0;
}

/* ----------------------- is_signed_number ------------------------- */

int is_signed_number (const char *text)
{
   int i;

   if ( text[0] != '+' && text[0] != '-' && !isdigit((int)text[0]) )
      return 0;
   for ( i=1; text[i] != '\0'; i++ )
      if ( !isdigit((int)text[i]) )
         return 0;
   return 1;
}

/* ---------------------- is_unsigned_number ------------------------ */

int is_unsigned_number (const char *text)
{
   int i;

   if ( !isdigit((int)text[0]) )
      return 0;
   for ( i=1; text[i] != '\0'; i++ )
      if ( !isdigit((int)text[i]) )
         return 0;
   return 1;
}

/* ------------------------ is_hex_number -------------------------- */

int is_hex_number (const char *text)
{
   int i;

   if ( !isxdigit((int)text[0]) )
      return 0;
   for ( i=1; text[i] != '\0'; i++ )
      if ( !isxdigit((int)text[i]) )
         return 0;
   return 1;
}

/* ------------------------ is_bin_number -------------------------- */

int is_bin_number (const char *text)
{
   int i;

   if ( text[0] != '0' && text[0] != '1' )
      return 0;
   for ( i=1; text[i] != '\0'; i++ )
      if ( text[i] != '0' && text[i] != '1' )
         return 0;
   return 1;
}

/* ---------------------- decode_bin_number ------------------------- */

unsigned long decode_bin_number (const char *text)
{
   unsigned int i, k;
   unsigned long l = 0;

   for ( i=0; text[i] != '\0'; i++ )
   {
      if ( text[i] == '1' )
      	 k = 1;
      else if ( text[i] == '0' )
      	 k = 0;
      else
      	 break;
      l = (l<<1) + k;
   }
   return l;
}

/* ------------------------ is_real_number -------------------------- */

int is_real_number (const char *text)
{
   int i, dot;

   dot = 0;
   if ( text[0] == '.' )
      dot = 1;
   else if ( text[0] != '+' && text[0] != '-' && !isdigit((int)text[0]) )
      return 0;
   for ( i=1; text[i] != '\0' && toupper((int)text[i]) != 'E'; i++ )
   {
      if ( text[i] == '.' )
      {
         if ( dot )
            return 0;
         else
            dot = 1;
      }
      else if ( !isdigit((int)text[i]) )
         return 0;
   }
   if ( toupper((int)text[i]) == 'E' )
      return(is_signed_number(text+i+1));
   return 1;
}

static char cfg_fname[1024];
static char preprocessor[4096] = ""; /* No default preprocessor! */
static char **cfg_stack;

#ifdef OS_UNIX
# define TMP_FORMAT "/tmp/cfg%d.tmp"
#else
# ifdef OS_OS9
#  define TMP_FORMAT "/r0/cfg%d.tmp"
# else
#  define TMP_FORMAT "cfg%d.tmp"
# endif
#endif

/* -------------------- set_config_filename ------------------ */
/**
    Set the name of the configuration file to be read by the
    function read_config_lines().

    @param  fname  Name of file to be used.

    @return (none)
*/

void set_config_filename (const char *fname)
{
#ifdef _REENTRANT
   pthread_mutex_lock(&mlock_config);
#endif
   strncpy(cfg_fname,fname,sizeof(cfg_fname)-1);
#ifdef _REENTRANT
   pthread_mutex_unlock(&mlock_config);
#endif
}

/* --------------------- get_config_filename ----------------- */
/**
    Return the current value of the configuration file name.

    @param -- (none)

    @return pointer to static file name string

 */

char *get_config_filename()
{
   char *val;
#ifdef _REENTRANT
   pthread_mutex_lock(&mlock_config);
#endif
   val = cfg_fname;
#ifdef _REENTRANT
   pthread_mutex_unlock(&mlock_config);
#endif
   return val;
}

/* -------------------- set_config_preprocessor ------------------ */
/**
    Set the command name and options of a preprocessor for
    configuration files to be read by function read_config_lines().
    The input and output file names will be appended to the
    command string set by this function.

    @param  preproc  Command string

    @return  (none)

 */

void set_config_preprocessor (char *preproc)
{
#ifdef _REENTRANT
   pthread_mutex_lock(&mlock_config);
#endif
   strncpy(preprocessor,preproc,sizeof(preprocessor)-1);
#ifdef _REENTRANT
   pthread_mutex_unlock(&mlock_config);
#endif
}

/* --------------------- get_config_preprocessor ----------------- */
/**
    Return the current value of the configuration preprocessor.

    @param --  (none)

    @return  pointer to static command string

 */

char *get_config_preprocessor()
{
   char *pp;
#ifdef _REENTRANT
   pthread_mutex_lock(&mlock_config);
#endif
   pp = preprocessor;
#ifdef _REENTRANT
   pthread_mutex_unlock(&mlock_config);
#endif
   return pp;
}

static int read_status;

/* ---------------------- set_config_stack ------------------- */
/**
    Set a list of configuration lines to be processed before
    any lines from a file are read by read_config_lines().

    @param  stack Pointer to NULL terminated vector of strings.

    @return  (none)

 */

void set_config_stack (char **stack)
{
#ifdef _REENTRANT
   pthread_mutex_lock(&mlock_config);
#endif
   cfg_stack = stack;
#ifdef _REENTRANT
   pthread_mutex_unlock(&mlock_config);
#endif
}

/* --------------------- read_config_lines ------------------ */
/**
    Read configuration data from a file and return it line
    by line to the calling function (one line per call).
    A NULL pointer is returned on end-of-file. This function
    is intended to be used as the usual 'fptr' argument for
    init_config().

    @param --  (none)

    @return  Pointer to character string or NULL.

 */

char *read_config_lines()
{
   static char tmp_fname[145];
   static char line[10241];
   static char *next_line, *this_line;
   int do_read;
   static FILE *cfg_file;

#ifdef _REENTRANT
   pthread_mutex_lock(&mlock_config);
#endif

   read_status = 0;

   if ( cfg_stack != (char **) NULL )
      if ( *cfg_stack != (char *) NULL )
      {
         strncpy(line,*cfg_stack,sizeof(line)-1);
         cfg_stack++;
#ifdef _REENTRANT
   pthread_mutex_unlock(&mlock_config);
#endif
         return(line);
      }

   /* Open the file if not open yet */
   if ( cfg_file == (FILE *) NULL )
   {
      char command[8000];
      char *fn;

      fn = cfg_fname;
      if ( *fn != '\0' )
      {
	 char message[1060];
         snprintf(message,sizeof(message),"Configuration file is '%s'.",fn);
         Information(message);
      }
      if ( *preprocessor != '\0' && *fn != '\0' )
      {
         if ( (cfg_file = fopen(fn,READ_TEXT)) == (FILE *) NULL )
         {
	    char message[1060];
            snprintf(message,sizeof(message),
	       "Cannot open configuration file '%s'",fn);
            Warning(message);
            Warning("Using default configuration");
            read_status = -2;
            goto null_return;
         }
         fclose(cfg_file);
         if ( strlen(preprocessor) + strlen(TMP_FORMAT) + 10 +
              strlen(cfg_fname) > sizeof(command)-6 )
            Warning("Preprocessor command string is too long.");
         else if ( strlen(TMP_FORMAT) > 128 )
            Warning("Preprocessor temporary file name is too long.");
         else
         {
	    char message[4200];
            sprintf(tmp_fname,TMP_FORMAT,getpid());
            snprintf(command,sizeof(command),"%s %s %s",preprocessor,
               cfg_fname,tmp_fname);
            snprintf(message,sizeof(message),
	       "Preprocessor is '%s'.",preprocessor);
            Information(message);
            if ( system(command) >= 0 )
               fn = tmp_fname;
            else
            {
               Warning("Configuration file preprocessing failed.");
               Warning("Using default configuration");
               read_status = -3;
               goto null_return;
            }
         }
      }
      if ( fn[0] == '\0' )
      {
         Warning("No configuration filename available");
         Warning("Using default configuration");
         read_status = -1;
         goto null_return;
      }
      if ( (cfg_file = fopen(fn,READ_TEXT)) == (FILE *) NULL )
      {
	 char message[1060];
         snprintf(message,sizeof(message),
	    "Cannot open configuration file '%s'",fn);
         Warning(message);
         Warning("Using default configuration");
         read_status = -2;
         goto null_return;
      }
   }
   do_read = 1;
   if ( next_line != (char *) NULL )
      if ( *next_line != '\0' )
         do_read = 0;
   if ( do_read )
   {
      fgets(line,sizeof(line)-1,cfg_file);
      if ( feof(cfg_file) || ferror(cfg_file) )
      {
         if ( ferror(cfg_file) )
         {
            Warning("Error reading configuration file");
            read_status = -4;
         }
         fclose(cfg_file);
         cfg_file = (FILE *) NULL;
         if ( *preprocessor != '\0' && *tmp_fname != '\0' )
            unlink(tmp_fname);
         goto null_return;
      }
      next_line = line;
   }
   this_line = next_line;
   /* Multi-lines separated by semicolons need special attention for interference with comments */
   if ( (next_line = strchr(this_line,';')) != (char *) NULL )
   {
      char *cmt = strchr(this_line,'%'), *s;
      int in_str = 0;
      if ( cmt != (char *) NULL && cmt < next_line )
      {
         Warning("Beware of side-effects with lines containing both percent and semicolon characters.");
         for ( s=this_line; *s != '\0'; s++ )
         {
            if ( *s == '\\' ) /* skip backslash-escaped next character */
            {
               s++;
               continue;
            }
            if ( *s == '\'' ) /* Start or end of string enclosed in single quotes */
            {
               if ( in_str == 0 )
                  in_str = 1;
               else if ( in_str == 1 )
                  in_str = 0;
            }
            if ( *s == '"' ) /* Start or end of string enclosed in double quotes */
            {
               if ( in_str == 0 )
                  in_str = 2;
               else if ( in_str == 2 )
                  in_str = 0;
            }
            if ( !in_str )
            {
               if ( *s == '%' ) /* Semicolon is in the string and will not split the line. */
               {
                  next_line = NULL;
                  break; /* Use the line as-is, not breaking it. */
               }
               else if ( *s == ';' ) /* Preceding percent seems to be escaped or in string */
               {
                  /* Note that other parts of hconfig might not respect that the percent
                     was enclosed in quotes and will consider it starting a comment anyway. */
                  *s = '\0';
                  next_line = s+1;
                  break;
               }
            }
         }
      }
      else
      {
         /* No %/; interference problem */
         *next_line = '\0';
         next_line++;
      }
   }

#ifdef _REENTRANT
   pthread_mutex_unlock(&mlock_config);
#endif

   return(this_line);
   
null_return:
#ifdef _REENTRANT
   pthread_mutex_unlock(&mlock_config);
#endif
   return((char *) NULL);
}

/* ----------------------- read_config_status ----------------------- */
/**
    Return the status of reading a configuration file with
    read_config_lines() in a preceding call to init_config().

    @param --  (none)

    @return         0 (o.k.),
		   -1 (no config file set),
		   -2 (config file open failed),
		   -3 (preprocessing failed),
		   -4 (read error).

 */

int read_config_status ()
{
#ifdef _REENTRANT
   pthread_mutex_lock(&mlock_config);
#endif
   return read_status;
#ifdef _REENTRANT
   pthread_mutex_unlock(&mlock_config);
#endif
}

/* ---------------------- config_syntax_error ----------------------- */

static void config_syntax_error (const char *name, const char *text)
{
   char message[1024];
   snprintf(message,sizeof(message),
      "Syntax error for configuration item '%s': %s",
      name, text);
   Warning(message);
}

/* ------------------------ config_info ---------------------------- */

static void config_info (const char *name, const char *text)
{
   char message[1024];
   snprintf(message,sizeof(message),
      "Configuration item '%s': %s", name, text);
   Information(message);
}

/* -------------------- config_binary_interface -------------------- */
/** 
    @short Define a binary interface for an I/O type.
 */
 
int define_config_binary_interface (int item_type, size_t elem_size,
   void *(*new_func)(int nelem, int item_type),
   int (*delete_func)(void *ptr, int nelem, int item_type),
   int (*read_func)(void *bin_item, IO_BUFFER *iobuf, int item_type),
   int (*write_func)(void *bin_item, IO_BUFFER *iobuf, int item_type),
   int (*readtext_func)(void *bin_item, char *text, int item_type),
   int (*list_func)(void *bin_item, int item_type),
   int (*copy_func)(void *bin_item_to, void *bin_item_from, int io_type))
{
   struct Binary_Interface_Chain *bc, *bc_last;
   struct Config_Binary_Item_Interface *ife;
   
   if ( item_type < 0 )
      item_type = -item_type;

   if ( readtext_func != NULL )
   {
      Warning("No text reading functions supported yet for binary-only items.");
   }

   for ( bc=bc_last=bin_chain_root; bc != NULL; bc++ )
   {
      if ( bc->interface == NULL )
      	 continue;
      if ( bc->interface->io_item_type == item_type )
      {
      	 char message[1024];
	 snprintf(message,sizeof(message),
	    "Item type %d already had a binary interface defined.\n",
	    item_type);
	 break;
      }
      bc_last = bc;
   }
   
   /* Existing interface to be replaced ? */
   if ( bc != NULL )
      ife = bc->interface;
   else if ( (ife = (struct Config_Binary_Item_Interface *) 
          calloc(1,sizeof(struct Config_Binary_Item_Interface))) == NULL ||
        (bc = (struct Binary_Interface_Chain *) 
          calloc(1,sizeof(struct Binary_Interface_Chain))) == NULL )
   {
      Error("Memory allocation failed");
      return -1;
   }
   else
   {
      bc->interface = ife;
      if ( bin_chain_root == NULL )
      	 bin_chain_root = bc;
      if ( bc_last != NULL )
      	 bc_last->next = bc;
   }
   
   ife->io_item_type = item_type;
   ife->elem_size = elem_size;
   ife->new_func = new_func;
   ife->delete_func = delete_func;
   ife->read_func = read_func;
   ife->write_func = write_func;
   ife->readtext_func = readtext_func;
   ife->list_func = list_func;
   ife->copy_func = copy_func;

   return 0;
}

/* ------------------ find_binary_interface ---------------------- */
/** 
    @short Find the matching binary interface for given item type.
 */

struct Config_Binary_Item_Interface *
   find_config_binary_interface (int item_type)
{
   struct Binary_Interface_Chain *bc, *bc_last;

   if ( item_type < 0 )
      item_type = -item_type;

   for ( bc=bc_last=bin_chain_root; bc != NULL; bc++ )
   {
      if ( bc->interface == NULL )
      	 continue;
      if ( bc->interface->io_item_type == item_type )
      	 return bc->interface;
   }

   return NULL;
}
