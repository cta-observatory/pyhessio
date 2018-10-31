/* ============================================================================

   Copyright (C) 1993, 2001, 2007, 2010  Konrad Bernloehr

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

/** @file hconfig.h
 *  @short Declare hconfig structures and functions
 *  @author  Konrad Bernloehr 
 *  @date    CVS $Date: 2017/08/02 13:29:04 $
 *  @version CVS $Revision: 1.8 $
 */

#ifndef HCONFIG_H__LOADED

#define HCONFIG_H__LOADED 1

#ifndef INITIAL_H__LOADED
# define NO_INITIAL_MACROS 1
# include "initial.h"
#endif
#ifndef IO_BASIC_H__LOADED
# include "io_basic.h"
#endif

#ifdef _REENTRANT
#include <pthread.h>
#endif

/* Two macros for expanding something into a string: */
#if ( defined __STDC__ && __STDC__ ) || defined __cplusplus
#define _XSTR_(s) _STR_(s) /**< Expand a macro first and then enclose in string */
#define _STR_(s) #s /**< Enclose in string without macro expansion. */
#else
#define _XSTR_(s) === /* No working conversion method known: let compilation fail. */
#define _STR_(s) "s" /**< This one also works with traditional C. */
#endif

#ifdef __cplusplus
// namespace hconfig {
extern "C" {
#endif

#ifndef CONST
#if defined (ANSI_C) || defined (__STDC__) || defined (__cplusplus)
# define CONST const
#else 
# define CONST
#endif
#endif

#define IO_TYPE_HCONFIG_ENVELOPE 900
#define IO_TYPE_HCONFIG_NAME     901
#define IO_TYPE_HCONFIG_TEXT     902
#define IO_TYPE_HCONFIG_INDEX    903
#define IO_TYPE_HCONFIG_NUMBERS  904

/** This union of pointers allows convenient access of various types of data. */

union ConfigDataPointer
{
   void *anything;
   char *cdata;
   unsigned char *ucdata;
   short *sdata;
   unsigned short *usdata;
   int *idata;
   unsigned int *uidata;
   long *ldata;
   unsigned long *uldata;
   float *fdata;
   double *ddata;
};

/** Configuration value may have optional lower and/or upper bounds. */

union ConfigBoundary
{
   long lval;
   unsigned long ulval;
   double *rval;
};

typedef void *(*PFVP)(char *, char *, int);
typedef int (*PFISI)(char *, int);
typedef int (*PFITI)(const char *, int);
typedef int (*PFISS)(char *, char *);

/** Configuration values and supporting data passed to user functions. */

struct ConfigValues
{
   void *data_changed;  /**< Pointer to the updated values. */
   void *data_saved;    /**< Pointer to the saved values. */
   int max_mod;         /**< How many elements can, at most, be modified. */
   int nmod;            /**< How many have been modified. */
   int *list_mod;       /**< List of indices to modified elements. */
   unsigned char *mod_flag; /**< Vector of size max_mod indicating modified elements. */
   int itype;           /**< Internal item type representation. */
   const char *name;    /**< The name of the element. */
   const char *section; /**< The section to which it belongs. */
   int elements;        /**< The number of elements it has. */
   int elem_size;       /**< The size of one element in bytes. */
   int binary_config;   /**< Set to one if binary configuration was used. */
};
typedef struct ConfigValues CONFIG_VALUES;

typedef int (*PFIX)(const char *name, CONFIG_VALUES *val);

/**
 *  @short Interface definitions for binary-only items.
 *
 *  Binary-only items are structures, classes, or unions which
 *  can only be filled via dedicated functions (methods) and not
 *  via the standard text-input.
 *
 *  This structure defines available interface methods.
 *  The item type is always passed to the functions, in case that
 *  a function can handle more than one type.
 */

struct Config_Binary_Item_Interface
{
   /** The eventio item type. */
   int io_item_type;
   /** The size of the elements. */
   int elem_size;
   /** The function to be called for allocating elements. */
   void *(*new_func)(int nelem, int item_type);
   /** The function to be called for deleting elements. */
   int (*delete_func)(void *ptr, int nelem, int item_type);
   /** The function to be called for reading elements from buffer. */
   int (*read_func)(void *bin_item, IO_BUFFER *iobuf, int item_type);
   /** The function to be called for writing elements to buffer. */
   int (*write_func)(void *bin_item, IO_BUFFER *iobuf, int item_type);
   /** The function to be called for reading elements from text line. */
   int (*readtext_func)(void *bin_item, char *text, int item_type);
   /** The optional function for listing element contents. */
   int (*list_func)(void *bin_item, int item_type);
   /** The optional function for copying elements.
    *  This is only needed if the element includes pointers to
    *  external or dynamically allocated material.
   */
   int (*copy_func)(void *bin_item_to, void *bin_item_from, int io_type);
};

/** Configuration elements used only internally. */

struct ConfigIntern
{
   int itype;        	/**< Parameter type code. */
   int elem_size;     	/**< Size of elements in bytes. */
   int locked;       	/**< Set to 1 if locked. */
   int bound;        	/**< Bits 0-3 set if lower soft, upper soft,  */
                     	/**<   lower hard, or upper hard bound present. */
   union ConfigBoundary lbound_soft; /**< Used for checking new values. */
   union ConfigBoundary ubound_soft; /**< Used for checking new values. */
   union ConfigBoundary lbound_hard; /**< Used for checking new values. */
   union ConfigBoundary ubound_hard; /**< Used for checking new values. */
   struct ConfigValues values;       /**< Passed to user function. */
   struct Config_Binary_Item_Interface *bin_interface;
   int bin_alloc_elements;
#ifdef _REENTRANT
   pthread_mutex_t mutex_setval;     /**< Mutex to control access. */
   pthread_mutex_t mutex_item;       /**< Mutex to control access. */
#endif
};

/** Configuration as used in definitions of configuration blocks. */

struct ConfigItemStruct
{
   const char *name;     /**< Parameter/function name. */
   const char *type;     /**< Data/function type. */
   int size;             /**< Number of elements. */
   void *data;           /**< Data pointer or NULL. */
   PFIX function;        /**< Associated function or NULL. */
   const char *initial;  /**< Initial values/argument or NULL. */
   const char *lbound;   /**< Lower bound (soft,hard) on values or NULL. */
   const char *ubound;   /**< Upper bound (soft,hard) on values or NULL. */
   int flags;            /**< Additional flag bits. */
   PFISS validate;       /**< Function to validate if change is possible or NULL. */
#ifdef _REENTRANT
   pthread_mutex_t *mutex_pointer; /**< Mutex to control concurrent access or NULL. */
#else
   void *res1;           /**< Placeholder to keep structure size the same. */
#endif
   void *res2;           /**< Not used. */
   struct ConfigIntern internal; /**< Internal data. */
};
typedef struct ConfigItemStruct CONFIG_ITEM;

/* You may set any of the following bits to flags: */
#define CFG_REQUIRE_DATA          1
#define CFG_REQUIRE_ALL_DATA      2
#define CFG_REJECT_MODIFICATION   4
#define CFG_HARD_BOUND            8
#define CFG_STRICT_BOUND         16
/* Do not set any of those bits: */
#define CFG_INITIALIZED          32
#define CFG_ALL_INITIALIZED      64
#define CFG_NOT_INITIAL         128

# define NULL_CONFIG_ITEM \
 (char *) NULL, (char *) NULL, 0, NULL, NULL, (char *) NULL, (char *) NULL, (char *) NULL, 0, NULL, NULL, NULL, {0}

/**
 * Mutexes are only inserted when pthreads are used.
 * In the multi-threaded variant: the address of the given mutex.
 * In the single-threaded variant: a null pointer.
*/
#ifdef _REENTRANT
#define CFG_MUTEX(mutex) (&mutex)
#else
#define CFG_MUTEX(mutex) (NULL)
#endif

/* config.c */
int build_config (CONFIG_ITEM *items, const char *section);
int init_config (char *(*fptr) (void));
void unhook_internal(void);
void rehook_internal(void);
int reload_config (char *(*fptr)(void));
void *config_alloc_data (char *name, char *type, int size);
int reconfig (char *text);
int verify_config_section (char *section);
int set_config_history (PFITI fptr);
void set_config_filename (const char *fname);
char *get_config_filename (void);
void set_config_preprocessor (char *preproc);
char *get_config_preprocessor (void);
void set_config_stack (char **stack);
char *read_config_lines (void);
int read_config_status (void);
CONFIG_ITEM *find_config_item (const char *name);
int define_config_binary_interface (int item_type, size_t elem_size,
   void *(*new_func)(int nelem, int item_type),
   int (*delete_func)(void *ptr, int nelem, int item_type),
   int (*read_func)(void *bin_item, IO_BUFFER *iobuf, int item_type),
   int (*write_func)(void *bin_item, IO_BUFFER *iobuf, int item_type),
   int (*readtext_func)(void *bin_item, char *text, int item_type),
   int (*list_func)(void *bin_item, int item_type),
   int (*copy_func)(void *bin_item_to, void *bin_item_from, int io_type));
struct Config_Binary_Item_Interface *
   find_config_binary_interface (int item_type);
int reconfig_binary (char *buffer, size_t buflen);
int config_binary_read_text (IO_BUFFER *iobuf, char *name, int maxlen);

int is_signed_number (const char *text);
int is_unsigned_number (const char *text);
int is_hex_number (const char *text);
int is_bin_number (const char *text);
int is_real_number (const char *text);
unsigned long decode_bin_number (const char *text);

/* straux.c */
int abbrev (CONST char *s, CONST char *t);
int getword (CONST char *s, int *spos, char *word,
   int maxlen, char blank, char endchar);
int config_binary_read_index (IO_BUFFER *iobuf, int *nidx, 
   int *idx_low, int *idx_high, int max_idx);

/* io_hconfig.c */
int config_binary_write_name(IO_BUFFER *iobuf, char *name);
int config_binary_write_text(IO_BUFFER *iobuf, char *text);
int config_binary_read_text(IO_BUFFER *iobuf, char *name, int maxlen);
int config_binary_text_length(IO_BUFFER *iobuf);
int config_binary_read_name(IO_BUFFER *iobuf, char *name, int maxlen);
int config_binary_write_index(IO_BUFFER *iobuf, int nidx, int *idx_low, int *idx_high);
int config_binary_read_index(IO_BUFFER *iobuf, int *nidx, int *idx_low, int *idx_high, int max_idx);
int config_binary_envelope_begin(IO_BUFFER *iobuf, IO_ITEM_HEADER *item_header);
int config_binary_envelope_end(IO_BUFFER *iobuf, IO_ITEM_HEADER *item_header);
int config_binary_inquire_numbers(IO_BUFFER *iobuf, int *ntype, int *nsize, int32_t *num, int *nopt);int config_binary_read_numbers(IO_BUFFER *iobuf, void *data, size_t max_size);
int config_binary_convert_data(void *out, int out_type, int out_size, void *in, int in_type, int in_size);

#ifdef __cplusplus
} /* extern "C" */
// } /* namespace hconfig */
#endif

#endif
