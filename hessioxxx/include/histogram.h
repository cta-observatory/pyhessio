/* ============================================================================

   Copyright (C) 1991, 2001, 2007, 2010  Konrad Bernloehr

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

/** @file histogram.h
 *  @short Declarations for handling one- and two-dimensional histograms.
 *
 *  The functions to work with these histograms is found in histogram.c .
 *  Eventio routines are available in
 *  io_histogram.c and conversion to HBOOK format is available
 *  through the 'cvt2' program. Handling of moments of a 1-D distribution
 *  is implemented in moments.c .
 *
 *  @author  Konrad Bernloehr 
 *  @date    1991 - 2010
 *  @date    CVS $Date: 2013/10/21 12:53:31 $
 *  @version CVS $Revision: 1.12 $
 */

#ifndef HISTOGRAM_H__LOADED            /* Ignore if included a second time */

#define HISTOGRAM_H__LOADED 1

#ifndef INITIAL_H__LOADED
#include "initial.h"
#endif

#ifdef _REENTRANT
#include <pthread.h>
#ifndef PTHREAD_ONCE_INIT
#define PTHREAD_ONCE_INIT pthread_once_init
#define PTHREAD_MUTEX_INITIALIZER 0
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
@verbatim
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   For compatibility reasons the following 'typedef's are kept, but   
   the defined types should not be used any more because all of them  
   were changed in histogram.c to 'long', 'double', etc.	      
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
@endverbatim

   HISTVALUE may be either an 'integer' type (recommended: long int)
   or a 'real' type (recommended: double). The method of calculating
   the array index corresponding to a given value is somewhat different 
   for these two alternatives. Using a float for the 'real' type 
   instead of a double would make no difference. However, a short int
   or an unsigned short int as 'integer' type  requires more care
   for the calculation of the array index compared to a long or a
   unsigned long (frequent overflows unless a type cast of intermediate 
   values to a long type is used).
*/

typedef double HISTVALUE_REAL;    /**< May be 'float' for ANSI C compiler */
typedef long HISTVALUE_INT;       /**< Short int is not recommended */

/**
 *  The histogram counts may be unsigned short or unsigned long.
 *  With a unsigned short the overflow of a bin might easily happen.
*/

typedef unsigned long  HISTCOUNT;
#ifdef ANSI_C
#define MAX_HISTCOUNT 4294967295UL  /* or ULONG_MAX from <limits.h> */
#else
#define MAX_HISTCOUNT 4294967295    /* or ULONG_MAX from <limits.h> */
#endif

/**
 *  To avoid loss of precision for adding many numbers, sums are of
 *  double type if 'real' type HISTVALUEs are used. 
*/

typedef double HISTSUM_REAL;
typedef long HISTSUM_INT;

typedef double HISTSTATVALUE;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/** Parameters defining the usable range of coordinates. */

union Histogram_Parameters
{
   /** Histogram parameters if it is some sort of 'F' or 'D' type. */
   struct
   {
      double lower_limit;  /**< Lower limit of histogram range */
      double upper_limit;  /**< Upper limit of histogram range */
      double sum;          /**< Sum of all values              */
      double tsum;         /**< Sum of values within range     */
      double inverse_binwidth;    /**< 1./(width_of_one_bin)   */
   } real;                 /**< Needed for real-type limits.   */
   /** Histogram parameters if it is some sort of 'I' (int) type. */
   struct
   {
      long lower_limit;    /**< Lower limit of histogram range */
      long upper_limit;    /**< Upper limit of histogram range */
      long sum;            /**< Sum of all values              */
      long tsum;           /**< Sum of values within range     */
      long width;          /**< Width of histogram range       */
   } integer;              /**< Needed for integer-type limits.*/
};

/** A histogram extension only allocated for weighted histograms. */

struct Histogram_Extension
{
   double content_all;           /**< Sum of all contents */
   double content_inside;        /**< Sum of contents within range */
   double content_outside[8];    /**< Contents outside range */
   float *fdata;                 /**< Data of each bin (ix+nx*iy) */
   double *ddata;                /**< in one of two precisions. */
};

/** A complete 1-D or 2-D histogram with control and data elements */

struct histogram
{
   char *title;                  /**< Histogram title (optional)     */
   long ident;                   /**< Histogram ID number (optional) */
   union Histogram_Parameters specific;
   union Histogram_Parameters specific_2d;
   int nbins;                    /**< Number of histogram bins       */
   int nbins_2d;                 /**< Same for 2nd coordinate of 2-D */
   unsigned long entries;        /**< No. of entries, incl. u.f./o.f.*/
   unsigned long tentries;       /**< No. of entries, without  """   */
   unsigned long underflow;      /**< No. of entries below range     */
   unsigned long underflow_2d;   /**< Same in 2nd coord of 2-D histo.*/
   unsigned long overflow;       /**< No. of entries above range     */
   unsigned long overflow_2d;    /**< Same in 2nd coord of 2-D histo.*/
   unsigned long *counts;        /**< Pointer to histogram data      */
   char type;                    /**< 'I' for integer histogram,     */
                                 /**< 'i' for int. lookup table,     */
                                 /**< 'R' for floating point histogr.*/
                                 /**< 'r' for fl. p. lookup table,   */
                                 /**< 'F'/'D' for single/double pre- */
                                 /**< cision weighted histograms.    */
   struct histogram *previous;   /**< References to neighbours in    */
   struct histogram *next;       /**< linked list of histograms.     */
   struct Histogram_Extension *extension;
                                 /**< Extension for weighted histos  */
#ifdef _REENTRANT
   pthread_mutex_t mlock_this;   /**< Mutex for locking concurrent access */
#endif
};

typedef struct histogram HISTOGRAM;

/** Statistics element for histogram analysis */

struct histstat
{
   /* For all entries */
   double mean, mean_2d;
   /* For entries in histogram range */
   double tmean, tmean_2d;
   /* Computed from histrogram entries */
   double hmean, hmean_2d;
   double sigma, sigma_2d;
   double median, median_2d;
};

/** First, second, and higher moments of a 1-D histogram. */

struct momstat
{
   /* For all entries */
   double mean;
   double sigma;
   double skewness;
   double kurtosis;
   /* For entries within range */
   double tmean;
   double tsigma;
   double tskewness;
   double tkurtosis;
};

/** Numbers to be summed up to obtain the moments */

struct moments
{
   double lower_limit, upper_limit;
   double sum,  tsum;
   double sum2, tsum2;
   double sum3, tsum3;
   double sum4, tsum4;
   unsigned long entries, tentries;
   int level;
};

typedef struct moments MOMENTS;

/* ----- Functions in histogram.c ----- */

void histogram_lock (HISTOGRAM *histo);
void histogram_unlock (HISTOGRAM *histo);
HISTOGRAM *get_first_histogram (void);
void set_first_histogram (HISTOGRAM * new_first_histogram);
HISTOGRAM *get_histogram_by_ident (long ident);
void list_histograms (long ident);
HISTOGRAM *book_histogram (long id, const char *title, const char *type,
   int dimension, double *low, double *high, int *nbins);
HISTOGRAM *book_int_histogram (long id, const char *title, int dimension, 
   long *low, long *high, int *nbins);
HISTOGRAM *book_1d_histogram (long id, const char *title, const char *type,
   double low, double high, int nbins);
HISTOGRAM *allocate_histogram (const char *type, int dimension,
   double *low, double *high, int *nbins);
HISTOGRAM *alloc_int_histogram (long low,
      long high, int nbins);
HISTOGRAM *alloc_real_histogram (double low,
      double high, int nbins);
HISTOGRAM *alloc_2d_int_histogram (long xlow,
      long xhigh, int nxbins, long ylow,
      long yhigh, int nybins);
HISTOGRAM *alloc_2d_real_histogram (double xlow,
      double xhigh, int nxbins, double ylow,
      double yhigh, int nybins);
void describe_histogram (HISTOGRAM *histo, const char *title, long ident);
void clear_histogram (HISTOGRAM *histo);
void free_histogram (HISTOGRAM *histo);
void free_all_histograms (void);
void unlink_histogram (HISTOGRAM *histo);
int fill_int_histogram (HISTOGRAM *histo, long value);
int fill_real_histogram (HISTOGRAM *histo, double value);
int fill_weighted_histogram (HISTOGRAM *histo, double value,
   double weight);
int fill_2d_int_histogram (HISTOGRAM *histo, long xvalue,
      long yvalue);
int fill_2d_real_histogram (HISTOGRAM *histo, double xvalue,
      double yvalue);
int fill_2d_weighted_histogram (HISTOGRAM *histo, double xvalue,
      double yvalue, double weight);
int fill_histogram (HISTOGRAM *histo, double xvalue,
      double yvalue, double weight);
int fill_histogram_by_ident (long id, double xvalue,
      double yvalue, double weight);
int stat_histogram (HISTOGRAM *histo, struct histstat *stbuf);
double locate_histogram_fraction (HISTOGRAM *histo, double fraction);
int fast_stat_histogram (HISTOGRAM *histo, struct histstat *stbuf);
int histogram_matching (HISTOGRAM *histo1, HISTOGRAM *histo2);
HISTOGRAM *add_histogram (HISTOGRAM *histo1, HISTOGRAM *histo2);
void print_histogram (HISTOGRAM *histo);
void display_histogram (HISTOGRAM *histo);
void display_all_histograms (void);
int histogram_to_lookup (HISTOGRAM *histo, HISTOGRAM *lookup);
long lookup_int (HISTOGRAM *lookup, long value, long factor);
double lookup_real (HISTOGRAM *lookup, double value, double factor);
int histogram_hashing (int tabsize);
void sort_histograms (void);
void release_histogram (HISTOGRAM *histo);

/* ---- Functions in moments.c ---- */

MOMENTS *alloc_moments (double low, double high);
void clear_moments (MOMENTS *mom);
void free_moments (MOMENTS *mom);
void fill_moments (MOMENTS *mom, double value);
void fill_mean (MOMENTS *mom, double value);
void fill_mean_and_sigma (MOMENTS *mom, double value);
void fill_real_moments (MOMENTS *mom, double value, double weight);
void fill_real_mean (MOMENTS *mom, double value, double weight);
void fill_real_mean_and_sigma (MOMENTS *mom, double value, 
     double weight);
int stat_moments (MOMENTS *mom, struct momstat *stmom);

/* ---- Functions in io_histogram.c (optional) ---- */

#ifdef EVENTIO_BASIC_H__LOADED
int write_histograms (HISTOGRAM **phisto, int nhisto,
    IO_BUFFER *iobuf);
int read_histograms (HISTOGRAM **phisto, int nhisto,
    IO_BUFFER *iobuf);
#endif

#ifdef __cplusplus
}
#endif

#endif
