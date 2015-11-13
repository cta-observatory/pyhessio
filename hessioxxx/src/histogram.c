/* ============================================================================

   Copyright (C) 1991, 2001, 2008, 2010  Konrad Bernloehr

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

/** @file histogram.c
 *  @short Manage, fill, and display one- and two-dimensional histograms.
 *
 *  Eventio routines for these types of histograms are available in
 *  io_histogram.c. Conversion to HBOOK format is available
 *  through the @c hdata2hbook (was @c cvt2) program.
 *  Conversion to ROOT format is available through the 
 *  @c hdata2root (was @c cvt3) program. 
 *
 *  Note: multi-threading safety of functions provided in this file
 *  has not been tested extensively. Threads must not delete histograms
 *  shared with other threads when referenced by pointers.
 *
 *  @author  Konrad Bernloehr 
 *  @date    1991 - 2010
 *  @date    CVS $Date: 2014/02/20 10:53:06 $
 *  @version CVS $Revision: 1.21 $
 */

#include "initial.h"
#include "histogram.h"
#include "warning.h"

#ifdef _REENTRANT
static pthread_mutex_t mlock_hist = PTHREAD_MUTEX_INITIALIZER;
#define _HLOCK_ pthread_mutex_lock(&mlock_hist);
#define _HUNLOCK_ pthread_mutex_unlock(&mlock_hist);
#define _WAIT_IF_BUSY_(histo) histogram_lock(histo);
#define _CLEAR_BUSY_(histo) histogram_unlock(histo);
#else
#define _HLOCK_
#define _HUNLOCK_
#define _WAIT_IF_BUSY_(histo)
#define _CLEAR_BUSY_(histo)
#endif

static void initialize_histogram (HISTOGRAM *histo);
static HISTOGRAM *aux_alloc_histogram (int nbins, const char *type);
static void free_histo_contents (HISTOGRAM *histo);
static void display_2d_histogram (HISTOGRAM *histo);

static HISTOGRAM *first_histogram = (HISTOGRAM *) NULL;
static HISTOGRAM *last_histogram = (HISTOGRAM *) NULL;

/* In multi-threaded programs, the following variable */
/* should on be set during program startup. */
FILE *histogram_file;  /* Global variable, initialized to NULL */

static HISTOGRAM **hash_table;
static long hash_size = 0;
static CONST_QUAL short primetab[] = 
     { 131, 233, 353, 541, 751, 1051, 1367, 1511, 1723,
       1931, 2393, 3163, 3907, 5261, 6143, 7187, 8623, 9749, 11321, 15031 };

/* --------------------- histo_lock ------------------------ */

void histogram_lock (HISTOGRAM *histo)
{
#ifdef _REENTRANT
   pthread_mutex_lock(&histo->mlock_this);
#endif
}

/* --------------------- histo_unlock ------------------------ */

void histogram_unlock (HISTOGRAM *histo)
{
#ifdef _REENTRANT
   pthread_mutex_unlock(&histo->mlock_this);
#endif
}


/* --------------------- get_first_histogram ------------------------ */
/**
 *  @short Get a pointer to the first histogram.
 *
 *  Get a pointer to the first histogram in the linked list
 *  of available histograms without making the corresponding
 *  variable global.
 *
 *  @return Pointer to the first histogram in the linked list.
 *
 */

HISTOGRAM *get_first_histogram ()
{
   HISTOGRAM *hptr;
_HLOCK_
   hptr = first_histogram;
_HUNLOCK_
   return(hptr);
}


#ifdef QSORT_FOR_HISTOGRAM_SORTING
static int compare_histograms (HISTOGRAM **h1, HISTOGRAM **h2)
{
   if ( h1 == (HISTOGRAM **) NULL || h2 == (HISTOGRAM **) NULL )
      return 0;
   if ( *h1 == (HISTOGRAM *) NULL || *h2 == (HISTOGRAM *) NULL )
      return 0;
   return ((*h1)->ident-(*h2)->ident);
}
#endif

/* ------------------------ sort_histograms ------------------------- */
/**
 *  Sort histograms in linked list by idents.
 *
 *  @return  (none)
 *
 */

void sort_histograms ()
{
   HISTOGRAM *histo, **phisto;
   size_t nhisto, ihisto, jhisto;

_HLOCK_
   if ( first_histogram == (HISTOGRAM *) NULL )
   {
_HUNLOCK_
      return;
   }

/*
   fprintf(stderr,"Before sorting:\n");
   for (histo=first_histogram; histo!=(HISTOGRAM *) NULL; histo=histo->next)
      fprintf(stderr," %d",histo->ident);
   fprintf(stderr,"\n");
*/

   nhisto = 0;
   for (histo=first_histogram; histo!=(HISTOGRAM *) NULL; histo=histo->next)
      nhisto++;
   if ( (phisto = (HISTOGRAM **) calloc(nhisto,sizeof(HISTOGRAM *))) ==
        (HISTOGRAM **) NULL )
   {
_HUNLOCK_
      return;
   }

   nhisto = 0;
   for (histo=first_histogram; histo!=(HISTOGRAM *) NULL; histo=histo->next)
      phisto[nhisto++] = histo;
#ifdef QSORT_FOR_HISTOGRAM_SORTING
   qsort ((void *)phisto,nhisto,sizeof(HISTOGRAM *),compare_histograms);
#else
   for (ihisto=0; ihisto<nhisto; ihisto++)
      for (jhisto=0; jhisto<ihisto; jhisto++)
         if ( phisto[ihisto]->ident < phisto[jhisto]->ident )
         {
            histo = phisto[ihisto];
            phisto[ihisto] = phisto[jhisto];
            phisto[jhisto] = histo;
         }
#endif
   histo = first_histogram = phisto[0];
   first_histogram->previous = (HISTOGRAM *) NULL;
   for ( ihisto=1; ihisto<nhisto-1; ihisto++ )
   {
      phisto[ihisto]->previous = histo;
      histo->next = phisto[ihisto];
      histo = histo->next;
   }
   histo->next = last_histogram = phisto[nhisto-1];
   last_histogram->previous = histo;
   last_histogram->next = (HISTOGRAM *) NULL;

/*
   fprintf(stderr,"After sorting:\n");
   for (histo=first_histogram; histo!=(HISTOGRAM *) NULL; histo=histo->next)
      fprintf(stderr," %d",histo->ident);
   fprintf(stderr,"\n");
*/

_HUNLOCK_
}

/* --------------------- set_first_histogram ------------------------ */
/**
 *  @short Set a new histogram as the first element (context switching).
 *
 *  To allow 'context switching' of histograms the first element
 *  of the linked list of histograms can be changed by this
 *  function. Before that, the old value should be obtained
 *  with get_first_histogram() and saved.
 *  Note: For context switching it is not necessary to specify the
 *  actually first member of a linked list but any member of a list
 *  can be specifed to activate that list.
 *
 *  @param new_first_histogram A histogram in the new list (may be NULL pointer).
 *
 *  @return none
 *
 */

void set_first_histogram (HISTOGRAM *new_first_histogram)
{
   HISTOGRAM *histo, *lhisto, *fhisto;

_HLOCK_
   histo = lhisto = new_first_histogram;
   while ( histo != (HISTOGRAM *) NULL )
   {
      if ( histo->next != (HISTOGRAM *) NULL )
         lhisto = histo->next;
      histo = histo->next;
   }
   histo = fhisto = new_first_histogram;
   while ( histo != (HISTOGRAM *) NULL )
   {
      if ( histo->previous != (HISTOGRAM *) NULL )
         fhisto = histo->previous;
      histo = histo->previous;
   }
   first_histogram = fhisto;
   last_histogram = lhisto;
_HUNLOCK_
}

/* ------------------- get_histogram_by_ident ----------------------- */
/**
 *  @short Get a histogram with the given ID.
 *
 *  Get the first histogram with a given ident (different from 0)
 *  or return NULL pointer if none exists.
 *
 *  @param  ident  --  The histogram ident to be searched for.
 *
 *  @return Histogram pointer or NULL
 *
 */

HISTOGRAM *get_histogram_by_ident (long ident)
{
   HISTOGRAM *histo;
   HISTOGRAM *hptr = NULL;

   if ( ident <= 0 )
      return ((HISTOGRAM *) NULL);
_HLOCK_
   if ( hash_table != (HISTOGRAM **) NULL )
   {
      if ( (histo = hash_table[ident % hash_size]) != (HISTOGRAM *) NULL )
         if ( histo->ident == ident )
	 {
	    hptr = histo;
            goto finished;
	 }
   }
   histo = first_histogram;
   while ( histo != (HISTOGRAM *) NULL )
   {
      if ( histo->ident == ident )
      {
         if ( hash_table != (HISTOGRAM **) NULL )
            hash_table[ident % hash_size] = histo;
	 {
	    hptr = histo;
	    goto finished;
         }
      }
      else
         histo = histo->next;
   }
   hptr = ((HISTOGRAM *) NULL);
finished:
_HUNLOCK_
   return hptr;
}

/* ------------------------- list_histograms ---------------------------- */
/**
 *  List all available histograms using the 'Output()' function.
 *
 *  @param  ident  --  histogram ident to search or 0
 *
 *  @return (none)
 *
 */

void list_histograms (long ident)
{
   HISTOGRAM *histo;
   int nfound = 0;
   char message[1024]; /* Should be long enough. */

_HLOCK_
   for (histo=first_histogram; histo!=(HISTOGRAM *) NULL; histo=histo->next)
   {
      if ( ident != 0 && ident != histo->ident )
         continue;
      nfound++;
      snprintf(message,sizeof(message),"Id %ld: ",(long)histo->ident);
      if ( histo->title != (char *) NULL &&  histo->title != NULL &&
      	   strlen(message) + strlen(histo->title) + 3 < sizeof(message))
      {
         strcat(message,"'");
         strcat(message,histo->title);
         strcat(message,"', ");
      }
      else
         strcat(message,"Unnamed, ");
      if ( strlen(message) + 50 < sizeof(message) )
      	 continue;
      if ( histo->type != 'I' )
         sprintf(message+strlen(message),"type '%c', ",histo->type);
      if ( histo->nbins_2d > 0 )
         sprintf(message+strlen(message),"%d*%d bins, ",
            histo->nbins,histo->nbins_2d);
      else
         sprintf(message+strlen(message),"%d bins, ",histo->nbins);
      if ( histo->entries == 0 )
         strcat(message,"emtpy.\n");
      else
         sprintf(message+strlen(message),"%ld/%ld entries.\n",
            histo->entries,histo->tentries);
      Output(message);
   }
   if ( nfound == 0 )
   {
      if ( ident == 0 )
         Output("No histograms found.\n");
      else
      {
         sprintf(message,"No histogram with ident %ld was found.\n",ident);
         Output(message);
      }
   }
_HLOCK_
}

/* ------------------------- book_histogram ---------------------- */
/**
 *  @short General histogram booking function, assigning ID and title.
 *
 *  Book a histogram of 1 or 2 dimensions, 'I', 'R', 'F', or 'D' type.
 *  The histogram is allocated (if possible) and the supplied ID number
 *  and title string are assigned.
 *
 *  @param  id     ID number
 *  @param  title  Histogram title string
 *  @param  type   "I" (int, no weights), "R" (real, no weights),
 *   	      "F" (float, with weights), "D" (double, w.w.)
 *  @param  dimension 1 or 2 for 1-D or 2-D histogram
 *  @param  low    Pointer to lower limits (x or x,y for 1-D or 2-D)
 *  @param  high   Pointer to upper limits
 *  @param  nbins  Pointer to no. of bins per dimension (nx or nx, ny)
 *
 *  @return Pointer to new histogram or NULL
 */

HISTOGRAM *book_histogram (long id, const char *title, const char *type, int dimension,
   double *low, double *high, int *nbins)
{
   HISTOGRAM *thisto;

   if ( (thisto = allocate_histogram(type,dimension,low,high,nbins)) !=
        (HISTOGRAM *) NULL )
      describe_histogram(thisto,title,id);
   else
   {
      char msg[200];
      char xbins[50];
      if ( dimension == 2 )
         sprintf(xbins,"%d*%d", nbins[0], nbins[1]);
      else
         sprintf(xbins,"%d", nbins[0]);
      sprintf(msg,"Allocation of %2d-D (%s bins) histogram ID %ld failed.\n", 
         dimension, xbins, id);
      Warning(msg);
   }
   return(thisto);
}

/* ------------------------- book_1d_histogram ---------------------- */
/**
 *  @short Simplified histogram booking function for one-dimensional histograms, assigning ID and title.
 *
 *  Book a histogram of one dimension, 'I', 'R', 'F', or 'D' type.
 *  The histogram is allocated (if possible) and the supplied ID number
 *  and title string are assigned.
 *
 *  @param  id     ID number
 *  @param  title  Histogram title string
 *  @param  type   "I" (int, no weights), "R" (real, no weights),
 *   	      "F" (float, with weights), "D" (double, w.w.)
 *  @param  low    Lower limit (x)
 *  @param  high   Upper limit (x)
 *  @param  nbins  No. of bins (nx)
 *
 *  @return Pointer to new histogram or NULL
 */

HISTOGRAM *book_1d_histogram (long id, const char *title, const char *type, 
   double low, double high, int nbins)
{
   HISTOGRAM *thisto;

   if ( (thisto = allocate_histogram(type,1,&low,&high,&nbins)) !=
        (HISTOGRAM *) NULL )
      describe_histogram(thisto,title,id);
   else
   {
      char msg[200];
      sprintf(msg,"Allocation of 1-D (%d bins) histogram ID %ld failed.\n", 
         nbins, id);
      Warning(msg);
   }
   return(thisto);
}

/* ------------------------ book_int_histogram --------------------- */
/**
 *  @short Book and integer-type histogram (content incremented by one per entry).
 *
 *  Like book_histogram() but for 'I' type histograms only (1-D or 2-D)
 *
 *  @param  id     ID number
 *  @param  title  Histogram title string
 *  @param  dimension 1 or 2 for 1-D or 2-D histogram
 *  @param  low    Pointer to lower limits (x or x,y for 1-D or 2-D)
 *  @param  high   Pointer to upper limits
 *  @param  nbins  Pointer to no. of bins per dimension (nx or nx, ny)
 *
 *  @return Pointer to new histogram or NULL
 */

HISTOGRAM *book_int_histogram (long id, const char *title, int dimension, 
   long *low, long *high, int *nbins)
{
   HISTOGRAM *thisto;

   if ( low == (long *) NULL || high == (long *) NULL || nbins == (int *) NULL )
   {
      Warning("Invalid parameters for histogram allocation.");
      return((HISTOGRAM *) NULL);
   }

   if ( dimension == 1 )
      thisto = alloc_int_histogram(low[0],high[0],nbins[0]);
   else if ( dimension == 2 )
      thisto = alloc_2d_int_histogram(low[0],high[0],nbins[0],
          low[1],high[1],nbins[1]);
   else
   {
      Warning("Invalid dimension for histogram allocation.");
      thisto = (HISTOGRAM *) NULL;
   }

   if ( thisto != (HISTOGRAM *) NULL )
      describe_histogram(thisto,title,id);
   return(thisto);
}

/* ------------------------ allocate_histogram --------------------- */
/**
 *  @short Allocate any histogram without ID and title.
 *
 *  Allocate a histogram of 1 or 2 dimensions, 'I', 'R', 'F' or 'D' type,
 *  without assigning an ID number and title string to it.
 *  To avoid the (long) <--> (double) typecasts, the direct calls to
 *  alloc_int_histogram() and alloc_2d_int_histogram() are recommended
 *  for integer-limits histograms (type 'I').
 *
 *  @param  type   "I" (int, no weights), "R" (real, no weights),
 *   	      "F" (float, with weights), "D" (double, w.w.)
 *  @param  dimension 1 or 2 for 1-D or 2-D histogram
 *  @param  low    Pointer to lower limits (x or x,y for 1-D or 2-D)
 *  @param  high   Pointer to upper limits
 *  @param  nbins  Pointer to no. of bins per dimension (nx or nx, ny)
 *
 *  @return  Pointer to new histogram or NULL
 */

HISTOGRAM *allocate_histogram (const char *type, int dimension, 
   double *low, double *high, int *nbins)
{
   HISTOGRAM *thisto;
   int idim;
   double size;
   int tbins;

   if ( dimension != 1 && dimension != 2 )
   {
      Warning("Invalid histogram dimension. Cannot allocate histogram.");
      return((HISTOGRAM *) NULL);
   }
   if ( type == (char *) NULL || low == (double *) NULL ||
        high == (double *) NULL || nbins == (int *) NULL )
   {
      Warning("Invalid parameters for histogram allocation.");
      return((HISTOGRAM *) NULL);
   }
   size = nbins[0];
   if ( dimension == 2 )
      size *= nbins[1];
#ifdef OS_MSDOS
   if ( nbins[0] < 0 || size > 16000. || size <= 0. )
#else
   if ( nbins[0] < 0 || size > 16000000. || size <= 0. )
#endif
   {
      Warning("Invalid number of bins for histogram allocation.");
      return((HISTOGRAM *) NULL);
   }
   tbins = 1;
   for ( idim=0; idim<dimension; idim++ )
   {
      if ( low[idim] >= high[idim] )
      {
         Warning("Invalid limits for allocation of histogram.");
         return((HISTOGRAM *) NULL);
      }
      tbins *= nbins[idim];
   }

   if ( *type == 'I' ) /* Supported but not recommended due to conversions */
   {
      for ( idim=0; idim<dimension; idim++ )
         if ( fabs(low[idim]) > 2.e9 || fabs(high[idim]) > 2.e9 )
         {
            Warning("Invalid limits for allocation of type 'I' histogram.");
            return((HISTOGRAM *) NULL);
         }
      if ( dimension == 1 )
         return(alloc_int_histogram((long)low[0],(long)high[0],nbins[0]));
      else
         return(alloc_2d_int_histogram((long)low[0],(long)high[0],nbins[0],
             (long)low[1],(long)high[1],nbins[1]));
   }
   else if ( *type != 'R' && *type != 'F' && *type != 'D' )
   {
      char message[1024];
      sprintf(message,"Don't know how to allocate a type '%c' histogram.",*type);
      Warning(message);
      return((HISTOGRAM *) NULL);
   }

   if ( (thisto = aux_alloc_histogram(tbins,type)) == (HISTOGRAM *) NULL )
   {
      Warning("Histogram allocation failed.");
      return ((HISTOGRAM *) NULL);
   }

   thisto->specific.real.lower_limit = low[0];
   thisto->specific.real.upper_limit = high[0];
   thisto->specific.real.inverse_binwidth =
          (double) nbins[0] / (high[0]-low[0]);
   thisto->nbins = nbins[0];
   if ( dimension < 2 )
      thisto->nbins_2d = 0;
   else
   {
      thisto->specific_2d.real.lower_limit = low[1];
      thisto->specific_2d.real.upper_limit = high[1];
      thisto->specific_2d.real.inverse_binwidth =
             (double) nbins[1] / (high[1]-low[1]);
      thisto->nbins_2d = nbins[1];
   }

   initialize_histogram(thisto);

   return (thisto);
}

/* -------------------------- alloc_int_histogram ---------------------- */
/**
 *  @short Allocate memory for a 1-D 'int' histogram and initialize it.
 *
 *  Resulting histogram has integer range limits and integer contents
 *  (incremented by one per entry).
 *
 *  @param  low  lower limit of values to be covered by histogram
 *  @param  high  upper limit ...
 *  @param  nbins  the number of bins to be allocated
 *
 *  @return pointer to allocated histogram or NULL
 *
 */

HISTOGRAM *alloc_int_histogram (long low, long high, int nbins)
{
   HISTOGRAM *thisto;

   if ( low >= high || nbins <= 0 )
      return ((HISTOGRAM *) NULL);
   if ( (thisto = aux_alloc_histogram(nbins,"I")) == (HISTOGRAM *) NULL )
   {
      Warning("Histogram allocation failed.");
      return ((HISTOGRAM *) NULL);
   }

   thisto->specific.integer.lower_limit = low;
   thisto->specific.integer.upper_limit = high;
   thisto->nbins = nbins;
   thisto->nbins_2d = 0;
   thisto->specific.integer.width = high-low;

   initialize_histogram(thisto);

   return (thisto);
}

/* ----------------------- alloc_real_histogram -------------------- */
/**
 *  @short Allocate memory for a 1-D 'real' histogram and initialize it.
 *
 *  Resulting histogram has floating point range limits and
 *  integer contents (incremented by one per entry).
 *
 *  @param  low  lower limit of values to be covered by histogram
 *  @param  high  upper limit ...
 *  @param  nbins  the number of bins to be allocated
 *
 *  @return pointer to allocated histogram or NULL
 *
 */

HISTOGRAM *alloc_real_histogram (double low, double high, int nbins)
{
   return(allocate_histogram("R",1,&low,&high,&nbins));
}

/* ------------------------- alloc_2d_int_histogram ---------------------- */
/**
 *  @short Allocate memory for a 2-D 'int' histogram and initialize it.
 *
 *  Resulting histogram has integer range limits and integer contents
 *  (incremented by one per entry).
 *
 *  @param  xlow    lower limit of values in X to be covered by histogram
 *  @param  xhigh   upper limit ...
 *  @param  nxbins  the number of bins to be allocated in X
 *  @param  ylow    lower limit of values in Y to be covered by histogram
 *  @param  yhigh   upper limit ...
 *  @param  nybins  the number of bins to be allocated in Y
 *
 *  @return pointer to allocated histogram or NULL
 */

HISTOGRAM *alloc_2d_int_histogram (long xlow, long xhigh,
   int nxbins, long ylow, long yhigh, int nybins)
{
   HISTOGRAM *thisto;

   if ( xlow >= xhigh || nxbins <= 0 || ylow >= yhigh || nybins <= 0 )
      return ((HISTOGRAM *) NULL);
   if ( (thisto = aux_alloc_histogram(nxbins*nybins,"I")) == (HISTOGRAM *) NULL )
   {
      Warning("Histogram allocation failed.");
      return ((HISTOGRAM *) NULL);
   }

   thisto->specific.integer.lower_limit = xlow;
   thisto->specific.integer.upper_limit = xhigh;
   thisto->specific.integer.width = xhigh-xlow;
   thisto->specific_2d.integer.lower_limit = ylow;
   thisto->specific_2d.integer.upper_limit = yhigh;
   thisto->specific_2d.integer.width = yhigh-ylow;
   thisto->nbins = nxbins;
   thisto->nbins_2d = nybins;

   initialize_histogram(thisto);

   return (thisto);
}

/* ------------------------- alloc_2d_real_histogram -------------------- */
/**
 *  @short Allocate memory for a 2-D 'int' histogram and initialize it.
 *
 *  Resulting histogram has floating point range limits and integer contents
 *  (incremented by one per entry).
 *
 *  @param  xlow    lower limit of values in X to be covered by histogram
 *  @param  xhigh   upper limit ...
 *  @param  nxbins  the number of bins to be allocated in X
 *  @param  ylow    lower limit of values in Y to be covered by histogram
 *  @param  yhigh   upper limit ...
 *  @param  nybins  the number of bins to be allocated in Y
 *
 *  @return pointer to allocated histogram or NULL
 */

HISTOGRAM *alloc_2d_real_histogram (double xlow, double xhigh,
   int nxbins, double ylow, double yhigh, int nybins)
{
   double low[2], high[2];
   int nbins[2];

   low[0] = xlow;
   low[1] = ylow;
   high[0] = xhigh;
   high[1] = yhigh;
   nbins[0] = nxbins;
   nbins[1] = nybins;

   return(allocate_histogram("R",2,low,high,nbins));
}

/* ----------------------- aux_alloc_histogram ------------------------ */
/** For internal purpose only */

static HISTOGRAM *aux_alloc_histogram (int ncounts, const char *type)
{
   HISTOGRAM *thisto;

   if ( (thisto = (HISTOGRAM *) calloc(1,(size_t)sizeof(HISTOGRAM))) ==
        (HISTOGRAM *) NULL )
      return ((HISTOGRAM *) NULL);

   if ( *type == 'I' || *type == 'R' )
   {
      if ( (thisto->counts = (unsigned long *)
           calloc(1,(size_t)(ncounts*sizeof(unsigned long)))) ==
           (unsigned long *) NULL )
      {
         free(thisto);
         return ((HISTOGRAM *) NULL);
      }
      thisto->extension = (struct Histogram_Extension *) NULL;
   }
   else if ( *type == 'F' || *type == 'D' )
   {
      struct Histogram_Extension *he;
      int err = 0;
      thisto->counts = (unsigned long *) NULL;
      if ( (he = thisto->extension = (struct Histogram_Extension *)
           calloc(1,sizeof(struct Histogram_Extension))) ==
           (struct Histogram_Extension *) NULL )
      {
         free(thisto);
         return ((HISTOGRAM *) NULL);
      }
      if ( *type == 'F' )
      {
         if ( (he->fdata = (float *) malloc((size_t)(ncounts*
               sizeof(float)))) == (float *) NULL )
            err = 1;
         he->ddata = NULL;
      }
      else
      {
         if ( (he->ddata = (double *) malloc((size_t)(ncounts*
               sizeof(double)))) == (double *) NULL )
            err = 1;
         he->fdata = NULL;
      }
      if ( err )
      {
         free(thisto->extension);
         thisto->extension = (struct Histogram_Extension *) NULL;
         free(thisto);
         return ((HISTOGRAM *) NULL);
      }
   }
   else
   {
      free(thisto);
      return ((HISTOGRAM *) NULL);
   }

   thisto->type = *type;

   return (thisto);
}

/* ---------------------- initialize_histogram ---------------------- */
/** For internal purpose only */

static void initialize_histogram (HISTOGRAM *histo)
{
   /* The histogram has no title so far */
   histo->title = (char *) NULL;
   histo->ident = 0;

_HLOCK_
   /* Add the new histogram to the double linked list of histograms */
   histo->previous = last_histogram;
   if ( last_histogram != (HISTOGRAM *) NULL )
      last_histogram->next = histo;
   last_histogram = histo;
   if ( first_histogram == (HISTOGRAM *) NULL )
      first_histogram = histo;
   histo->next = (HISTOGRAM *) NULL;
_HUNLOCK_

   /* Initialize the histogram contents */
   clear_histogram(histo);
}

/* ---------------------- describe_histogram ------------------------ */
/**
 *  Add a describing title to a histogram previously allocated.
 *
 *  @param  histo  Histogram to which the title should be added
 *  @param  title  The title string. This is ignored if the histogram
 *                 already has a title.
 *  @param  ident  Identification number, must be unique (or 0)
 *                 if any I/O is intended, because read_histogram()
 *                 deletes a pre-existing histogram with the same ID.
 *
 *  @return none
 */

void describe_histogram ( HISTOGRAM *histo, const char *title, long ident )
{
   char *temp;
   HISTOGRAM *thisto;

   if ( histo == (HISTOGRAM *) NULL )
      return;
_WAIT_IF_BUSY_(histo)
   if ( title != (const char *) NULL && histo->title == (char *) NULL )
      if ( *title )
         if ( (temp = (char *)malloc(strlen(title)+1)) != (char *) NULL )
         {
            strcpy(temp,title);
            histo->title = temp;
         }
   if ( ident < 0 )
      ident = 0;
   else if ( (thisto = get_histogram_by_ident(ident)) != (HISTOGRAM *) NULL &&
              thisto != histo )
   {
      char message[1024];
      sprintf(message,
      "Histogram id %ld used for more than one histogram ('%s' and '%s')",
           ident,(title==(char*)NULL)?"unnamed":title,
           (thisto->title==(char*)NULL)?"unnamed":thisto->title);
      Information(message);
   }
   histo->ident = ident;
_CLEAR_BUSY_(histo)

_HLOCK_
   if ( hash_table != (HISTOGRAM **) NULL && ident > 0 )
      hash_table[ident % hash_size] = histo;
_HUNLOCK_
}

/* ---------------------- clear_histogram --------------------------- */
/**
 *  Initialize an existing histogram.
 *
 *  @param histo -- pointer to histogram
 *
 *  @return (none)
 *
 */

void clear_histogram (HISTOGRAM *histo)
{
   REGISTER int i, ncounts;
   struct Histogram_Extension *he;

   if ( histo == (HISTOGRAM *) NULL )
      return;
_WAIT_IF_BUSY_(histo)
   histo->underflow = histo->overflow =
      histo->underflow_2d = histo->overflow_2d = 0;
   if ( histo->type == 'I' || histo->type == 'i' )
      histo->specific.integer.sum = histo->specific.integer.tsum =
         histo->specific_2d.integer.sum = histo->specific_2d.integer.tsum = 0;
   else
      histo->specific.real.sum = histo->specific.real.tsum =
         histo->specific_2d.real.sum = histo->specific_2d.real.tsum = 0.;
   histo->entries = histo->tentries = 0;
   if ( histo->nbins_2d > 0 )
      ncounts = histo->nbins * histo->nbins_2d;
   else
      ncounts = histo->nbins;
   if ( (he = histo->extension) != (struct Histogram_Extension *) NULL )
   {
      if ( he->fdata != (float *) NULL )
         for ( i=0; i<ncounts; i++ )
            he->fdata[i] = (float) 0.;
      else if ( he->ddata != (double *) NULL )
         for ( i=0; i<ncounts; i++ )
            he->ddata[i] = 0.;
      he->content_all = he->content_inside = 0.;
      for (i=0; i<8; i++)
         he->content_outside[i] = 0.;
   }
   if ( histo->counts != (unsigned long *) NULL )
      for ( i=0; i<ncounts; i++ )
         histo->counts[i] = 0;
_CLEAR_BUSY_(histo)
}

/* ------------------------ free_histogram ------------------------- */
/**
 *  Free a histogram completely (both data and control structure).
 *
 *  Deallocates memory previously allocated to a histogram.
 *  If release_histogram was applied to that histogram before,
 *  it cannot be reallocated.
 *
 *  @param histo -- pointer to previously allocated histogram
 *
 *  @return (none)
 */

void free_histogram (HISTOGRAM *histo)
{
   /* Always check for a NULL pointer before freeing memory. */
   if ( histo == (HISTOGRAM *) NULL )
      return;

_WAIT_IF_BUSY_(histo)

   /* Free all pointers inside the histogram structure. */
   free_histo_contents(histo);

   /* Remove the histogram from the linked histogram list. */
   unlink_histogram(histo);

_CLEAR_BUSY_(histo)

   /* Free the entire histogram structure. */
   free(histo);
}

/* --------------------- free_histo_contents ----------------------- */
/**
 *  Free the contents (data pointers) of a histogram to be released or removed.
 *
 *  @param  Pointer to histogram that should be 'cleaned'.
 *
 *  @return  (none)
 *
 */

static void free_histo_contents (HISTOGRAM *histo)
{
   if ( histo->counts != (unsigned long *) NULL )
   {
      free((void*)histo->counts);
      histo->counts = (unsigned long *) NULL;
   }
   if ( histo->title != (char *) NULL )
   {
      free((void*)histo->title);
      histo->title = (char *) NULL;
   }
   if ( histo->extension != (struct Histogram_Extension *) NULL )
   {
      struct Histogram_Extension *he = histo->extension;
      if ( he->fdata != (float *) NULL )
      {
         free((void*)he->fdata);
         he->fdata = (float *) NULL;
      }
      if ( he->ddata != (double *) NULL )
      {
         free((void*)he->ddata);
         he->ddata = (double *) NULL;
      }
      free((void*)he);
      histo->extension = (struct Histogram_Extension *) NULL;
   }
}

/* --------------------- free_all_histograms ----------------------- */
/**
 *  Deletes all histograms which are included in the linked list of histograms.
 *
 *  @return (none)
 *
 */

void free_all_histograms ()
{
   HISTOGRAM *histo, *nhisto;
   histo = first_histogram;
   while ( histo != (HISTOGRAM *) NULL )
   {
      nhisto = histo->next;
      free_histogram(histo);
      histo = nhisto;
   }
}

/* --------------------- unlink_histogram --------------------- */
/**
 *  @short Remove a histogram from the list without destroying it.
 *
 *  Remove a histogram from the linked list of histograms.
 *  That histogram will therefore not be found by any subsequent
 *  call to 'free_all_histograms()', display_all_histograms()',
 *  and 'get_histogram_by_ident()'.
 *
 *  @param  histo   Pointer to histogram.
 *
 *  @return (none)
 */

void unlink_histogram (HISTOGRAM *histo)
{
   if ( histo == (HISTOGRAM *) NULL )
      return;

_HLOCK_
   /* Remove the histogram from the double linked list. */
   if ( histo->previous != (HISTOGRAM *) NULL )
      (histo->previous)->next = histo->next;
   if ( histo->next != (HISTOGRAM *) NULL )
      (histo->next)->previous = histo->previous;
   if ( last_histogram == histo )
      last_histogram = histo->previous;
   if ( first_histogram == histo )
      first_histogram = histo->next;
   if ( hash_table != (HISTOGRAM **) NULL && histo->ident > 0 )
      hash_table[histo->ident % hash_size] = (HISTOGRAM *) NULL;
_HUNLOCK_
}

/* ------------------------ fill_int_histogram -------------------- */
/**
 *  @short Increment a bin of a 1-D 'int' histogram by one.
 *
 *  Either a count for one of the bins in the histogram range
 *  is incremented or an underflow or overflow count. For the
 *  calculation of the mean value and truncated mean value sums
 *  of values and number of histogram entries are updated as well.
 *
 *  @param  histo  Pointer to histogram
 *  @param  value  Position where an entry is to be added
 *                 (may be outside the given range)
 *
 *  @return 0 (o.k.), -1 (no histogram that can be filled)
 *
 */

int fill_int_histogram (HISTOGRAM *histo, long value)
{
   int indx;

   if ( histo == (HISTOGRAM *) NULL )
      return -1;
   if ( histo->type != 'I' )
   {
      if ( histo->type == 'R' || histo->type == 'F' || histo->type == 'D' )
         return(fill_real_histogram(histo,(double)value));
      return -1;
   }

_WAIT_IF_BUSY_(histo)
   /* Sum of values and no. of entries (for mean value) */
   histo->specific.integer.sum += (long) value;
   histo->entries++;

   if ( value < histo->specific.integer.lower_limit )
   { 
      histo->underflow++; 
_CLEAR_BUSY_(histo)
      return 0;
   }
   if ( value >= histo->specific.integer.upper_limit )
   { 
      histo->overflow++; 
_CLEAR_BUSY_(histo)
      return 0; 
   }

   indx = (int) (((value-histo->specific.integer.lower_limit) *
         histo->nbins) / histo->specific.integer.width);

   if ( indx>=histo->nbins )  /* check for rounding errors */
      indx = histo->nbins - 1;
   else if ( indx<0 )  /* check for rounding errors */
      indx = 0;

   /* Sum of values and no. of entries (for truncated mean value) */
   histo->specific.integer.tsum += (long) value;
   histo->tentries++;
   /* Increment histogram bin but avoid count overflow (wrap around) */
   if ( histo->counts[indx] < MAX_HISTCOUNT )
      histo->counts[indx]++;

_CLEAR_BUSY_(histo)
   return 0;
}

/* ----------------------- fill_real_histogram ----------------------- */
/**
 *  @short Increment a bin of a 1-D 'real' histogram by one.
 *
 *  Either a count for one of the bins in the histogram range
 *  is incremented or an underflow or overflow count. For the
 *  calculation of the mean value and truncated mean value sums
 *  of values and number of histogram entries are updated as well.
 *
 *  @param  histo  Pointer to histogram
 *  @param  value  Position where an entry is to be added
 *                 (may be outside the given range)
 *
 *  @return 0 (o.k.), -1 (no histogram that can be filled)
 *
 */

int fill_real_histogram (HISTOGRAM *histo, double value)
{
   int indx;

   if ( histo == (HISTOGRAM *) NULL )
      return -1;
   if ( histo->type != 'R' )
   {
      if ( histo->type == 'F' || histo->type == 'D' )
         return(fill_weighted_histogram(histo,value,1.0));
      else
         return -1;
   }

_WAIT_IF_BUSY_(histo)
   /* Sum of values and no. of entries (for mean value) */
   histo->specific.real.sum += (double) value;
   histo->entries++;

   if ( value < histo->specific.real.lower_limit )
   { 
      histo->underflow++; 
_CLEAR_BUSY_(histo)
      return 0;
   }
   if ( value >= histo->specific.real.upper_limit )
   { 
      histo->overflow++; 
_CLEAR_BUSY_(histo)
      return 0;
   }

   indx = (int) ((value-histo->specific.real.lower_limit) *
         histo->specific.real.inverse_binwidth);

   if ( indx>=histo->nbins )  /* check for rounding errors */
      indx = histo->nbins - 1;
   else if ( indx<0 )  /* check for rounding errors */
      indx = 0;

   /* Sum of values and no. of entries (for truncated mean value) */
   histo->specific.real.tsum += (double) value;
   histo->tentries++;
   /* Increment histogram bin but avoid count overflow (wrap around) */
   if ( histo->counts[indx] < MAX_HISTCOUNT )
      histo->counts[indx]++;

_CLEAR_BUSY_(histo)
   return 0;
}

/* ---------------------- fill_weighted_histogram ---------------------- */
/**
 *  @short Add an entry to a weighted 1-D histogram.
 *
 *  Increment a bin of a histogram by a given weight rather than by 1.
 *  This requires a suitable histogram type 'F' or 'D'.
 *
 *  @param  histo   Pointer to histogram.
 *  @param  value   Position where an entry is to be added.
 *  @param  weight  The weight of that entry.
 *
 *  @return 0 (o.k.),  -1 (no histogram that can be filled with weights)
 */

int fill_weighted_histogram (HISTOGRAM *histo, double value,
   double weight)
{
   int indx;
   struct Histogram_Extension *he;

   if ( histo == (HISTOGRAM *) NULL )
      return -1;
   if ( histo->type != 'F' && histo->type != 'D' )
   {
      char message[1024];
      sprintf(message,"Weighted filling invalid for type '%c' histogram %ld",
         histo->type,histo->ident);
      Warning(message);
      return -1;
   }

_WAIT_IF_BUSY_(histo)
   /* Sum of values and no. of entries (for mean value) */
   histo->specific.real.sum += weight * value;
   histo->entries++;

   if ( (he = histo->extension) == (struct Histogram_Extension *) NULL )
   {
_CLEAR_BUSY_(histo)
      return -1;
   }

   he->content_all += weight;

   if ( value < histo->specific.real.lower_limit )
   { 
      histo->underflow++; 
      he->content_outside[0] += weight; 
_CLEAR_BUSY_(histo)
      return 0;
   }
   else if ( value >= histo->specific.real.upper_limit )
   { 
      histo->overflow++; 
      he->content_outside[1] += weight; 
_CLEAR_BUSY_(histo)
      return 0;
   }

   he->content_inside += weight;

   indx = (int) ((value-histo->specific.real.lower_limit) *
         histo->specific.real.inverse_binwidth);

   if ( indx>=histo->nbins )  /* check for rounding errors */
      indx = histo->nbins - 1;
   else if ( indx<0 )  /* check for rounding errors */
      indx = 0;

   /* Sum of values and no. of entries (for truncated mean value) */
   histo->specific.real.tsum += weight * value;
   histo->tentries++;
   if ( he->fdata != (float *) NULL )   /* 'F' type */
      he->fdata[indx] += (float) weight;
   else if ( he->ddata != (double *) NULL )  /* 'D' type */
      he->ddata[indx] += weight;

_CLEAR_BUSY_(histo)
   return 0;
}

/* --------------------- fill_2d_int_histogram --------------------- */
/**
 *  @short Increment a bin of a 2-D 'int' histogram by one.
 *
 *  Increment a bin of a 2-D histogram by one.
 *  Either a count for one of the bins in the histogram range
 *  is incremented or an underflow or overflow count. For the
 *  calculation of the mean value and truncated mean value sums
 *  of values and number of histogram entries are updated as well.
 *
 *  Arguments: histo -- pointer to histogram
 *	       xvalue, yvalue -- X and Y positions where an entry
 *		       is to be to the histogram (they may be outside
 *		       the given ranges)
 *
 *  Return value:  0 (o.k.),  -1 (no histogram that can be filled)
 */

int fill_2d_int_histogram (HISTOGRAM *histo, long xvalue,
     long yvalue)
{
   int indx, indy;

   if ( histo == (HISTOGRAM *) NULL )
      return  -1;
   if ( histo->type != 'I' )
      if ( histo->type == 'R' || histo->type == 'F' || histo->type == 'D' )
         return(fill_2d_real_histogram(histo,(double)xvalue,(double)yvalue));
   if ( histo->nbins_2d <= 0 )
      return(fill_int_histogram(histo,xvalue));

_WAIT_IF_BUSY_(histo)
   histo->specific.integer.sum += xvalue;
   histo->specific_2d.integer.sum += yvalue;
   histo->entries++;

   if ( xvalue < histo->specific.integer.lower_limit )
   { 
      histo->underflow++; 
_CLEAR_BUSY_(histo)
      return 0;
   }
   if ( xvalue >= histo->specific.integer.upper_limit )
   {
      histo->overflow++; 
_CLEAR_BUSY_(histo)
      return 0;
   }
   if ( yvalue < histo->specific_2d.integer.lower_limit )
   {
      histo->underflow_2d++; 
_CLEAR_BUSY_(histo)
      return 0;
   }
   if ( yvalue >= histo->specific_2d.integer.upper_limit )
   {
      histo->overflow_2d++;
_CLEAR_BUSY_(histo)
      return 0;
   }

   indx = (int) (((xvalue-histo->specific.integer.lower_limit) *
         histo->nbins) / histo->specific.integer.width);
   indy = (int) (((yvalue-histo->specific_2d.integer.lower_limit) *
         histo->nbins_2d) / histo->specific_2d.integer.width);
   /* Check for rounding errors */
   if ( indx >= histo->nbins )
      indx = histo->nbins - 1;
   else if ( indx < 0 )
      indx = 0;
   else if ( indy >= histo->nbins_2d )
      indy = histo->nbins_2d - 1;
   else if ( indy < 0 )
      indy = 0;

   /* Sum of values and no. of entries (for truncated mean value) */
   histo->specific.integer.tsum += xvalue;
   histo->specific_2d.integer.tsum += yvalue;
   histo->tentries++;
   histo->counts[indy*histo->nbins + indx]++;

_CLEAR_BUSY_(histo)
   return 0;
}

/* ----------------------- fill_2d_real_histogram --------------------- */
/**
 *  @short Increment a bin of a 2-D 'real' histogram by one.
 *
 *  Increment a bin of a 2-D histogram by one.
 *  Either a count for one of the bins in the histogram range
 *  is incremented or an underflow or overflow count. For the
 *  calculation of the mean value and truncated mean value sums
 *  of values and number of histogram entries are updated as well.
 *
 *  @param histo Pointer to histogram
 *  @param xvalue X position where an entry
 *	is to be to the histogram (may be outside the given ranges)
 *  @param yvalue Y position where an entry
 *	is to be to the histogram (may be outside the given ranges)
 *
 *  @return  0 (o.k.),  -1 (no histogram that can be filled)
 */

int fill_2d_real_histogram (HISTOGRAM *histo, double xvalue,
     double yvalue)
{
   int indx, indy;

   if ( histo == (HISTOGRAM *) NULL )
      return -1;
   if ( histo->type != 'R' )
      if ( histo->type == 'F' || histo->type == 'D' )
         return(fill_2d_weighted_histogram(histo,xvalue,yvalue,1.));
   if ( histo->nbins_2d <= 0 )
      return(fill_real_histogram(histo,xvalue));

_WAIT_IF_BUSY_(histo)
   histo->specific.real.sum += xvalue;
   histo->specific_2d.real.sum += yvalue;
   histo->entries++;

   if ( xvalue < histo->specific.real.lower_limit )
   {
      histo->underflow++;
_CLEAR_BUSY_(histo)
      return 0;
   }
   if ( xvalue >= histo->specific.real.upper_limit )
   {
      histo->overflow++;
_CLEAR_BUSY_(histo)
      return 0;
   }
   if ( yvalue < histo->specific_2d.real.lower_limit )
   {
      histo->underflow_2d++;
_CLEAR_BUSY_(histo)
      return 0;
   }
   if ( yvalue >= histo->specific_2d.real.upper_limit )
   {
      histo->overflow_2d++;
_CLEAR_BUSY_(histo)
      return 0;
   }

   indx = (int) ((xvalue-histo->specific.real.lower_limit) *
         histo->specific.real.inverse_binwidth);
   indy = (int) ((yvalue-histo->specific_2d.real.lower_limit) *
         histo->specific_2d.real.inverse_binwidth);
   /* Check for rounding errors */
   if ( indx >= histo->nbins )
      indx = histo->nbins - 1;
   else if ( indx < 0 )
      indx = 0;
   else if ( indy >= histo->nbins_2d )
      indy = histo->nbins_2d - 1;
   else if ( indy < 0 )
      indy = 0;

   /* Sum of values and no. of entries (for truncated mean value) */
   histo->specific.real.tsum += xvalue;
   histo->specific_2d.real.tsum += yvalue;
   histo->tentries++;
   histo->counts[indy*histo->nbins + indx]++;

_CLEAR_BUSY_(histo)
   return 0;
}

/* --------------------- fill_2d_weighted_histogram ------------------- */
/**
 *  @short Add an entry to a weighted 2-D histogram.
 *
 *  Increment a bin of a 2-D histogram by a given weight rather than by 1.
 *  This requires a suitable histogram type 'F' or 'D'.
 *
 *  @param  histo   Pointer to histogram.
 *  @param  xvalue  X posistion where an entry is to be added.
 *  @param  yvalue  Y posistion where an entry is to be added.
 *  @param  weight  The weight of that entry.
 *
 *  @return 0 (o.k.),  -1 (no histogram that can be filled with weights)
 */

int fill_2d_weighted_histogram (HISTOGRAM *histo, double xvalue,
     double yvalue, double weight)
{
   int indx, indy, izone;

   if ( histo == (HISTOGRAM *) NULL )
      return -1;
   if ( histo->type != 'F' && histo->type != 'D' )
   {
      char message[1024];
      sprintf(message,"Weighted filling invalid for type '%c' histogram %ld",
         histo->type,histo->ident);
      Warning(message);
      return -1;
   }
   if ( histo->nbins_2d <= 0 )
      return(fill_weighted_histogram(histo,xvalue,weight));

_WAIT_IF_BUSY_(histo)
   histo->specific.real.sum += weight * xvalue;
   histo->specific_2d.real.sum += weight * yvalue;
   histo->entries++;
   (histo->extension)->content_all += weight;

   izone = 8;
   if ( xvalue < histo->specific.real.lower_limit )
      { izone -= 2; histo->underflow++; }
   else if ( xvalue >= histo->specific.real.upper_limit )
      { izone -= 1; histo->overflow++; }
   if ( yvalue < histo->specific_2d.real.lower_limit )
      { if ( (izone -= 6) == 2 ) histo->underflow_2d++; }
   else if ( yvalue >= histo->specific_2d.real.upper_limit )
      { if ( (izone -= 3) == 5 ) histo->overflow_2d++; }

   if ( izone < 8 )
   {
      (histo->extension)->content_outside[izone] += weight;
_CLEAR_BUSY_(histo)
      return 0;
   }
   (histo->extension)->content_inside += weight;

   indx = (int) ((xvalue-histo->specific.real.lower_limit) *
         histo->specific.real.inverse_binwidth);
   indy = (int) ((yvalue-histo->specific_2d.real.lower_limit) *
         histo->specific_2d.real.inverse_binwidth);
   /* Check for rounding errors */
   if ( indx >= histo->nbins )
      indx = histo->nbins - 1;
   else if ( indx < 0 )
      indx = 0;
   if ( indy >= histo->nbins_2d )
      indy = histo->nbins_2d - 1;
   else if ( indy < 0 )
      indy = 0;

   /* Sum of values and no. of entries (for truncated mean value) */
   histo->specific.real.tsum += weight * xvalue;
   histo->specific_2d.real.tsum += weight * yvalue;
   histo->tentries++;
   if ( histo->type == 'F' )
      (histo->extension)->fdata[indy*histo->nbins + indx] += (float) weight;
   else
      (histo->extension)->ddata[indy*histo->nbins + indx] += weight;

_CLEAR_BUSY_(histo)
   return 0;
}

/* --------------------------- fill_histogram ------------------------- */
/**
 *  @short Fill any type of 1-D or 2-D histogram known by its pointer.
 *
 *  Generic histogram fill function that can be used for type 'I', 'R',
 *  'F', and 'D' histograms, although it is not recommended for type 'I'
 *  histograms, due to type conversions.
 *
 *  @param  histo   Pointer to histogram.
 *  @param  xvalue  X posistion where an entry is to be added.
 *  @param  yvalue  Y posistion (ignored for 1-D histograms)
 *  @param  weight  The weight of that entry (must be 1.0 for
 *			 'I' and 'R' type histograms).
 *
 *  @return 0 (o.k.),  -1 (no histogram that can be filled)
 */

int fill_histogram (HISTOGRAM *histo, double xvalue, double yvalue,
   double weight)
{
   if ( histo == (HISTOGRAM *) NULL )
      return -1;
   if ( histo->type == 'I' || histo->type == 'R' )
   {
      if ( weight != 1. )
      {
      	 char message[1024];
         sprintf(message,"Weighted filling invalid for type '%c' histogram %ld",
            histo->type,histo->ident);
         Warning(message);
         return -1;
      }
   }
   else if ( histo->type != 'F' && histo->type != 'D' )
   {
      char message[1024];
      sprintf(message,"Don't know how to fill a type '%c' histogram.",
         histo->type);
      Warning(message);
      return -1;
   }

   if ( histo->nbins_2d <= 0 )
   {
      switch ( histo->type )
      {
         case 'I':
            return(fill_int_histogram(histo,(long)xvalue));
         case 'R':
            return(fill_real_histogram(histo,xvalue));
         case 'F':
         case 'D':
            return(fill_weighted_histogram(histo,xvalue,weight));
      }
   }
   else
   {
      switch ( histo->type )
      {
         case 'I':
            return(fill_2d_int_histogram(histo,(long)xvalue,(long)yvalue));
         case 'R':
            return(fill_2d_real_histogram(histo,xvalue,yvalue));
         case 'F':
         case 'D':
            return(fill_2d_weighted_histogram(histo,xvalue,yvalue,weight));
      }
   }
   return 0;
}

/* -------------------- fill_histogram_by_ident ---------------------- */
/**
 *  @short Fill any type of 1-D or 2-D histogram known by its ID number.
 *
 *  Generic histogram fill function that can be used for type 'I', 'R',
 *  'F', and 'D' histograms, although it is not recommended for type 'I'
 *  histograms, due to type conversions.
 *
 *  @param  id      Identifier number of the histogram.
 *  @param  xvalue  X posistion where an entry is to be added.
 *  @param  yvalue  Y posistion (ignored for 1-D histograms)
 *  @param  weight  The weight of that entry (must be 1.0 for
 *			 'I' and 'R' type histograms).
 *
 *  @return 0 (o.k.),  -1 (no histogram that can be filled)
 *
 */

int fill_histogram_by_ident (long id, double xvalue, double yvalue,
   double weight)
{
   /* This static variable is not considered problematic with multi-threading */
   static int warned = 0;

   if ( hash_table == (HISTOGRAM **) NULL )
      if ( !warned ) /* Warning is displayed only once */
      {
         Warning("Inefficient histogram filling due to missing hashing table.");
         warned = 1;
      }

   return(fill_histogram(get_histogram_by_ident(id),xvalue,yvalue,weight));
}

/* ------------------------ histogram_matching ----------------------- */
/**
 *  @short Check if two histograms have exactly matching definitions
 *         (same type, dimension, size, ranges).
 *
 *  @param  histo1  pointer to first histogram
 *  @param  histo2  pointer to second histogram
 *
 *  @return 0 (not matching) or 1 (matching)
 */

int histogram_matching (HISTOGRAM *histo1, HISTOGRAM *histo2)
{
   if ( histo1->type != histo2->type || 
        histo1->nbins != histo2->nbins ||
        histo1->nbins_2d != histo2->nbins_2d ||
        (histo1->counts == NULL && histo2->counts != NULL) ||
        (histo1->counts != NULL && histo2->counts != NULL) ||
        (histo1->extension == NULL && histo2->extension !=NULL) ||
        (histo1->extension != NULL && histo2->extension == NULL) )
   {
      return 0;
   }

   if ( histo1->type == 'I' || histo1->type == 'i' )
   {
      if ( histo1->specific.integer.lower_limit != histo2->specific.integer.lower_limit ||
           histo1->specific.integer.upper_limit != histo2->specific.integer.upper_limit ||
           ( histo1->nbins_2d != 0 && 
              (histo1->specific_2d.integer.lower_limit != 
               histo2->specific_2d.integer.lower_limit ||
               histo1->specific_2d.integer.upper_limit != 
               histo2->specific_2d.integer.upper_limit) ) )
      {
         return 0;
      }
   }
   else
   {
      if ( histo1->specific.real.lower_limit != histo2->specific.real.lower_limit ||
           histo1->specific.real.upper_limit != histo2->specific.real.upper_limit ||
           ( histo1->nbins_2d != 0 && 
              (histo1->specific_2d.real.lower_limit != 
               histo2->specific_2d.real.lower_limit ||
               histo1->specific_2d.real.upper_limit != 
               histo2->specific_2d.real.upper_limit) ) )
      {
         return 0;
      }
   }

   return 1;
}

/* --------------------------- add_histogram ------------------------- */
/**
 *  @short Add a second histogram to a first one.
 *
 *  The histograms must exactly match in their definitions.
 *  The first histogram will be modified, the second is unchanged.
 *
 *  @param  histo1  pointer to first histogram
 *  @param  histo2  pointer to second histogram
 *
 *  @return NULL pointer indicates failure.
 *
 */

HISTOGRAM *add_histogram (HISTOGRAM *histo1, HISTOGRAM *histo2)
{
   int nbins, ibin;
   
   if ( histo1 == NULL || histo2 == NULL )
      return NULL;
_WAIT_IF_BUSY_(histo1)
_WAIT_IF_BUSY_(histo2)

   if ( ! histogram_matching(histo1,histo2) )
   {
      char msg[1024];
      sprintf(msg,"Histograms %ld and %ld not matching.",
         histo1->ident,histo2->ident);
      Warning(msg);
_CLEAR_BUSY_(histo2);
_CLEAR_BUSY_(histo1);
      return NULL;
   }

   if ( histo1->type == 'I' || histo1->type == 'i' )
   {
      histo1->specific.integer.sum += histo2->specific.integer.sum;
      histo1->specific.integer.tsum += histo2->specific.integer.tsum;
      if ( histo1->nbins_2d > 0 )
      {
         histo1->specific_2d.integer.sum += histo2->specific_2d.integer.sum;
         histo1->specific_2d.integer.tsum += histo2->specific_2d.integer.tsum;
      }
   }
   else
   {
      histo1->specific.real.sum += histo2->specific.real.sum;
      histo1->specific.real.tsum += histo2->specific.real.tsum;
      if ( histo1->nbins_2d > 0 )
      {
         histo1->specific_2d.real.sum += histo2->specific_2d.real.sum;
         histo1->specific_2d.real.tsum += histo2->specific_2d.real.tsum;
      }
   }
   nbins = histo1->nbins;
   if ( histo1->nbins_2d > 0 )
      nbins *= histo1->nbins_2d;
   histo1->underflow += histo2->underflow;
   histo1->overflow += histo2->overflow;
   histo1->underflow_2d += histo2->underflow_2d;
   histo1->overflow_2d += histo2->overflow_2d;
   if ( histo1->counts != NULL )
      for ( ibin=0; ibin<nbins; ibin++ )
         histo1->counts[ibin] += histo2->counts[ibin];
   if ( histo1->extension != NULL )
   {
      int j;
      histo1->extension->content_all += histo2->extension->content_all;
      histo1->extension->content_inside += histo2->extension->content_inside;
      for (j=0; j<8; j++)
         histo1->extension->content_outside[j] += histo2->extension->content_outside[j];
      if ( histo1->extension->fdata != NULL &&
           histo2->extension->fdata != NULL )
      {
         float *f1 = histo1->extension->fdata;
         float *f2 = histo2->extension->fdata;
         for ( ibin=0; ibin<nbins; ibin++ )
            f1[ibin] += f2[ibin];
      }
      if ( histo1->extension->ddata != NULL &&
           histo2->extension->ddata != NULL )
      {
         double *d1 = histo1->extension->ddata;
         double *d2 = histo2->extension->ddata;
         for ( ibin=0; ibin<nbins; ibin++ )
            d1[ibin] += d2[ibin];
      }
   }
   
_CLEAR_BUSY_(histo2);
_CLEAR_BUSY_(histo1);

   return histo1;
}

/* --------------------------- stat_histogram ------------------------- */
/**
 *  @short Statistical analysis of a histogram.
 *
 *  The median calculation is implemented for 1-D 'I' and 'R' types histograms only.
 *
 *  @param  histo  pointer to histogram
 *  @param  stbuf  pointer to histogram statistics structure
 *
 *  @return Nonzero result indicates failure
 *
 */

int stat_histogram (HISTOGRAM *histo, struct histstat *stbuf)
{
   int i, j;
   unsigned long hentries = 0;
   double dhentries = 0.;
   double sum = 0., sum2 = 0., hmean = 0.;
   double sum_2d = 0., sum2_2d = 0., hmean_2d = 0.;
   double x, lower_limit, upper_limit, step;
   double y, lower_limit_2d, upper_limit_2d, step_2d;
   int is_2d = 0;

   stbuf->mean = stbuf->tmean = stbuf->hmean =
      stbuf->median = stbuf->sigma = (double) 0.;

   if ( histo == (HISTOGRAM *) NULL )
      return(-1); /* Histogram has not been allocated */
_WAIT_IF_BUSY_(histo)
   if ( histo->nbins_2d > 0 )
      is_2d = 1;
   if ( histo->counts == (unsigned long *) NULL )
   {
_CLEAR_BUSY_(histo)
      return(-1); /* This can't happen for dynamically allocated histograms */
   }
   if ( histo->tentries == 0 )
   {
_CLEAR_BUSY_(histo)
      return(1);  /* No entries in histogram range */
   }
   if ( hentries == 0 )
   {
_CLEAR_BUSY_(histo)
      return(1);  /* That shouldn't happen (histogram corrupted) */
   }

   if ( histo->type == 'I' )
   {
      lower_limit = (double) histo->specific.integer.lower_limit;
      upper_limit = (double) histo->specific.integer.upper_limit;
   }
   else
   {
      lower_limit = (double) histo->specific.real.lower_limit;
      upper_limit = (double) histo->specific.real.upper_limit;
   }
   step = (upper_limit-lower_limit) / (double)histo->nbins;
   
   if ( is_2d )
   {
      sum_2d = sum2_2d = 0.;
      if ( histo->type == 'I' )
      {
         lower_limit_2d = (double) histo->specific.integer.lower_limit;
         upper_limit_2d = (double) histo->specific.integer.upper_limit;
      }
      else
      {
         lower_limit_2d = (double) histo->specific.real.lower_limit;
         upper_limit_2d = (double) histo->specific.real.upper_limit;
      }
      step_2d = (upper_limit_2d-lower_limit_2d) / (double)histo->nbins_2d;

      if ( histo->type == 'F' )
      {
         for ( j=0; j<histo->nbins_2d; j++ )
         {
            int io = j*histo->nbins;
            y = lower_limit_2d + (j+0.5)*step_2d;
            for ( i=0; i<histo->nbins; i++ )
            {
               int ie = io+i;
               dhentries += histo->extension->fdata[ie];
               x = lower_limit + (i+0.5)*step;
               sum += x*histo->extension->fdata[ie];
               sum2 += x*x*histo->extension->fdata[ie];
               sum_2d += y*histo->extension->fdata[ie];
               sum2_2d += y*y*histo->extension->fdata[ie];
            }
         }
      }
      else if ( histo->type == 'D' )
      {
         for ( j=0; j<histo->nbins_2d; j++ )
         {
            int io = j*histo->nbins;
            y = lower_limit_2d + (j+0.5)*step_2d;
            for ( i=0; i<histo->nbins; i++ )
            {
               int ie = io+i;
               dhentries += histo->extension->ddata[ie];
               x = lower_limit + (i+0.5)*step;
               sum += x*histo->extension->ddata[ie];
               sum2 += x*x*histo->extension->ddata[ie];
               sum_2d += y*histo->extension->ddata[ie];
               sum2_2d += y*y*histo->extension->ddata[ie];
            }
         }
      }
      else
      {
         for ( j=0; j<histo->nbins_2d; j++ )
         {
            int io = j*histo->nbins;
            y = lower_limit_2d + (j+0.5)*step_2d;
            for ( i=0; i<histo->nbins; i++ )
            {
               int ie = io+i;
               hentries += histo->counts[ie];
               x = lower_limit + (i+0.5)*step;
               sum += x*(double)histo->counts[ie];
               sum2 += x*x*(double)histo->counts[ie];
               sum_2d += y*(double)histo->counts[ie];
               sum2_2d += y*y*(double)histo->counts[ie];
            }
         }
      }
   }
   else
   {
      if ( histo->type == 'F' )
      {
         for ( i=0; i<histo->nbins; i++ )
         {
            dhentries += histo->extension->fdata[i];
            x = lower_limit + (i+0.5)*step;
            sum += x*histo->extension->fdata[i];
            sum2 += x*x*histo->extension->fdata[i];
         }
      }
      else if ( histo->type == 'D' )
      {
         for ( i=0; i<histo->nbins; i++ )
         {
            dhentries += histo->extension->ddata[i];
            x = lower_limit + (i+0.5)*step;
            sum += x*histo->extension->ddata[i];
            sum2 += x*x*histo->extension->ddata[i];
         }
      }
      else
      {
         for ( i=0; i<histo->nbins; i++ )
         {
            hentries += histo->counts[i];
            x = lower_limit + (i+0.5)*step;
            sum += x*(double)histo->counts[i];
            sum2 += x*x*(double)histo->counts[i];
         }
      }
   }
   
   if ( histo->type == 'F' || histo->type == 'D' )
   {
      if ( dhentries != 0. )
      {
         stbuf->hmean = (hmean=sum/dhentries);
         stbuf->hmean_2d = (hmean_2d=sum_2d/dhentries);
         if ( hentries > 1 )
         {
            double dhc = dhentries * (double) (hentries-1) / (double) hentries; 
            stbuf->sigma = sqrt((sum2-sum*hmean)/dhc);
            stbuf->sigma_2d = sqrt((sum2_2d-sum_2d*hmean_2d)/dhc);
         }
         else
            stbuf->sigma = stbuf->sigma_2d = 0.;
      }
      else
         stbuf->hmean = stbuf->sigma = stbuf->hmean_2d = stbuf->sigma_2d = 0.;
   }
   else
   {
      if ( hentries > 0 )
      {
         stbuf->hmean = (hmean=sum/(double)hentries);
         stbuf->hmean_2d = (hmean_2d=sum_2d/(double)hentries);
      }
      else
         stbuf->hmean = stbuf->hmean_2d = 0.;
      if ( hentries > 1 )
      {
         stbuf->sigma = sqrt((sum2-sum*hmean)/(double)(hentries-1));
         stbuf->sigma_2d = sqrt((sum2_2d-sum_2d*hmean_2d)/(double)(hentries-1));
      }
      else
         stbuf->sigma = stbuf->sigma_2d = 0.;
   }

   stbuf->mean = stbuf->tmean = stbuf->mean_2d = stbuf->tmean_2d = 0.;

   if ( histo->type == 'I' )
   {
      if ( histo->entries > 0 )
         stbuf->mean = (double)histo->specific.integer.sum /
            (double)histo->entries;
      if ( histo->tentries > 0 )
         stbuf->tmean = (double)histo->specific.integer.tsum /
            (double)histo->tentries;
      if ( is_2d )
      {
         if ( histo->entries > 0 )
            stbuf->mean_2d = (double)histo->specific_2d.integer.sum /
               (double)histo->entries;
         if ( histo->tentries > 0 )
            stbuf->tmean_2d = (double)histo->specific_2d.integer.tsum /
               (double)histo->tentries;
      }
   }
   else if ( histo->type == 'F' )
   {
      if ( histo->entries > 0 )
         stbuf->mean = histo->specific.real.sum /
            (double)histo->entries;
      if ( histo->tentries > 0 )
         stbuf->tmean = histo->specific.real.tsum /
            (double)histo->tentries;
      if ( is_2d )
      {
         if ( histo->entries > 0 )
            stbuf->mean_2d = histo->specific_2d.real.sum /
               (double)histo->entries;
         if ( histo->tentries > 0 )
            stbuf->tmean_2d = histo->specific_2d.real.tsum /
               (double)histo->tentries;
      }
   }
   else
   {
      if ( histo->extension->content_all != 0. )
         stbuf->mean = histo->specific.real.sum /
            histo->extension->content_all;
      if ( histo->extension->content_inside != 0. )
         stbuf->tmean = histo->specific.real.tsum /
            histo->extension->content_inside;
      if ( is_2d )
      {
         if ( histo->extension->content_all != 0. )
            stbuf->mean_2d = histo->specific_2d.real.sum /
               histo->extension->content_all;
         if ( histo->extension->content_inside != 0. )
            stbuf->tmean_2d = histo->specific_2d.real.tsum /
               histo->extension->content_inside;
      }
   }

_CLEAR_BUSY_(histo)

   /* After the unlocking because locate_histogram_fraction has its own locking */

   if ( is_2d || (histo->type != 'I' && histo->type != 'R') )
      stbuf->median_2d = stbuf->median = 0.; // Not implemented for 2-D or 'F'/'D' type.
   else
      stbuf->median = locate_histogram_fraction(histo,0.5);

   return(0);
}

/* ------------------- locate_histogram_fraction ---------------------- */
/**
 *  Locate point of arbitrary fraction of entries (quantile).
 *
 *  Locate the place in a 1-D histogram where a given fraction of the
 *  entries is to the 'left' of this place ('I' and 'R' type only).
 *
 *  @param  histo     Pointer to histogram
 *  @param  fraction  Fraction of entries to the left.
 *
 *  @return x-coordinate of given fraction or 0. for error.
 */

double locate_histogram_fraction (HISTOGRAM *histo, double fraction)
{
   int i, last;
   unsigned long hentries, mentries, lastcount;
   double lower_limit, upper_limit, step, centries, location;
   double a, b, c;

   if ( histo == (HISTOGRAM *) NULL )
      return 0.;

_WAIT_IF_BUSY_(histo)
   if ( histo->nbins_2d > 0 || fraction < 0. || fraction > 1. )
   {
_CLEAR_BUSY_(histo)
      return 0.;
   }

   hentries = 0;
   if ( histo->type == 'I' )
   {
      lower_limit = (double) histo->specific.integer.lower_limit;
      upper_limit = (double) histo->specific.integer.upper_limit;
      step = (histo->specific.integer.upper_limit-
              histo->specific.integer.lower_limit) / (double)histo->nbins;
   }
   else
   {
      lower_limit = (double) histo->specific.real.lower_limit;
      upper_limit = (double) histo->specific.real.upper_limit;
      step = (histo->specific.real.upper_limit-
              histo->specific.real.lower_limit) / (double)histo->nbins;
   }

   for ( i=0; i<histo->nbins; i++ )
      hentries += histo->counts[i];

   hentries += histo->overflow + histo->underflow;
   if ( (centries = hentries*fraction) == 0. )
   {
_CLEAR_BUSY_(histo)
      return(lower_limit);
   }
   lastcount = 0;
   if ( (double) histo->underflow >= centries )
      location = lower_limit;
   else
   {
      mentries = lastcount = histo->underflow;
      last = -1;
      location = upper_limit;
      for ( i=0; i<histo->nbins; i++ )
      {
         if ( (double)mentries+histo->counts[i] >= centries )
         {
#ifdef SIMPLE_MEDIAN_METHOD
            location = (double) (lower_limit +
                step * (0.5 + last +
                  (i-last)*(double)(hentries-2.0*mentries)/
                  (double)(lastcount+histo->counts[i])));
#else
            a = ((double)histo->counts[i]-lastcount)/(double)(i-last);
            b = lastcount;
            c = (double) mentries - (double) centries;
            if ( a == 0. || (i==0 && histo->underflow==0) )
               location = lower_limit + step * (last + 1.0 +
                  (centries-mentries)/(double)histo->counts[i]);
            else
               location = lower_limit + step * (last + 1.0 +
                  (sqrt(fabs(b*b-4*a*c))-b)/(2*a));
#endif
            break;
         }
         if ( histo->counts[i] != (size_t) -1 )
         {
            last = i;
            lastcount = histo->counts[i];
            mentries += histo->counts[i];
         }
      }
   }

_CLEAR_BUSY_(histo)
   return(location);
}

/* ------------------------ fast_stat_histogram ----------------------- */
/**
 *  @short Fast and basic histogram statistics.
 *
 *  Compute mean and truncated mean for histogram.
 *  For this kind of histogram analysis actually no histogram
 *  is required. A 'moments' structure would be sufficient.
 *
 *  @param  histo  pointer to histogram (1-D)
 *  @param  stbuf  pointer to histogram statistics structure
 *
 *  @return Nonzero result indicates failure
 */

int fast_stat_histogram (HISTOGRAM *histo, struct histstat *stbuf)
{
   stbuf->mean = stbuf->tmean = stbuf->hmean =
      stbuf->median = stbuf->sigma = (double) 0.;

   if ( histo == (HISTOGRAM *) NULL )
      return(-1); /* Histogram has not been allocated */
   if ( histo->tentries == 0 )
      return(1);  /* No entries in histogram range */
   if ( histo->type != 'I' && histo->type != 'R' &&
        histo->type != 'F' && histo->type != 'D' )
      return(-1); /* Histogram has neither integer nor real type */
   if ( histo->nbins_2d > 0 )
      return(-1); /* No statistics accumulated for 2-D histograms. */

_WAIT_IF_BUSY_(histo)
   if ( histo->type == 'I' )
   {
      if ( histo->entries > 0 )
         stbuf->mean = (double) (histo->specific.integer.sum /
            (double)histo->entries);
      if ( histo->tentries > 0 )
         stbuf->tmean = (double) (histo->specific.integer.tsum /
               (double)histo->tentries);
   }
   else if ( histo->type == 'R' )
   {
      if ( histo->entries > 0 )
         stbuf->mean = (double) (histo->specific.real.sum /
            (double)histo->entries);
      if ( histo->tentries > 0 )
         stbuf->tmean = (double) (histo->specific.real.tsum /
               (double)histo->tentries);
   }
   else if (histo->extension != (struct Histogram_Extension *) NULL)
   {
      if ( (histo->extension)->content_all != 0. )
         stbuf->mean = histo->specific.real.sum /
            (histo->extension)->content_all;
      if ( (histo->extension)->content_inside != 0. )
         stbuf->mean = histo->specific.real.tsum /
            (histo->extension)->content_inside;
   }

_CLEAR_BUSY_(histo)
   return(0);
}

static CONST_QUAL int zero = 0;
#define HistOutput(a) do { if ( histogram_file == (FILE *) NULL ) \
        Output(a); \
     else \
        fputs(a,histogram_file); } while(zero)

/* --------------------- print_histogram -------------------- */
/**
 *  @short Print contents of a histogram on the terminal.
 *
 *  Showing the actual content of each bin.
 *
 *  @param  histo  Pointer to histogram
 *
 *  @return  (none)
 */

void print_histogram (HISTOGRAM *histo)
{
   double xmin, xmax;
   double ymin, ymax;
   double dx, dy;
   int nx=0, ny=1, iy, ix;

   if ( histo == (HISTOGRAM *) NULL )
      return; /* Histogram has not been allocated */
   if ( histo->nbins <= 0 )
      return;
   if ( histo->type == 'F' || histo->type == 'D' )
   {
      if (histo->extension == (struct Histogram_Extension *) NULL)
         return;
      if ( (histo->type=='F' && (histo->extension)->fdata==(float*)NULL) ||
           (histo->type=='D' && (histo->extension)->ddata==(double*)NULL) )
         return;
   }
   else if ( histo->counts == (unsigned long *) NULL )
      return; /* This can't happen for dynamically allocated histograms */
   if ( histo->tentries == 0 )
      return;  /* No entries in histogram range */
   if ( histo->type != 'I' && histo->type != 'R' &&
        histo->type != 'i' && histo->type != 'r' &&
        histo->type != 'F' && histo->type != 'D' )
      return; /* Histogram has neither integer nor real type */

   if ( histo->type == 'I' || histo->type == 'i' )
   {
      xmin = (double) histo->specific.integer.lower_limit;
      xmax = (double) histo->specific.integer.upper_limit;
   }
   else
   {
      xmin = (double) histo->specific.real.lower_limit;
      xmax = (double) histo->specific.real.upper_limit;
   }
   nx = histo->nbins;
   dx = 0.;
   if ( nx > 0 )
      dx = (xmax-xmin)/nx;

_HLOCK_
   if ( histo->ident > 0 )
      printf("\n\n# ID=%ld: ",histo->ident);
   else
      printf("\n\n");
   if ( histo->title != (char *) NULL )
   {
      printf("%s\n",histo->title);
   }

   if ( (ny = histo->nbins_2d) > 0 )
   {
      if ( histo->type == 'I' || histo->type == 'i' )
      {
         ymin = (double) histo->specific_2d.integer.lower_limit;
         ymax = (double) histo->specific_2d.integer.upper_limit;
      }
      else
      {
         ymin = (double) histo->specific_2d.real.lower_limit;
         ymax = (double) histo->specific_2d.real.upper_limit;
      }
      dy = 0.;
      if ( ny > 0 )
         dy = (ymax-ymin)/ny;
      printf("# X = %g to %g, Y = %g to %g, E=%lu, U=%lu, O=%lu\n",
         xmin,xmax,ymin,ymax,histo->entries,histo->underflow,histo->overflow);
      if ( histo->type == 'F' || histo->type == 'D' )
      {
         printf("# Content: %12.7g total (%12.7g inside range)\n",
            (histo->extension)->content_all, (histo->extension)->content_inside);
      }
      printf("\n");
      for ( iy = 0; iy < ny; iy++ )
      {
         printf("\n");
         for ( ix = 0; ix < nx; ix++ )
         {
            if ( histo->type != 'F' && histo->type != 'D' )
               printf("   %d %d %ld\n", ix, iy, histo->counts[iy*nx+ix]);
            else
            {
               double z = 0.;
               if ( histo->type == 'F' )
                  z = (histo->extension)->fdata[iy*nx+ix];
               else
                  z = (histo->extension)->ddata[iy*nx+ix];
               printf("   %g\t%g\t%g\n", xmin+dx*(ix+0.5), ymin+dy*(iy+0.5), z);
            }
         }
      }
   }
   else
   {
      printf("# X = %g to %g, E=%lu, U=%lu, O=%lu\n",
         xmin,xmax,histo->entries,histo->underflow,histo->overflow);
      if ( histo->type == 'F' || histo->type == 'D' )
      {
         printf("# Content: %12.7g total (%12.7g inside range)\n",
            (histo->extension)->content_all, (histo->extension)->content_inside);
      }
      printf("\n");
      for ( ix = 0; ix < nx; ix++ )
      {
         if ( histo->type != 'F' && histo->type != 'D' )
            printf("   %d %ld\n", ix, histo->counts[ix]);
         else
         {
            double z = 0.;
            if ( histo->type == 'F' )
               z = (histo->extension)->fdata[ix];
            else
               z = (histo->extension)->ddata[ix];
            printf("   %g\t%g\n", xmin+dx*(ix+0.5), z);
         }
      }
   }

_HUNLOCK_
}

/* --------------------- display_histogram -------------------- */
/**
 *  @short Display contents of a histogram on the terminal.
 *
 *  This is a simple 'HPRINT' type display on one screen.
 *
 *  @param  histo  Pointer to histogram
 *
 *  @return  (none)
 */

void display_histogram (HISTOGRAM *histo)
{
   double xmin, xmax;
   double ymin, ymax;
   int kbins, lbins;
   int i, j;
   double *array;
   char message[1024];

   if ( histo == (HISTOGRAM *) NULL )
      return; /* Histogram has not been allocated */
   if ( histo->nbins <= 0 )
      return;
   if ( histo->type == 'F' || histo->type == 'D' )
   {
      if (histo->extension == (struct Histogram_Extension *) NULL)
         return;
      if ( (histo->type=='F' && (histo->extension)->fdata==(float*)NULL) ||
           (histo->type=='D' && (histo->extension)->ddata==(double*)NULL) )
         return;
   }
   else if ( histo->counts == (unsigned long *) NULL )
      return; /* This can't happen for dynamically allocated histograms */
   if ( histo->tentries == 0 )
      return;  /* No entries in histogram range */
   if ( histo->type != 'I' && histo->type != 'R' &&
        histo->type != 'i' && histo->type != 'r' &&
        histo->type != 'F' && histo->type != 'D' )
      return; /* Histogram has neither integer nor real type */

_HLOCK_
   if ( histo->ident > 0 )
      sprintf(message,"\n\nID=%ld: ",histo->ident);
   else
      sprintf(message,"\n\n");
   if ( histo->title != (char *) NULL )
   {
      strcat(message,histo->title);
      strcat(message,"\n");
   }
   HistOutput(message);

   if ( histo->nbins_2d > 0 )
   {
      display_2d_histogram(histo);
_HUNLOCK_
      return;
   }

   if ( histo->type == 'I' || histo->type == 'i' )
   {
      xmin = (double) histo->specific.integer.lower_limit;
      xmax = (double) histo->specific.integer.upper_limit;
   }
   else
   {
      xmin = (double) histo->specific.real.lower_limit;
      xmax = (double) histo->specific.real.upper_limit;
   }

   kbins = histo->nbins;
   lbins = 1;
   while ( kbins>78 )
   {
      kbins = (histo->nbins+lbins-1) / lbins;
      lbins++;
   }
      /* { kbins = (kbins+1)/2; lbins *= 2; } */
   if ( (array = (double *) malloc((size_t)(sizeof(double)*kbins))) ==
        (double *) NULL )
   {
_HUNLOCK_
      return;
   }
   for (i=0; i<kbins; i++)
      array[i] = 0;
   if ( histo->type == 'F' )
      for (i=0; i<histo->nbins; i++)
         array[i/lbins] += (histo->extension)->fdata[i];
   else if ( histo->type == 'D' )
      for (i=0; i<histo->nbins; i++)
         array[i/lbins] += (histo->extension)->ddata[i];
   else
      for (i=0; i<histo->nbins; i++)
         array[i/lbins] += histo->counts[i];
   ymin = ymax = array[0];
   for (i=1; i<kbins; i++)
   {
      if ( array[i] < ymin )
         ymin = array[i];
      else if ( array[i] > ymax )
         ymax = array[i];
   }
   for (i=0; i<kbins; i++)
      array[i] = 19.99*(array[i]-ymin)/(ymax-ymin+1.e-20);

   sprintf(message,"X = %g to %g, Y = %g to %g, E=%lu, U=%lu, O=%lu\n",
      xmin,xmax,ymin,ymax,histo->entries,histo->underflow,histo->overflow);
   HistOutput(message);
   if ( histo->type == 'F' || histo->type == 'D' )
   {
      sprintf(message,"Content: %12.7g total (%12.7g inside range)\n",
         (histo->extension)->content_all, (histo->extension)->content_inside);
      HistOutput(message);
   }
   for (j=20; j>=0; j--)
   {
      if ( j==20 )
         HistOutput("^");
      else
         HistOutput("|");
      for (i=0; i<kbins; i++)
         if ( array[i] > j )
         {
            if ( array[i]-j >= 0.5 )
               HistOutput("X");
            else if ( array[i]-j >= 0.2 )
               HistOutput("x");
            else
               HistOutput(".");
         }
         else
            HistOutput(" ");
      HistOutput("\n");
   }
   for (i=0; i<kbins; i++)
      HistOutput("-");
   HistOutput(">\n");
   free((void*)array);
_HUNLOCK_
}

/* ------------------- display_2d_histogram ------------------- */
/**
 *  @short Display contents of a 2D histogram. Called by display_histogram().
 *
 *  The histogram has already been checked by display_histogram()
 *  and its title has been printed.
 *
 *  @param  histo -- Pointer to histogram
 *
 *  @return (none)
 */

static void display_2d_histogram (HISTOGRAM *histo)
{
   double xmin, xmax, ymin, ymax;
   double zmin, zmax, z;
   int icount, ncounts, ix, iy, iz;
   char message[1024];

   if ( histo == (HISTOGRAM *) NULL )
      return;
   if ( histo->nbins <= 0 || histo->nbins_2d <= 0 )
      return;
   if ( histo->type == 'F' || histo->type == 'D' )
   {
      struct Histogram_Extension *he;
      if ( (he = histo->extension) == (struct Histogram_Extension *) NULL )
         return;
      if ( (histo->type == 'F' && he->fdata == (float *) NULL) ||
           (histo->type == 'D' && he->ddata == (double *) NULL) )
         return;
   }
   else if ( histo->counts == (unsigned long *) NULL )
      return;

_HLOCK_
   if ( histo->type == 'I' || histo->type == 'i' )
   {
      xmin = (double) histo->specific.integer.lower_limit;
      xmax = (double) histo->specific.integer.upper_limit;
      ymin = (double) histo->specific_2d.integer.lower_limit;
      ymax = (double) histo->specific_2d.integer.upper_limit;
   }
   else
   {
      xmin = (double) histo->specific.real.lower_limit;
      xmax = (double) histo->specific.real.upper_limit;
      ymin = (double) histo->specific_2d.real.lower_limit;
      ymax = (double) histo->specific_2d.real.upper_limit;
   }

   ncounts = histo->nbins*histo->nbins_2d;
   if ( histo->type == 'F' )
   {
      zmin = zmax = (histo->extension)->fdata[0];
      for ( icount=1; icount<ncounts; icount++ )
      {
         if ( (double)(histo->extension)->fdata[icount] < zmin )
            zmin = (histo->extension)->fdata[icount];
         else if ( (double)(histo->extension)->fdata[icount] > zmax )
            zmax = (histo->extension)->fdata[icount];
      }
   }
   else if ( histo->type == 'D' )
   {
      zmin = zmax = (histo->extension)->ddata[0];
      for ( icount=1; icount<ncounts; icount++ )
      {
         if ( (histo->extension)->ddata[icount] < zmin )
            zmin = (histo->extension)->ddata[icount];
         else if ( (histo->extension)->ddata[icount] > zmax )
            zmax = (histo->extension)->ddata[icount];
      }
   }
   else
   {
      zmin = zmax = histo->counts[0];
      for ( icount=1; icount<ncounts; icount++ )
      {
         if ( (double)histo->counts[icount] < zmin )
            zmin = histo->counts[icount];
         else if ( (double)histo->counts[icount] > zmax )
            zmax = histo->counts[icount];
      }
   }
   sprintf(message,"\nX = %g to %g, Y = %g to %g, Z = %g to %g\n",
      xmin,xmax,ymin,ymax,zmin,zmax);
   HistOutput(message);
   sprintf(message,"Nx*Ny=%d*%d, E=%lu, U=%lu/%lu, O=%lu/%lu\n",
      histo->nbins, histo->nbins_2d,
      histo->entries,histo->underflow,histo->underflow_2d,
      histo->overflow,histo->overflow_2d);
   HistOutput(message);
   if ( histo->type == 'F' || histo->type == 'D' )
   {
      sprintf(message,"Content: %12.7g total (%12.7g inside range)\n",
         (histo->extension)->content_all, (histo->extension)->content_inside);
      HistOutput(message);
   }
   if ( zmin == zmax )
   {
      HistOutput("All bins have the same contents.\n\n");
_HUNLOCK_
      return;
   }
   HistOutput("^\n");
   for (iy=histo->nbins_2d-1; iy>=0; iy--)
   {
      HistOutput("|");
      for (ix=0; ix<histo->nbins; ix++)
      {
         if ( histo->type == 'F' )
            z = (histo->extension)->fdata[iy*histo->nbins+ix];
         else if ( histo->type == 'D' )
            z = (histo->extension)->ddata[iy*histo->nbins+ix];
         else
            z = histo->counts[iy*histo->nbins+ix];
         iz = (int) ((27*(z-zmin))/(zmax-zmin));
         if ( iz > 0 )
         {
            sprintf(message,"%c",'A'-1+ /*min(iz,26)*/ (iz<26?iz:26));
            HistOutput(message);
         }
         else if ( z>zmin )
            HistOutput(".");
         else
            HistOutput(" ");
      }
      HistOutput("\n");
   }
   for (ix=0; ix<histo->nbins; ix++)
      HistOutput("-");
   HistOutput(">\n");
_HUNLOCK_
}

#undef HistOutput

/* ------------------ display_all_histograms ------------------------ */
/**
 *  Display all histograms in list of histograms
 *
 *  Arguments: none
 *
 *  Return value: none
 */

void display_all_histograms ()
{
   HISTOGRAM *thisto;

_HLOCK_
   for ( thisto = first_histogram;
         thisto != (HISTOGRAM *) NULL;
         thisto = thisto->next )
   {
_HUNLOCK_ /* Must not nest locking */
      display_histogram(thisto);
_HLOCK_
   }
_HUNLOCK_
}

/* ---------------------- histogram_to_lookup ----------------------- */
/**
 *  Convert a histogram to a lookup table by integrating the histogram.
 *
 *  @param  histo   input histogram
 *  @param  lookup  output lookup table
 *
 *  @return 0 if ok or -1 for failure
 */

int histogram_to_lookup (HISTOGRAM *histo, HISTOGRAM *lookup)
{
   int i;

   if ( histo == (HISTOGRAM *) NULL || lookup == (HISTOGRAM *) NULL ||
        histo == lookup )
      return(-1);
_WAIT_IF_BUSY_(histo)
_WAIT_IF_BUSY_(lookup)
   if ( histo->nbins != lookup->nbins || histo->nbins <= 0 || histo->nbins_2d > 0 )
   {
_CLEAR_BUSY_(lookup)
_CLEAR_BUSY_(histo)
      return(-1);
   }
   if ( lookup->type == 'I' )
      lookup->type = 'i';
   else if ( lookup->type == 'R' )
      lookup->type = 'r';
   if ( (lookup->type == 'i' && histo->type != 'I') ||
        (lookup->type == 'r' && histo->type != 'R') )
   {
_CLEAR_BUSY_(lookup)
_CLEAR_BUSY_(histo)
      return(-1);
   }
   clear_histogram(lookup);
   lookup->counts[0] = 0;
   for ( i=1; i<lookup->nbins; i++ )
      lookup->counts[i] = lookup->counts[i-1] + histo->counts[i-1];
   lookup->entries = histo->entries;
   lookup->tentries = histo->tentries;
   lookup->underflow = histo->underflow;
   lookup->overflow = histo->overflow;

_CLEAR_BUSY_(lookup)
_CLEAR_BUSY_(histo)
   return(0);
}

/* ---------------------------- lookup_int -------------------------- */
/**
 *  Look up a table created from an integer histogram.
 *
 *  @param  lookup  the lookup table
 *  @param  value   the value at which to look up
 *  @param  factor  the scaling factor of the lookup result or 0
 *
 *  @return If 'value' is inside the range of the lookup table
 *    (that is the range of the histogram from which the lookup table was
 *    created), a value between 0 and 'factor' (or the number of entries
 *    in the range, if factor==0) is returned.
 */

long lookup_int (HISTOGRAM *lookup, long value, long factor)
{
   long indx, interpol;
   long look1, look2;

   if ( lookup == (HISTOGRAM *) NULL )
      return(-1);
   if ( lookup->type != 'i' )
      return -1;
   if ( lookup->nbins <= 0 || lookup->tentries <= 0 || lookup->nbins_2d > 0 )
      return(-1);
   if ( factor == 0L )
      factor = lookup->tentries;
   if ( value < lookup->specific.integer.lower_limit )
      return(-1);
   if ( value >= lookup->specific.integer.upper_limit )
      return(factor+1);

   indx = (1000*lookup->nbins*(value-lookup->specific.integer.lower_limit)) /
           lookup->specific.integer.width;
   interpol = indx - (indx/1000)*1000;
   indx /= 1000;

   if ( indx <0 || indx >= lookup->nbins )
      return(-1);

   look1 = lookup->counts[indx];
   if ( indx < lookup->nbins-1 )
      look2 = lookup->counts[indx+1];
   else
      look2 = lookup->tentries;

   look1 = (look1*(1000-interpol)+look2*interpol)/1000;
   return(look1*factor/lookup->tentries);
}

/* --------------------------- lookup_real -------------------------- */
/**
 *  Look up a table created from an 'real' histogram.
 *
 *  @param  lookup  the lookup table
 *  @param  value   the value at which to look up
 *  @param  factor  the scaling factor of the lookup result or 0
 *
 *  @return If 'value' is inside the range of the lookup table
 *    (that is the range of the histogram from which the lookup table was
 *    created), a value between 0 and 'factor' (or the number of entries
 *    in the range, if factor==0) is returned.
 */

double lookup_real (HISTOGRAM *lookup, double value, double factor)
{
   long indx;
   double indxr, interpol;
   long look1, look2;

   if ( lookup == (HISTOGRAM *) NULL )
      return(-1.);
   if ( lookup->type != 'r' )
      return -1.;
   if ( lookup->nbins <= 0 || lookup->tentries <= 0 || lookup->nbins_2d > 0 )
      return(-1.);
   if ( factor == 0L )
      factor = (double) lookup->tentries;
   if ( value < lookup->specific.real.lower_limit )
      return(0.);
   if ( value >= lookup->specific.real.upper_limit )
      return(factor);

   indxr = (value-lookup->specific.real.lower_limit) *
         lookup->specific.real.inverse_binwidth;
   indx = (long) indxr;

   interpol = indxr - indx;

   if ( indx <0 || indx >= lookup->nbins )
      return(-1.);

   look1 = lookup->counts[indx];
   if ( indx < lookup->nbins-1 )
      look2 = lookup->counts[indx+1];
   else
      look2 = lookup->tentries;

   look1 = look1*(1.-interpol) + look2*interpol;
   return(look1*factor/lookup->tentries);
}

/* --------------------- histogram_hashing --------------------- */
/**
 *  Turn hashing of histograms (using their ident as key) on or off.
 *
 *  @param  tabsize   Minimum number of elements in
 *		      hashing table or 0 if hash table
 *		      should be released (max: 15000).
 *
 *  @return  0 (o.k.),  -1 (error)
 */

int histogram_hashing (int tabsize)
{
   size_t i;
   HISTOGRAM *histo;

_HLOCK_
   if ( tabsize != 0 )
   {
      /* No dynamical change of hash table size: do nothing if */
      /* hashing is already active. */
      if ( hash_table != (HISTOGRAM **) NULL )
      {
_HUNLOCK_
         return 0;
      }
      for (i=0; i < sizeof(primetab)/sizeof(primetab[0])-1; i++ )
         if ( (int)primetab[i] >= tabsize )
            break;
      hash_size = primetab[i];
      if ( (hash_table = (HISTOGRAM **) 
             calloc((size_t)hash_size,sizeof(HISTOGRAM *)))
           == (HISTOGRAM **) NULL )
      {
         Warning("Too little memory for histogram hashing");
         hash_size = 0;
_HUNLOCK_
         return -1;
      }
      /* Initialize the hashing table with all known histograms. */
      for ( histo = first_histogram; histo != (HISTOGRAM *) NULL;
            histo=histo->next )
         if ( histo->ident > 0 )
            hash_table[histo->ident % hash_size] = histo;
   }
   else
   {
      if ( hash_table != (HISTOGRAM **) NULL )
         free((void *) hash_table);
      hash_size = 0;
   }
_HUNLOCK_
   return 0;
}
