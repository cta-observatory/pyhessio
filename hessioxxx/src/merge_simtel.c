/* ============================================================================

Copyright (C) 2013, 2014, 2015  Konrad Bernloehr

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

/** @file merge_simtel.c
 *  @short A program for merging events from separate telescope simulations
 *     of the same showers.
 *
 *  The program will read sim_telarray raw or DST data on two input files,
 *  map telescope ID according to a mapping file and write the merged
 *  blocks to an output file.
 *
 *  Inputs expected - and the action to be performed:
 *     Type
 *   Once per run:
 *       70 (history)    -  Write as-is, impossible to merge
 *     2000 (run_header) -  Merging needed for telescope list and positions
 *     2001 (MC run header)  -  Only one of two MC run-headers needed (should be identical)
 *     1212 (input config = CORSIKA inputs) -  Only one needed (should be identical, duplicate)
 *   Once per telescope (and per run for raw & DST levels 0-2; just once for DST level 3):
 *     2002 (camera settings) - Write after mapping of telescope ID (if mapped)
 *     2003 (camera organization)  - Write after mapping of telescope ID (if mapped)
 *     2004 (pixel settings)  - Write after mapping of telescope ID (if mapped)
 *     2005 (pixel disable)  - Write after mapping of telescope ID (if mapped)
 *     2006 (camera software settings)  - Write after mapping of telescope ID (if mapped)
 *     2008 (tracking settings)   - Write after mapping of telescope ID (if mapped)
 *     2007 (pointing corrections)  - Write after mapping of telescope ID (if mapped)
 *     2022 (telescope monitoring)   - Write after mapping of telescope ID (if mapped)
 *     2023 (Laser calibration)  - Write after mapping of telescope ID (if mapped)
 *   Per shower:
 *    once:
 *     2020 (MC shower) -  Only one of two MC run-headers needed (should be identical)
 *    per array:
 *     2021 (MC event) - Only one of two blocks needed (anything to get merged?)
 *   Optional per event; not immediately written but delayed until next MC etc. block:
 *     2026 (MC pe sum) - ???
 *     1204 (photo-electrons individually) - ???
 *     2010 (event) - Needs remapping and merging at all levels
 *   At end of run:
 *     2024 (run statistics - usually not present)
 *     2025 (MC run statistics - usually not present)
 *      100 (histograms) - Cannot be merged properly. Histograms of generated showers
 *            should agree, but for triggered showers we cannot tell how many are common.
 *
 *   FIXME: Ignoring 'trgmask' files initially - include them later on.
 *
@verbatim
Syntax:  merge_simtel [ options ] map-file input1 input2 output
Options:
     --auto-trgmask  : Load trgmask.gz files for each input file where available.
     --min-trg-tel n : Require at least n telescopes in merged event (default: 2).
     --verbose       : Show events being merged.
@endverbatim
 *
 *  @author  Konrad Bernloehr
 *  @date    @verbatim CVS $Date: 2015/05/31 13:02:40 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.6 $ @endverbatim
 */

/** @defgroup merge_simtel_c The merge_simtel program */
/** @{ */

#include "initial.h"      /* This file includes others as required. */
#include "io_basic.h"     /* This file includes others as required. */
#include "mc_tel.h"
#include "history.h"
#include "io_hess.h"
#include "histogram.h"
#include "io_histogram.h"
#include "fileopen.h"
#include "straux.h"
#include "warning.h"
#include "io_trgmask.h"
#include "eventio_version.h"

#ifndef _UNUSED_
# ifdef __GNUC__
#  define _UNUSED_ __attribute__((unused))
# else
#  define _UNUSED_
# endif
#endif
#include <signal.h>

void stop_signal_function (int isig);

static int interrupted;
static int verbose = 0;

/* ---------------------- stop_signal_function -------------------- */
/**
 *  Stop the program gracefully when it catches an INT or TERM signal.
 *
 *  @param isig  Signal number.
 *
 *  @return (none)
 */

void stop_signal_function (int isig)
{
   if ( isig >= 0 )
   {
      fprintf(stderr,"Received signal %d\n",isig);
   }
   if ( !interrupted )
      fprintf(stderr,
      "Program stop signal. Stand by until current data block is finished.\n");

   interrupted = 1;
   
   /* Reset the signal handler: pressing Ctrl-C again will terminate the program */
   signal(SIGINT,SIG_DFL);
   signal(SIGTERM,SIG_DFL);
}


/* --------------------------- syntax -------------------------------- */

/** Show program syntax */

static void syntax (const char *program);

static void syntax (const char *program)
{
   printf("Syntax: %s [ options ] map-file input1 input2 output\n",program);
   printf("Options:\n");
   printf("     --auto-trgmask  : Load trgmask.gz files for each input file where available.\n");
   printf("     --min-trg-tel n : Require at least n telescopes in merged event (default: 2).\n");
   printf("     --verbose       : Show events being merged.\n");
   printf("\nCompiled for a maximum of %d telescopes before and after merging.\n", H_MAX_TEL);
   printf("Linked against eventIO/hessio library\n");
#ifdef EVENTIO_VERSION
         printf("   version %s\n", EVENTIO_VERSION);
#endif
#ifdef EVENTIO_RELEASE
         printf("   release %s\n", EVENTIO_RELEASE);
#else
# ifdef EVENTIO_BASE_RELEASE
         printf("   based on release %s\n", EVENTIO_BASE_RELEASE);
# endif
#endif
   show_hessio_max();
   H_CHECK_MAX();
   exit(1);
}

/**
 *  Structure with per output telescope information keeping track of prerequisites.
 */

struct map_tel_struct
{
   int tel_id;              ///< Telescope ID on output
   int ifn;                 ///< Input file number (1 or 2)
   int inp_id;              ///< Telescope ID on input
   int inp_itel;            ///< Sequential telescope count on input
   int have_camset;         ///< Have camera_settings for this telescope
   int have_camorg;         ///< Have camera organisation for this telescope
   int have_pixset;         ///< Have pixel settings for this telescope
   int have_pixdis;         ///< Have pixels disabled for this telescope (optional)
   int have_camsoft;        ///< Have camera software settings for this telescope
   int have_pointcor;       ///< Have pointing correction for this telescope
   int have_trackset;       ///< Have tracking settings for this telescope
};

struct map_tel_struct map_tel[H_MAX_TEL];


/** Mapping structures from input telescope ID to output telescope ID.
    Not mapped telescopes are defined by output telescope ID of -1. */
int map_to[2][H_MAX_TEL+1];   ///< The telescope ID to which a given input telescope ID should get mapped.

/** Mapping from telescope IDs to offsets in the data structures, first for input telescope IDs.
 *  We restrict the ID/index mapping here to well behaved cases (0<ID<=H_MAX_TEL).
 *  An index value of -1 indicates a non-existant/ignored telescope. */
int tel_idx[2][H_MAX_TEL+1];  ///< Where is a telescope of given ID in the input data structures?
/** Mapping from output telescope ID to offset in output data structures. */
int tel_idx_out[H_MAX_TEL+1]; ///< Where is a telescope of given ID in the output data structures?

int find_in_tel_idx(int tel_id, int ifile);
int find_out_tel_idx(int tel_id, int ifile);
int find_mapped_telescope (int tel_id, int ifile);
int ntel1, ntel2, ntel;
int nrtel1, nrtel2;
long event1 = -1, event2 = 0;
long ev_hess_event = 0, ev_pe_sum = 0; /**< For delayed writing */
int run1 = -1, run2 = -1;
int min_trg = 2;

/**
 *   Offset of an input telescope of given ID within the input structures.
 */

int find_in_tel_idx(int tel_id, int ifile)
{
   if ( tel_id >= 0 && tel_id <= H_MAX_TEL && ifile > 0 && ifile <= 2 )
      return tel_idx[ifile-1][tel_id];
   else
      return -1;
}

/**
 *   Offset of an input telescope of given ID within the output structures.
 */

int find_out_tel_idx(int tel_id, int ifile)
{
   int itel_in, itel_out, tel_id_out;
   if ( tel_id >= 0 && tel_id <= H_MAX_TEL && ifile > 0 && ifile <= 2 )
   {
      /* First check that the telescope is present in the input */
      itel_in = tel_idx[ifile-1][tel_id];
      if ( itel_in >= 0 && itel_in < H_MAX_TEL )
      {
         /* Then we must have a valid mapped telescope ID. */
         tel_id_out = map_to[ifile-1][tel_id];
         if ( tel_id_out >= 0 && tel_id_out <= H_MAX_TEL )
         {
            itel_out = tel_idx_out[tel_id_out];
            return itel_out;
         }
      }
   }
   return -1;
}

/**
 *   Mapping from telescope ID on input to telescope ID on output, with check.
 */

int find_mapped_telescope (int tel_id, int ifile)
{
   if ( tel_id >= 0 && tel_id <= H_MAX_TEL && ifile > 0 && ifile <= 2 )
      return map_to[ifile-1][tel_id];
   else
      return -1;
}

/* ----------------------------------------------------------------------- */
/** 
 *  Write an I/O block as-is to another file than foreseen for the I/O buffer.
 */

int write_io_block_to_file (IO_BUFFER *iobuf, FILE *f);

int write_io_block_to_file (IO_BUFFER *iobuf, FILE *f)
{
   int rc = 0;
   if ( iobuf != NULL && f != NULL )
   {
      FILE *t = iobuf->output_file;
      iobuf->output_file = f;
      rc = write_io_block(iobuf);
      iobuf->output_file = t;
   }
   return rc;
}

/* ----------------------------------------------------------------------- */

int check_for_delayed_write(IO_ITEM_HEADER *item_header, int ifile, AllHessData *hsdata_out, IO_BUFFER *iobuf_out);

int check_for_delayed_write(IO_ITEM_HEADER *item_header, int ifile, AllHessData *hsdata_out, IO_BUFFER *iobuf_out)
{
   int rc = 0;
   int type = item_header->type;
   long id  = item_header->ident;
   long event = 0;

   if ( type == 1204 ) /* ignore */
      return 0;
   if ( type == 2020 )
      event = id * 100; /* shower number */
   else if ( type == 2021 || type == 2026 /* || type == 1204 */ || type == 2010 )
      event = id; /* event number */
   else if ( type == 2000 || type == 2001 || type == 100 )
      event = -1;
   else
   {
#ifdef DEBUG_MERGE
printf("Unexpected check_for_delayed_write(): type: %d, id: %ld)\n", type, id);
#endif
      return 0;
   }
   
   if ( type != IO_TYPE_HESS_EVENT && type != IO_TYPE_HESS_MC_PE_SUM )
   {
      if ( ev_pe_sum > 0 && ev_pe_sum < event )
      {
#ifdef DEBUG_MERGE
printf("Delayed writing of pe sums for event %ld (current data block type: %d, id: %ld)\n", ev_pe_sum, type, id);
#endif
         rc = write_hess_mc_pe_sum(iobuf_out,&hsdata_out->mc_event.mc_pesum);
         ev_pe_sum = 0;
         if ( rc != 0 )
            return rc;
      }
      if ( ev_hess_event > 0 && ev_hess_event < event )
      {
#ifdef DEBUG_MERGE
printf("Delayed writing of data event %ld (current data block type: %d, id: %ld)\n", ev_hess_event, type, id);
#endif
         if ( hsdata_out->event.num_teldata >= min_trg ||
              hsdata_out->event.central.num_teltrg >= min_trg ||
              hsdata_out->event.central.num_teldata >= min_trg )
         {
            set_tel_idx_ref(0);
            rc = write_hess_event(iobuf_out,&hsdata_out->event,-1);
         }
         ev_hess_event = 0;
      }
   }
   else if ( type == IO_TYPE_HESS_EVENT )
   {
      if ( ev_pe_sum > 0 && ev_pe_sum <= event ) /* p.e. sums for this or any preceding event */
      {
#ifdef DEBUG_MERGE
printf("Delayed writing of pe sums for event %ld (current data block type: %d, id: %ld)\n", ev_pe_sum, type, id);
#endif
         rc = write_hess_mc_pe_sum(iobuf_out,&hsdata_out->mc_event.mc_pesum);
         ev_pe_sum = 0;
         if ( rc != 0 )
            return rc;
      }
      if ( ev_hess_event > 0 && ev_hess_event < event ) /* event data from preceding event ? */
      {
#ifdef DEBUG_MERGE
printf("Delayed writing of data event %ld (current data block type: %d, id: %ld)\n", ev_hess_event, type, id);
#endif
         if ( hsdata_out->event.num_teldata >= min_trg ||
              hsdata_out->event.central.num_teltrg >= min_trg ||
              hsdata_out->event.central.num_teldata >= min_trg )
         {
            set_tel_idx_ref(0);
            rc = write_hess_event(iobuf_out,&hsdata_out->event,-1);
         }
         ev_hess_event = 0;
      }
   }
   else if ( type == IO_TYPE_HESS_MC_PE_SUM )
   {
      if ( ev_pe_sum > 0 && ev_pe_sum < event ) /* p.e. sums from preceding event ? */
      {
#ifdef DEBUG_MERGE
printf("Delayed writing of pe sums for event %ld (current data block type: %d, id: %ld)\n", ev_pe_sum, type, id);
#endif
         rc = write_hess_mc_pe_sum(iobuf_out,&hsdata_out->mc_event.mc_pesum);
         ev_pe_sum = 0;
         if ( rc != 0 )
            return rc;
      }
      if ( ev_hess_event > 0 && ev_hess_event < event ) /* event data from preceding event ? */
      {
#ifdef DEBUG_MERGE
printf("Delayed writing of data event %ld (current data block type: %d, id: %ld)\n", ev_hess_event, type, id);
#endif
         if ( hsdata_out->event.num_teldata >= min_trg ||
              hsdata_out->event.central.num_teltrg >= min_trg ||
              hsdata_out->event.central.num_teldata >= min_trg )
         {
            set_tel_idx_ref(0);
            rc = write_hess_event(iobuf_out,&hsdata_out->event,-1);
         }
         ev_hess_event = 0;
      }
   }
   return rc;
}

/* ----------------------------------------------------------------------- */

static struct trgmask_set *tms[2] = { NULL, NULL };
static struct trgmask_hash_set *ths[2] = { NULL, NULL };
static int events[2] = { 0, 0 };
static int mcshowers[2] = { 0, 0 };
static int mcevents[2] = { 0, 0 };
static int max_list = 999;

/**
 *  Processing and merging of I/O blocks from the two input files,
 *  hopefully presented in the right order.
 */

int merge_data_from_io_block (IO_BUFFER *iobuf, IO_ITEM_HEADER *item_header, 
   int ifile, AllHessData *hsdata, 
   AllHessData *hsdata_out, IO_BUFFER *iobuf_out );

int merge_data_from_io_block (IO_BUFFER *iobuf, IO_ITEM_HEADER *item_header, 
   int ifile, AllHessData *hsdata, 
   AllHessData *hsdata_out, IO_BUFFER *iobuf_out )
{
   int tel_id, itel;
   int tel_id3, itel3;
   int rc = -9;
   static int nlist = 0;

   if ( iobuf == NULL || item_header == NULL || hsdata == NULL || 
         hsdata_out == NULL || iobuf_out == NULL || ifile < 1 || ifile > 2 )
   {
      fprintf(stderr,"Invalid parameter in merge_data_from_io_block()\n");
      exit(1);
   }

#ifdef DEBUG_MERGE
printf("Now processing item type %d (ID=%ld) from file %d ...\n", 
  (int) item_header->type, item_header->ident, ifile);
#endif

   set_tel_idx_ref(ifile);

   switch ( (int) item_header->type )
   {
      /* =================================================== */
      case IO_TYPE_HESS_RUNHEADER: /* 2000 */
         check_for_delayed_write(item_header, ifile, hsdata_out, iobuf_out);
         /* Free memory allocated inside ... */
         for (itel=0; itel<hsdata->run_header.ntel; itel++)
         {
            if ( hsdata->event.teldata[itel].raw != NULL )
	    {
               free(hsdata->event.teldata[itel].raw);
	       hsdata->event.teldata[itel].raw = NULL;
	    }
            if ( hsdata->event.teldata[itel].pixtm != NULL )
	    {
               free(hsdata->event.teldata[itel].pixtm);
	       hsdata->event.teldata[itel].pixtm = NULL;
	    }
            if ( hsdata->event.teldata[itel].img != NULL )
	    {
               free(hsdata->event.teldata[itel].img);
	       hsdata->event.teldata[itel].img = NULL;
               hsdata->event.teldata[itel].num_image_sets = 0;
	    }
         }
         hsdata->run_header.ntel = 0;
         /* Clear absolutely everything else (Warning: problem with DST level 3!!!) */
         // memset(hsdata,0,sizeof(AllHessData));
         // memset(&hsdata->event,0,sizeof(FullEvent));
         memset(&hsdata->mc_event,0,sizeof(MCEvent));
         /* Clearing MCShower is a bit more tricky because of allocated shower profiles */
         if ( ifile == 1 )
         {
            for (itel3=0; itel3<hsdata_out->run_header.ntel; itel3++)
            {
               if ( hsdata_out->event.teldata[itel3].raw != NULL )
	       {
                  free(hsdata_out->event.teldata[itel3].raw);
	          hsdata_out->event.teldata[itel3].raw = NULL;
	       }
               if ( hsdata_out->event.teldata[itel3].pixtm != NULL )
	       {
                  free(hsdata_out->event.teldata[itel3].pixtm);
	          hsdata_out->event.teldata[itel3].pixtm = NULL;
	       }
               if ( hsdata_out->event.teldata[itel3].img != NULL )
	       {
                  free(hsdata_out->event.teldata[itel3].img);
	          hsdata_out->event.teldata[itel3].img = NULL;
                  hsdata_out->event.teldata[itel3].num_image_sets = 0;
	       }
               hsdata_out->run_header.tel_id[itel3] = -1;
               hsdata_out->run_header.tel_pos[itel3][0] = 0.;
               hsdata_out->run_header.tel_pos[itel3][1] = 0.;
               hsdata_out->run_header.tel_pos[itel3][2] = 0.;
            }
            // memset(&hsdata_out->event,0,sizeof(FullEvent)); /* ???? */
            memset(&hsdata_out->mc_event,0,sizeof(MCEvent)); /* ???? */
            hsdata_out->run_header.ntel = 0;
         }
         /* Read the run header data from the buffer */
         fflush(stdout);
         if ( (rc = read_hess_runheader(iobuf,&hsdata->run_header)) < 0 )
         {
            Warning("Reading run header failed.");
            exit(1);
         }
         else
         {
            printf("Input file %d starting now run %d with %d telescopes produced at time %ld\n", 
               ifile,
               hsdata->run_header.run, hsdata->run_header.ntel,
               (long) hsdata->run_header.time);
         }
         if ( ifile == 2 && hsdata->run_header.run != hsdata_out->run_header.run )
         {
            fflush(stdout);
            fprintf(stderr,"Run numbers do not match!\n");
            exit(1);
         }
         else if ( ifile == 1 )
            hsdata_out->run_header.run = hsdata->run_header.run;
         if ( hsdata->run_header.ntel > H_MAX_TEL )
         {
            fflush(stdout);
            fprintf(stderr,"Too many telescopes (%d, supporting only up to %d)!\n", 
               hsdata->run_header.ntel, H_MAX_TEL);
            exit(1);
         }
         if ( ifile == 1 )
         {
            hsdata_out->run_header.time = hsdata->run_header.time;
            hsdata_out->run_header.run_type = hsdata->run_header.run_type;
            hsdata_out->run_header.tracking_mode = hsdata->run_header.tracking_mode;
            hsdata_out->run_header.reverse_flag = hsdata->run_header.reverse_flag;
            hsdata_out->run_header.direction[0] = hsdata->run_header.direction[0];
            hsdata_out->run_header.direction[1] = hsdata->run_header.direction[1];
            hsdata_out->run_header.offset_fov[0] = hsdata->run_header.offset_fov[0];
            hsdata_out->run_header.offset_fov[1] = hsdata->run_header.offset_fov[1];
            hsdata_out->run_header.conv_depth = hsdata->run_header.conv_depth;
            hsdata_out->run_header.conv_ref_pos[0] = hsdata->run_header.conv_ref_pos[0];
            hsdata_out->run_header.conv_ref_pos[1] = hsdata->run_header.conv_ref_pos[1];
            nrtel1 = hsdata->run_header.ntel;
            hsdata_out->event.num_tel = hsdata_out->run_header.ntel = ntel;
            hsdata_out->run_header.min_tel_trig = hsdata->run_header.min_tel_trig;
            hsdata_out->run_header.duration = hsdata->run_header.duration;
            if ( hsdata_out->run_header.target != NULL )
               free(hsdata_out->run_header.target);
            hsdata_out->run_header.target = strdup(hsdata->run_header.target);
            if ( hsdata_out->run_header.observer != NULL )
               free(hsdata_out->run_header.observer);
            hsdata_out->run_header.observer = strdup(hsdata->run_header.observer);
            hsdata_out->run_header.max_len_target = hsdata->run_header.max_len_target;
            hsdata_out->run_header.max_len_observer = hsdata->run_header.max_len_observer;
            hsdata_out->event.num_tel = hsdata->run_header.ntel;
         }
         else
         {
            nrtel2 = hsdata->run_header.ntel;
            hsdata_out->event.num_tel = hsdata_out->run_header.ntel = ntel;
            rc = 0;
            if ( hsdata->run_header.time != hsdata_out->run_header.time )
               printf("Simulation in file at time %ld, in file at time %ld.\n",
                  hsdata_out->run_header.time, hsdata->run_header.time);
            if ( hsdata->run_header.run_type != hsdata_out->run_header.run_type )
            {
               printf("File 1 has run type %d, file 2 has run type %d.\n",
                  hsdata_out->run_header.run_type, hsdata->run_header.run_type);
               rc = 1;
            }
            if ( hsdata->run_header.tracking_mode != hsdata_out->run_header.tracking_mode )
            {
               printf("File 1 has tracking mode %d, file 2 has tracking mode %d.\n",
                  hsdata_out->run_header.tracking_mode, hsdata->run_header.tracking_mode);
               rc = 1;
            }
            if ( hsdata->run_header.reverse_flag != hsdata_out->run_header.reverse_flag )
            {
               printf("File 1 has reverse flag %d, file 2 has reverse flag %d.\n",
                  hsdata_out->run_header.reverse_flag, hsdata->run_header.reverse_flag);
               rc = 1;
            }
            if ( hsdata->run_header.direction[0] != hsdata_out->run_header.direction[0] ||
                 hsdata->run_header.direction[1] != hsdata_out->run_header.direction[1] )
            {
               printf("File 1 has direction %f, %f, file 2 has %f, %f.\n",
                  hsdata_out->run_header.direction[0], hsdata_out->run_header.direction[1], 
                  hsdata->run_header.direction[0], hsdata->run_header.direction[1]);
               rc = 1;
            }
            if ( hsdata->run_header.offset_fov[0] != hsdata_out->run_header.offset_fov[0] ||
                 hsdata->run_header.offset_fov[1] != hsdata_out->run_header.offset_fov[1] )
            {
               printf("File 1 has FoV offset %f, %f, file 2 has %f, %f.\n",
                  hsdata_out->run_header.offset_fov[0], hsdata_out->run_header.offset_fov[1], 
                  hsdata->run_header.offset_fov[0], hsdata->run_header.offset_fov[1]);
               rc = 1;
            }
            if ( hsdata->run_header.conv_depth != hsdata_out->run_header.conv_depth )
            {
               printf("File 1 has convergence depth %f, file 2 has %f.\n",
                  hsdata_out->run_header.conv_depth, hsdata->run_header.conv_depth);
               rc = 1;
            }
            if ( hsdata->run_header.conv_ref_pos[0] != hsdata_out->run_header.conv_ref_pos[0] ||
                 hsdata->run_header.conv_ref_pos[1] != hsdata_out->run_header.conv_ref_pos[1] )
            {
               printf("File 1 has convergence reference position %f, %f, file 2 has %f, %f.\n",
                  hsdata_out->run_header.conv_ref_pos[0], hsdata_out->run_header.conv_ref_pos[1], 
                  hsdata->run_header.conv_ref_pos[0], hsdata->run_header.conv_ref_pos[1]);
               rc = 1;
            }
            if ( hsdata->run_header.min_tel_trig != hsdata_out->run_header.min_tel_trig )
            {
               printf("File 1 requires by default %d for an array trigger, file 2 by default %d.\n",
                  hsdata_out->run_header.min_tel_trig, hsdata->run_header.min_tel_trig);
               /* No error */
            }
            if ( rc != 0 )
            {
               fflush(stdout);
               fprintf(stderr,"Run headers do not match. Giving up now.\n");
               exit(1);
            }
         }
         /* Reset telescope ID to index number lookup */
         for (itel=0; itel<=H_MAX_TEL; itel++)
            tel_idx[ifile-1][itel] = -1;
         /* Now initialize the telescope IDs from the new run */
         for (itel=0; itel<hsdata->run_header.ntel; itel++)
         {
            tel_id = hsdata->run_header.tel_id[itel];
            tel_id3 = find_mapped_telescope(tel_id, ifile);
            if ( tel_id >= 0 && tel_id <= H_MAX_TEL )
            {
               tel_idx[ifile-1][tel_id] = itel;
            }
            else
               fprintf(stderr,"Telescope ID %d in file %d is outside of valid range (0 to %d)\n",
                  tel_id, ifile, H_MAX_TEL);

            hsdata->camera_set[itel].tel_id = tel_id;
            hsdata->camera_org[itel].tel_id = tel_id;
            hsdata->pixel_set[itel].tel_id = tel_id;
            hsdata->pixel_disabled[itel].tel_id = tel_id;
            hsdata->cam_soft_set[itel].tel_id = tel_id;
            hsdata->tracking_set[itel].tel_id = tel_id;
            hsdata->point_cor[itel].tel_id = tel_id;
            hsdata->event.teldata[itel].tel_id = tel_id;
            hsdata->event.trackdata[itel].tel_id = tel_id;
            if ( (hsdata->event.teldata[itel].raw = 
                   (AdcData *) calloc(1,sizeof(AdcData))) == NULL )
            {
               Warning("Not enough memory");
               exit(1);
            }
            hsdata->event.teldata[itel].raw->tel_id = tel_id;
            if ( (hsdata->event.teldata[itel].pixtm =
                  (PixelTiming *) calloc(1,sizeof(PixelTiming))) == NULL )
            {
               Warning("Not enough memory");
               exit(1);
            }
            hsdata->event.teldata[itel].pixtm->tel_id = tel_id;
            if ( (hsdata->event.teldata[itel].img = 
                   (ImgData *) calloc(2,sizeof(ImgData))) == NULL )
            {
               Warning("Not enough memory");
               exit(1);
            }
            hsdata->event.teldata[itel].max_image_sets = 2;
            hsdata->event.teldata[itel].img[0].tel_id = tel_id;
            hsdata->event.teldata[itel].img[1].tel_id = tel_id;
            hsdata->tel_moni[itel].tel_id = tel_id;
            hsdata->tel_lascal[itel].tel_id = tel_id;


            if ( tel_id3 > 0 && tel_id3 <= H_MAX_TEL )
            {
               itel3 = tel_idx_out[tel_id3];
            }
            else
            {
               itel3 = -1;
            }
            if ( itel3 >= 0 && itel3 < H_MAX_TEL )
            {
               if ( map_tel[itel3].tel_id != tel_id3 ||
                    map_tel[itel3].ifn != ifile ||
                    map_tel[itel3].inp_id != tel_id )
               {
                  fflush(stdout);
                  fprintf(stderr,"Telescope ID lookup problem:\n"
                     "  At offset %d expected ID %d (input as ID %d from file %d)\n"
                     "  but instead found ID %d (ID %d from file %d at offset %d)\n",
                     itel3, map_tel[itel3].tel_id, map_tel[itel3].inp_id, map_tel[itel3].ifn,
                     tel_id3, tel_id, ifile, itel);
                  exit(1);
               }
               map_tel[itel3].inp_itel = itel;
               map_tel[itel3].have_camset = 0;
               map_tel[itel3].have_camorg = 0;
               map_tel[itel3].have_pixset = 0;
               map_tel[itel3].have_pixdis = 0;
               map_tel[itel3].have_camsoft = 0;
               map_tel[itel3].have_pointcor = 0;
               map_tel[itel3].have_trackset = 0;

               hsdata_out->camera_set[itel3].tel_id = tel_id3;
               hsdata_out->camera_org[itel3].tel_id = tel_id3;
               hsdata_out->pixel_set[itel3].tel_id = tel_id3;
               hsdata_out->pixel_disabled[itel3].tel_id = tel_id3;
               hsdata_out->cam_soft_set[itel3].tel_id = tel_id3;
               hsdata_out->tracking_set[itel3].tel_id = tel_id3;
               hsdata_out->point_cor[itel3].tel_id = tel_id3;
               hsdata_out->event.teldata[itel3].tel_id = tel_id3;
               hsdata_out->event.trackdata[itel3].tel_id = tel_id3;
               if ( hsdata_out->event.teldata[itel3].raw == NULL )
                  if ( (hsdata_out->event.teldata[itel3].raw = 
                         (AdcData *) calloc(1,sizeof(AdcData))) == NULL )
                  {
                     Warning("Not enough memory");
                     exit(1);
                  }
               hsdata_out->event.teldata[itel3].raw->tel_id = tel_id3;
               if ( hsdata_out->event.teldata[itel3].pixtm == NULL )
                  if ( (hsdata_out->event.teldata[itel3].pixtm =
                        (PixelTiming *) calloc(1,sizeof(PixelTiming))) == NULL )
                  {
                     Warning("Not enough memory");
                     exit(1);
                  }
               hsdata_out->event.teldata[itel3].pixtm->tel_id = tel_id3;
               if ( hsdata_out->event.teldata[itel3].img == NULL )
                  if ( (hsdata_out->event.teldata[itel3].img = 
                         (ImgData *) calloc(2,sizeof(ImgData))) == NULL )
                  {
                     Warning("Not enough memory");
                     exit(1);
                  }
               hsdata_out->event.teldata[itel3].max_image_sets = 2;
               hsdata_out->event.teldata[itel3].img[0].tel_id = tel_id3;
               hsdata_out->event.teldata[itel3].img[1].tel_id = tel_id3;
               hsdata_out->tel_moni[itel3].tel_id = tel_id3;
               hsdata_out->tel_lascal[itel3].tel_id = tel_id3;
               
               hsdata_out->run_header.tel_id[itel3] = tel_id3;
               hsdata_out->run_header.tel_pos[itel3][0] = hsdata->run_header.tel_pos[itel][0];
               hsdata_out->run_header.tel_pos[itel3][1] = hsdata->run_header.tel_pos[itel][1];
               hsdata_out->run_header.tel_pos[itel3][2] = hsdata->run_header.tel_pos[itel][2];
            }

            events[ifile-1] = mcshowers[ifile-1] = mcevents[ifile-1] = 0;
         }
         if ( ifile == 2 )
         {
            set_tel_idx_ref(0);
            write_hess_runheader(iobuf_out,&hsdata_out->run_header);
         }
         break;

      /* =================================================== */
      case IO_TYPE_HESS_MCRUNHEADER: /* 2001 */
         check_for_delayed_write(item_header, ifile, hsdata_out, iobuf_out);
         rc = read_hess_mcrunheader(iobuf,&hsdata->mc_run_header);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_mcrunheader(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         if ( ifile == 1 )
            memcpy(&hsdata_out->mc_run_header,&hsdata->mc_run_header,sizeof(hsdata->mc_run_header));
         /* Anything different between headers / to be merged ? */
         /* Except for the detector simulation version, which could differ, */
         /* pretty much everything should be identical. */
         /* Nevertheless, record both MC run headers as they came in. */
         rc = write_io_block_to_file (iobuf,iobuf_out->output_file);
         break;

      /* =================================================== */
      case IO_TYPE_MC_INPUTCFG: /* 1212 */
         /* Copy to output buffer without unpacking it, if possible ... */
         rc = write_io_block_to_file(iobuf, iobuf_out->output_file);
	 if ( rc != 0 || verbose )
            printf("writing input cfg to output, rc=%d\n", rc);
         break;

      /* =================================================== */
      case 70: /* How sim_hessarray was run and how it was configured. */
         /* Copy to output buffer without unpacking it, if possible ... */
         rc = write_io_block_to_file(iobuf, iobuf_out->output_file);
	 if ( rc != 0 || verbose )
            printf("writing history to output, rc=%d\n", rc);
         break;

      /* =================================================== */
      case IO_TYPE_HESS_CAMSETTINGS:  /* 2002: Per-telescope block */
         tel_id = item_header->ident; // Telescope ID is in the header
         itel = find_in_tel_idx(tel_id,ifile);
         itel3 = find_out_tel_idx(tel_id,ifile);
         if ( itel3 >= 0 && itel3 < H_MAX_TEL )
            tel_id3 = map_tel[itel3].tel_id;
         else
            tel_id3 = -1;
#ifdef DEBUG_MERGE
printf("itel = %d, tel_id = %d;  itel3 = %d, tel_id3 = %d\n", itel, tel_id, itel3, tel_id3);
#endif
         if ( itel < 0 || tel_id3 < 0 )
            return 0;
         rc = read_hess_camsettings(iobuf,&hsdata->camera_set[itel]);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_camsettings(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         /* Use memcpy for copying simple struct without any pointers in it. */
         memcpy(&hsdata_out->camera_set[itel3],&hsdata->camera_set[itel],sizeof(hsdata->camera_set[0]));
         /* All we need to change is the telescope ID */
         hsdata_out->camera_set[itel3].tel_id = tel_id3;
         /* Write the modified camera settings to output */
         rc = write_hess_camsettings(iobuf_out, &hsdata_out->camera_set[itel3]);
         break;

      /* =================================================== */
      case IO_TYPE_HESS_CAMORGAN:  /* 2003: Per-telescope block */
         tel_id = item_header->ident; // Telescope ID is in the header
         itel = find_in_tel_idx(tel_id,ifile);
         itel3 = find_out_tel_idx(tel_id,ifile);
         if ( itel3 >= 0 && itel3 < H_MAX_TEL )
            tel_id3 = map_tel[itel3].tel_id;
         else
            tel_id3 = -1;
#ifdef DEBUG_MERGE
printf("itel = %d, tel_id = %d;  itel3 = %d, tel_id3 = %d\n", itel, tel_id, itel3, tel_id3);
#endif
         if ( itel < 0 || tel_id3 < 0 )
            return 0;
         rc = read_hess_camorgan(iobuf,&hsdata->camera_org[itel]);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_camorgan(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         /* Use memcpy for copying simple struct without any pointers in it. */
         memcpy(&hsdata_out->camera_org[itel3],&hsdata->camera_org[itel],sizeof(hsdata->camera_org[0]));
         /* All we need to change is the telescope ID */
         hsdata_out->camera_org[itel3].tel_id = tel_id3;
         /* Write the modified camera organisation to output */
         rc = write_hess_camorgan(iobuf_out, &hsdata_out->camera_org[itel3]);
         break;

      /* =================================================== */
      case IO_TYPE_HESS_PIXELSET:  /* 2004: Per-telescope block */
         tel_id = item_header->ident; // Telescope ID is in the header
         itel = find_in_tel_idx(tel_id,ifile);
         itel3 = find_out_tel_idx(tel_id,ifile);
         if ( itel3 >= 0 && itel3 < H_MAX_TEL )
            tel_id3 = map_tel[itel3].tel_id;
         else
            tel_id3 = -1;
#ifdef DEBUG_MERGE
printf("itel = %d, tel_id = %d;  itel3 = %d, tel_id3 = %d\n", itel, tel_id, itel3, tel_id3);
#endif
         if ( itel < 0 || tel_id3 < 0 )
            return 0;
         rc = read_hess_pixelset(iobuf,&hsdata->pixel_set[itel]);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_pixelset(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         /* Use memcpy for copying simple struct without any pointers in it. */
         memcpy(&hsdata_out->pixel_set[itel3],&hsdata->pixel_set[itel],sizeof(hsdata->pixel_set[0]));
         /* All we need to change is the telescope ID */
         hsdata_out->pixel_set[itel3].tel_id = tel_id3;
         /* Write the modified camera organisation to output */
         rc = write_hess_pixelset(iobuf_out, &hsdata_out->pixel_set[itel3]);
         break;

      /* =================================================== */
      case IO_TYPE_HESS_PIXELDISABLE:  /* 2005: Per-telescope block */
         tel_id = item_header->ident; // Telescope ID is in the header
         itel = find_in_tel_idx(tel_id,ifile);
         itel3 = find_out_tel_idx(tel_id,ifile);
         if ( itel3 >= 0 && itel3 < H_MAX_TEL )
            tel_id3 = map_tel[itel3].tel_id;
         else
            tel_id3 = -1;
#ifdef DEBUG_MERGE
printf("itel = %d, tel_id = %d;  itel3 = %d, tel_id3 = %d\n", itel, tel_id, itel3, tel_id3);
#endif
         if ( itel < 0 || tel_id3 < 0 )
            return 0;
         rc = read_hess_pixeldis(iobuf,&hsdata->pixel_disabled[itel]);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_pixeldis(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         /* Use memcpy for copying simple struct without any pointers in it. */
         memcpy(&hsdata_out->pixel_disabled[itel3],&hsdata->pixel_disabled[itel],sizeof(hsdata->pixel_disabled[0]));
         /* All we need to change is the telescope ID */
         hsdata_out->pixel_disabled[itel3].tel_id = tel_id3;
         /* Write the modified camera organisation to output */
         rc = write_hess_pixeldis(iobuf_out, &hsdata_out->pixel_disabled[itel3]);
         break;

      /* =================================================== */
      case IO_TYPE_HESS_CAMSOFTSET:  /* 2006: Per-telescope block */
         tel_id = item_header->ident; // Telescope ID is in the header
         itel = find_in_tel_idx(tel_id,ifile);
         itel3 = find_out_tel_idx(tel_id,ifile);
         if ( itel3 >= 0 && itel3 < H_MAX_TEL )
            tel_id3 = map_tel[itel3].tel_id;
         else
            tel_id3 = -1;
#ifdef DEBUG_MERGE
printf("itel = %d, tel_id = %d;  itel3 = %d, tel_id3 = %d\n", itel, tel_id, itel3, tel_id3);
#endif
         if ( itel < 0 || tel_id3 < 0 )
            return 0;
         rc = read_hess_camsoftset(iobuf,&hsdata->cam_soft_set[itel]);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_camsoftset(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         /* Use memcpy for copying simple struct without any pointers in it. */
         memcpy(&hsdata_out->cam_soft_set[itel3],&hsdata->cam_soft_set[itel],sizeof(hsdata->cam_soft_set[0]));
         /* All we need to change is the telescope ID */
         hsdata_out->cam_soft_set[itel3].tel_id = tel_id3;
         /* Write the modified camera organisation to output */
         rc = write_hess_camsoftset(iobuf_out, &hsdata_out->cam_soft_set[itel3]);
         break;

      /* =================================================== */
      case IO_TYPE_HESS_POINTINGCOR:  /* 2007: Per-telescope block */
         tel_id = item_header->ident; // Telescope ID is in the header
         itel = find_in_tel_idx(tel_id,ifile);
         itel3 = find_out_tel_idx(tel_id,ifile);
         if ( itel3 >= 0 && itel3 < H_MAX_TEL )
            tel_id3 = map_tel[itel3].tel_id;
         else
            tel_id3 = -1;
#ifdef DEBUG_MERGE
printf("itel = %d, tel_id = %d;  itel3 = %d, tel_id3 = %d\n", itel, tel_id, itel3, tel_id3);
#endif
         if ( itel < 0 || tel_id3 < 0 )
            return 0;
         rc = read_hess_pointingcor(iobuf,&hsdata->point_cor[itel]);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_pointingco(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         /* Use memcpy for copying simple struct without any pointers in it. */
         memcpy(&hsdata_out->point_cor[itel3],&hsdata->point_cor[itel],sizeof(hsdata->point_cor[0]));
         /* All we need to change is the telescope ID */
         hsdata_out->point_cor[itel3].tel_id = tel_id3;
         /* Write the modified camera organisation to output */
         rc = write_hess_pointingcor(iobuf_out, &hsdata_out->point_cor[itel3]);
         break;

      /* =================================================== */
      case IO_TYPE_HESS_TRACKSET:  /* 2008: Per-telescope block */
         tel_id = item_header->ident; // Telescope ID is in the header
         itel = find_in_tel_idx(tel_id,ifile);
         itel3 = find_out_tel_idx(tel_id,ifile);
         if ( itel3 >= 0 && itel3 < H_MAX_TEL )
            tel_id3 = map_tel[itel3].tel_id;
         else
            tel_id3 = -1;
#ifdef DEBUG_MERGE
printf("itel = %d, tel_id = %d;  itel3 = %d, tel_id3 = %d\n", itel, tel_id, itel3, tel_id3);
#endif
         if ( itel < 0 || tel_id3 < 0 )
            return 0;
         rc = read_hess_trackset(iobuf,&hsdata->tracking_set[itel]);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_trackset(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         /* Use memcpy for copying simple struct without any pointers in it. */
         memcpy(&hsdata_out->tracking_set[itel3],&hsdata->tracking_set[itel],sizeof(hsdata->tracking_set[0]));
         /* All we need to change is the telescope ID */
         hsdata_out->tracking_set[itel3].tel_id = tel_id3;
         /* Write the modified camera organisation to output */
         rc = write_hess_trackset(iobuf_out, &hsdata_out->tracking_set[itel3]);
        break;

      /* =================================================== */
      /* Normal data event: to be merged */
      case IO_TYPE_HESS_EVENT: /* 2010 */
if ( ev_hess_event > 0 && ev_hess_event < item_header->ident )
{
   fflush(stdout);
   fprintf(stderr,"Unprocessed event data block for event %ld gets lost.\n",
      ev_hess_event);
}
      {
         FullEvent *evi = &hsdata->event;
         FullEvent *evo = &hsdata_out->event;
         int j, ntrg=0, ndata=0, ndt=0;
         rc = read_hess_event(iobuf,evi,-1);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_event(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         events[ifile-1]++;

         if ( ev_hess_event > 0 && ev_hess_event != item_header->ident )
         {
            fflush(stdout);
            fprintf(stderr,
               "Encountered full event data block for event %ld before %ld was written.\n",
               item_header->ident, ev_hess_event);
         }
         if ( ev_hess_event == 0 )
         {
            evo->central.glob_count = evo->central.teltrg_pattern = 
               evo->central.teldata_pattern = 
               evo->central.num_teltrg = evo->central.num_teldata = 0;
            evo->central.num_teldata = evo->central.num_teltrg = 0;
            /* evo->num_tel was set based on run header data */
            for ( j=0; j<evo->num_tel; j++ )
            {
               evo->teldata[j].known = 0;
               evo->trackdata[j].raw_known = evo->trackdata[j].cor_known = 0;
            }
            evo->shower.known = 0;
            evo->num_teldata = 0;
         }

         /* First of all fix trigger type bit masks, if necessary */
         if ( ths[ifile-1] != NULL && ths[ifile-1]->run == hsdata->run_header.run )
         {
            struct trgmask_entry *tme;
            for ( j=0; j<evi->central.num_teltrg; j++ )
            {
               tel_id = evi->central.teltrg_list[j];
               tme = find_trgmask(ths[ifile-1],item_header->ident,tel_id);
               if ( tme == NULL )
               {
                  printf("No extra trigger type mask found for telescope %d, setting from %d to 0.\n",
                     tel_id, evi->central.teltrg_type_mask[j]);
                  evi->central.teltrg_type_mask[j] = 0;
                  evi->central.teltrg_time_by_type[j][0] = 
                     evi->central.teltrg_time_by_type[j][1] = 
                     evi->central.teltrg_time_by_type[j][2] = 9999.;
               }
               else
               {
                  /* Note: Any trigger time by type corresponding to a missing bit was lost. */
                  if ( !(tme->trg_mask & 0x01) )
                     evi->central.teltrg_time_by_type[j][0] = 9999.;
                  if ( !(tme->trg_mask & 0x02) )
                     evi->central.teltrg_time_by_type[j][1] = 9999.;
                  if ( !(tme->trg_mask & 0x04) )
                     evi->central.teltrg_time_by_type[j][2] = 9999.;
#ifdef DEBUG_MERGE
   printf("Fixing trigger type bits for telescope ID %d in file %d from %d to %d.\n",
      tel_id, ifile, evi->central.teltrg_type_mask[j], tme->trg_mask);
#endif
                  evi->central.teltrg_type_mask[j] = tme->trg_mask;
               }
            }
         }
#ifdef DEBUG_MERGE
         else
            printf("Not fixing trigger type bits in file %d\n", ifile);
#endif

         /* Central trigger event data */
         evo->central.glob_count = evi->central.glob_count;
         memcpy(&evo->central.cpu_time,&evi->central.cpu_time,sizeof(HTime));
         memcpy(&evo->central.gps_time,&evi->central.gps_time,sizeof(HTime));
         /* Clear old-style bit patterns. Use lists only. */
         evo->central.teltrg_pattern = 0;
         evo->central.teldata_pattern = 0;
         ntrg = evo->central.num_teltrg;
         ndata = evo->central.num_teldata;
         for ( j=0; j<evi->central.num_teltrg; j++ )
         {
            tel_id = evi->central.teltrg_list[j];
            itel3 = find_out_tel_idx(tel_id,ifile);
            if ( itel3 >= 0 && itel3 < H_MAX_TEL && ntrg < H_MAX_TEL )
            {
               tel_id3 = map_tel[itel3].tel_id;
               evo->central.teltrg_list[ntrg] = tel_id3;
               evo->central.teltrg_time[ntrg] = evi->central.teltrg_time[j];
               evo->central.teltrg_type_mask[ntrg] = evi->central.teltrg_type_mask[j];
               evo->central.teltrg_time_by_type[ntrg][0] = evi->central.teltrg_time_by_type[j][0];
               evo->central.teltrg_time_by_type[ntrg][1] = evi->central.teltrg_time_by_type[j][1];
               evo->central.teltrg_time_by_type[ntrg][2] = evi->central.teltrg_time_by_type[j][2];
               ntrg++;
            }
         }
         evo->central.num_teltrg = ntrg;
         for ( j=0; j<evi->central.num_teldata; j++ )
         {
            tel_id = evi->central.teldata_list[j];
            itel3 = find_out_tel_idx(tel_id,ifile);
            if ( itel3 >= 0 && itel3 < H_MAX_TEL && ndata < H_MAX_TEL )
            {
               tel_id3 = map_tel[itel3].tel_id;
               evo->central.teldata_list[ndata] = tel_id3;
               ndata++;
            }
         }
         evo->central.num_teldata = ndata;

         /* Telescope event data */
         ndt = evo->num_teldata; /* May differ from central trigger value */
         for ( j=0; j<evi->num_teldata; j++ )
         {
            tel_id = evi->teldata_list[j];
            itel = find_in_tel_idx(tel_id,ifile);
            itel3 = find_out_tel_idx(tel_id,ifile);
            if ( itel3 >= 0 && itel3 < H_MAX_TEL )
            {
               int k;
               TelEvent *tei = &evi->teldata[itel];
               TelEvent *teo = &evo->teldata[itel3];
               tel_id3 = map_tel[itel3].tel_id;

               /* Tracking data */
               memcpy(&evo->trackdata[itel3],&evi->trackdata[itel],sizeof(TrackEvent));
               evo->trackdata[itel3].tel_id = tel_id3;
 
               if ( !tei->known )
               {
                  printf("Warning: telescope with event data is not marked as known (file %d, tel. ID %d -> %d).\n",
                     ifile, tel_id, tel_id3);
                  continue;
               }
               teo->loc_count = tei->loc_count;
               teo->glob_count = tei->glob_count;
               memcpy(&teo->cpu_time,&tei->cpu_time,sizeof(HTime));
               memcpy(&teo->gps_time,&tei->gps_time,sizeof(HTime));
               teo->trg_source = tei->trg_source;
               teo->num_list_trgsect = tei->num_list_trgsect;
               teo->known_time_trgsect = tei->known_time_trgsect;
               for ( k=0; k<teo->num_list_trgsect; k++ )
               {
                  teo->list_trgsect[k] = tei->list_trgsect[k];
                  teo->time_trgsect[k] = tei->time_trgsect[k];
               }
               teo->readout_mode = tei->readout_mode;
               teo->num_image_sets = tei->num_image_sets;
               if ( teo->num_image_sets > teo->max_image_sets )
                  teo->num_image_sets = teo->max_image_sets;

               /* Copy raw data */
               if ( teo->raw != NULL && tei->raw != NULL )
               {
                  AdcData *adi = tei->raw;
                  AdcData *ado = teo->raw;
                  if ( adi != NULL && ado != NULL )
                  {
                   ado->known = adi->known;
                   if ( ado->known )
                   {
                     int kg, kp, ks;
                     ado->tel_id = tel_id3; /* Not: adi->tel_id */
                     ado->num_pixels = adi->num_pixels;
                     ado->num_gains = adi->num_gains;
                     ado->num_samples = adi->num_samples;
                     ado->zero_sup_mode = adi->zero_sup_mode;
                     ado->data_red_mode = adi->data_red_mode;
                     ado->offset_hg8 = adi->offset_hg8;
                     ado->scale_hg8 = adi->scale_hg8;
                     ado->threshold = adi->threshold;
                     ado->list_known = adi->list_known;
                     ado->list_size = adi->list_size;
                     for ( k=0; k<ado->list_size; k++ )
                        ado->adc_list[k] = adi->adc_list[k];
                     for ( kp=0; kp<ado->num_pixels; kp++ )
                        ado->significant[kp] = adi->significant[kp];
                     for ( kg=0; kg<ado->num_gains; kg++ )
                     {
                        for ( kp=0; kp<ado->num_pixels; kp++ )
                        {
                           ado->adc_known[kg][kp] = adi->adc_known[kg][kp];
                           ado->adc_sum[kg][kp] = adi->adc_sum[kg][kp];
                           for ( ks=0; ks<ado->num_samples; ks++ )
                              ado->adc_sample[kg][kp][ks] = adi->adc_sample[kg][kp][ks];
                        }
                     }
                   }
                  }
               }
               else
               {
                  fprintf(stderr,"Missing raw data pointer for file %d, tel. ID %d -> %d\n",
                     ifile, tel_id, tel_id3);
               }
               /* Copy pixel timing */
               if ( teo->pixtm != NULL && tei->pixtm != NULL )
               {
                teo->pixtm->known = tei->pixtm->known;
                if ( teo->pixtm->known )
                {
                  int kg, kp, kt;
                  PixelTiming *pi = tei->pixtm;
                  PixelTiming *po = teo->pixtm;
                  po->tel_id = tel_id3; /* Not: pi->tel_id */
                  po->num_pixels = pi->num_pixels;
                  po->num_gains = pi->num_gains;
                  po->list_type = pi->list_type;
                  po->list_size = pi->list_size;
                  for ( k=0; k<po->list_size; k++ )
                  {
                     if ( pi->list_type == 1 )
                        po->pixel_list[k] = pi->pixel_list[k];
                     else
                     {
                        po->pixel_list[2*k] = pi->pixel_list[2*k];
                        po->pixel_list[2*k+1] = pi->pixel_list[2*k+1];
                     }
                  }
                  po->threshold = pi->threshold;
                  po->before_peak = pi->before_peak;
                  po->after_peak = pi->after_peak;
                  po->num_types = pi->num_types;
                  for ( k=0; k<po->num_types; k++ )
                  {
                     po->time_type[k] = pi->time_type[k];
                     po->time_level[k] = pi->time_level[k];
                  }
                  po->granularity = pi->granularity;
                  po->peak_global = pi->peak_global;
                  for ( kp=0; kp<po->num_pixels; kp++ )
                     for ( kt=0; kt<po->num_types; kt++ )
                        po->timval[kp][kt] = pi->timval[kp][kt];
                  for ( kg=0; kg<po->num_gains; kg++ )
                     for ( kp=0; kp<po->num_pixels; kp++ )
                     {
                        po->pulse_sum_loc[kg][kp] = pi->pulse_sum_loc[kg][kp];
                        po->pulse_sum_glob[kg][kp] = pi->pulse_sum_glob[kg][kp];
                     }
                }
               }
               else
               {
                  fprintf(stderr,"Missing pixel timing pointer for file %d, tel. ID %d -> %d\n",
                     ifile, tel_id, tel_id3);
               }
               
               /* Copy image data */
               for ( k=0; k<teo->num_image_sets; k++ )
               {
                  if ( teo->img != NULL && tei->img != NULL )
                  {
                     memcpy(&teo->img[k], &tei->img[k], sizeof(ImgData));
                     teo->img[k].tel_id = tel_id3;
                  }
                  else if ( k == 0 )
                  {
                     fprintf(stderr,"Missing image pointer for file %d, tel. ID %d -> %d\n",
                        ifile, tel_id, tel_id3);
                  }
               }
               // teo->num_phys_addr = tei->num_phys_addr; /* not used */
               /* Copy trigger pixels list */
               teo->trigger_pixels.code = tei->trigger_pixels.code;
               teo->trigger_pixels.pixels = tei->trigger_pixels.pixels;
               for ( k=0; k<tei->trigger_pixels.pixels; k++ )
               {
                  teo->trigger_pixels.pixel_list[k] = tei->trigger_pixels.pixel_list[k];
               }
               /* Copy image pixels list */
               teo->image_pixels.code = tei->image_pixels.code;
               teo->image_pixels.pixels = tei->image_pixels.pixels;
               for ( k=0; k<tei->image_pixels.pixels; k++ )
               {
                  teo->image_pixels.pixel_list[k] = tei->image_pixels.pixel_list[k];
               }
               teo->known = tei->known;
               evo->teldata_list[ndt] = tel_id3;
               ndt++;
            }
         }
         ndt = evo->num_teldata = ndt;
         ev_hess_event = item_header->ident;
if ( nlist < max_list || nlist%1000 ==0 )
printf("Event data for event %ld from file %d: ntrg=%d, ndata=%d, ndt=%d\n",
   ev_hess_event, ifile, ntrg, ndata, ndt);
else if ( nlist == max_list )
   printf("...\n");
nlist++;
      }
         break;

      /* =================================================== */
      /* Calibration data event: Not merging for now. They should be */
      /* split off with the extract_hess tool before merging and */
      /* then handled separately. Not used in normal production. */
      case IO_TYPE_HESS_CALIBEVENT: /* 2028 */
      {
         int type = -1;
         rc = read_hess_calib_event(iobuf,&hsdata->event,-1,&type);
	 if ( rc != 0 || verbose )
            printf("read_hess_calib_event(), rc = %d, type=%d\n",rc,type);
printf("Not merging calibration event data yet. Stopping here.\n");
exit(1);
      }
         break;

      /* =================================================== */
      case IO_TYPE_HESS_MC_SHOWER: /* 2020: Expected to be identical */
         check_for_delayed_write(item_header, ifile, hsdata_out, iobuf_out);
         rc = read_hess_mc_shower(iobuf,&hsdata->mc_shower);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_mc_shower(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         mcshowers[ifile-1]++;
         /* If we see this MC sower block for the first time, write it out as-is. */
         /* If we got it first from the other input, check that it matches - */
         /* and then ignore it. */
         if ( hsdata->mc_shower.shower_num != hsdata_out->mc_shower.shower_num )
         {
            if ( hsdata->mc_shower.shower_num < hsdata_out->mc_shower.shower_num )
            {
               fflush(stdout);
               fprintf(stderr,"Wrong order of processing, MC shower %d after shower %d\n",
                  hsdata->mc_shower.shower_num, hsdata_out->mc_shower.shower_num);
               exit(1);
            }
            memcpy(&hsdata_out->mc_shower,&hsdata->mc_shower,sizeof(hsdata->mc_shower));
            /* Write the original data block to the output file. */
            rc = write_io_block_to_file(iobuf, iobuf_out->output_file);
         }
         else
         {
            if ( hsdata->mc_shower.primary_id != hsdata_out->mc_shower.primary_id ||
                 hsdata->mc_shower.energy != hsdata_out->mc_shower.energy ||
                 hsdata->mc_shower.azimuth != hsdata_out->mc_shower.azimuth ||
                 hsdata->mc_shower.altitude != hsdata_out->mc_shower.altitude )
            {
               fflush(stdout);
               fprintf(stderr,"MC shower blocks do not match.\n");
               exit(1);
            }
         }
         break;

      /* =================================================== */
      case IO_TYPE_HESS_MC_EVENT: /* 2021 */
         check_for_delayed_write(item_header, ifile, hsdata_out, iobuf_out);
         rc = read_hess_mc_event(iobuf,&hsdata->mc_event);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_mc_event(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         mcevents[ifile-1]++;
      {
         MCEvent *mce  = &hsdata->mc_event;
         MCEvent *mceo = &hsdata_out->mc_event;
         /* If we see this MC event block for the first time, write it out as-is. */
         /* If we got it first from the other input, just ignore this data block. */
         if ( hsdata->mc_event.event != hsdata_out->mc_event.event ||
              hsdata->mc_event.shower_num != hsdata_out->mc_event.shower_num )
         {
            if ( mce->event < mceo->event ||
                 mce->shower_num < mceo->shower_num )
            {
               fflush(stdout);
               fprintf(stderr,"Wrong order of processing, MC event %d after event %d\n",
                  mce->event, mceo->event);
               exit(1);
            }
            /* Only copy what is really needed obtained with this data block. */
            /* The other parts of this structure can be filled elsewhere. */
            mceo->event = mce->event;
            mceo->shower_num = mce->shower_num;
            mceo->xcore = mce->xcore;
            mceo->ycore = mce->ycore;
            mceo->aweight = mce->aweight;
            /* Reset some of the sub-structure contents. */
            mceo->mc_pesum.num_tel = 0;
            mceo->mc_pesum.event = -1;
            for ( itel=0; itel< hsdata_out->run_header.ntel; itel++ )
            {
               mceo->mc_pesum.num_pe[itel] = -1;
               mceo->mc_photons[itel].nbunches = -1;
               mceo->mc_pe_list[itel].npe = -1;
               mceo->mc_pesum.photons[itel] = 
               mceo->mc_pesum.photons_atm[itel] = 
               mceo->mc_pesum.photons_atm_3_6[itel] = 
               mceo->mc_pesum.photons_atm_400[itel] = 
               mceo->mc_pesum.photons_atm_qe[itel] = 0.;
            }
            /* Write the original data block to the output file. */
            /* The other data reset above is not part of the I/O block. */
            rc = write_io_block_to_file(iobuf, iobuf_out->output_file);
         }
         else
         {
            if ( mce->xcore != mceo->xcore || mce->ycore != mceo->ycore )
            {
               fflush(stdout);
               fprintf(stderr,"MC event blocks do not match.\n");
               exit(1);
            }
         }
      }
         break;

      /* =================================================== */
      /* Photons/photo-electrons for redoing simulations: ignored. */
      case IO_TYPE_MC_TELARRAY: /* 1204 */
         if ( hsdata && hsdata->run_header.ntel > 0 )
         {
            rc = read_hess_mc_phot(iobuf,&hsdata->mc_event);
	    if ( rc != 0 || verbose )
            {
               printf("read_hess_mc_phot(), rc = %d\n",rc);
               if ( rc != 0 )
                  return rc;
            }
         }
#ifdef DEBUG_MERGE
printf("Ignoring photons/photo-electrons data block\n");
#endif
         break;

      /* =================================================== */
      case IO_TYPE_HESS_MC_PE_SUM: /* 2026 */
if ( ev_pe_sum > 0 && ev_pe_sum < item_header->ident )
{
   fflush(stdout);
   fprintf(stderr,"Unprocessed p.e. sum block for event %ld gets lost.\n",
      ev_pe_sum);
}
         rc = read_hess_mc_pe_sum(iobuf,&hsdata->mc_event.mc_pesum);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_mc_pe_sum(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         hsdata_out->mc_event.mc_pesum.event = hsdata->mc_event.event;
         hsdata_out->mc_event.mc_pesum.shower_num = hsdata->mc_event.shower_num;
         hsdata_out->mc_event.mc_pesum.num_tel = hsdata_out->run_header.ntel;
         for ( itel=0; itel<hsdata->run_header.ntel; itel++ )
         {
            tel_id = hsdata->run_header.tel_id[itel];
            itel3 = find_out_tel_idx(tel_id,ifile);
            if ( itel3 >= 0 && itel3 < H_MAX_TEL )
               tel_id3 = map_tel[itel3].tel_id;
            else
               tel_id3 = -1;
            if ( itel < 0 || tel_id3 < 0 )
               continue;
            hsdata_out->mc_event.mc_pesum.num_pe[itel3] = 
               hsdata->mc_event.mc_pesum.num_pe[itel];
            hsdata_out->mc_event.mc_pesum.num_pixels[itel3] = 
               hsdata->mc_event.mc_pesum.num_pixels[itel];
            if ( hsdata->mc_event.mc_pesum.num_pe[itel] > 0 && 
                 hsdata->mc_event.mc_pesum.num_pixels[itel] > 0 )
            {
               int ipix, npix = hsdata->mc_event.mc_pesum.num_pixels[itel];
               for ( ipix=0; ipix<npix; ipix++ )
               {
                  hsdata_out->mc_event.mc_pesum.pix_pe[itel3][ipix] = 
                     hsdata->mc_event.mc_pesum.pix_pe[itel][ipix];
               }
            }
            hsdata_out->mc_event.mc_pesum.photons[itel3] = 
               hsdata->mc_event.mc_pesum.photons[itel];
            hsdata_out->mc_event.mc_pesum.photons_atm[itel3] = 
               hsdata->mc_event.mc_pesum.photons_atm[itel];
            hsdata_out->mc_event.mc_pesum.photons_atm_3_6[itel3] = 
               hsdata->mc_event.mc_pesum.photons_atm_3_6[itel];
            hsdata_out->mc_event.mc_pesum.photons_atm_400[itel3] = 
               hsdata->mc_event.mc_pesum.photons_atm_400[itel];
            hsdata_out->mc_event.mc_pesum.photons_atm_qe[itel3] = 
               hsdata->mc_event.mc_pesum.photons_atm_qe[itel];
         }
         /* Writing of the data block is delayed until we have merged data. */
         ev_pe_sum = item_header->ident;
         break;

      /* =================================================== */
      case IO_TYPE_HESS_TEL_MONI: /* 2022: Per-telescope block */
         // Telescope ID among others in the header
         tel_id = (item_header->ident & 0xff) | 
                  ((item_header->ident & 0x3f000000) >> 16); 
         itel = find_in_tel_idx(tel_id,ifile);
         itel3 = find_out_tel_idx(tel_id,ifile);
         if ( itel3 >= 0 && itel3 < H_MAX_TEL )
            tel_id3 = map_tel[itel3].tel_id;
         else
            tel_id3 = -1;
         if ( itel < 0 || tel_id3 < 0 )
            return 0;
         rc = read_hess_tel_monitor(iobuf,&hsdata->tel_moni[itel]);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_tel_monitor(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         /* Use memcpy for copying simple struct without any pointers in it. */
         memcpy(&hsdata_out->tel_moni[itel3],&hsdata->tel_moni[itel],sizeof(hsdata->tel_moni[0]));
         /* All we need to change is the telescope ID */
         hsdata_out->tel_moni[itel3].tel_id = tel_id3;
         /* Write the modified monitor data for this telescope to output. */
         rc = write_hess_tel_monitor(iobuf_out,&hsdata_out->tel_moni[itel3],
            hsdata_out->tel_moni[itel3].new_parts);
         break;

      /* =================================================== */
      case IO_TYPE_HESS_LASCAL: /* 2023: Per-telescope block */
         tel_id = item_header->ident; // Telescope ID is in the header
         itel = find_in_tel_idx(tel_id,ifile);
         itel3 = find_out_tel_idx(tel_id,ifile);
         if ( itel3 >= 0 && itel3 < H_MAX_TEL )
            tel_id3 = map_tel[itel3].tel_id;
         else
            tel_id3 = -1;
         if ( itel < 0 || tel_id3 < 0 )
            return 0;
         rc = read_hess_laser_calib(iobuf,&hsdata->tel_lascal[itel]);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_laser_calib(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         /* Use memcpy for copying simple struct without any pointers in it. */
         memcpy(&hsdata_out->tel_lascal[itel3],&hsdata->tel_lascal[itel],sizeof(hsdata->tel_lascal[0]));
         /* All we need to change is the telescope ID */
         hsdata_out->tel_lascal[itel3].tel_id = tel_id3;
         /* Write the modified laser calibration data for this telescope to output. */
         rc = write_hess_laser_calib(iobuf_out,&hsdata_out->tel_lascal[itel3]);        
         break;

      /* =================================================== */
      case IO_TYPE_HESS_RUNSTAT: /* 2024 */
         check_for_delayed_write(item_header, ifile, hsdata_out, iobuf_out);
         /* Is that type actually being used ? */
         rc = read_hess_run_stat(iobuf,&hsdata->run_stat);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_run_stat(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         /* Use memcpy for copying simple struct without any pointers in it. Actually not needed. */
         memcpy(&hsdata_out->run_stat,&hsdata->run_stat,sizeof(hsdata->run_stat));
         /* Copy the data block from input 1 and ignore that from input 2. */
         if ( ifile == 1 )
            rc = write_io_block_to_file(iobuf, iobuf_out->output_file);
         break;

      /* =================================================== */
      case IO_TYPE_HESS_MC_RUNSTAT: /* 2025 */
         check_for_delayed_write(item_header, ifile, hsdata_out, iobuf_out);
         /* Is that type actually being used ? */
         rc = read_hess_mc_run_stat(iobuf,&hsdata->mc_run_stat);
	 if ( rc != 0 || verbose )
         {
            printf("read_hess_mc_run_stat(), rc = %d\n",rc);
            if ( rc != 0 )
               return rc;
         }
         /* Use memcpy for copying simple struct without any pointers in it. Actually not needed. */
         memcpy(&hsdata_out->mc_run_stat,&hsdata->mc_run_stat,sizeof(hsdata->mc_run_stat));
          /* Copy the data block from input 1 and ignore that from input 2. */
         if ( ifile == 1 )
            rc = write_io_block_to_file(iobuf, iobuf_out->output_file);
         break;

      /* =================================================== */
      /* In-data trigger pattern block, or in extra file (auto-loaded or explicit) */
      case IO_TYPE_HESS_XTRGMASK: /* 2090 */
         if ( tms[ifile-1] == NULL )
            tms[ifile-1] = calloc(1,sizeof(struct trgmask_set));
         if ( ths[ifile-1] == NULL )
            ths[ifile-1] = calloc(1,sizeof(struct trgmask_hash_set));
         rc = read_trgmask(iobuf, tms[ifile-1]);
	 if ( rc != 0 || verbose )
            printf("read_trgmask(), rc = %d\n",rc);
         trgmask_fill_hashed(tms[ifile-1],ths[ifile-1]);
         /* No need to merge or remap since applied to inputs. */
         /* No output needed since the trigger patterns will be fixed. */
         break;

      /* =================================================== */
      /* (End-of-job or DST) histograms */
      case 100:
         check_for_delayed_write(item_header, ifile, hsdata_out, iobuf_out);
         /* The processing loop should just offer the histogram block from one file, */
         /* typically from file 2, if available. But if file 2 has no histogram */
         /* block at a matching position, it could be from file 1.  */
         /* Thus no extra check here. */
         rc = write_io_block_to_file(iobuf, iobuf_out->output_file);
         break;

      /* =================================================== */
      default:
         fflush(stdout);
         fprintf(stderr,"Ignoring unknown data block type %ld from input file %d\n",
            item_header->type, ifile);

   }
#ifdef DEBUG_MERGE
printf("... done.\n");
#endif
   return rc;
}

/* --------------------- check_autoload_trgmask -------------------------- */
/**
 *   Check for a 'trgmask.gz' file matching the given input data file name and,
 *   if it exists, extract the corrected trigger bit patterns from it.
 *   (Note: this is only relevant for multi-trigger data produced with a
 *   bug in recording the trigger bit pattern.)
 *
 *   We do not need to merge the contents of this file since the
 *   trigger bit patterns are corrected after reading the data.
 */

int check_autoload_trgmask(const char *input_fname, IO_BUFFER *iobuf, int ifile);

int check_autoload_trgmask(const char *input_fname, IO_BUFFER *iobuf, int ifile)
{
   char fname[10240];
   char *s;
   int rc = 0;
   FILE *f1 = iobuf->input_file;
   IO_ITEM_HEADER item_header;
   
   if ( input_fname == NULL || iobuf == NULL || ifile < 1 || ifile > 2 )
      return -99;
 
   if ( tms[ifile-1] != NULL )
   {
      free(tms[ifile-1]);
      tms[ifile-1] = NULL;
   }
   if ( ths[ifile-1] != NULL )
   {
      free(ths[ifile-1]);
      ths[ifile-1] = NULL;
   }

   strncpy(fname,input_fname,sizeof(fname)-30);
   if ( (s = strstr(fname,".sim")) != NULL )     /* Not universal but covering the usual cases */
      *s = '\0';
   strcat(fname,".trgmask.gz");
   iobuf->input_file = fileopen(fname,READ_BINARY);
   if ( iobuf->input_file == NULL )
   {
      printf("No such file: %s\n", fname);
      return 0;
   }

   printf("Using extra trigger type bit patterns from %s\n",fname);
   if ( find_io_block(iobuf,&item_header) == 0 )
   {
      printf("Found I/O block of type %ld\n",item_header.type); // Expecting: 2090
      if ( read_io_block(iobuf,&item_header) == 0 )
      {
         if ( tms[ifile-1] == NULL )
            tms[ifile-1] = calloc(1,sizeof(struct trgmask_set));
         if ( ths[ifile-1] == NULL )
            ths[ifile-1] = calloc(1,sizeof(struct trgmask_hash_set));
         rc = read_trgmask(iobuf, tms[ifile-1]);
	 if ( rc != 0 )
            printf("read_trgmask(), rc = %d\n",rc);
         if ( rc == 0 )
            trgmask_fill_hashed(tms[ifile-1],ths[ifile-1]);
      }
   }
   fileclose(iobuf->input_file);
   iobuf->input_file = f1; /* Restore original file pointer */
   if ( rc < 0 )
      return rc;
   else
      return 1; /* Success, trgmask values loaded */
}

/* ---------------------- print_process_status --------------------------- */

void print_process_status (int prev_type1, int this_type1, int prev_type2, int this_type2);

void print_process_status (int prev_type1, int this_type1, int prev_type2, int this_type2)
{
   printf("Last processed data type from file 1: %d\n", prev_type1);
   if ( this_type1 != 0 )
      printf("Not yet processed data type from file 1: %d\n", this_type1);
   printf("Last processed data type from file 2: %d\n", prev_type2);
   if ( this_type2 != 0 )
      printf("Not yet processed data type from file 2: %d\n", this_type2);
}

/* ----------------------------- read_map -------------------------------- */

int read_map(const char *map_fname);

int read_map(const char *map_fname)
{
   int itel=0;
   FILE *map_file = NULL;
   char line[1000];
   char hl[120];

   /* Set defaults */

   for ( itel=0; itel<H_MAX_TEL; itel++ )
   {
      map_tel[itel].tel_id = -1;
      map_tel[itel].ifn = 0;
      map_tel[itel].inp_id = -1;
      map_tel[itel].inp_itel = -1;
      map_tel[itel].have_camset = 0;
      map_tel[itel].have_camorg = 0;
      map_tel[itel].have_pixset = 0;
      map_tel[itel].have_pixdis = 0;
      map_tel[itel].have_camsoft = 0;
      map_tel[itel].have_pointcor = 0;
      map_tel[itel].have_trackset = 0;

      map_to[0][itel] = map_to[1][itel] = -1;
   }
   for ( itel=0; itel<=H_MAX_TEL; itel++ )
      tel_idx[1][itel] = tel_idx[0][itel] = -1;

   /* Open map file */

   if ( (map_file = fileopen(map_fname,"r")) == NULL )
   {
      perror(map_fname);
      exit(1);
   }

   /* Read from map file */

   while ( fgets(line,sizeof(line)-1,map_file) != NULL )
   {
      char word[20];
      int ipos=0;
      char *s = strchr(line,'#');
      int to, ni, ti;
      if ( s != NULL )
         *s = '\0';
      while ( (s=strchr(line,'\t')) != NULL )
         *s = ' ';
      // push_config_history(line,1);

      if ( getword(line,&ipos,word,sizeof(word)-1,' ','\n') <= 0 )
         continue;
      ni = atoi(word);
      
      while ( getword(line,&ipos,word,sizeof(word)-1,',','\n') > 0 )
      {
         int tel_idx1 = atoi(word), tel_idx2;
         char *sl;
         if ( tel_idx1 > 0 )
         {
            if ( (sl = strchr(word,'-')) != NULL )
               tel_idx2 = atoi(sl+1);
            else
               tel_idx2 = tel_idx1;
            
            if ( tel_idx2 < tel_idx1 || tel_idx2 > H_MAX_TEL )
            {
               fprintf(stderr,"Invalid mapping line: %s\n", line);
               exit(1);
            }
            for (ti=tel_idx1; ti<=tel_idx2; ti++)
            {
               if ( ntel >= H_MAX_TEL )
               {
                  fprintf(stderr,"Too many mappings. Can only handle %d telescopes.\n",H_MAX_TEL);
                  exit(1);
               }
               to = ntel+1; /* Force incremental telescope IDs on output */
               if ( (ni != 1 && ni != 2) || 
                    to < 1 || to > H_MAX_TEL ||
                    ti < 1 || ti > H_MAX_TEL )
               {
                  fprintf(stderr,"Invalid mapping line: %s\n", line);
                  exit(1);
               }
               for (itel=0; itel<ntel; itel++)
               {
                  if ( map_tel[itel].tel_id == to )
                  {
                     fprintf(stderr,"Duplicated output telescope ID %d\n", to);
                     exit(1);
                  }
                  if ( map_tel[itel].inp_id == ti && map_tel[itel].ifn == ni )
                  {
                     fprintf(stderr,"Duplicated input telescope ID %d from file %d\n", ti, ni);
                     exit(1);
                  }
               }
               sprintf(hl,"   Telescope ID %3d from input file %d mapped to telescope ID %3d.\n",
                  ti, ni, to);
               push_config_history(hl,1);
               map_tel[ntel].tel_id = to;
               map_tel[ntel].ifn = ni;
               map_tel[ntel].inp_id = ti;
               map_to[ni-1][ti] = to;
               ntel++;
               if ( ni == 1 )
                  ntel1++;
               else
                  ntel2++;
            }
         }
      }

   }

   fileclose(map_file);

   printf("Mapping for %d telescopes (%d from input 1, %d from input 2)\n", ntel, ntel1, ntel2);
   if ( ntel1 == 0 || ntel2 == 0 )
   {
      fprintf(stderr,"No event merging needed.\n");
      exit(1);
   }
   
   for (itel=0; itel<H_MAX_TEL; itel++)
      tel_idx_out[itel] = -1;
   for (itel=0; itel<ntel; itel++)
   {
      if ( map_tel[itel].tel_id >= 0 && map_tel[itel].tel_id <= H_MAX_TEL )
         tel_idx_out[map_tel[itel].tel_id] = itel;
      printf("   Tel. no. %d: ID %d in input %d mapped to ID %d\n", 
         itel, map_tel[itel].inp_id, map_tel[itel].ifn, map_tel[itel].tel_id);
   }
   printf("\n");
   
   return 0;
}

/* -------------------- main program ---------------------- */
/** 
 *  @short Main program 
 */

int main (int argc, char **argv)
{
   const char *prg = argv[0];
   IO_BUFFER *iobuf1 = NULL, *iobuf2 = NULL, *iobuf3 = NULL;
   IO_ITEM_HEADER item_header1, item_header2;
   const char *map_fname = NULL;
   const char *input_fname1 = NULL, *input_fname2 = NULL;
   const char *output_fname = NULL;
   int max_events = -1;
   int auto_trgmask = 0;
   char hl[512];

   int iarg=0, jarg=0;
   char *program = argv[0];
   int done1 = 0, done2 = 0;
   int read_from = 1;
   int prev_type1 = 0, prev_type2 = 0;
   int this_type1 = 0, this_type2 = 0;
   
   H_CHECK_MAX();

   push_command_history(argc, argv);

   AllHessData *hsdata1 = calloc(1,sizeof(AllHessData)); /* First input */
   AllHessData *hsdata2 = calloc(1,sizeof(AllHessData)); /* Second input */
   AllHessData *hsdata3 = calloc(1,sizeof(AllHessData)); /* Merged, output */

   if ( hsdata1 == NULL || hsdata2 == NULL || hsdata3 == NULL )
   {
      perror("Allocating data structures");
      exit(1);
   }

   /* Show command line on output */
   if ( getenv("SHOWCOMMAND") != NULL )
   {
      for (iarg=0; iarg<argc; iarg++)
         printf("%s ",argv[iarg]);
      printf("\n");
   }
   if ( argc < 2 )
      syntax(prg);
   for (iarg=1; iarg<argc; iarg++)
   {
      if ( argv[iarg][0] == '-' && argv[iarg][1] != '\0' )
      {
         if ( strcmp(argv[iarg],"--auto-trgmask") == 0 )
         {
       	    auto_trgmask = 1; /* Find matching files with extra trigger mask patterns */
	    continue;
         }
         else if ( strcmp(argv[iarg],"--verbose") == 0 || strcmp(argv[iarg],"-V") == 0 ||
                   strcmp(argv[iarg],"-v") == 0 )
         {
            verbose = 1;
	    continue;
         }
         else if ( strcmp(argv[iarg],"--min-trg-tel") == 0 && iarg+1<argc )
         {
            min_trg = atoi(argv[iarg+1]);
            iarg++;
            continue;
         }
         else if ( strcmp(argv[iarg],"--max-list") == 0 && iarg+1<argc )
         {
            max_list = atoi(argv[iarg+1]);
            iarg++;
            continue;
         }
         else
         {
            if ( strcmp(argv[iarg],"--help") != 0 && strcmp(argv[iarg],"--version") != 0 )
               fprintf(stderr,"Syntax error at '%s'\n", argv[iarg]);
            syntax(program);
            continue;
         }
      }
      else if ( jarg == 0 )
         map_fname = argv[iarg];
      else if ( jarg == 1 )
         input_fname1 = argv[iarg];
      else if ( jarg == 2 )
         input_fname2 = argv[iarg];
      else if ( jarg == 3 )
         output_fname = argv[iarg];
      else
         syntax(program);
      jarg++;
   }

   if ( strcmp(output_fname,"-") == 0 )
   {
      fprintf(stderr,"Output cannot be standard output ('-') because of interfering printf calls.\n");
      fprintf(stderr,"A possible alternative solution may be an output pipe ('|command')\n");
      exit(1);
   }

   snprintf(hl,sizeof(hl),"Mapping file is %s\n", map_fname);
   push_config_history(hl,1);
   snprintf(hl,sizeof(hl),"Input file 1 is %s\n", input_fname1);
   push_config_history(hl,1);
   snprintf(hl,sizeof(hl),"Input file 2 is %s\n", input_fname2);
   push_config_history(hl,1);
   snprintf(hl,sizeof(hl),"Output file is %s\n", output_fname);
   push_config_history(hl,1);

   if ( read_map(map_fname) != 0 )
      exit(1);

   /* Catch INTerrupt and TERMinate signals to stop program */
   signal(SIGINT,stop_signal_function);
   signal(SIGTERM,stop_signal_function);
   interrupted = 0;

   /* Check assumed limits with the ones compiled into the library. */
   H_CHECK_MAX();

   /* Start with quite large buffers; don't worry about wasting memory */
   if ( (iobuf1 = allocate_io_buffer(40000000L)) == NULL ||
        (iobuf2 = allocate_io_buffer(40000000L)) == NULL ||
        (iobuf3 = allocate_io_buffer(40000000L)) == NULL )
   {
      Error("Cannot allocate I/O buffers");
      exit(1);
   }
   /* Allow buffers much larger than during simulation: most likely never exceeded */
   iobuf1->max_length = 400000000L;
   iobuf2->max_length = 400000000L;
   iobuf3->max_length = 800000000L;

   if ( auto_trgmask )
      check_autoload_trgmask(input_fname1, iobuf1, 1);
   if ( (iobuf1->input_file = fileopen(input_fname1,"r")) == NULL )
   {
      perror(input_fname1);
      exit(1);
   }
   if ( auto_trgmask )
      check_autoload_trgmask(input_fname2, iobuf2, 2);
   if ( (iobuf2->input_file = fileopen(input_fname2,"r")) == NULL )
   {
      perror(input_fname2);
      exit(1);
   }
   if ( (iobuf3->output_file = fileopen(output_fname,"w")) == NULL )
   {
      perror(output_fname);
      exit(1);
   }
   write_history(9,iobuf3);

   for (;;) /* Loop over all data in both input files */
   {
      int merge_now = 0, merge_next = 0;

      int rc = 0;
      if ( interrupted )
      {
         rc = 1;
         break;
      }

      if ( max_events > 0 && (events[0] >= max_events || events[1] >= max_events) )
         break;

#ifdef DEBUG_MERGE
printf("Processing loop: r=%d, t1=%d, t2=%d, p1=%d, p2=%d, e1=%ld, e2=%ld, d1=%d, d2=%d, mnow=%d, mnext=%d\n",
   read_from, this_type1, this_type2, prev_type1, prev_type2, event1, event2, done1, done2, merge_now, merge_next);
#endif

      /* Can the request to read the next data block be satisfied? */
      if ( done1 || done2 )
      {
         if ( done1 && done2 )
         {
            printf("Reading from both input files completed.\n");
            if ( this_type1 == 100 && this_type2 == 100 )
            {
               merge_data_from_io_block(iobuf1, &item_header1, 1, hsdata1, hsdata3, iobuf3);
               /* Ignore the histogram block from the second file */
               this_type1 = this_type2 = 0;
            }
            if ( this_type1 != 0 || this_type2 != 0 )
            {
               printf("There is still some unprocessed data:\n");
               print_process_status(prev_type1,this_type1,prev_type2,this_type2);
            }
            break;
         }
         if ( (read_from == 1 && done1) || (read_from == 2 && done2) )
         {
            printf("Should read from input file %d but no more data available.\n", read_from);
            print_process_status(prev_type1,this_type1,prev_type2,this_type2);
            break;
         }
      }

      if ( (read_from == 1 && this_type1 != 0) || (read_from == 2 && this_type2 != 0) )
      {
         fprintf(stderr,"Should read from input file %d but have not processed previous block yet.\n", read_from);
         print_process_status(prev_type1,this_type1,prev_type2,this_type2);
         break;
      }

      /* Should we read a data block from the first file? */
      if ( read_from == 1 && !done1 )
      {
         /* Find and read the next block of data. */
         /* In case of problems with the data, just give up. */
         if ( find_io_block(iobuf1,&item_header1) != 0 )
         {
            rc = 2;
            done1 = 1; /* No more data to come from the first file. */
            read_from = 2;
            continue;
         }
         else if ( read_io_block(iobuf1,&item_header1) != 0 )
         {
            rc = -1; /* An unexpected failure to read the data. We better give up. */
            break;
         }
         this_type1 = item_header1.type;
#ifdef DEBUG_MERGE
printf("Type1 = %d, ID = %ld\n", this_type1, item_header1.ident);
#endif
         if ( (this_type1 == 70 || this_type1 == 2000 || 
               this_type1 == 2001 || this_type1 == 1212) && 
              this_type2 == 0 ) /* switch to second file */
         {
            read_from = 2;
            continue;
         }
         switch ( this_type1 )
         {
            case 70:   /* Immediate output without other actions */
            case 2000: /* IO_TYPE_HESS_RUNHEADER: Merging without immediate output (waits for file 2) */
            case 2001: /* IO_TYPE_HESS_MCRUNHEADER: Merging without immediate output (waits for file 2) */
            case 1212: /* IO_TYPE_MC_INPUTCFG: Immediate output (once) without other actions */
               if ( this_type1 == 2000 )
                  run1 = item_header1.ident;
               merge_now = 1;
               if ( this_type1 == this_type2 )
                  merge_next = 2;
               read_from = 1;
               if ( this_type1 == 2001 )
               {
                  merge_next = 2;
                  read_from = 2;
               }
               event1 = -2;
               break;

            /* These can all be merged with immediate output, unless file 2 is behind: */
            case 2002: /* IO_TYPE_HESS_CAMSETTINGS */
            case 2003: /* IO_TYPE_HESS_CAMORGAN */
            case 2004: /* IO_TYPE_HESS_PIXELSET */
            case 2005: /* IO_TYPE_HESS_PIXELDISABLE */
            case 2006: /* IO_TYPE_HESS_CAMSOFTSET */
            case 2007: /* IO_TYPE_HESS_POINTINGCOR */
            case 2008: /* IO_TYPE_HESS_TRACKSET */
            case 2022: /* IO_TYPE_HESS_TEL_MONI */
            case 2023: /* IO_TYPE_HESS_LASCAL */
               if ( this_type2 == 70 || this_type2 == 2000 || 
                    this_type2 == 2001 || this_type2 == 1212 )
               {
                  merge_now = 2;
                  read_from = 2;
               }
               else
               {
                  merge_now = 1;
                  // merge_next = 2;
                  read_from = 1;
               }
               if ( this_type1 < 2020 )
                  event1 = -1;
               break;

            case 2020: /* IO_TYPE_HESS_MC_SHOWER */
            case 2021: /* IO_TYPE_HESS_MC_EVENT */
            case 2026: /* IO_TYPE_HESS_MC_PE_SUM */
            case 1204: /* IO_TYPE_MC_TELARRAY */
            case 2010: /* IO_TYPE_HESS_EVENT */
               if ( this_type1 == 2020 )
                  event1 = item_header1.ident * 100; /* shower number */
               else
                  event1 = item_header1.ident; /* event number */
               if ( this_type2 == 70 || this_type2 == 1212 || 
                    (this_type2 >= 2000 && this_type2 <= 2008 ) ||
                    this_type2 == 2022 || this_type2 == 2023 )
               {
                  merge_now = 2; /* File 2 is far behind and needs to catch up */
                  read_from = 2;
               }
               else if ( event2 < event1 ) /* File 2 may be behind */
               {
                  merge_now = 2;
                  read_from = 2;
               }
               else if ( event2 > event1 ) /* File 1 is behind */
               {
                  merge_now = 1;
                  read_from = 1;
               }
               else /* event2 == event1 */
               {
                  if ( this_type1 == this_type2 ) /* We caught up with file 2. Merge both. */
                  {
                     merge_now = 1;
                     merge_next = 2;
                     read_from = 1;
                  }
                  else if ( this_type1 == 2020 ) /* First the MC shower block */
                  {
                     merge_now = 1;
                     read_from = 1;
                  }
                  else if ( this_type2 == 2020 ) /* First the MC shower block */
                  {
                     merge_now = 2;
                     read_from = 2;
                  }
                  else if ( this_type1 == 2021 ) /* Next the MC event block */
                  {
                     merge_now = 1;
                     read_from = 1;
                  }
                  else if ( this_type2 == 2021 ) /* Next the MC event block */
                  {
                     merge_now = 2;
                     read_from = 2;
                  }
                  else if ( this_type1 == 2026 || this_type1 == 1204 ) /* then anything else than the event */
                  {
                     merge_now = 1;
                     read_from = 1;
                  }
                  else if ( this_type2 == 2026 || this_type2 == 1204 ) /* then anything else than the event */
                  {
                     merge_now = 2;
                     read_from = 2;
                  }
                  else
                  {
                     merge_now = 1;
                     read_from = 1;
                  }
               }
               break;

            case 100: /* Probably the last block of file 1 */
               merge_now = 0;
               if ( this_type2 == 100 )
               {
                  merge_now = 2;
                  this_type1 = 0;
               }
               else if ( done2 )
                  merge_now = 1;
               read_from = 2;
               break;

            default: /* Probably complaining about unknown block type but try anyway */
               merge_now = 1;
               read_from = 2;
               break;
         }
      }
      /* Or should we read a data block from the second file? */
      else if ( read_from == 2 && !done2 && rc >= 0 )
      {
         /* Find and read the next block of data. */
         /* In case of problems with the data, just give up. */
         if ( find_io_block(iobuf2,&item_header2) != 0 )
         {
            rc = 2;
            done2 = 1; /* No more data to come from the first file. */
            continue;
         }
         else if ( read_io_block(iobuf2,&item_header2) != 0 )
         {
            rc = -1; /* An unexpected failure to read the data. We better give up. */
            break;
         }
         this_type2 = item_header2.type;
#ifdef DEBUG_MERGE
printf("Type2 = %d, ID = %ld\n", this_type2, item_header2.ident);
#endif
         switch ( this_type2 )
         {
            case 70:   /* Immediate output without other actions */
            case 2000: /* IO_TYPE_HESS_RUNHEADER: Merging with immediate output (must have merged run header of file 1) */
            case 2001: /* IO_TYPE_HESS_MCRUNHEADER: Merging without immediate output (must have MC merged run header of file 1) */
            case 1212: /* Immediate output (once) without other actions */
               if ( this_type2 == 2000 )
                  run2 = item_header2.ident;
               if ( this_type1 == this_type2 ) /* Caught up with file 1 */
               {
                  merge_now = 1;
                  merge_next = 2;
                  read_from = 1;
               }
               else if ( this_type1 == 70 || this_type1 == 2000 || this_type1 == 2001 ) /* First process data from file 1 */
               {
                  merge_now = 1;
                  read_from = 1;
               }
               else
               {
                  merge_now = 2;
                  read_from = 2;
               }
               event2 = -2;
               break;

            /* These can all be merged with immediate output, unless file 1 is behind: */
            case 2002: /* IO_TYPE_HESS_CAMSETTINGS */
            case 2003: /* IO_TYPE_HESS_CAMORGAN */
            case 2004: /* IO_TYPE_HESS_PIXELSET */
            case 2005: /* IO_TYPE_HESS_PIXELDISABLE */
            case 2006: /* IO_TYPE_HESS_CAMSOFTSET */
            case 2007: /* IO_TYPE_HESS_POINTINGCOR */
            case 2008: /* IO_TYPE_HESS_TRACKSET */
            case 2022: /* IO_TYPE_HESS_TEL_MONI */
            case 2023: /* IO_TYPE_HESS_LASCAL */
               if ( this_type1 == 70 || this_type1 == 2000 || 
                    this_type1 == 2001 || this_type1 == 1212 ||
                    this_type1 == 2002 || this_type1 == 2003 || this_type1 == 2004 || this_type1 == 2005 ||
                    this_type1 == 2006 || this_type1 == 2007 || this_type1 == 2008 ||
                    ((this_type2 == 2022 || this_type2 == 2023) && (this_type1 == 2022 || this_type1 == 2023)) )
               {
                  merge_now = 1; /* Process file one as far as we get */
                  read_from = 1;
               }
               else
               {
                  merge_now = 2;
                  read_from = 2;
               }
               if ( this_type2 < 2020 )
                  event2 = -1;
               break;

            case 2020: /* IO_TYPE_HESS_MC_SHOWER */
            case 2021: /* IO_TYPE_HESS_MC_EVENT */
            case 2026: /* IO_TYPE_HESS_MC_PE_SUM */
            case 1204:
            case 2010: /* IO_TYPE_HESS_EVENT */
               if ( this_type2 == 2020 )
                  event2 = item_header2.ident * 100; /* shower number */
               else
                  event2 = item_header2.ident; /* event number */
               if ( this_type1 == 70 || this_type1 == 1212 || 
                    (this_type1 >= 2000 && this_type1 <= 2008 ) ||
                    this_type1 == 2022 || this_type1 == 2023 )
               {
                  merge_now = 1; /* File 1 is far behind and needs to catch up */
                  read_from = 1;
               }
               else if ( event1 < event2 ) /* File 1 may be behind */
               {
                  merge_now = 1;
                  read_from = 1;
               }
               else if ( event1 > event2 ) /* File 2 is behind */
               {
                  merge_now = 2;
                  read_from = 2;
               }
               else /* event2 == event1 */
               {
                  if ( this_type1 == this_type2 ) /* We caught up with file 1. Merge both. */
                  {
                     merge_now = 1;
                     merge_next = 2;
                     read_from = 1;
                  }
                  else if ( this_type1 == 2020 ) /* First the MC shower block */
                  {
                     merge_now = 1;
                     read_from = 1;
                  }
                  else if ( this_type2 == 2020 ) /* First the MC shower block */
                  {
                     merge_now = 2;
                     read_from = 2;
                  }
                  else if ( this_type1 == 2021 ) /* Next the MC event block */
                  {
                     merge_now = 1;
                     read_from = 1;
                  }
                  else if ( this_type2 == 2021 ) /* Next the MC event block */
                  {
                     merge_now = 2;
                     read_from = 2;
                  }
                  else if ( this_type1 == 2026 || this_type1 == 1204 ) /* then anything else than the event */
                  {
                     merge_now = 1;
                     read_from = 1;
                  }
                  else if ( this_type2 == 2026 || this_type2 == 1204 ) /* then anything else than the event */
                  {
                     merge_now = 2;
                     read_from = 2;
                  }
                  else
                  {
                     merge_now = 2;
                     read_from = 2;
                  }
               }
               break;

            case 100: /* Probably the last block of file 2 */
               merge_now = 2;
               if ( this_type1 == 100 )
                  this_type1 = 0;
               read_from = 1;
               break;

            default: /* Probably complaining about unknown block type but try anyway */
               merge_now = 2;
               read_from = 1;
               break;
         }
      }

      if ( merge_now == 1 && this_type1 != 0 )
      {
         merge_data_from_io_block(iobuf1, &item_header1, 1, hsdata1, hsdata3, iobuf3);
         prev_type1 = this_type1;
         this_type1 = 0;
         merge_now = 0;
      }
      else if ( merge_now == 2 && this_type2 != 0 )
      {
         merge_data_from_io_block(iobuf2, &item_header2, 2, hsdata2, hsdata3, iobuf3);
         prev_type2 = this_type2;
         this_type2 = 0;
         merge_now = 0;
      }
      if ( merge_next == 1 && this_type1 != 0 )
      {
         merge_data_from_io_block(iobuf1, &item_header1, 1, hsdata1, hsdata3, iobuf3);
         prev_type1 = this_type1;
         this_type1 = 0;
         merge_next = 0;
      }
      else if ( merge_next == 2 && this_type2 != 0 )
      {
         merge_data_from_io_block(iobuf2, &item_header2, 2, hsdata2, hsdata3, iobuf3);
         prev_type2 = this_type2;
         this_type2 = 0;
         merge_next = 0;
      }

   }
   
   fileclose(iobuf1->input_file);
   fileclose(iobuf2->input_file);
   fileclose(iobuf3->output_file);
   /* Avoid cppcheck false positives: */
   free(hsdata1);
   free(hsdata2);
   free(hsdata3);

   return 0;
}

/** @} */
