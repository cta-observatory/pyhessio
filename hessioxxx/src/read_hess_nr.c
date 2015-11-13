/* ============================================================================

Copyright (C) 2008, 2009, 2010  Konrad Bernloehr

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

/** @file read_hess_nr.c
 *  @short A skeleton program reading H.E.S.S. data.
 *
 *  As a skeleton for programs reading H.E.S.S. data in eventio format,
 *  this program reads the whole range of hessio item types into a
 *  single tree of data structures but normally does nothing with the data.
 *
 *  It can be instructed, though, to create nice camera images similar to
 *  those generated in sim_hessarray.
 *
 *
@verbatim
Syntax: read_hess_nr [ options ] [ - | input_fname ... ]
Options:
   -p ps_filename  (Write a PostScript file with camera images.)
   -r level        (Reconstruction level not fully used in this program version.)
                   level >= 1: show parameters from sim_hessarray.
   -v              (More verbose output)
   -q              (Much more quiet output)
   -s              (Show data explained)
   -S              (Show data explained, including raw data)
   --history (-h)  (Show contents of history data block)
   -i              (Ignore unknown data block types)
   -u              (Call user-defined analysis function)
   --powerlaw x    (Use this spectral index for events weights in output.)
                   (Default spectral index is -2.7)
   --max-events n  (Skip remaining data after so many triggered events.)
@endverbatim
 *
 *  @author  Konrad Bernloehr
 *  @date    @verbatim CVS $Date: 2011/07/21 16:07:26 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.16 $ @endverbatim
 */

/** @defgroup read_hess_nr_c The read_hess_nr program */
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
#include "rec_tools.h"
#include "warning.h"
#include "camera_image.h"

#ifndef _UNUSED_
# ifdef __GNUC__
#  define _UNUSED_ __attribute__((unused))
# else
#  define _UNUSED_
# endif
#endif

/** The factor needed to transform from mean p.e. units to units of the single-p.e. peak:
    Depends on the collection efficiency, the asymmetry of the single p.e. amplitude 
    distribution and the electronic noise added to the signals. */
#define CALIB_SCALE 0.92

/* ---------------------- calibrate_pixel_amplitude ----------------------- */

/** Calibrate a single pixel amplitude, for cameras with two gains per pixel.
 *  This version does not include amplitude clipping nor obtaining amplitudes
 *  from the pixel timing data structure.
 *
 * @return Pixel amplitude in peak p.e. units.
 */

double calibrate_pixel_amplitude(AllHessData *hsdata, int itel, int ipix, 
   int dummy, double cdummy);

double calibrate_pixel_amplitude(AllHessData *hsdata, int itel, int ipix, 
   _UNUSED_ int dummy, _UNUSED_ double cdummy)
{
   int i = ipix, npix, significant, hg_known;
   double npe, sig_hg, npe_hg;
#if (H_MAX_GAINS >= 2)
   double sig_lg, npe_lg;
   int lg_known;
#endif
   AdcData *raw;

   if ( hsdata == NULL || itel < 0 || itel >= H_MAX_TEL )
      return 0.;
   npix = hsdata->camera_set[itel].num_pixels;
   if ( ipix < 0 || ipix >= npix )
      return 0.;
   raw = hsdata->event.teldata[itel].raw;
   if ( raw == NULL )
      return 0.;
   if ( ! raw->known )
      return 0.;
 
   significant = hsdata->event.teldata[itel].raw->significant[i];

   hg_known = hsdata->event.teldata[itel].raw->adc_known[HI_GAIN][i];
   sig_hg = hg_known ? (hsdata->event.teldata[itel].raw->adc_sum[HI_GAIN][i] -
        hsdata->tel_moni[itel].pedestal[HI_GAIN][i]) : 0.;
   npe_hg = sig_hg * hsdata->tel_lascal[itel].calib[HI_GAIN][i];

#if (H_MAX_GAINS >= 2 )
   lg_known = hsdata->event.teldata[itel].raw->adc_known[LO_GAIN][i];
   sig_lg = lg_known ? (hsdata->event.teldata[itel].raw->adc_sum[LO_GAIN][i] -
        hsdata->tel_moni[itel].pedestal[LO_GAIN][i]) : 0.;
   npe_lg = sig_lg * hsdata->tel_lascal[itel].calib[LO_GAIN][i];
#endif

   if ( !significant ) 
      npe = 0.;
#if (H_MAX_GAINS >= 2 )
   /* FIXME: we should be more flexible here: */
   else if ( hg_known && sig_hg < 10000 && sig_hg > -1000 )
      npe = npe_hg;
   else if ( raw->num_gains >= 2 )
      npe = npe_lg;
#endif
   else
      npe = npe_hg;

   /* npe is in units of 'mean photo-electrons'. */
   /* We convert to experimentalist's */
   /* 'peak photo-electrons' now. */
   return CALIB_SCALE * npe;
}

#include <signal.h>

void stop_signal_function (int isig);

static int interrupted;

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
   
   signal(SIGINT,SIG_DFL);
   signal(SIGTERM,SIG_DFL);
}

static void show_run_summary(AllHessData *hsdata, int nev, int ntrg, double plidx,
   double wsum_all, double wsum_trg, double rmax_x, double rmax_y, double rmax_r);

static void show_run_summary(AllHessData *hsdata, int nev, int ntrg, double plidx,
   double wsum_all, double wsum_trg, double rmax_x, double rmax_y, double rmax_r)
{
   static int explained = 0;

   if ( ! explained )
   {
      printf("\n#@; Column 1: Run number\n"
             "#@;        2: ID of primary particle\n"
             "#@;        3: Number of events, total\n"
             "#@;        4: Number of events, triggered\n"
             "#@;        5: Altitude (mean) [deg.]\n"
             "#@;        6: Azimuth (mean) [deg.]\n"
             "#@;        7: Cone (max) [deg.]\n"
             "#@;        8: Lower limit of energy range [TeV]\n"
             "#@;        9: Upper limit of energy range [TeV]\n"
             "#@;       10: Spectral index in simulation\n"
             "#@;       11: Spectral index in weighting\n"
             "#@;       12: Weighted sum of events, total\n"
             "#@;       13: Weighted sum of events, triggered\n"
             "#@;       14: Maximum horizontal core distance in X\n"
             "#@;       15: Maximum horizontal core distance in Y\n"
             "#@;       16: Maximum core distance in shower plane\n"
             "#@;       17: Supposed maximum core distance\n");
      explained = 1;
   }
   printf("\n@; %d %d %d %d   %5.2f %5.2f %4.2f    %6.4f %6.4f %5.3f %5.3f   %f %f   %3.1f %3.1f %3.1f %3.1f\n", 
     hsdata->run_header.run, hsdata->mc_shower.primary_id, nev, ntrg,
     0.5*(180./M_PI)*
      (hsdata->mc_run_header.alt_range[0]+hsdata->mc_run_header.alt_range[1]),
     0.5*(180./M_PI)*
      (hsdata->mc_run_header.az_range[0]+hsdata->mc_run_header.az_range[1]),
     (180./M_PI)*hsdata->mc_run_header.viewcone[1],
     hsdata->mc_run_header.E_range[0], hsdata->mc_run_header.E_range[1],
     hsdata->mc_run_header.spectral_index, plidx,
     wsum_all, wsum_trg, rmax_x, rmax_y, rmax_r,
     hsdata->mc_run_header.core_range[1]);
}

/** Show program syntax */

static void syntax (char *program);

static void syntax (char *program)
{
   printf("Syntax: %s [ options ] [ - | input_fname ... ]\n",program);
   printf("Options:\n");
   printf("   -p ps_filename  (Write a PostScript file with camera images.)\n");
   printf("   -r level        (Reconstruction level not fully used in this program version.)\n");
   printf("                   level >= 1: show parameters from sim_hessarray.\n");
   printf("   -v              (More verbose output)\n");
   printf("   -q              (Much more quiet output)\n");
   printf("   -s              (Show data explained)\n");
   printf("   -S              (Show data explained, including raw data)\n");
   printf("   --history (-h)  (Show contents of history data block)\n");
   printf("   -i              (Ignore unknown data block types)\n");
   printf("   --powerlaw x    (Use this spectral index for events weights in output.)\n");
   printf("                   (Default spectral index is -2.7)\n");
   printf("   --max-events n  (Skip remaining data after so many triggered events.)\n");

   exit(1);
}

/* -------------------- main program ---------------------- */
/** 
 *  @short Main program 
 *
 *  Main program function of read_hess.c program.
 */

int main (int argc, char **argv)
{
   IO_BUFFER *iobuf = NULL;
   IO_ITEM_HEADER item_header;
   const char *input_fname = NULL;
   int itel, rc = 0;
   int tel_id;
   const char *ps_fname = "none";
   int verbose = 0, ignore = 0, quiet = 0;
   int reco_flag = 0;
   int iprint_mc = 0, iprint_sim = 0;
   double plidx = -2.7;
   double wsum_all = 0., wsum_trg = 0.;
   double rmax_x = 0., rmax_y = 0., rmax_r = 0., rs;
   int ntel_trg = 0, min_tel_trg = 0;
   int nev = 0, ntrg = 0;
   char *program = argv[0];
   int showdata = 0, showhistory = 0;
   size_t events = 0, max_events = 0;
   int iarg;
   
   static AllHessData *hsdata;
   
   /* Show command line on output */
   if ( getenv("SHOWCOMMAND") != NULL )
   {
      for (iarg=0; iarg<argc; iarg++)
         printf("%s ",argv[iarg]);
      printf("\n");
   }

   /* Catch INTerrupt and TERMinate signals to stop program */
   signal(SIGINT,stop_signal_function);
   signal(SIGTERM,stop_signal_function);
   interrupted = 0;

   /* Check assumed limits with the ones compiled into the library. */
   H_CHECK_MAX();

   if ( argc < 2 )
      input_fname = "iact.out";

   if ( (iobuf = allocate_io_buffer(1000000L)) == NULL )
   {
      Error("Cannot allocate I/O buffer");
      exit(1);
   }
   iobuf->max_length = 100000000L;

   /* Command line options */
   while ( argc > 1 )
   {
      if ( strcmp(argv[1],"-p") == 0 && argc > 2 )
      {
	 ps_fname = argv[2];
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strncmp(argv[1],"-p",2) == 0 && strlen(argv[1]) > 2 )
      {
	 ps_fname = argv[1]+2;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"-r") == 0 && argc > 2 )
      {
	 reco_flag = atoi(argv[2]);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strncmp(argv[1],"-r",2) == 0 && strlen(argv[1]) > 2 )
      {
	 reco_flag = atoi(argv[1]+2);
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"-i") == 0 )
      {
       	 ignore = 1;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"-v") == 0 )
      {
       	 verbose = 1;
         quiet = 0;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"-q") == 0 )
      {
       	 quiet = 1;
         verbose = 0;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--history") == 0 )
      {
       	 showhistory = 1;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"-s") == 0 )
      {
       	 showdata = 1;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"-S") == 0 )
      {
       	 showdata = 1;
#ifndef __cplusplus
         putenv((char *)"PRINT_VERBOSE=1");
#else
         putenv(const_cast<char *>("PRINT_VERBOSE=1"));
#endif
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--powerlaw") == 0 && argc > 2 )
      {
       	 plidx = atof(argv[2]);
         if ( plidx > 0. )
            plidx = -plidx;
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--max-events") == 0 && argc > 2 )
      {
       	 max_events = atol(argv[2]);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--help") == 0 )
      {
        printf("\nread_hess: A program for viewing sim_hessarray data.\n\n");
        syntax(program);
      }
      else if ( argv[1][0] == '-' && argv[1][1] != '\0' )
      {
        printf("Syntax error at '%s'\n", argv[1]);
        syntax(program);
      }
      else
        break;
   }
   
   if ( verbose && !quiet )
      showhistory = 1;
    
   /* Now go over rest of the command line */
   while ( argc > 1 || input_fname != NULL )
   {
    if ( interrupted )
      break;
    if ( argc > 1 )
    {
      if ( argv[1][0] == '-' && argv[1][1] != '\0' )
         syntax(program);
      else
      {
	 input_fname = argv[1];
	 argc--;
	 argv++;
      }
    }
    if ( strcmp(input_fname ,"-") == 0 )
      iobuf->input_file = stdin;
    else if ( (iobuf->input_file = fileopen(input_fname,READ_BINARY)) == NULL )
    {
      perror(input_fname);
      Error("Cannot open input file.");
      break;
    }

    fflush(stdout);
    fprintf(stderr,"%s\n",input_fname);
    printf("\nInput file '%s' has been opened.\n",input_fname);
    input_fname = NULL;

#ifdef READ_NUM_SLICES_FIRST
    {
      TelMoniData tmp_moni;
      item_header.type = 0;
      do
      {
	 if ( find_io_block(iobuf,&item_header) != 0 )
            break;
	 if ( read_io_block(iobuf,&item_header) != 0 )
            break;
      } while ( item_header.type != IO_TYPE_HESS_TEL_MONI );
      if ( item_header.type == IO_TYPE_HESS_TEL_MONI )
      {
         tmp_moni.tel_id = 
            (item_header.ident & 0xff) |
            ((item_header.ident & 0x3f000000) >> 16);
	 printf("Found monitor block for telescope %d\n",tmp_moni.tel_id);
      	 read_hess_tel_monitor(iobuf,&tmp_moni);
	 printf("%d slices!\n",tmp_moni.num_ped_slices);
      }
      else
         printf("No monitor block found!\n");
      // exit(0);
      if ( iobuf->input_file != stdin )
         rewind(iobuf->input_file);
      else
      {
      	 fprintf(stderr,"This will never work with anything else but plain disk files!\n");
	 exit(1);
      }
    }
#endif

    for (;;) /* Loop over all data in the input file */
    {
      if ( interrupted )
         break;

      /* Find and read the next block of data. */
      /* In case of problems with the data, just give up. */
      if ( find_io_block(iobuf,&item_header) != 0 )
         break;
      if ( max_events > 0 && events >= max_events )
      {
         if ( iobuf->input_file != stdin )
            break;
         if ( skip_io_block(iobuf,&item_header) != 0 )
            break;
         continue;
      }
      if ( read_io_block(iobuf,&item_header) != 0 )
         break;

      if ( hsdata == NULL && 
           item_header.type > IO_TYPE_HESS_RUNHEADER &&
           item_header.type < IO_TYPE_HESS_RUNHEADER + 200)
      {
         fprintf(stderr,"Trying to read event data before run header.\n");
         fprintf(stderr,"Skipping this data block.\n");
         continue;
      }
      
      /* What did we actually get? */
      switch ( (int) item_header.type )
      {
         /* =================================================== */
         case IO_TYPE_HESS_RUNHEADER:
            /* Summary of a preceding run in the same file ? */
            if ( !quiet && hsdata != NULL && hsdata->run_header.run > 0 )
               show_run_summary(hsdata,nev,ntrg,plidx,wsum_all,wsum_trg,
                  rmax_x,rmax_y,rmax_r);
            else if ( nev > 0 )
	       printf("%d of %d events triggered.\n", ntrg, nev);
            nev = ntrg = 0;
            wsum_all = wsum_trg = 0.;
            /* Structures might be allocated from previous run */
            if ( hsdata != NULL )
            {
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
		  }
               }
               /* Free main structure */
               free(hsdata);
               hsdata = NULL;
               
               /* Perhaps some cleaning needed in ROOT as well ... */
               
            }
            hsdata = (AllHessData *) calloc(1,sizeof(AllHessData));
            if ( (rc = read_hess_runheader(iobuf,&hsdata->run_header)) < 0 )
            {
               Warning("Reading run header failed.");
               exit(1);
            }
            if ( !quiet )
               printf("Reading simulated data for %d telescope(s)\n",hsdata->run_header.ntel);
	    if ( verbose || rc != 0 )
               printf("read_hess_runheader(), rc = %d\n",rc);
            fprintf(stderr,"\nStarting run %d\n",hsdata->run_header.run);
            if ( showdata )
               print_hess_runheader(iobuf);

            for (itel=0; itel<hsdata->run_header.ntel; itel++)
            {
               tel_id = hsdata->run_header.tel_id[itel];
               hsdata->camera_set[itel].tel_id = tel_id;
               hsdata->camera_org[itel].tel_id = tel_id;
               hsdata->pixel_set[itel].tel_id = tel_id;
               hsdata->pixel_disabled[itel].tel_id = tel_id;
               hsdata->cam_soft_set[itel].tel_id = tel_id;
               hsdata->tracking_set[itel].tel_id = tel_id;
               hsdata->point_cor[itel].tel_id = tel_id;
               hsdata->event.num_tel = hsdata->run_header.ntel;
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
            }
            break;

         /* =================================================== */
         case IO_TYPE_HESS_MCRUNHEADER:
            rc = read_hess_mcrunheader(iobuf,&hsdata->mc_run_header);
	    if ( verbose || rc != 0 )
               printf("read_hess_mcrunheader(), rc = %d\n",rc);
            if ( showdata )
               print_hess_mcrunheader(iobuf);
            break;

         /* =================================================== */
	 case IO_TYPE_MC_INPUTCFG:
            {
               struct linked_string corsika_inputs;
               corsika_inputs.text = NULL;
               corsika_inputs.next = NULL;
	       read_input_lines(iobuf,&corsika_inputs);
	       if ( corsika_inputs.text != NULL )
	       {
	          struct linked_string *xl = NULL, *xln = NULL;
                  if ( ! quiet )
	             printf("\nCORSIKA was run with the following input lines:\n");
	          for (xl = &corsika_inputs; xl!=NULL; xl=xln)
                  {
                     if ( ! quiet )
		        printf("   %s\n",xl->text);
		     free(xl->text);
		     xl->text = NULL;
		     xln = xl->next;
		     xl->next = NULL;
		     if ( xl != &corsika_inputs )
		        free(xl);
                  }
               }
            }
            break;

         /* =================================================== */
         case 70: /* How sim_hessarray was run and how it was configured. */
            if ( showhistory )
               list_history(iobuf,NULL);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_CAMSETTINGS:
            tel_id = item_header.ident; // Telescope ID is in the header
            if ( (itel = find_tel_idx(tel_id)) < 0 )
            {
               char msg[256];
               snprintf(msg,sizeof(msg)-1,
                  "Camera settings for unknown telescope %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_camsettings(iobuf,&hsdata->camera_set[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_camsettings(), rc = %d\n",rc);
            if ( showdata )
               print_hess_camsettings(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_CAMORGAN:
            tel_id = item_header.ident; // Telescope ID is in the header
            if ( (itel = find_tel_idx(tel_id)) < 0 )
            {
               char msg[256];
               snprintf(msg,sizeof(msg)-1,
                  "Camera organisation for unknown telescope %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_camorgan(iobuf,&hsdata->camera_org[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_camorgan(), rc = %d\n",rc);
            if ( showdata )
               print_hess_camorgan(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_PIXELSET:
            tel_id = item_header.ident; // Telescope ID is in the header
            if ( (itel = find_tel_idx(tel_id)) < 0 )
            {
               char msg[256];
               snprintf(msg,sizeof(msg)-1,
                  "Pixel settings for unknown telescope %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_pixelset(iobuf,&hsdata->pixel_set[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_pixelset(), rc = %d\n",rc);
            if ( showdata )
               print_hess_pixelset(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_PIXELDISABLE:
            tel_id = item_header.ident; // Telescope ID is in the header
            if ( (itel = find_tel_idx(tel_id)) < 0 )
            {
               char msg[256];
               snprintf(msg,sizeof(msg)-1,
                  "Pixel disable block for unknown telescope %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_pixeldis(iobuf,&hsdata->pixel_disabled[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_pixeldis(), rc = %d\n",rc);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_CAMSOFTSET:
            tel_id = item_header.ident; // Telescope ID is in the header
            if ( (itel = find_tel_idx(tel_id)) < 0 )
            {
               char msg[256];
               snprintf(msg,sizeof(msg)-1,
                  "Camera software settings for unknown telescope %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_camsoftset(iobuf,&hsdata->cam_soft_set[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_camsoftset(), rc = %d\n",rc);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_POINTINGCOR:
            tel_id = item_header.ident; // Telescope ID is in the header
            if ( (itel = find_tel_idx(tel_id)) < 0 )
            {
               char msg[256];
               snprintf(msg,sizeof(msg)-1,
                  "Pointing correction for unknown telescope %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_pointingcor(iobuf,&hsdata->point_cor[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_pointingco(), rc = %d\n",rc);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_TRACKSET:
            tel_id = item_header.ident; // Telescope ID is in the header
            if ( (itel = find_tel_idx(tel_id)) < 0 )
            {
               char msg[256];
               snprintf(msg,sizeof(msg)-1,
                  "Tracking settings for unknown telescope %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_trackset(iobuf,&hsdata->tracking_set[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_trackset(), rc = %d\n",rc);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_EVENT:
            rc = read_hess_event(iobuf,&hsdata->event,-1);
	    if ( verbose || rc != 0 )
               printf("read_hess_event(), rc = %d\n",rc);
            events++;
            if ( showdata )
               print_hess_event(iobuf);
            /* Count number of telescopes (still) present in data and triggered */
            ntel_trg = 0;
            for (itel=0; itel<hsdata->run_header.ntel; itel++)
               if ( hsdata->event.teldata[itel].known )
               {
                  /* If non-triggered telescopes record data (like HEGRA),
                     we may have to check the central trigger bit as well,
                     but ignore this for now. */
                  ntel_trg++;
               }
	    if ( hsdata->event.shower.known )
	       hsdata->event.shower.num_trg = ntel_trg;
            if ( ntel_trg < min_tel_trg )
               continue;
            wsum_trg += pow(hsdata->mc_shower.energy,
               plidx-hsdata->mc_run_header.spectral_index);
            ntrg++;

            if ( reco_flag && !quiet )
            {
               if ( iprint_sim++ == 0 )
               {
                  printf("#@: Lines starting with '@:' contain the following columns from sim_hessarray:\n"
                         "#@:  (1): event\n"
                         "#@:  (2): number of telescopes triggered\n"
                         "#@:  (3): number of images used in shower reconstruction\n"
                         "#@:  (4): results bit pattern\n"
                         "#@:  (5): shower azimuth [deg] (with bit 0)\n"
                         "#@:  (6): shower altitude [deg] (with bit 0)\n"
                         "#@:  (7): angle between reconstructed and true direction [deg]\n"
                         "#@:  (8): core position x [m] (with bit 2)\n"
                         "#@:  (9): core position y [m] (with bit 2)\n"
                         "#@: (10): horizontal displacement between reconstructed and true core [m]\n"
                         "#@: (11): MSCL [deg] (with bit 4)\n"
                         "#@: (12): MSCW [deg] (with bit 4)\n"
                         "#@: (13): energy [TeV] (with bit 6)\n"
                         "#@: (14): Xmax [g/cm^2] (with bit 8)\n");
                  printf("#@+ Lines starting with '@+' contain the following columns from sim_hessarray:\n"
                         "#@+  (1): event\n"
                         "#@+  (2): telescope\n"
                         "#@+  (3): energy\n"
                         "#@+  (4): core distance to telescope\n"
                         "#@+  (5): image size (amplitude) [p.e.]\n"
                         "#@+  (6): number of pixels in image\n"
                         "#@+  (7): width [deg.]\n"
                         "#@+  (8): length [deg.]\n"
                         "#@+  (9): distance [deg.]\n"
                         "#@+ (10): miss [deg.]\n"
                         "#@+ (11): alpha [deg.]\n"
                         "#@+ (12): orientation [deg.]\n"
                         "#@+ (13): direction [deg.]\n"
                         "#@+ (14): image c.o.g. x [deg.]\n"
                         "#@+ (15): image c.o.g. y [deg.]\n"
                         "#@+ (16): Xmax [g/cm^2]\n"
                         "#@+ (17): Hmax [m]\n"
                         "#@+ (18): Npe (true number of photo-electrons)\n"
                         "#@+ (19-23): Hottest pixel amplitudes)\n");
               }
               
               if ( hsdata->event.shower.known )
                  printf("@: %d %d %d %d   %f %f %f   %f %f %f   %f %f  %f %f\n",
                     hsdata->mc_event.event, 
                     hsdata->event.shower.num_trg, hsdata->event.shower.num_img,
                     hsdata->event.shower.result_bits,
                     hsdata->event.shower.Az*(180./M_PI), hsdata->event.shower.Alt*(180./M_PI),
                     angle_between(hsdata->event.shower.Az, hsdata->event.shower.Alt,
                        hsdata->mc_shower.azimuth, hsdata->mc_shower.altitude) * (180./M_PI),
                     hsdata->event.shower.xc, hsdata->event.shower.yc,
                     sqrt((hsdata->event.shower.xc-hsdata->mc_event.xcore)*
                          (hsdata->event.shower.xc-hsdata->mc_event.xcore) +
                          (hsdata->event.shower.yc-hsdata->mc_event.ycore)*
                          (hsdata->event.shower.yc-hsdata->mc_event.ycore)),
                     hsdata->event.shower.mscl*(180./M_PI), 
                     hsdata->event.shower.mscw*(180./M_PI),
                     hsdata->event.shower.energy, hsdata->event.shower.xmax);

               for (itel=0; itel<hsdata->run_header.ntel; itel++)
                  if ( hsdata->event.teldata[itel].known &&
                       hsdata->event.teldata[itel].num_image_sets > 0 )
                  {
                     double direction = (180./M_PI) * hsdata->event.teldata[itel].img[0].phi;
                     double orientation = (180./M_PI) * atan2(hsdata->event.teldata[itel].img[0].y,
                                                         hsdata->event.teldata[itel].img[0].x);
                     double alpha = (180./M_PI) * (hsdata->event.teldata[itel].img[0].phi -
                        atan2(hsdata->event.teldata[itel].img[0].y,
                              hsdata->event.teldata[itel].img[0].x));
                     if ( orientation < 0. )
                        orientation += 360.;
                     if ( direction < 0. )
                        direction += 360.;
                     while ( alpha < 0. )
                        alpha += 180.;
                     while ( alpha > 180. )
                        alpha -= 180.;
                     if ( alpha > 90. )
                        alpha = 180.-alpha;
                     /* Number of pixels not directly in image struct for older format. Fix it. */
                     if ( hsdata->event.teldata[itel].img[0].pixels == 0 &&
                          hsdata->event.teldata[itel].image_pixels.pixels > 0 )
                        hsdata->event.teldata[itel].img[0].pixels =
                           hsdata->event.teldata[itel].image_pixels.pixels;
                     printf("@+ %d %d %6.4f %7.2f %7.2f %d %7.5f %7.5f %7.5f %7.5f %7.5f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f   %d  %4.2f %4.2f %4.2f %4.2f %4.2f\n",
                        hsdata->mc_event.event, 
                        hsdata->camera_set[itel].tel_id,
                        hsdata->mc_shower.energy, 
                        line_point_distance(hsdata->mc_event.xcore, hsdata->mc_event.ycore,0.,
                                       cos(hsdata->mc_shower.altitude)*cos(hsdata->mc_shower.azimuth),
                                       cos(hsdata->mc_shower.altitude)*sin(-hsdata->mc_shower.azimuth),
                                       sin(hsdata->mc_shower.altitude), 
                                       hsdata->run_header.tel_pos[itel][0],
                                       hsdata->run_header.tel_pos[itel][1], 
                                       hsdata->run_header.tel_pos[itel][2]),
                        CALIB_SCALE*hsdata->event.teldata[itel].img[0].amplitude,
                        hsdata->event.teldata[itel].img[0].pixels,
                        hsdata->event.teldata[itel].img[0].w*(180./M_PI),
                        hsdata->event.teldata[itel].img[0].l*(180./M_PI),
                        sqrt(hsdata->event.teldata[itel].img[0].x*hsdata->event.teldata[itel].img[0].x +
                          hsdata->event.teldata[itel].img[0].y*hsdata->event.teldata[itel].img[0].y) *
                          (180./M_PI),
                        -1., alpha, orientation, direction,
                        hsdata->event.teldata[itel].img[0].x*(180./M_PI),
                        hsdata->event.teldata[itel].img[0].y*(180./M_PI),
                        hsdata->mc_shower.xmax, hsdata->mc_shower.hmax,
                        hsdata->mc_event.mc_pesum.num_pe[itel],
                        hsdata->event.teldata[itel].img[0].num_hot>0 ? CALIB_SCALE*hsdata->event.teldata[itel].img[0].hot_amp[0] : -1.,
                        hsdata->event.teldata[itel].img[0].num_hot>1 ? CALIB_SCALE*hsdata->event.teldata[itel].img[0].hot_amp[1] : -1.,
                        hsdata->event.teldata[itel].img[0].num_hot>2 ? CALIB_SCALE*hsdata->event.teldata[itel].img[0].hot_amp[2] : -1.,
                        hsdata->event.teldata[itel].img[0].num_hot>3 ? CALIB_SCALE*hsdata->event.teldata[itel].img[0].hot_amp[3] : -1.,
                        hsdata->event.teldata[itel].img[0].num_hot>4 ? CALIB_SCALE*hsdata->event.teldata[itel].img[0].hot_amp[4] : -1.);
                     if ( hsdata->event.teldata[itel].img[0].num_hot>0 &&
                          hsdata->event.teldata[itel].raw != NULL )
                     {
                        int ihp = hsdata->event.teldata[itel].img[0].hot_pixel[0];
                        printf("Calibration cross check for hottest pixel of telescope %d: "
                           "%4.2f versus %4.2f peak p.e. in pixel %d.\n",
                           hsdata->camera_set[itel].tel_id,
                           CALIB_SCALE*hsdata->event.teldata[itel].img[0].hot_amp[0],
                           calibrate_pixel_amplitude(hsdata,itel,ihp,0,0.), ihp);
                     }
                  }
            }

	    for (itel=0; itel<hsdata->run_header.ntel; itel++)
	       if ( hsdata->event.teldata[itel].known )
	          hesscam_ps_plot(ps_fname, hsdata, itel, -1, 0, 0.);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_CALIBEVENT:
         {
            int type = -1;
            rc = read_hess_calib_event(iobuf,&hsdata->event,-1,&type);
	    if ( verbose || rc != 0 )
               printf("read_hess_calib_event(), rc = %d, type=%d\n",rc,type);
	    for (itel=0; itel<hsdata->run_header.ntel; itel++)
	       if ( hsdata->event.teldata[itel].known )
	          hesscam_ps_plot(ps_fname, hsdata, itel, type, 0, 0.);
         }
            break;

         /* =================================================== */
         case IO_TYPE_HESS_MC_SHOWER:
            rc = read_hess_mc_shower(iobuf,&hsdata->mc_shower);
	    if ( verbose || rc != 0 )
               printf("read_hess_mc_shower(), rc = %d\n",rc);
            if ( showdata )
               print_hess_mc_shower(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_MC_EVENT:
            rc = read_hess_mc_event(iobuf,&hsdata->mc_event);
	    if ( verbose || rc != 0 )
               printf("read_hess_mc_event(), rc = %d\n",rc);
            if ( showdata )
               print_hess_mc_event(iobuf);
            wsum_all += pow(hsdata->mc_shower.energy,
               plidx-hsdata->mc_run_header.spectral_index);
            if ( fabs(hsdata->mc_event.xcore) > rmax_x )
               rmax_x = fabs(hsdata->mc_event.xcore);
            if ( fabs(hsdata->mc_event.ycore) > rmax_y )
               rmax_y = fabs(hsdata->mc_event.ycore);
            rs = line_point_distance(hsdata->mc_event.xcore, 
                     hsdata->mc_event.ycore, 0.,
                     cos(hsdata->mc_shower.altitude) *
                        cos(hsdata->mc_shower.azimuth),
                     cos(hsdata->mc_shower.altitude) *
                        sin(-hsdata->mc_shower.azimuth),
                     sin(hsdata->mc_shower.altitude),
                     0., 0., 0.);
            if ( rs > rmax_r )
               rmax_r = rs;
            nev ++;
            if ( ( reco_flag>3 || 
                   (reco_flag==3 && hsdata->mc_shower.primary_id == 0) ) && 
                 rc == 0 && !quiet )
            {
               if ( iprint_mc++ == 0 )
                  printf("#@! Lines starting with '@!' contain the following columns:\n"
                         "#@!  (1): Event number (Shower*100+Array)\n"
                         "#@!  (2): Simulated shower azimuth [deg.]\n"
                         "#@!  (3): Simulated shower altitude [deg.]\n"
                         "#@!  (4): Simulated X core [m]\n"
                         "#@!  (5): Simulated Y core [m]\n"
                         "#@!  (6): Simulated core offset in shower plane [m]\n"
                         "#@!  (7): Simulated primary ID\n"
                         "#@!  (8): Simulated energy [TeV]\n"
                         "#@!  (9): Simulated Xmax [g/cm^2]\n"
                         "#@! (10): Simulated Xmax(e) [g/cm^2]\n"
                         "#@! (11): Simulated Xmax(C) [g/cm^2]\n"
                         "#@! (12): Simulated Hmax [km]\n");
	       printf("@! %d %7.5f %7.5f  %5.2f %5.2f %5.2f   %d %7.4f   "
                  "%5.2f %5.2f %5.2f %5.3f\n",
                  hsdata->mc_event.event,
	          (180./M_PI)*hsdata->mc_shower.azimuth, 
                  (180./M_PI)*hsdata->mc_shower.altitude, 
                  hsdata->mc_event.xcore, 
                  hsdata->mc_event.ycore,
                  rs,
                  hsdata->mc_shower.primary_id,
	          hsdata->mc_shower.energy,
                  hsdata->mc_shower.xmax,
                  hsdata->mc_shower.emax,
                  hsdata->mc_shower.cmax,
                  hsdata->mc_shower.hmax*1e-3);
            }
            break;

         /* =================================================== */
         case IO_TYPE_MC_TELARRAY:
            if ( hsdata && hsdata->run_header.ntel > 0 )
            {
               rc = read_hess_mc_phot(iobuf,&hsdata->mc_event);
	       if ( verbose || rc != 0 )
                  printf("read_hess_mc_phot(), rc = %d\n",rc);
            }
            break;

         /* =================================================== */
         /* With extended output option activated, the particles
            arriving at ground level would be stored as seemingly
            stray photon bunch block. */
         case IO_TYPE_MC_PHOTONS:
            if ( showdata )
               print_tel_photons(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_MC_RUNH:
         case IO_TYPE_MC_EVTH:
         case IO_TYPE_MC_EVTE:
         case IO_TYPE_MC_RUNE:
            if ( showdata )
               print_tel_block(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_MC_PE_SUM:
            rc = read_hess_mc_pe_sum(iobuf,&hsdata->mc_event.mc_pesum);
	    if ( verbose || rc != 0 )
               printf("read_hess_mc_pe_sum(), rc = %d\n",rc);
            if ( showdata )
               print_hess_mc_pe_sum(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_TEL_MONI:
            // Telescope ID among others in the header
            tel_id = (item_header.ident & 0xff) | 
                     ((item_header.ident & 0x3f000000) >> 16); 
            if ( (itel = find_tel_idx(tel_id)) < 0 )
            {
               char msg[256];
               snprintf(msg,sizeof(msg)-1,
                  "Telescope monitor block for unknown telescope %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_tel_monitor(iobuf,&hsdata->tel_moni[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_tel_monitor(), rc = %d\n",rc);
            if ( showdata )
               print_hess_tel_monitor(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_LASCAL:
            tel_id = item_header.ident; // Telescope ID is in the header
            if ( (itel = find_tel_idx(tel_id)) < 0 )
            {
               char msg[256];
               snprintf(msg,sizeof(msg)-1,
                  "Laser/LED calibration for unknown telescope %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_laser_calib(iobuf,&hsdata->tel_lascal[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_laser_calib(), rc = %d\n",rc);
            if ( showdata )
               print_hess_laser_calib(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_RUNSTAT:
            rc = read_hess_run_stat(iobuf,&hsdata->run_stat);
	    if ( verbose || rc != 0 )
               printf("read_hess_run_stat(), rc = %d\n",rc);
            if ( showdata )
               print_hess_run_stat(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_MC_RUNSTAT:
            rc = read_hess_mc_run_stat(iobuf,&hsdata->mc_run_stat);
	    if ( verbose || rc != 0 )
               printf("read_hess_mc_run_stat(), rc = %d\n",rc);
            if ( showdata )
               print_hess_mc_run_stat(iobuf);
            break;

         /* (End-of-job or DST) histograms */
         case 100:
            {
               int nhist;
               HISTOGRAM *h_list[1] = { NULL };
               nhist = read_histograms(h_list,-1,iobuf);

               if ( !quiet )
                  printf("read_histograms() got %d histograms\n",nhist);
               if ( nhist == 1 && h_list[0] != NULL && h_list[0]->ident == 99 )
               {
                  double plidx_new = plidx;
                  if ( h_list[0]->extension != NULL )
                     if ( h_list[0]->extension->ddata != NULL )
                        plidx_new = h_list[0]->extension->ddata[0];
                  if ( fabs(plidx_new - plidx) > 0.001 )
                  {
                     fprintf(stderr,"\nWeighting of spectrum in MC DST data (for E^%5.3f)"
                        " and in analysis (for E^%5.3f) does not match!\n\n",
                        plidx_new, plidx);
                     /* plidx = plidx_new; */
                  }
               }
            }
            break;

         default:
            if ( !ignore )
               fprintf(stderr,"Ignoring unknown data block type %ld\n",item_header.type);
      }
    }

    if ( iobuf->input_file != NULL && iobuf->input_file != stdin )
      fileclose(iobuf->input_file);
    iobuf->input_file = NULL;
    reset_io_block(iobuf);

    if ( hsdata != NULL && !quiet )
      show_run_summary(hsdata,nev,ntrg,plidx,wsum_all,wsum_trg,rmax_x,rmax_y,rmax_r);
    else if ( nev > 0 )
      printf("%d of %d events triggered\n", ntrg, nev);

    if ( hsdata != NULL )
       hsdata->run_header.run = 0;
   }

   return 0;
}

/** @} */
