/* ============================================================================

Copyright (C) 2001, 2003, 2005, 2006, 2008, 2009, 2010, 2011, 2013, 2014, 2015  Konrad Bernloehr

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

/** @file read_hess.c
 *  @short A program reading simulated data, optionally analysing the data,
 *      and also optionally also writing summary ("DST") data.
 *
 *  This program started as a skeleton for reading H.E.S.S. data in eventio
 *  format (which is what the read_hess_nr program is now intended for).
 *  The read_hess program reads the whole range of hessio item types into a
 *  single tree of data structures but normally does nothing with the data.
 *
 *  It can be instructed to create nice camera images similar to
 *  those generated in sim_hessarray.
 *
 *  It can also be instructed to redo the image cleaning (with the simple
 *  10/5 tail-cut algorithm) and the shower reconstruction, writing
 *  ASCII output of the results.
 *
 *  In addition, it includes an interface for a full-scale analysis which
 *  can optionally be activated.
 *
 *  And finally, it can be instructed to extract DST-level data in order
 *  to reduce the amount of data by a large factor. This depends on the
 *  dst-level flag: 1) Remove all raw data (you cannot redo image cleaning)
 *  afterwards. 2) Remove also all MC data from non-triggered event (you
 *  should better stay with the spectral index used for DST extraction
 *  because you have to rely on its histograms for MC energy distribution).
 *  3) and 4) Keep only user-defined events (with or without raw data).
 *
@verbatim
read_hess: A program for viewing and analyzing sim_telarray (sim_hessarray) data.

Syntax: read_hess [ options ] [ - | input_fname ... ]
Options:
   -p ps_filename  (Write a PostScript file with camera images.)
   --plot-with-true-pe (If data available, include true p.e. plot in PS file.)
   --plot-with-sum-only (Show only sum image even if we have traces.)
   --plot-with-pixel-id (Show pixel ID number on top of pixel.)
   --plot-with-pixel-amp (Show pixel amplitude value on top of pixel.)
   --plot-with-pixel-pe (Show count of true Cherenkov p.e. on top of pixel.)
   --plot-without-reco (Do not show reconstructed image/shower parameters.)
   --plot-with-type-sum (Plot sum of pixel intensities over telescope types.)
   --plot-with-title text (User-defined title on top of page.)
   -r level        (Use 10/5 tail-cut image cleaning and redo reconstruction.)
                   level >= 1: show parameters from sim_hessarray.
                   level >= 2: redo shower reconstruction
                   level >= 3: redo image cleaning (and shower reconstruction
                               with new image parameters)
                   level >= 4: redo amplitude summation
                   level >= 5: PostScript file includes original and
                               new shower reconstruction.
   -v              (More verbose output)
   -q              (Much more quiet output)
   -s              (Show data explained)
   -S              (Show data explained, including raw data)
   --history (-h)  (Show contents of history data block)
   --clean-history (Drop previous history data blocks)
   -i              (Ignore unknown data block types)
   -u              (Call user-defined analysis function)
   --global-peak   (For image analysis use amplitude sums around global peak
                    in 'on-line' pulse shape analysis.)
   --local-peak    (For image analysis use amplitude sums around local peaks
                    in 'on-line' pulse shape analysis.)
   --powerlaw x    (Use this spectral index for events weights in output.)
                   (Default spectral index is -2.7)
   --only-run run1[,run2-run3[,...]] (Select runs being processed.)
   --not-run run1[,run2-run3[,...]]
   --only-telescope id1[,id2-i3[,...]] (Select telescopes being used.)
   --not-telescope id1[,id2-id3[,...]]
   --auto-trgmask  (Automatically load matching .trgmask.gz files.)
   --trgmask-path dir (Search the trgmask files in this path first.)
   --trg-required b *(Required trigger bits, e.g. 5=1|4 -> majo or asum)
   --type nt[,id1,id2,A,f,npix] (Set [requirements for] telescope type nt.)
   --focal-length f *(Set telescope imaging effective focal length [m].)
   --min-tel tmn   *(The minimum number of tel. images required in analysis.)
   --max-tel tmx   (The maximum number of tel. images required in analysis.)
   --min-trg-tel n (Minimum number of telescopes in system trigger.)
   --hard-stereo id1,id2,.. (Telescope of ID id1 etc. only use if stereo.)
   --min-amp npe   *(Minimum image amplitude for shower reconstruction.)
   --min-pix npix  *(Minimum number of pixels for shower reconstruction.)
   --max-events n  (Skip remaining data after so many triggered events.)
   --max-theta d   (Maximum angle between source and shower direction [deg].)
   --min-theta d   (Where cut angle is multiplicity dependent, use this
                    as the lower limit [deg].)
   --theta-scale f (Scale fixed and optimized theta cut by this factor.)
   --theta-E-scale t0,ts,min,max (Energy-dependent scaling beyond multiplicity.)
   --tail-cuts l,h[,n,f] *(Low and high level tail cuts to be applied in analysis.)
   --nb-radius r1[,r2[,r3]] *(Maximum distance of neighbour pixels [px diam.])
   --ext-radius r  *(Radius to extend preserved pixels beyond cleaning [px diam.])
   --dE2-cut c     (Cut parameter for dE2 cut.)
   --hess-standard-cuts (Apply HESS-style selection with standard cuts.)
   --hess-hard-cuts (Apply HESS-style selection with hard cuts.)
   --hess-loose-cuts (Apply HESS-style selection with loose cuts.)
   --hess-style-cuts (No shape parameter rescaling as HESS-style.)
   --shape-cuts wmn,wmx,lmn,lmx (Shape cut parameters: mscrw/l min/max).
   --dE-cut c      (Scale parameter for dE cut strictness, def=1.0).
   --hmax-cut c    (Scale parameter for hmax cut strictness, def=1.0).
   --min-img-angle a (Only use image pairs intersecting at angle > a deg, def=0).
   --min-disp d    *(Do not use round images with disp = (1-w/l) < d, def=0).
   --max-core-distance r *(Only use images from telescope not further from core).
   --impact-range r,x,y (Accept only events with reconstructed core in range).
   --true-impact-range r,x,y (Accept only events with true core in range).
                    Note that r is in shower plane but x,y ranges are on surface.
   --min-true-energy e (Completely skip events below given true energy.
   --clip-camera-radius r *(In image reconstruction clip camera at radius r deg.)
   --clip-camera-diameter d *(Same as before but with diameter d deg.)
   --clip-pixel-amplitude a *(Calibrated pixel ampl. does not exceed a mean p.e.)
   --only-high-gain (Use only high-gain channel and ignore low gain.)
   --only-low-gain (Use only low-gain channel and ignore high gain.)
   --max-events    (Stop after having processed this many events.)
   --pure-raw      (Discard any sub-items of TelescopeEvent which are not raw data.)
   --no-mc-data    (Discard MC shower and MC event data.)
   --broken-pixels-fraction (Add random broken/dead pixels on run-by-run basis.)
   --dead-time-fraction (Set telescopes randomly as dead from prior triggers.)
   --integration-scheme n *(Set the integration scheme for sample-mode data.
                   Use '--integration-scheme help' to show available schemes.)
   --integration-window w,o[,ps] *(Set integration window width and offset.)
                   For some integration schemes there is a pulse shaping option.
   --integration-treshold h[,l] *(Set significance thresholds for integration.)
   --integration-no-rescale *(Don't rescale pulse sum for integration with
                   windows narrower than a single-p.e. pulse.)
   --integration-rescale *(Rescale for single-p.e. fraction in window; default)
   --calib-scale f *(Rescale from mean p.e. to experiment units. Default: 0.92)
   --calib-error f (Random pixel relative calibration error. Default: 0.)
   --calibrate     (Store calibrated pixel intensities to DST file, if possible.)
   --only-calibrated (Like '--calibrate' but omit raw data from DST.)
   --diffuse-mode  (True shower position assumed as source position.)
   --random-seed n|auto (Initialize random number generator.)
   --off-axis-range a1,a2 (Only for diffuse mode, restricting range in deg.)
   --auto-lookup   (Automatically generate lookup table (gammas only).)
   --lookup-file name (Override automatic naming of lookup files.)
   --cleaning n    (Imaging cleaning setting: 0=no, 1-5=yes, see '--cleaning help')
   --zero-suppression n (Zero suppression scheme; 0: off, 3=auto)
   -z              (Equivalent to '--zero-suppression auto')
   --dst-level n   (Level of data reduction when writing DST-type output.)
                   Valid levels: 0, 1, 2, 3, 10, 11, 12, 13.
                   Raw data is stripped off at all levels except 0 and 10.
                   Level 0 has any sample mode data reduced to sums,
                   Level 1 includes all MC shower/event blocks,
                   level 2 only for triggered events,
                   level 3 has many config/calib blocks only once, not per run.
                   Levels 10-13 include only selected gamma-like events.
   --raw-level n   (Re-write original raw data or processed data, with possible
                   selection or reduction of other data according to level.)
                   Level 0 has all data written as available.
                   Level 1 has MC data only for triggered events.
                   Level 2 has no MC data (--no-mc-data).
                   Level 3 has only raw data for telescopes and nothing else (--pure-raw).
                   Level 4 also cleans past history data (--clean-history).
   --dst-file name (Name of output file for DST-type output.)
                   A DST file is needed for cleaning > 0 or DST level >= 0.
   --output-file   (Synonym to --dst-file)
   --histogram-file name (Name of histogram file.)
   -f fname        (Get list of input file names from fname.)

Parameters followed by a '*' can be type-specific if preceded by a
'--type' option. Their interpretation is thus position-dependent.

@endverbatim
 *
 *  @author Konrad Bernloehr
 *
 *  @date    @verbatim CVS $Date: 2018/09/18 15:10:17 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.143 $ @endverbatim
 */

/** @defgroup read_hess_c The read_hess (aka read_simtel, read_cta) program */
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
#include "reconstruct.h"
#include "user_analysis.h"
#include "warning.h"
#include "camera_image.h"
#include "basic_ntuple.h"
#include "io_trgmask.h"
#include "eventio_version.h"
#ifdef WITH_RANDFLAT
#include "rndm2.h"
#endif
#include <sys/time.h>
#include <strings.h>

struct basic_ntuple bnt;

/** The factor needed to transform from mean p.e. units to units of the single-p.e. peak:
    Depends on the collection efficiency, the asymmetry of the single p.e. amplitude 
    distribution and the electronic noise added to the signals. */
#define CALIB_SCALE 0.92

#include <signal.h>

void stop_signal_function (int isig);

static int interrupted;

static int dst_processing;

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

/* ----------------------- init_rand ---------------------- */
/* @short For a number of options (broken pixels and dead time)
 *        we need random numbers. These can be initialized by
 *        a predetermined seed or a random seed.
 */

static void init_rand(int is);

static void init_rand(int is)
{
   struct timeval tv;
//   struct timezone tz;
   unsigned long rnd_seed;

   if ( is == 0 ) /* random seed */
   {
//      gettimeofday(&tv,&tz);
      gettimeofday(&tv,NULL);
      rnd_seed = (tv.tv_sec%147483647)+2000*tv.tv_usec+getpid()*12345;
      // Even though we don't strictly need /dev/urandom, we still use it if available.
      FILE *frn = fopen("/dev/urandom","r");
      if ( frn != 0 )
      {
         unsigned int r;
         if ( fread(&r,sizeof(unsigned int),1,frn) > 0 )
            rnd_seed = (rnd_seed ^ r) % 2000000000;
         fclose(frn);
      }
   }
   else
      rnd_seed = (unsigned long) is;

#ifdef WITH_RANDFLAT
   Ranlux_setSeed((long)rnd_seed,5);
#else
   srand48((long)rnd_seed);
#endif
}

#ifndef WITH_RANDFLAT
static int g48_set;
static double g48_next;

double grand48(double mean, double sigma);

/** Like RandFlat() from rndm2.c but using the drand48 engine */

double grand48(double mean, double sigma)
{
  register double r;
  double v1,v2,fac,val;

  if ( g48_set ) 
  {
    g48_set = 0;
    return mean+sigma*g48_next;
  }

  do 
  {
    v1 = 2.0 * drand48() - 1.0;
    v2 = 2.0 * drand48() - 1.0;
    r = v1*v1 + v2*v2;
  } while ( r > 1.0 );

  fac = sqrt(-2.0*log(r)/r);
  val = v1*fac;
  g48_next = val;
  g48_set = 1;
  return mean+sigma*v2*fac;
}
#endif

/* ---------------------- mc_event_fill ------------------------- */
/** 
 *  @short Fill histogram(s) for DST writing which require all MC shower
 *         and event data and which cannot be filled from DST level >= 2 data.
 */

static void mc_event_fill (AllHessData *hsdata, double d_sp_idx)
{
   static HISTOGRAM *h12000 = NULL, *h11000 = NULL, *h11100 = NULL, *h22000 = NULL, *h22100 = NULL;
   static HISTOGRAM *h11010 = NULL, *h11110 = NULL, *h11020 = NULL, *h11120 = NULL;
   double ewt = pow(hsdata->mc_shower.energy, d_sp_idx);
   double cx = cos(hsdata->mc_shower.altitude) * cos(hsdata->mc_shower.azimuth);
   double cy = cos(hsdata->mc_shower.altitude) * sin(-hsdata->mc_shower.azimuth);
   double cz = sin(hsdata->mc_shower.altitude);

   double rs = line_point_distance(hsdata->mc_event.xcore, 
                     hsdata->mc_event.ycore, 0.,
                     cos(hsdata->mc_shower.altitude) *
                        cos(hsdata->mc_shower.azimuth),
                     cos(hsdata->mc_shower.altitude) *
                        sin(-hsdata->mc_shower.azimuth),
                     sin(hsdata->mc_shower.altitude),
                     0., 0., 0.);

   double sys_az = hsdata->run_header.direction[0], sys_alt = hsdata->run_header.direction[1];
   double csx = cos(sys_alt) * cos(sys_az);
   double csy = cos(sys_alt) * sin(-sys_az);
   double csz = sin(sys_alt);
   double coa = cx*csx + cy*csy + cz*csz;
   double oa = (coa >= -1. && coa <= 1. ) ? acos(coa)*(180./M_PI) : 999.;

   if ( h12000 == NULL )
   {
      double xylow[2], xyhigh[2];
      int nbins[2];

      fprintf(stderr,"DST production with event weight E^%5.3f\n",d_sp_idx);
      xylow[0]  = 0.;    xylow[1]  = -3.;
      xyhigh[0] = 4000.; xyhigh[1] = 4.;
      nbins[0]  = 400;   nbins[1]  = 140;
      histogram_hashing(10000);
      h12000 = book_histogram(12000, 
         "lg(true E) versus array core distance, all, before DST", 
         "D", 2, xylow, xyhigh, nbins);

      h22000 = book_1d_histogram(22000,  "lg(true E), all, before DSTs - no weights",
         "D", xylow[1], xyhigh[1], nbins[1]*5);
      h22100 = book_1d_histogram(22100,  "lg(true E), all, before DSTs - with weights",
         "D", xylow[1], xyhigh[1], nbins[1]*5);
   }
   if ( h11000 == NULL )
   {
      double xl = -3., xh = 4.;
      int nx = 280;

      double max_vc = 180.0;
      double max_r = 4000.;

      h11000 = book_histogram(11000,"lg(true E), all, before DST, weights",
         "D", 1, &xl, &xh, &nx);
      book_histogram(11001,"lg(true E), triggered, before DST, weights",
         "D", 1, &xl, &xh, &nx);
      h11100 = book_histogram(11100,"lg(true E), all, before DST, no weights",
         "D", 1, &xl, &xh, &nx);
      book_histogram(11101,"lg(true E), triggered, before DST, no weights",
         "D", 1, &xl, &xh, &nx);

      if ( hsdata->mc_run_header.core_pos_mode == 1 )
         max_r = hsdata->mc_run_header.core_range[1];
      if ( max_r < 1000. )
         max_r = 1000.;

      h11010 = book_1d_histogram(11010,  "impact parameter, all, before DSTs - no weights",
         "D", 0., max_r*1.25, 1250);
      h11110 = book_1d_histogram(11110,  "impact parameter, all, before DSTs - with weights",
         "D", 0., max_r*1.25, 1250);

      if ( hsdata->run_header.tracking_mode == 0 )
         max_vc =  hsdata->mc_run_header.viewcone[1]*(180./M_PI);
      if ( max_vc < 1.0 )
         max_vc = 1.0;
      h11020 = book_1d_histogram(11020,  "system off-axis angle, all, before DSTs - no weights",
         "D", 0., 1.5*max_vc, 750);
      h11120 = book_1d_histogram(11120,  "system off-axis angle, all, before DSTs - with weights",
         "D", 0., 1.5*max_vc, 750);
   }

   fill_histogram(h11000,log10(hsdata->mc_shower.energy),0.,ewt);
   fill_histogram(h11100,log10(hsdata->mc_shower.energy),0.,1.);
   fill_histogram(h12000,rs,log10(hsdata->mc_shower.energy),ewt);
   fill_histogram(h22000,log10(hsdata->mc_shower.energy),0.,1.);
   fill_histogram(h22100,log10(hsdata->mc_shower.energy),0.,ewt);
   fill_histogram(h11010,rs,0.,1.);
   fill_histogram(h11110,rs,0.,ewt);
   fill_histogram(h11020,oa,0.,1.);
   fill_histogram(h11120,oa,0.,ewt);
}

/* ---------------------- write_dst_histos --------------------------- */
/** 
 * @short Write histograms for DST book-keeping and clear them afterwards. 
 */

static int write_dst_histos (IO_BUFFER *iobuf2);
static int write_dst_histos (IO_BUFFER *iobuf2)
{
   int hlist[] = { 12000, 11000, 11001, 11100, 11101, 22000, 22100, 11010, 11110, 11020, 11120 };
   size_t i, nh = sizeof(hlist)/sizeof(hlist[0]);
   HISTOGRAM *h;
   for (i=0; i<nh; i++)
   {
      h = get_histogram_by_ident(hlist[i]);
      if ( h != NULL && h->entries > 0 )
      {
         write_histograms(&h,1,iobuf2);
         clear_histogram(h);
      }
   }
   return 0;
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
             "#@;       14: Maximum horizontal core distance in X [m]\n"
             "#@;       15: Maximum horizontal core distance in Y [m]\n"
             "#@;       16: Maximum core distance in shower plane [m]\n"
             "#@;       17: Supposed maximum core distance [m]\n");
      explained = 1;
   }
   printf("\n@; %d %d %d %d   %5.2f %5.2f %4.2f    %6.4f %6.4f %5.3f %5.3f   %g %g   %3.1f %3.1f %3.1f %3.1f\n", 
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
   printf("   --plot-with-true-pe (If data available, include true p.e. plot in PS file.)\n");
   printf("   --plot-with-sum-only (Show only sum image even if we have traces.)\n");
   printf("   --plot-with-pixel-id (Show pixel ID number on top of pixel.)\n");
   printf("   --plot-with-pixel-amp (Show pixel amplitude value on top of pixel.)\n");
   printf("   --plot-with-pixel-pe (Show count of true Cherenkov p.e. on top of pixel.)\n");
   printf("   --plot-without-reco (Do not show reconstructed image/shower parameters.)\n");
   printf("   --plot-with-type-sum (Plot sum of pixel intensities over telescope types.)\n");
   printf("   --plot-with-title text (User-defined title on top of page.)\n");
   printf("   -r level        (Use 10/5 tail-cut image cleaning and redo reconstruction.)\n");
   printf("                   level >= 1: show parameters from sim_hessarray.\n");
   printf("                   level >= 2: redo shower reconstruction\n");
   printf("                   level >= 3: redo image cleaning (and shower reconstruction\n");
   printf("                               with new image parameters)\n");
   printf("                   level >= 4: redo amplitude summation\n");
   printf("                   level >= 5: PostScript file includes original and\n");
   printf("                               new shower reconstruction.\n");
   printf("   -v              (More verbose output)\n");
   printf("   -q              (Much more quiet output)\n");
   printf("   -s              (Show data explained)\n");
   printf("   -S              (Show data explained, including raw data)\n");
   printf("   --history (-h)  (Show contents of history data block)\n");
   printf("   --clean-history (Drop previous history data blocks)\n");
   printf("   -i              (Ignore unknown data block types)\n");
   printf("   -u              (Call user-defined analysis function)\n");
   printf("   --global-peak   (For image analysis use amplitude sums around global peak\n");
   printf("                    in 'on-line' pulse shape analysis.)\n");
   printf("   --local-peak    (For image analysis use amplitude sums around local peaks\n");
   printf("                    in 'on-line' pulse shape analysis.)\n");
   printf("   --powerlaw x    (Use this spectral index for events weights in output.)\n");
   printf("                   (Default spectral index is -2.7)\n");
   printf("   --only-run run1[,run2-run3[,...]] (Select runs being processed.)\n");
   printf("   --not-run run1[,run2-run3[,...]]\n");
   printf("   --only-telescope id1[,id2-i3[,...]] (Select telescopes being used.)\n");
   printf("   --not-telescope id1[,id2-id3[,...]]\n");
   printf("   --auto-trgmask  (Automatically load matching .trgmask.gz files.)\n");
   printf("   --trgmask-path dir (Search the trgmask files in this path first.)\n");
   printf("   --trg-required b *(Required trigger bits, e.g. 5=1|4 -> majo or asum)\n");
   printf("   --type nt[,id1,id2,A,f,npix] (Set [requirements for] telescope type nt.)\n");
   printf("   --focal-length f *(Set telescope imaging effective focal length [m].)\n");
   printf("   --min-tel tmn   *(The minimum number of tel. images required in analysis.)\n");
   printf("   --max-tel tmx   (The maximum number of tel. images required in analysis.)\n");
   printf("   --min-trg-tel n (Minimum number of telescopes in system trigger.)\n");
   printf("   --hard-stereo id1,id2,.. (Telescope of ID id1 etc. only use if stereo.)\n");
   printf("   --min-amp npe   *(Minimum image amplitude for shower reconstruction.)\n");
   printf("   --min-pix npix  *(Minimum number of pixels for shower reconstruction.)\n");
   printf("   --max-events n  (Skip remaining data after so many triggered events.)\n");
   printf("   --max-theta d   (Maximum angle between source and shower direction [deg].)\n");
   printf("   --min-theta d   (Where cut angle is multiplicity dependent, use this\n");
   printf("                    as the lower limit [deg].)\n");
   printf("   --theta-scale f (Scale fixed and optimized theta cut by this factor.)\n");
   printf("   --theta-E-scale t0,ts,min,max (Energy-dependent scaling beyond multiplicity.)\n");
   printf("   --tail-cuts l,h[,n,f] *(Low and high level tail cuts to be applied in analysis.)\n");
   printf("   --nb-radius r1[,r2[,r3]] *(Maximum distance of neighbour pixels [px diam.])\n");
   printf("   --ext-radius r  *(Radius to extend preserved pixels beyond cleaning [px diam.])\n");
   printf("   --dE2-cut c     (Cut parameter for dE2 cut.)\n");
   printf("   --hess-standard-cuts (Apply HESS-style selection with standard cuts.)\n");
   printf("   --hess-hard-cuts (Apply HESS-style selection with hard cuts.)\n");
   printf("   --hess-loose-cuts (Apply HESS-style selection with loose cuts.)\n");
   printf("   --hess-style-cuts (No shape parameter rescaling as HESS-style.)\n");
   printf("   --shape-cuts wmn,wmx,lmn,lmx (Shape cut parameters: mscrw/l min/max).\n");
   printf("   --dE-cut c      (Scale parameter for dE cut strictness, def=1.0).\n");
   printf("   --hmax-cut c    (Scale parameter for hmax cut strictness, def=1.0).\n");
   printf("   --min-img-angle a (Only use image pairs intersecting at angle > a deg, def=0).\n");
   printf("   --min-disp d    *(Do not use round images with disp = (1-w/l) < d, def=0).\n");
   printf("   --max-core-distance r *(Only use images from telescope not further from core).\n");
   printf("   --impact-range r,x,y (Accept only events with reconstructed core in range).\n");
   printf("   --true-impact-range r,x,y (Accept only events with true core in range).\n");
   printf("                    Note that r is in shower plane but x,y ranges are on surface.\n");
   printf("   --min-true-energy e (Completely skip events below given true energy.\n");
   printf("   --clip-camera-radius r *(In image reconstruction clip camera at radius r deg.)\n");
   printf("   --clip-camera-diameter d *(Same as before but with diameter d deg.)\n");
   printf("   --clip-pixel-amplitude a *(Calibrated pixel ampl. does not exceed a mean p.e.)\n");
   printf("   --only-high-gain (Use only high-gain channel and ignore low gain.)\n");
   printf("   --only-low-gain (Use only low-gain channel and ignore high gain.)\n");
   printf("   --max-events    (Stop after having processed this many events.)\n");
   printf("   --pure-raw      (Discard any sub-items of TelescopeEvent which are not raw data.)\n");
   printf("   --no-mc-data    (Discard MC shower and MC event data.)\n");
   printf("   --broken-pixels-fraction (Add random broken/dead pixels on run-by-run basis.)\n");
   printf("   --dead-time-fraction (Set telescopes randomly as dead from prior triggers.)\n");
   printf("   --integration-scheme n *(Set the integration scheme for sample-mode data.\n");
   printf("                   Use '--integration-scheme help' to show available schemes.)\n");
   printf("   --integration-window w,o[,ps] *(Set integration window width and offset.)\n");
   printf("                   For some integration schemes there is a pulse shaping option.\n");
   printf("   --integration-treshold h[,l] *(Set significance thresholds for integration.)\n");
   printf("   --integration-no-rescale *(Don't rescale pulse sum for integration with\n");
   printf("                   windows narrower than a single-p.e. pulse.)\n");
   printf("   --integration-rescale *(Rescale for single-p.e. fraction in window; default)\n");
   printf("   --calib-scale f *(Rescale from mean p.e. to experiment units. Default: 0.92)\n");
   printf("   --calib-error f (Random pixel relative calibration error. Default: 0.)\n");
   printf("   --calibrate     (Store calibrated pixel intensities to DST file, if possible.)\n");
   printf("   --only-calibrated (Like '--calibrate' but omit raw data from DST.)\n");
   printf("   --diffuse-mode  (True shower position assumed as source position.)\n");
   printf("   --random-seed n|auto (Initialize random number generator.)\n");
   printf("   --off-axis-range a1,a2 (Only for diffuse mode, restricting range in deg.)\n");
   printf("   --auto-lookup   (Automatically generate lookup table (gammas only).)\n");
   printf("   --lookup-file name (Override automatic naming of lookup files.)\n");
   printf("   --cleaning n    (Imaging cleaning setting: 0=no, 1-5=yes, see '--cleaning help')\n");
   printf("   --zero-suppression n (Zero suppression scheme; 0: off, 3=auto)\n");
   printf("   -z              (Equivalent to '--zero-suppression auto')\n");
   printf("   --dst-level n   (Level of data reduction when writing DST-type output.)\n");
   printf("                   Valid levels: 0, 1, 2, 3, 10, 11, 12, 13.\n");
   printf("                   Raw data is stripped off at all levels except 0 and 10.\n");
   printf("                   Level 0 has any sample mode data reduced to sums,\n");
   printf("                   Level 1 includes all MC shower/event blocks,\n");
   printf("                   level 2 only for triggered events,\n");
   printf("                   level 3 has many config/calib blocks only once, not per run.\n");
   printf("                   Levels 10-13 include only selected gamma-like events.\n");
   printf("   --raw-level n   (Re-write original raw data or processed data, with possible\n");
   printf("                   selection or reduction of other data according to level.)\n");
   printf("                   Level 0 has all data written as available.\n");
   printf("                   Level 1 has MC data only for triggered events.\n");
   printf("                   Level 2 has no MC data (--no-mc-data).\n");
   printf("                   Level 3 has only raw data for telescopes and nothing else (--pure-raw).\n");
   printf("                   Level 4 also cleans past history data (--clean-history).\n");
   printf("   --dst-file name (Name of output file for DST-type output.)\n");
   printf("                   A DST file is needed for cleaning > 0 or DST level >= 0.\n");
   printf("   --output-file   (Synonym to --dst-file)\n");
   printf("   --histogram-file name (Name of histogram file.)\n");
#ifdef CHECK_MISSING_PE_LIST
   printf("   --check-missing-pe-list (Check if any p.e. lists are missing.)\n");
#endif
   printf("   -f fname        (Get list of input file names from fname.)\n");

   printf("\nParameters followed by a '*' can be type-specific if preceded by a\n"
          "'--type' option. Their interpretation is thus position-dependent.\n");

   exit(1);
}

struct next_file_struct
{
   char *fname;
   struct next_file_struct *next;
};
typedef struct next_file_struct NextFile;

NextFile *add_next_file (const char *fn, NextFile *nxt);

NextFile *add_next_file (const char *fn, NextFile *nxt)
{
   char *s;

   if ( fn == NULL || nxt == NULL )
      return NULL;

   /* If not already at the end of the chain, then find the end. */
   while ( nxt->next != NULL )
      nxt = nxt->next;

   /* In particular the first node may not have any name filled in. */
   if ( nxt->fname == NULL )
   {
      nxt->fname = strdup(fn);
      if ( (s = strchr(nxt->fname,'\n')) != NULL )
         *s = '\0';
      return nxt;
   }

   /* Later we have to allocate a new node before we can store the name. */
   nxt->next = (struct next_file_struct*) malloc(sizeof(NextFile));
   if ( nxt->next == NULL )
      return NULL;
   nxt = nxt->next;
   nxt->next = NULL;
   nxt->fname = strdup(fn);
   if ( (s = strchr(nxt->fname,'\n')) != NULL )
      *s = '\0';

   return nxt;   
}

struct range_list_struct
{
   long from, to;
   struct range_list_struct *next;
};
typedef struct range_list_struct RangeList;

RangeList *add_range(long f, long t, RangeList *rl);
int is_in_range(long n, RangeList *rl);

RangeList *add_range(long f, long t, RangeList *rl)
{
   if ( rl == NULL || (f == 0 && t == 0) )
      return NULL;

   /* If not already at the end of the chain, then find the end. */
   while ( rl->next != NULL )
      rl = rl->next;

   if ( rl->from == 0 && rl->to == 0 ) /* Typical for first node */
   {
      rl->from = f;
      rl->to = t;
      return rl;
   }

   /* Later we have to allocate a new node. */
   rl->next = (struct range_list_struct*) malloc(sizeof(RangeList));
   if ( rl->next == NULL )
      return NULL;
   rl = rl->next;
   rl->from = f;
   rl->to = t;
   rl->next = NULL;
   
   return rl;
}

int is_in_range(long n, RangeList *rl)
{
   if ( rl == NULL )
      return 0;

   while ( rl != NULL )
   {
      if ( n >= rl->from && n <= rl->to )
         return 1;
      rl = rl->next;
   }
   
   return 0;
}

/* -------------------- main program ---------------------- */
/** 
 *  @short Main program 
 *
 *  Main program function of read_hess.c program.
 */

int main (int argc, char **argv)
{
   IO_BUFFER *iobuf = NULL, *iobuf2 = NULL;
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
   double last_clip_amp = 0.;
   int ntel_trg = 0, min_tel_trg = 0, max_tel_trg = 0;
   int nev = 0, ntrg = 0, nrun = 0;
   char *program = argv[0];
   char base_program[1024];
   int showdata = 0, showhistory = 0;
   int clean_history = 0;
   int flag_amp_tm = 0;
   size_t num_only = 0, num_not = 0;
   FILE *ntuple_file = NULL;
   double oa_range[2] = { 0., 90. };
   int diffuse_mode = 0;
#ifdef CHECK_MISSING_PE_LIST
   int got_pe_list = 0, check_missing_pe_list = 0;
#endif

   int only_telescope[H_MAX_TEL], not_telescope[H_MAX_TEL], hard_stereo[H_MAX_TEL];
   int nhard_st = 0;
   int pure_raw = 0;
   int no_mc_data = 0;
   int show_true_pe = 0;

   static double min_amp_tel[H_MAX_TEL];
   static double tailcut_low_tel[H_MAX_TEL];
   static double tailcut_high_tel[H_MAX_TEL];
   static size_t min_pix_tel[H_MAX_TEL];
   static int lref_tel[H_MAX_TEL];
   static double minfrac_tel[H_MAX_TEL];

   double min_amp = 80., tailcut_low = 5., tailcut_high = 10.;
   double max_theta = 0., min_theta = 0., theta_scale = 1.0;
   double de_cut[] = { 1., 0., 1., 1. };
   double de2_cut[] = { 0.5, 0.0, 0.5, 0.5};
   double hmax_cut = 1.;
   size_t min_pix = 2, min_tel_img = 2, max_tel_img = 999;
   int user_ana = 0;
   double impact_range[3] = { 0., 0., 0. };
   double true_impact_range[3] = { 0., 0., 0. };
   double min_true_energy = 0.;
   size_t events = 0, max_events = 0;
   int dst_level = -1; /* <0: No data summary processing; 0: samples -> sums; ...; 3: Hillas parameters only; ...  */
   int cleaning = 0; /* 0: no cleaning, 1: clean + store sums, 2: clean + store samples, 3: clean + store both */
   int zero_suppression = -1;
   char *dst_fname = NULL;
   int mc_shower_stored = -1, mc_event_stored = -1, mc_pe_sum_stored = -1;
   int last_event_selected = 0, last_mc_pe_sum_id = 0;
   int user_flags = 0, auto_lookup = 0;
   int lref = 0;
   double minfrac = 0.0;
   int iarg;
   int with_clipping = 0;
   int do_calibrate = 0, only_calibrated = 0;
   double broken_pixels_fraction = 0., dead_time_fraction = 0.;
   double calerror = 0.;

   static AllHessData *hsdata;
   static NextFile first_file;
   NextFile *next_file = &first_file;
   static RangeList only_runs = { 0, 0, NULL }, not_runs = { -1, -1, NULL };
   int skip_run = 0;
   int show_type_sum = 0, ntypes = 0;
   int auto_trgmask = 0;
   struct trgmask_set *tms = NULL;
   struct trgmask_hash_set *ths = NULL;
   char trgmask_path[1024];
   int trg_req = 0, type_set = 0;
   trgmask_path[0] = '\0';

   /* In case we write a DST file, we want to record how we did it. */
   push_command_history(argc,argv);

   /* Show command line on output */
   if ( getenv("SHOWCOMMAND") != NULL )
   {
      for (iarg=0; iarg<argc; iarg++)
         printf("%s ",argv[iarg]);
      printf("\n");
   }
   {
      char *s = strrchr(program,'/');
      if ( s == NULL )
         strncpy(base_program, program, sizeof(base_program)-1);
      else
         strncpy(base_program, s+1, sizeof(base_program)-1);
   }

   /* Initialize them no matter if we need them */
   user_init_parameters();

   /* Catch INTerrupt and TERMinate signals to stop program */
   signal(SIGINT,stop_signal_function);
   signal(SIGTERM,stop_signal_function);
   interrupted = 0;

   /* First pass over command line for switching to quiet mode */
   for ( iarg=1; iarg < argc; iarg++ )
   {
      if ( strcmp(argv[iarg],"-q") == 0 || strcmp(argv[iarg],"--quiet") == 0 )
         quiet = 1;
      else if ( strcmp(argv[iarg],"-v") == 0 || strcmp(argv[iarg],"--verbose") == 0 )
         quiet = 0;
   }
   if ( ! quiet )
      show_hessio_max();
   /* Check assumed limits with the ones compiled into the library. */
   H_CHECK_MAX();

   if ( argc < 2 )
      input_fname = "iact.out";

   if ( (iobuf = allocate_io_buffer(5000000L)) == NULL )
   {
      Error("Cannot allocate I/O buffer");
      exit(1);
   }
   iobuf->max_length = 200000000L;

   for ( iarg=1; iarg<argc; iarg++ )
   {
      if ( strcmp(argv[iarg],"-u") == 0 )
      {
       	 user_ana = 1;
         break;
      }
   }

   if ( user_ana )
   {
      user_set_flags(user_flags);
      user_set_min_amp(min_amp);
      user_set_tail_cuts(tailcut_low,tailcut_high,lref,minfrac);
      user_set_min_pix(min_pix);
      user_set_reco_flag(reco_flag);
      user_set_tel_img(min_tel_img,max_tel_img);
      user_set_max_theta(max_theta,theta_scale,min_theta);
      user_set_de_cut(de_cut);
      user_set_de2_cut(de2_cut);
      user_set_hmax_cut(hmax_cut);
      user_set_auto_lookup(auto_lookup);
   }

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
      else if ( strcmp(argv[1],"--plot-with-true-pe") == 0 )
      {
	 show_true_pe = 1;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--plot-with-sum-only") == 0 )
      {
	 setenv("PLOT_WITH_SUM_ONLY","1",1);
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--plot-with-pixel-id") == 0 )
      {
	 setenv("PLOT_WITH_PIXEL_ID","1",1);
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--plot-with-pixel-amp") == 0 )
      {
	 setenv("PLOT_WITH_PIXEL_AMP","1",1);
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--plot-with-pixel-pe") == 0 )
      {
	 setenv("PLOT_WITH_PIXEL_PE","1",1);
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--plot-without-reco") == 0 )
      {
	 setenv("PLOT_WITHOUT_RECO","1",1);
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--plot-with-type-sum") == 0 )
      {
	 show_type_sum = 1;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--plot-with-title") == 0 && argc > 2 )
      {
	 setenv("PLOT_WITH_TITLE",argv[2],1);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"-r") == 0 && argc > 2 )
      {
	 reco_flag = atoi(argv[2]);
         user_set_reco_flag(reco_flag);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strncmp(argv[1],"-r",2) == 0 && strlen(argv[1]) > 2 )
      {
	 reco_flag = atoi(argv[1]+2);
         user_set_reco_flag(reco_flag);
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
      else if ( strcmp(argv[1],"-v") == 0 || strcmp(argv[1],"--verbose") == 0 )
      {
       	 verbose = 1;
         quiet = 0;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"-V") == 0  ||  strcmp(argv[1],"--version") == 0 )
      {
         printf("Program %s linked against eventIO/hessio library\n", base_program);
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

#if ( defined(__STDC__) && __STDC__ ) || defined(__cplusplus)
#ifdef SHOW
# undef SHOW
#endif
#ifdef _STR_
# undef _STR_
#endif
#define _XSTR_(s) _STR_(s) /**< Expand a macro first and then enclose in string */
#define _STR_(s) #s /**< Enclose in string without macro expansion. */

#define SHOW(s) if ( strcmp(#s,_XSTR_(s)) != 0 ) printf("   " #s " = " _XSTR_(s) "\n" )

         printf("Preprocessor definitions include:\n");
         SHOW(CTA);
         SHOW(CTA_KIND);
         SHOW(CTA_ULTRA);
         SHOW(CTA_ULTRA3);
         SHOW(CTA_ULTRA5);
         SHOW(CTA_MAX);
         SHOW(CTA_MAX_SC);
         SHOW(CTA_SC);
         SHOW(CTA_PROD1);
         SHOW(CTA_PROD2);
         SHOW(CTA_PROD2_SC);
         SHOW(CTA_PROD3);
         SHOW(CTA_PROD3_DEMO);
         SHOW(CTA_PROD3_SC);
         SHOW(CTA_PROD3_MERGE);
         SHOW(CTA_PROD4);
         SHOW(CTA_MINI);
         SHOW(CTA_MINI2);
         SHOW(CTA_MINI3);
         SHOW(CTA_MINI4);
         SHOW(HESS_PHASE_3);
         SHOW(HESS_PHASE_2);
         SHOW(HESS_PHASE_1);
         SHOW(LARGE_TELESCOPE);
         SHOW(NO_LARGE_TELESCOPE);
         SHOW(ULTRA_FINE_PIXELS);
         SHOW(MEGAPIX);
         SHOW(SMARTPIXEL);
         SHOW(NO_SMARTPIXEL);
         SHOW(MAXIMUM_TELESCOPES);
         SHOW(MAXIMUM_PIXELS);
         SHOW(MAXIMUM_SECTORS);
         SHOW(MAXIMUM_DRAWERS);
         SHOW(MAXIMUM_SLICES);
         SHOW(H_SAVE_MEMORY);
         SHOW(STORE_PIX_PHOTONS);
#undef SHOW
#else
         printf("Other preprocessor definitions not shown.\n");
#endif
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"-q") == 0 || strcmp(argv[1],"--quiet") == 0 )
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
      else if ( strcmp(argv[1],"--clean-history") == 0 ||
                strcmp(argv[1],"--clear-history") == 0 )
      {
       	 clean_history = 1;
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
#ifdef __cplusplus
         putenv(const_cast<char *>("PRINT_VERBOSE=1"));
#else
         putenv((char *)"PRINT_VERBOSE=1");
#endif
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"-u") == 0 )
      {
       	 user_ana = 1;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--auto-trgmask") == 0 )
      {
       	 auto_trgmask = 1; /* Find matching files with extra trigger mask patterns */
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--trgmask-path") == 0 && argc > 2 )
      {
       	 strncpy(trgmask_path,argv[2],sizeof(trgmask_path)-1);
	 argc-=2;
	 argv+=2;
	 continue;
      }
      else if ( strcmp(argv[1],"--trg-required") == 0 && argc > 2  )
      {
       	 trg_req = atoi(argv[2]);
         user_set_trg_req(trg_req);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--calibrate") == 0 )
      {
       	 do_calibrate = 1; /* If writing DST files, add calibrated pixel intensities */
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--only-calibrated") == 0 )
      {
       	 do_calibrate = only_calibrated = 1; /* Only calibrated pixel intensities, no raw data */
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--global-peak") == 0 )
      {
       	 flag_amp_tm = 1; /* Use amplitude sums around global peak if available. */
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--local-peak") == 0 )
      {
       	 flag_amp_tm = 2; /* Use amplitude sums around local peaks if available. */
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--dst-level") == 0 && argc > 2 )
      {
       	 dst_level = atoi(argv[2]);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--raw-level") == 0 && argc > 2 )
      {
       	 int lvl = atoi(argv[2]);
         dst_level = 100;
         if ( lvl >= 1 )
            dst_level = 101;
         if ( lvl >= 2 )
            no_mc_data = 1;
         if ( lvl >= 3 )
            pure_raw = 1;
         if ( lvl >= 4 )
            clean_history = 1;
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( (strcmp(argv[1],"--dst-file") == 0 || 
         strcmp(argv[1],"--output-file") == 0) && argc > 2 )
      {
       	 dst_fname = argv[2];
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--dst-process") == 0 )
      {
         /* In normal data, the configuration for every telescope
            is found in every run. In DST-level data, depending
            on the details, only the first set of configuration
            blocks may be included -- assuming that the configuration
            will not change from one run to another. For processing
            DST-level data we should thus not reset the configuration
            from one run to another. */
       	 dst_processing = 1;
	 argc -= 1;
	 argv += 1;
	 continue;
      }
      else if ( strcmp(argv[1],"--ntuple-file") == 0 && argc > 2 )
      {
       	 ntuple_file = fileopen(argv[2],"w");
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--histogram-file") == 0 && argc > 2 )
      {
       	 user_set_histogram_file(argv[2]);
	 argc -= 2;
	 argv += 2;
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
      else if ( strcmp(argv[1],"--type" ) == 0 && argc > 2 )
      {
         user_set_tel_type_param_by_str(argv[2]);
         type_set = 1;
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--min-amp") == 0 && argc > 2 )
      {  // Type-specific parameter (set immediately)
       	 min_amp = atof(argv[2]);
         user_set_min_amp(min_amp);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--min-pix") == 0 && argc > 2 )
      {  // Type-specific parameter (set immediately)
       	 min_pix = atoi(argv[2]);
         user_set_min_pix(min_pix);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--max-theta") == 0 && argc > 2 )
      {  // Global parameter (setting is delayed)
       	 max_theta = atof(argv[2]);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--min-theta") == 0 && argc > 2 )
      {  // Global parameter (setting is delayed)
       	 min_theta = atof(argv[2]);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--theta-scale") == 0 && argc > 2 )
      {  // Global parameter (setting is delayed)
       	 theta_scale = atoi(argv[2]);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--theta-E-scale") == 0 && argc > 2 )
      {  // Global parameter only (set immediately)
         double the[4] = { 1., 0., 0., 10. };
         int nt = sscanf(argv[2],"%lf,%lf,%lf,%lf", &the[0], &the[1], &the[2], &the[3]);
         if ( nt > 0 )
       	    user_set_theta_escale(the);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( (strcmp(argv[1],"--tail-cuts") == 0 ||
                 strcmp(argv[1],"--tailcuts") == 0) && argc > 2 )
      {  // Type-specific parameter (set immediately)
         int lr = 0;
         double mf = 0.0;
         int na = 0;
         if ( (na=sscanf(argv[2],"%lf,%lf,%d,%lf",&tailcut_low,&tailcut_high,&lr,&mf)) < 2 || na == 3 )
	 {
	    fprintf(stderr,"Syntax error in tail-cuts parameters.\n");
            exit(1);
	 }
         if ( na == 4 )
         {
            lref = lr;
            minfrac = mf; // FIXME: not used type-specific yet
         }
         user_set_tail_cuts(tailcut_low,tailcut_high,lr,mf);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( (strcmp(argv[1],"--nb-radius") == 0 ||
                 strcmp(argv[1],"--nb-distance") == 0) && argc > 2 )
      {  // Type-specific parameter (set immediately)
         double rnb[3] = { 0., 0., 0. };
         int na = 0;
         if ( (na=sscanf(argv[2],"%lf,%lf,%lf",&rnb[0],&rnb[1],&rnb[2])) < 1 )
	 {
	    fprintf(stderr,"Syntax error in neighbourhood radii.\n");
            exit(1);
	 }
         user_set_nb_radius(rnb);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--ext-radius") == 0 && argc > 2 )
      {  // Type-specific parameter (set immediately)
         double rxt = atof(argv[2]);
         user_set_nxt_radius(rxt);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--cleaning") == 0 && argc > 2 )
      {
         if ( strcmp(argv[2],"help") == 0 )
         {
            printf("Known data cleaning values:\n"
                   "   0: No cleaning, leave data as-is\n"
                   "   1: Sums for clean+extension region but no traces\n"
                   "   2: Sums for all pixels and traces only for clean+extension region\n"
                   "   3: Sums and traces only for clean+extension region\n"
                   "   4: Sums for clean+extension region, traces only for clean+nb region\n"
                   "   5: Sums for clean+extension region, traces only for clean region\n"
                   "Clean region includes pixels passing image cleaning, clean+nb also\n"
                   "their neighbours, clean+extension the neighbours within --ext-radius\n"
                   "rather than the normal neighbours.\n");
            exit(1);
         }
       	 cleaning = atoi(argv[2]);
         if ( cleaning > 0 )
         {
            if ( reco_flag < 3 )
               reco_flag = 3; /* Implied image cleaning needs reconstruction flag of at least 3. */
            if ( zero_suppression <= 0 )
               zero_suppression = 3; /* Implied auto zero suppression */
         }
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--zero-suppression") == 0 && argc > 2 )
      {
         if ( strcasecmp(argv[2],"auto") == 0 )
            zero_suppression = 3;
         else
         {
            zero_suppression = -1;
            sscanf(argv[2], "%d", &zero_suppression);
         }
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"-z") == 0 )
      {
         zero_suppression = 3;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( (strcmp(argv[1],"--dE-cut") == 0 ||
                 strcmp(argv[1],"--eres-cut") == 0 ||
                 strcmp(argv[1],"--de-cut") == 0) && argc > 2 )
      {  // Global parameter (setting is delayed)
      	 int i, nv = 0;
	 double pv[] = { 0., 0., 0., 0. };
         if ( (nv = sscanf(argv[2],"%lf,%lf,%lf,%lf",pv+0,pv+1,pv+2,pv+3)) < 1 )
	 {
	    fprintf(stderr,"Syntax error in dE cut parameter.\n");
            exit(1);
	 }
         for (i=0; i<nv; i++)
            de_cut[i] = pv[i];
	 if ( nv < 2 )
	    de_cut[1] = 0.;
	 if ( de_cut[1] == 0. )
	    de_cut[2] = de_cut[3] = de_cut[0];
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( (strcmp(argv[1],"--dE2-cut") == 0 ||
                 strcmp(argv[1],"--eres2-cut") == 0 ||
                 strcmp(argv[1],"--de2-cut") == 0) && argc > 2 )
      {  // Global parameter (setting is delayed)
      	 int i, nv = 0;
	 double pv[] = { 0., 0., 0., 0. };
         if ( (nv = sscanf(argv[2],"%lf,%lf,%lf,%lf",pv+0,pv+1,pv+2,pv+3)) < 1 )
	 {
	    fprintf(stderr,"Syntax error in dE2 cut parameter.\n");
            exit(1);
	 }
         for (i=0; i<nv; i++)
            de2_cut[i] = pv[i];
	 if ( nv < 2 )
	    de2_cut[1] = 0.;
	 if ( de2_cut[1] == 0. )
	    de2_cut[2] = de2_cut[3] = de2_cut[0];
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--hmax-cut") == 0 && argc > 2 )
      {  // Global parameter (setting is delayed)
         if ( sscanf(argv[2],"%lf",&hmax_cut) != 1 )
	 {
	    fprintf(stderr,"Syntax error in hmax cut parameter.\n");
            exit(1);
	 }
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( (strcmp(argv[1],"--shape-cuts") == 0 ||
                 strcmp(argv[1],"--shapecuts") == 0) && argc > 2 )
      {
         double wmn = -2.0, wmx = 0.6, lmn = -2.0, lmx = 1.2;
         size_t n = sscanf(argv[2],"%lf,%lf,%lf,%lf",&wmn,&wmx,&lmn,&lmx);
         if ( n < 4 )
         {
            fprintf(stderr,"Syntax error: You must specify all four parameters for '--shape-cuts'.\n");
            exit(1);
         }
         user_set_shape_cuts(wmn,wmx,lmn,lmx);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( (strcmp(argv[1],"--width-cut") == 0) && argc > 2 )
      {
      	 int nv = 0;
	 double pv[] = { 0.6, 0., 0.6, 0.6 };
         if ( (nv = sscanf(argv[2],"%lf,%lf,%lf,%lf",pv+0,pv+1,pv+2,pv+3)) < 1 )
	 {
	    fprintf(stderr,"Syntax error in mscrw_max cut parameter.\n");
            exit(1);
	 }
	 if ( nv < 2 )
	    pv[1] = 0.;
	 if ( pv[1] == 0. )
	    pv[2] = pv[3] = pv[0];
         user_set_width_max_cut(pv);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( (strcmp(argv[1],"--min-true-energy") == 0) && argc > 2 )
      {
         /* Mainly useful to pick nice images for plots */
         min_true_energy = atof(argv[2]);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( (strcmp(argv[1],"--min-img-angle") == 0) && argc > 2 )
      {
         /* One global parameter, needs to be passed to the image intersections */
         fprintf(stderr,"Not yet implemented: --min-img-angle\n");
         exit(1);
      }
      else if ( (strcmp(argv[1],"--min-disp") == 0) && argc > 2 )
      {
         /* Telescope-type specific, to be used for image selection */
         fprintf(stderr,"Not yet implemented: --min-disp\n");
         exit(1);
      }
      else if ( (strcmp(argv[1],"--length-cut") == 0) && argc > 2 )
      {
      	 int nv = 0;
	 double pv[] = { 1.2, 0., 0., 0. };
         if ( (nv = sscanf(argv[2],"%lf,%lf,%lf,%lf",pv+0,pv+1,pv+2,pv+3)) < 1 )
	 {
	    fprintf(stderr,"Syntax error in mscrl_max cut parameter.\n");
            exit(1);
	 }
	 if ( nv < 2 )
	    pv[1] = 0.;
	 if ( pv[1] == 0. )
	    pv[2] = pv[3] = pv[0];
         user_set_length_max_cut(pv);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--max-core-distance") == 0 && argc > 2 )
      {
         double rt = atoi(argv[2]);
         user_set_max_core_distance(rt);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--impact-range") == 0 && argc > 2 )
      {
         int nv = 0;
         double rr=0., rx=0., ry=0.;
         impact_range[0] = impact_range[1] = impact_range[2] = 0.;
         nv = sscanf(argv[2],"%lf,%lf,%lf",&rr,&rx,&ry);
         if ( nv > 0 )
            impact_range[0] = rr;
         if ( nv > 1 )
            impact_range[1] = rx;
         if ( nv > 2 )
            impact_range[2] = ry;
         user_set_impact_range(impact_range);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--true-impact-range") == 0 && argc > 2 )
      {
         int nv = 0;
         double rr=0., rx=0., ry=0.;
         true_impact_range[0] = true_impact_range[1] = true_impact_range[2] = 0.;
         nv = sscanf(argv[2],"%lf,%lf,%lf",&rr,&rx,&ry);
         if ( nv > 0 )
            true_impact_range[0] = rr;
         if ( nv > 1 )
            true_impact_range[1] = rx;
         if ( nv > 2 )
            true_impact_range[2] = ry;
         user_set_true_impact_range(true_impact_range);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--clip-camera-radius") == 0 && argc > 2 )
      {
       	 double dc = atof(argv[2]);
         user_set_clipping(dc);
         with_clipping = 1;
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--clip-camera-diameter") == 0 && argc > 2 )
      {
       	 double dc = atof(argv[2]) * 0.5;
         user_set_clipping(dc);
         with_clipping = 1;
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--clip-pixel-amplitude") == 0 && argc > 2 )
      {
       	 double cpa = atof(argv[2]);
         user_set_clipamp(cpa);
         last_clip_amp = cpa; /* Used type-independent for PS plots */
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--only-high-gain") == 0 )
      {
         select_calibration_channel(1); // Has global impact, for all telescope types.
         argc--;
         argv++;
         continue;
      }
      else if ( strcmp(argv[1],"--only-low-gain") == 0 )
      {
         select_calibration_channel(2); // Has global impact, for all telescope types.
         argc--;
         argv++;
         continue;
      }
      else if ( strcmp(argv[1],"--max-events") == 0 && argc > 2 )
      {
       	 max_events = atol(argv[2]);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--broken-pixels-fraction") == 0 && argc > 2 )
      {
         broken_pixels_fraction = atof(argv[2]);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--random-seed") == 0 && argc > 2 )
      {
         if ( strcmp(argv[2],"auto") == 0 )
            init_rand(0); // random seed will be generated
         else
            init_rand(atoi(argv[2]));
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      
      else if ( strcmp(argv[1],"--dead-time-fraction") == 0 && argc > 2 )
      {
         dead_time_fraction = atof(argv[2]);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--pure-raw") == 0 )
      {
       	 pure_raw = 1;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--no-mc-data") == 0 )
      {
       	 no_mc_data = 1;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( ( strcmp(argv[1],"--only-telescope") == 0 ||
                  strcmp(argv[1],"--only-telescopes") == 0 ) &&
                argc > 2 )
      {
       	 char word[20];
         int ipos=0;
         while ( getword(argv[2],&ipos,word,sizeof(word)-1,',','\n') > 0 &&
                 num_only < H_MAX_TEL )
         {
            int tel_idx = atoi(word), tel_idx2;
            char *sl;
            if ( tel_idx > 0 )
            {
               if ( (sl = strchr(word,'-')) != NULL )
               {
                  if ( ( tel_idx2 = atoi(sl+1) ) >= tel_idx )
                  {
                     int ii;
                     for (ii=tel_idx; ii<=tel_idx2 && num_only<H_MAX_TEL; ii++)
                        only_telescope[num_only++] = ii;
                     printf("Only telescopes in the range %d to %d\n",
                        tel_idx, tel_idx2);
                  }
                  else
                     fprintf(stderr,"Bad telescope range: %s\n",word);
               }
               else
               {
                  only_telescope[num_only++] = tel_idx;
                  printf("Only telescope %d\n",tel_idx);
               }
            }
         }
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( ( strcmp(argv[1],"--not-telescope") == 0 ||
                  strcmp(argv[1],"--not-telescopes") == 0 ) &&
                argc > 2 )
      {
       	 char word[20];
         int ipos=0;
         while ( getword(argv[2],&ipos,word,sizeof(word)-1,',','\n') > 0 &&
                 num_not < H_MAX_TEL )
         {
            int tel_idx = atoi(word), tel_idx2;
            char *sl;
            if ( tel_idx > 0 )
            {
               if ( (sl = strchr(word,'-')) != NULL )
               {
                  if ( ( tel_idx2 = atoi(sl+1) ) >= tel_idx )
                  {
                     int ii;
                     for (ii=tel_idx; ii<=tel_idx2 && num_not<H_MAX_TEL; ii++)
                        not_telescope[num_not++] = ii;
                     printf("Not telescopes in the range %d to %d\n",
                        tel_idx, tel_idx2);
                  }
                  else
                     fprintf(stderr,"Bad telescope range: %s\n",word);
               }
               else
               {
                  not_telescope[num_not++] = tel_idx;
                  printf("Not telescope %d\n",tel_idx);
               }
            }
         }
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--hard-stereo") == 0 && argc > 2 )
      {
       	 char word[20];
         int ipos=0;
         int j;
         while ( getword(argv[2],&ipos,word,sizeof(word)-1,',','\n') > 0 )
         {
            int dupl = 0;
            tel_id = atoi(word);
            for ( j=0; j<nhard_st; j++ )
               if ( tel_id == hard_stereo[j] )
                  dupl = 1;
            if ( dupl )
               continue;
            if ( nhard_st >= H_MAX_TEL )
            {
               fprintf(stderr,"Too many telescopes in hardware stereo list.\n");
               exit(1);
            }
            hard_stereo[nhard_st++] = tel_id;
         }
         printf("Hardware stereo required for a subset of %d telescopes.\n", nhard_st);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( ( strcmp(argv[1],"--only-run") == 0 ||
                  strcmp(argv[1],"--only-runs") == 0 ) &&
                argc > 2 )
      {
       	 char word[20];
         int ipos=0;
         while ( getword(argv[2],&ipos,word,sizeof(word)-1,',','\n') > 0 )
         {
            long run_idx = atol(word), run_idx2;
            char *sl;
            if ( run_idx > 0 )
            {
               run_idx2 = run_idx;
               if ( (sl = strchr(word,'-')) != NULL )
               {
                  run_idx2 = atol(sl+1);
               }
               add_range(run_idx, run_idx2, &only_runs);
            }
         }
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( ( strcmp(argv[1],"--not-run") == 0 ||
                  strcmp(argv[1],"--not-runs") == 0 ) &&
                argc > 2 )
      {
       	 char word[20];
         int ipos=0;
         while ( getword(argv[2],&ipos,word,sizeof(word)-1,',','\n') > 0 )
         {
            long run_idx = atol(word), run_idx2;
            char *sl;
            if ( run_idx > 0 )
            {
               run_idx2 = run_idx;
               if ( (sl = strchr(word,'-')) != NULL )
               {
                  run_idx2 = atol(sl+1);
               }
               add_range(run_idx, run_idx2, &not_runs);
            }
         }
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--focal-length") == 0 && argc > 2 )
      {
       	 double flen = atof(argv[2]);
         user_set_focal_length(flen);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--min-tel") == 0 && argc > 2 )
      {
       	 min_tel_img = atoi(argv[2]);
         user_set_tel_img(min_tel_img,max_tel_img);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--max-tel") == 0 && argc > 2 )
      {
       	 max_tel_img = atoi(argv[2]);
         user_set_tel_img(min_tel_img,max_tel_img);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--min-trg-tel") == 0 && argc > 2 )
      {
       	 min_tel_trg = atoi(argv[2]);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--max-trg-tel") == 0 && argc > 2 )
      {
       	 max_tel_trg = atoi(argv[2]);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--select") == 0 && argc > 2 )
      {
       	 int list[H_MAX_TEL];
         char *s = strchr(argv[2],'(');
         size_t ntel = 0, mintel = 0;
         char word[10];
         int jpos = 0;
         if ( s == NULL )
            continue;
         *s++ = '\0';
         mintel = atoi(argv[2]);
         while (ntel < H_MAX_TEL)
         {
            if ( getword(s,&jpos,word,9,',',')') <= 0 )
               break;
            if ( (list[ntel] = atoi(word)) > 0 )
               ntel++;
         }
         if ( ntel > 0 )
            user_set_tel_list(mintel,ntel,list);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--integration-scheme") == 0 && argc > 2 )
      {
         int scheme = 0;
         if ( strcmp(argv[2],"help") == 0 )
         {
            printf("Supported signal integration schemes:\n"
                   "   1 = simple = fixed (o: after start of readout)\n"
                   "   2 = peak = global (o: before global peak bin from significant pixels)\n"
                   "   3 = local (o: before pixel local peak bin)\n"
                   "   4 = nb = neighbour = neighbor (o: before peak bin in sum of neighbours)\n"
                   "   5 = nb+local (o: before peak bin in sum of neighbours + 3 * local pixel)\n"
                   "   6 = gradient (first fitting time gradient along major axis)\n"
                   "   7 = nb-fc (like no. 6 = 'nb' but after FlashCam-style pulse shaping)\n"
                   "   where the offset 'o' is given by the '--integration-window' option.\n"
                   "   The width parameter 'w' of the '--integration-window' option is used\n"
                   "   in the same way by all integration schemes.\n"
                   "   The high-gain (and low-gain) threshold(s) of the '--integration-threshold'\n"
                   "   option are only used with integration schemes no. 2, 3, and 6.\n"
                   "   In scheme no. 6 it defines pixels to be used for the gradient fit.\n"
                   "   In scheme no. 7 the h-g threshold gets used to re-evaluate pixel timing.\n"
                   "   Scheme no. 7 also uses a third parameter of '--integration-window'\n"
                   "   for defining which pulse shaping and differencing option to use\n"
                   "   (0: original pzpsa shaping, 1: equivalent to pzpsa with 4 ns differencing,\n" 
                   "   or 2: with 8 ns differencing)\n"
                  );
            exit(1);
         }
         if ( strcasecmp(argv[2],"simple") == 0 || strcasecmp(argv[2],"fixed") == 0 )
            scheme = 1;
         else if ( strcasecmp(argv[2],"peak") == 0 || strcasecmp(argv[2],"global") == 0 )
            scheme = 2;
         else if ( strcasecmp(argv[2],"local") == 0 )
            scheme = 3;
         else if ( strcasecmp(argv[2],"nb") == 0 || strcasecmp(argv[2],"neighbour") == 0 ||
               strcasecmp(argv[2],"neighbor") == 0 )
            scheme = 4;
         else if ( strcasecmp(argv[2],"nb+local") == 0 )
            scheme = 5;
         else if ( strcasecmp(argv[2],"gradient") == 0 )
            scheme = 6;
         else if ( strcasecmp(argv[2],"nb-fc") == 0 )
            scheme = 7;
         else
            scheme = atoi(argv[2]);
         if ( scheme < 1 || scheme > 7 )
         {
            fprintf(stderr, "Invalid signal integration scheme %s\n", argv[2]);
            exit(1);
         }
         /* Without sufficient reconstruction level the integration is not actually
            carried out; instead the plain sum over all samples gets used.
            Avoid such a surprising result and force reconstruction level 4. */
         if ( reco_flag < 3 )
         {
	    reco_flag = 4;
            user_set_reco_flag(reco_flag);
            fprintf(stderr,"Integration scheme forcing image reconstruction.\n");
         }
         user_set_integrator(scheme);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--integration-window") == 0 && argc > 2 )
      {
         int nsum = 0, noff = 0, ps_opt = 0;
         sscanf(argv[2], "%d,%d,%d", &nsum, &noff, &ps_opt);
         user_set_integ_window(nsum, noff, ps_opt);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--integration-threshold") == 0 && argc > 2 )
      {
         int ithg = 0, itlg = 0;
         sscanf(argv[2], "%d,%d", &ithg, &itlg);
         user_set_integ_threshold(ithg, itlg);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--integration-no-rescale") == 0 )
      {
         user_set_integ_no_rescale(1);
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--integration-rescale") == 0 )
      {
         user_set_integ_no_rescale(0);
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--calib-scale") == 0 && argc > 2 )
      {
         double s = 0.0;
         sscanf(argv[2],"%lf",&s);
         user_set_calib_scale(s);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--calib-error") == 0 && argc > 2 )
      {
         double s = 0.0;
         sscanf(argv[2],"%lf",&s);
         if ( s >= 0. && s < 1. )
            calerror = s;
         else
            fprintf(stderr,"Calibration scale error %f is out of range and ignored.\n", s);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--hess-standard-cuts") == 0 )
      {
         user_flags = 1;
         user_set_flags(user_flags);
         min_tel_img = 2;
         user_set_tel_img(min_tel_img,max_tel_img);
         min_pix = 2;
         user_set_min_pix(min_pix);
         min_amp = 80.;
         user_set_min_amp(min_amp);
         max_theta = min_theta = sqrt(0.0125);
         user_set_max_theta(max_theta,theta_scale,min_theta);
         tailcut_low = 5.;
         tailcut_high = 10.;
         user_set_tail_cuts(tailcut_low,tailcut_high,0,0.);
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--hess-hard-cuts") == 0 )
      {
         user_flags = 2;
         user_set_flags(user_flags);
         min_tel_img = 2;
         user_set_tel_img(min_tel_img,max_tel_img);
         min_pix = 2;
         user_set_min_pix(min_pix);
         min_amp = 200.;
         user_set_min_amp(min_amp);
         max_theta = min_theta = 0.1;
         user_set_max_theta(max_theta,theta_scale,min_theta);
         tailcut_low = 5.;
         tailcut_high = 10.;
         user_set_tail_cuts(tailcut_low,tailcut_high,0,0.);
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--hess-loose-cuts") == 0 )
      {
         user_flags = 3;
         user_set_flags(user_flags);
         min_tel_img = 2;
         user_set_tel_img(min_tel_img,max_tel_img);
         min_pix = 2;
         user_set_min_pix(min_pix);
         min_amp = 40.;
         user_set_min_amp(min_amp);
         max_theta = min_theta = 0.2;
         user_set_max_theta(max_theta,theta_scale,min_theta);
         tailcut_low = 5.;
         tailcut_high = 10.;
         user_set_tail_cuts(tailcut_low,tailcut_high,0,0.);
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--hess-style-cuts") == 0 )
      {
         user_flags = 4;
         user_set_flags(user_flags);
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--diffuse-mode") == 0 )
      {
         diffuse_mode = 1;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--off-axis-range") == 0 && argc > 2 )
      {
         int n = sscanf(argv[2],"%lf,%lf", &oa_range[0], &oa_range[1]);
         if ( n == 0 )
            oa_range[0] = oa_range[1] = 0.;
         else if ( n == 1 )
         {
            oa_range[1] = oa_range[0];
            oa_range[0] = 0.;
         }
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( strcmp(argv[1],"--auto-lookup") == 0 ||
                strcmp(argv[1],"--auto-lookups") == 0 )
      {
         auto_lookup = 1;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--lookup-file") == 0 && argc > 2 )
      {
         user_set_lookup_file(argv[2]);
	 argc -= 2;
	 argv += 2;
	 continue;
      }
      else if ( ( strcmp(argv[1],"-f") == 0 ||
                  strcmp(argv[1],"--files-from") == 0 ) && argc > 2 )
      {
         FILE *f = fileopen(argv[2],"r");
         char line[1024];
         if ( f == NULL )
         {
            fprintf(stderr, "Cannot read list of input files from '%s'.\n", argv[2]);
            exit(1);
         }
         while ( fgets(line,sizeof(line)-1,f) != NULL )
         {
            NextFile *nxt = add_next_file(line,next_file);
            if ( nxt == NULL )
               break;
            next_file = nxt;
         }
         fileclose(f);
	 argc -= 2;
	 argv += 2;
      }
#ifdef CHECK_MISSING_PE_LIST
      else if ( strcmp(argv[1],"--check-missing-pe-list") == 0 )
      {
         printf("Checking for event data not preceded by a p.e. list\n");
         check_missing_pe_list = 1;
         argc--;
         argv++;
      }
#endif
      else if ( strcmp(argv[1],"--help") == 0 )
      {
        printf("\n%s: A program for viewing and analyzing sim_telarray (sim_hessarray) data.\n\n", base_program);
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

   if ( with_clipping )
      printf("Some telescope will effectively get their field of view clipped.\n");
   if ( diffuse_mode )
   {
      if ( oa_range[1] <= oa_range[0] )
      {
         fprintf(stderr,"Invalid off-axis angle range for  diffuse mode.\n");
         exit(1);
      }
      if ( theta_scale != 1. )
      {
         fprintf(stderr,"Forcing theta_scale=1 in diffuse mode.\n");
         theta_scale = 1.0;
      }
   }
    
   set_reco_verbosity(verbose-quiet);
   user_set_verbosity(verbose-quiet);
   /* If really verbose output wanted, we don't require '--history' to show also the history. */
   if ( verbose && !quiet )
      showhistory = 1;

   if ( dst_fname == NULL )
   {
      if ( dst_level >= 0 )
         fprintf(stderr,"\nDST level reset due to missing DST/output file name.\n\n");
      if ( cleaning > 0 )
         fprintf(stderr,"\nCleaning option reset due to missing DST/output file name.\n\n");
      dst_level = -1;
      cleaning = 0;
   }
   else if ( (dst_level >= 0 || cleaning > 0) && iobuf->output_file == NULL )
   {
      iobuf->output_file = fileopen(dst_fname,"w");
      if ( iobuf->output_file == NULL )
      {
         perror(dst_fname);
         exit(1);
      }
      else
         printf("DST output file opened: %s\n", dst_fname);
      if ( (dst_level%10) >= 2 && dst_level < 100 )
      {
         double xl = -0.5, xh = 6.5;
         int nx = 7;
         HISTOGRAM *h;
         /* For DST level 2 we need a second buffer but writing */
         /* into the same file/stream as the first buffer. */
         if ( (iobuf2 = allocate_io_buffer(100000L)) == NULL )
         {
            Error("Cannot allocate 2nd I/O buffer");
            exit(1);
         }
         iobuf2->max_length = 200000000L;
         iobuf2->output_file = iobuf->output_file;
         if ( (h = get_histogram_by_ident(99)) != NULL )
            free_histogram(h);
         h = book_histogram(99,"Spectral index in DST production",
               "D", 1, &xl, &xh, &nx);
         fill_histogram(h,0.,0.,plidx);
         fill_histogram(h,1.,0.,plidx+2.0); /* this is just a guess */
         fill_histogram(h,2.,0.,min_amp);
         fill_histogram(h,3.,0.,tailcut_low);
         fill_histogram(h,4.,0.,tailcut_high);
         fill_histogram(h,5.,0.,min_pix);
         fill_histogram(h,6.,0.,reco_flag);
         write_histograms(&h,1,iobuf2);
      }
      /* Save the command line history */
      write_history(0,iobuf);
   }

   if ( user_ana )
   {
      user_set_telescope_type(0);
      user_set_flags(user_flags);
      // user_set_min_amp(min_amp);
      // user_set_tail_cuts(tailcut_low,tailcut_high);
      // user_set_min_pix(min_pix);
      user_set_reco_flag(reco_flag);
      // user_set_tel_img(min_tel_img,max_tel_img);
      user_set_max_theta(max_theta,theta_scale,min_theta);
      user_set_de_cut(de_cut);
      user_set_de2_cut(de2_cut);
      user_set_hmax_cut(hmax_cut);
      user_set_auto_lookup(auto_lookup);
      user_set_diffuse_mode(diffuse_mode,oa_range);
   }
   
   /* Now go over rest of the command line */
   while ( argc > 1 )
   {
      if ( argv[1][0] == '-' && argv[1][1] != '\0' )
         syntax(program);
      else
      {
         NextFile *nxt = add_next_file(argv[1],next_file);
         if ( nxt == NULL )
            break;
         next_file = nxt;
	 argc--;
	 argv++;
      }
   }
   if ( first_file.fname == NULL && input_fname != NULL )
      first_file.fname = strdup(input_fname);

   /* ================= Big loop over all input files =============== */


   for ( next_file = &first_file; 
         next_file != NULL && next_file->fname != NULL; 
         next_file = next_file->next )
   {
    if ( interrupted )
      break;
    input_fname = next_file->fname;
    if ( strcmp(input_fname ,"-") == 0 )
    {
      iobuf->input_file = stdin;
      if ( auto_trgmask )
         fprintf(stderr,"Cannot auto-load trgmask file matching standard input.\n");
    }
    else 
    {
      /* Before opening the actual input file we may want to get hold of a companion 
         file with trigger type bit patterns which may live in one of several places,
         a special directory specified with the --trgmask-path option, the same
         directory as the data, a companion directory for Log files, or from the
         current working directory or any input directories specified for fileopen().
      */
      if ( auto_trgmask ) /* Fix-up for bug in storing storing multi-trigger bit pattern */
      {
         if ( strstr(input_fname,".trgmask.gz") == NULL ) 
         /* Expected to be a regular input file and not an explicit trigger bit pattern file */
         {
            char fname1[10240], fname2[10240];
            char *bsn = fname1;
            char *s;
            strncpy(fname1,input_fname,sizeof(fname1)-20);
            if ( (s = strstr(fname1,".sim")) != NULL )
               *s = '\0';
            strcat(fname1,".trgmask.gz");
            if ( strrchr(fname1,'/') != NULL )
               bsn = strrchr(fname1,'/')+1;
            /* First possibility: in a special directory for these files */
            if ( trgmask_path[0] != '\0' && strlen(trgmask_path)+strlen(bsn)+2 < sizeof(fname2) )
            {
               strcpy(fname2,trgmask_path);
               strcat(fname2,"\n");
               strcat(fname2,bsn);
               iobuf->input_file = fileopen(fname2,READ_BINARY);
               if ( verbose && iobuf->input_file == NULL )
                  printf("No such file: %s\n", fname2);
            }
            /* Next possibility: in same directory as the sim_telarray data file */
            if ( iobuf->input_file == NULL )
             if ( (iobuf->input_file = fileopen(fname1,READ_BINARY)) == NULL )
            {
               if ( verbose )
                  printf("No such file: %s\n", fname1);
               if ( (bsn = strrchr(fname1,'/')) == NULL )
               {
                  strcpy(fname2,"../Log/");
                  bsn = fname1;
               }
               else
               {
                  *bsn = '\0';
                  bsn++;
                  s = strrchr(fname1,'/');
                  if ( s != NULL )
                  {
                     *s = '\0';
                     strcpy(fname2,fname1);
                     strcat(fname2,"/Log/");
                  }
                  else
                  {
                     strcpy(fname2,fname1);
                     strcat(fname2,"/Log/");
                  }
               }
               strcat(fname2,bsn);
               /* Second possibility: in the log file directory */
               if ( (iobuf->input_file = fileopen(fname2,READ_BINARY)) == NULL )
               {
                  if ( verbose )
                     printf("No such file: %s\n", fname2);
                  /* Third possibility: in the current directory or standard search paths */
                  iobuf->input_file = fileopen(bsn,READ_BINARY);
               }
            }
            if ( iobuf->input_file != NULL )
            {
               printf("Using extra trigger type bit patterns from %s\n",bsn);
	       if ( find_io_block(iobuf,&item_header) == 0 )
               {
                  printf("Found I/O block of type %ld\n",item_header.type); // Expecting: 2090
	          if ( read_io_block(iobuf,&item_header) == 0 )
                  {
                     if ( tms == NULL )
                        tms = (struct trgmask_set *) calloc(1,sizeof(struct trgmask_set));
                     if ( ths == NULL )
                        ths = (struct trgmask_hash_set *) calloc(1,sizeof(struct trgmask_hash_set));
                     rc = read_trgmask(iobuf, tms);
	             if ( verbose || rc != 0 )
                        printf("read_trgmask(), rc = %d\n",rc);
                     if ( showdata )
                        print_trgmask(iobuf);
                     if ( rc == 0 )
                        trgmask_fill_hashed(tms,ths);
                  }
               }
               fileclose(iobuf->input_file);
               iobuf->input_file = NULL;
            }
            else
               fprintf(stderr,"No such trgmask file: %s\n",bsn);
         }
      }

      /* Now open the actual input file */
      if ( (iobuf->input_file = fileopen(input_fname,READ_BINARY)) == NULL )
      {
         perror(input_fname);
         Error("Cannot open input file.");
         break;
      }
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
	 printf("Found monitor block for telescope ID %d\n",tmp_moni.tel_id);
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

      if ( hsdata == NULL && !skip_run &&
           item_header.type > IO_TYPE_HESS_RUNHEADER &&
           item_header.type < IO_TYPE_HESS_RUNHEADER + 200 &&
           item_header.type != IO_TYPE_HESS_XTRGMASK )
      {
         if ( item_header.type == IO_TYPE_HESS_MCRUNHEADER )
         {
            MCRunHeader tmp_mc_run_header;
            rc = read_hess_mcrunheader(iobuf,&tmp_mc_run_header);
	    if ( verbose || rc != 0 )
               printf("read_hess_mcrunheader(), rc = %d\n",rc);
            if ( showdata )
               print_hess_mcrunheader(iobuf);
	       
            if ( user_ana )
	       user_set_spectrum(plidx - tmp_mc_run_header.spectral_index);
            fprintf(stderr,"Stray MC run header block before run header (possible in merged files)\n"
                 "is basically ignored. There should be more such blocks later-on.\n");
         }
         else
         {
            fprintf(stderr,"Trying to read data block of type %ld before run header.\n", item_header.type);
            fprintf(stderr,"Skipping this data block.\n");
         }
         continue;
      }
      
      if ( clean_history && (int) item_header.type == 70 )
         continue;
      if ( cleaning > 0 || (dst_level >= 0 && dst_level <= 3) ||
           (dst_level >= 10 && dst_level <= 13) || dst_level >= 100 ) /* DST-type data extraction */
      {
         if ( dst_level == 0 || dst_level == 1 || dst_level == 100 ) /* Everything copied verbatim except event data */
         {
            if ( (int) item_header.type != IO_TYPE_HESS_EVENT &&
                 (int) item_header.type != IO_TYPE_MC_TELARRAY )
            {
               if ( no_mc_data && 
                    ( (int) item_header.type == IO_TYPE_HESS_MC_EVENT ||
                      (int) item_header.type == IO_TYPE_HESS_MC_SHOWER ||
                      (int) item_header.type != IO_TYPE_HESS_MC_PE_SUM ) )
                  ; /* Don't write these block types if asking for no mc data */
               else if ( write_io_block(iobuf) )
               {
                  fprintf(stderr,"Writing DST output failed on item type %lu, id %ld.\n",
                     item_header.type, item_header.ident);
                  exit(1);
               }
            }
            else if ( dst_level == 100 && (int) item_header.type == IO_TYPE_MC_TELARRAY && ! no_mc_data )
            {
               if ( write_io_block(iobuf) )
               {
                  fprintf(stderr,"Writing DST output failed on item type %lu, id %ld.\n",
                     item_header.type, item_header.ident);
                  exit(1);
               }
            }
         }
         else if ( (dst_level%10) >= 2 || dst_level == 11 || dst_level == 10 || dst_level >= 101 ) /* MC data only for triggered/selected events */
         {
            if ( (int) item_header.type != IO_TYPE_HESS_EVENT &&
                 (int) item_header.type != IO_TYPE_HESS_MC_EVENT &&
                 (int) item_header.type != IO_TYPE_HESS_MC_SHOWER &&
                 (int) item_header.type != IO_TYPE_HESS_MC_PE_SUM &&
                 (int) item_header.type != IO_TYPE_MC_TELARRAY )
            {
               if ( (dst_level%10) >= 3 && nrun > 1 && 
                    ( (int) item_header.type == IO_TYPE_HESS_CAMSETTINGS ||
                      (int) item_header.type == IO_TYPE_HESS_CAMORGAN ||
                      (int) item_header.type == IO_TYPE_HESS_PIXELSET ||
                      (int) item_header.type == IO_TYPE_HESS_TEL_MONI || 
                      (int) item_header.type == IO_TYPE_HESS_LASCAL) )
                  ; /* Copy these blocks for levels 3,13 to DSTs only for the first run. */
               else if ( write_io_block(iobuf) )
               {
                  fprintf(stderr,"Writing DST output failed on item type %lu, id %ld.\n",
                     item_header.type, item_header.ident);
                  exit(1);
               }
            }
            else if ( (int) item_header.type == IO_TYPE_HESS_MC_EVENT )
            {
               last_event_selected = 0;
               mc_event_stored = 0; /* Delay until there is a trigger */
            }
            else if ( (int) item_header.type == IO_TYPE_HESS_MC_SHOWER )
            {
               mc_shower_stored = 0; /* Delay until there is a trigger */
               mc_event_stored = -1; /* Ignore until next event gets read */
            }
            else if ( (int) item_header.type == IO_TYPE_HESS_MC_PE_SUM )
            {
               mc_pe_sum_stored = 0; /* Delay until there is a trigger */
               last_mc_pe_sum_id = item_header.ident;
            }
         }
         else if ( cleaning > 0 && dst_level < 0 )
         {
            if ( (int) item_header.type != IO_TYPE_HESS_EVENT )
               if ( write_io_block(iobuf) )
               {
                  fprintf(stderr,"Writing DST output failed on item type %lu, id %ld.\n",
                     item_header.type, item_header.ident);
                  exit(1);
               }
         }
         else
         {
            fprintf(stderr,"Unsupported DST level %d\n",dst_level);
            exit(1);
         }
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
                  if ( hsdata->event.teldata[itel].pixcal != NULL )
		  {
                     free(hsdata->event.teldata[itel].pixcal);
		     hsdata->event.teldata[itel].pixcal = NULL;
		  }
               }
               /* Free main structure */
               if ( !dst_processing )
               {
                  free(hsdata);
                  hsdata = NULL;
               }
            }
            
            nev = ntrg = 0;
            wsum_all = wsum_trg = 0.;

            nrun++;

            if ( hsdata == NULL )
               hsdata = (AllHessData *) calloc(1,sizeof(AllHessData));

            fflush(stdout);
            if ( (rc = read_hess_runheader(iobuf,&hsdata->run_header)) < 0 )
            {
               Warning("Reading run header failed.");
               exit(1);
            }
            if ( !quiet )
               printf("Reading simulated data for %d telescope(s)\n",hsdata->run_header.ntel);
	    if ( verbose || rc != 0 )
               printf("read_hess_runheader(), rc = %d (run %d)\n",rc,hsdata->run_header.run);
            fprintf(stderr,"\nStarting run %d\n",hsdata->run_header.run);
            if ( showdata )
               print_hess_runheader(iobuf);
            if ( user_ana )
               do_user_ana(hsdata,item_header.type,0);
            
            /* Initialize MC data structures in case file contains no MC data */
            
            hsdata->mc_shower.shower_num = -1;
            hsdata->mc_shower.primary_id = 0;
            hsdata->mc_shower.energy = 0.;
            hsdata->mc_event.event = hsdata->mc_event.shower_num = -1;
            hsdata->mc_event.xcore = hsdata->mc_event.ycore = 0.;
            hsdata->mc_event.aweight = 0.;

            /* Allocate dynamic sub-structures and set up telescope ID in all sub-structures */

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
                  Warning("Not enough memory for AdcData");
                  exit(1);
               }
               hsdata->event.teldata[itel].raw->tel_id = tel_id;

               if ( (hsdata->event.teldata[itel].pixtm =
                     (PixelTiming *) calloc(1,sizeof(PixelTiming))) == NULL )
               {
                  Warning("Not enough memory for PixelTiming");
                  exit(1);
               }
               hsdata->event.teldata[itel].pixtm->tel_id = tel_id;

               if ( do_calibrate && dst_level >= 0 ) /* Only when needed */
               {
                  if ( (hsdata->event.teldata[itel].pixcal = 
                         (PixelCalibrated *) calloc(1,sizeof(PixelCalibrated))) == NULL )
                  {
                     Warning("Not enough memory for PixelCalibrated");
                     exit(1);
                  }
                  hsdata->event.teldata[itel].pixcal->tel_id = tel_id;
               }

               if ( (hsdata->event.teldata[itel].img = 
                      (ImgData *) calloc(2,sizeof(ImgData))) == NULL )
               {
                  Warning("Not enough memory for ImgData");
                  exit(1);
               }
               hsdata->event.teldata[itel].max_image_sets = 2;
               hsdata->event.teldata[itel].img[0].tel_id = tel_id;
               hsdata->event.teldata[itel].img[1].tel_id = tel_id;

               hsdata->tel_moni[itel].tel_id = tel_id;
               hsdata->tel_lascal[itel].tel_id = tel_id;
            }

            skip_run = 0;

            if ( only_runs.from != 0 || only_runs.to != 0 )
            {
               if ( !is_in_range(item_header.ident, &only_runs) )
               {
                  skip_run = 1;
                  printf("Ignoring data of run %ld\n", item_header.ident);
                  if ( nrun > 0 )
                     continue;
               }
            }
            if ( is_in_range(item_header.ident, &not_runs) )
            {
               skip_run = 1;
               printf("Ignoring data of run %ld\n", item_header.ident);
               if ( nrun > 0 )
                  continue;
            }
            break;

         /* =================================================== */
         case IO_TYPE_HESS_MCRUNHEADER:
            if ( skip_run )
               continue;
            rc = read_hess_mcrunheader(iobuf,&hsdata->mc_run_header);
	    if ( verbose || rc != 0 )
               printf("read_hess_mcrunheader(), rc = %d\n",rc);
            if ( showdata )
               print_hess_mcrunheader(iobuf);
	       
            if ( user_ana )
	    {
	       user_set_spectrum(plidx - hsdata->mc_run_header.spectral_index);
               do_user_ana(hsdata,item_header.type,0);
	    }
            ntypes = 0;
            break;

         /* =================================================== */
	 case IO_TYPE_MC_INPUTCFG:
            if ( skip_run )
               continue;
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
                  "Camera settings for unknown telescope ID %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_camsettings(iobuf,&hsdata->camera_set[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_camsettings(), rc = %d (tel. ID=%d, itel=%d)\n",rc,tel_id,itel);
            if ( showdata )
               print_hess_camsettings(iobuf);
            if ( reco_flag > 0 )
            {
               int ttype = which_telescope_type(&hsdata->camera_set[itel]);
               struct user_parameters* up = user_get_parameters(ttype);
               min_amp_tel[itel] = up->d.min_amp;
               tailcut_low_tel[itel] = up->d.tailcut_low;
               tailcut_high_tel[itel] = up->d.tailcut_high;
               min_pix_tel[itel] = up->i.min_pix;
               lref_tel[itel] = up->i.lref;
               minfrac_tel[itel] = up->d.minfrac;
               set_disabled_pixels(hsdata,itel,broken_pixels_fraction);
            }
            if ( num_only > 0 )
            {
               size_t jtel, keep_known = 0;
               for (jtel=0; jtel<num_only; jtel++)
               {
                  if ( hsdata->run_header.tel_id[itel] == only_telescope[jtel] )
                     keep_known = 1;
               }
               if ( !keep_known )
                  hsdata->camera_set[itel].num_mirrors = -1;
            }
            if ( num_not > 0 )
            {
               size_t jtel;
               for (jtel=0; jtel<num_not; jtel++)
               {
                  if ( hsdata->run_header.tel_id[itel] == not_telescope[jtel] )
                     hsdata->camera_set[itel].num_mirrors = -2;
               }
            }
            if ( user_ana )
               do_user_ana(hsdata,item_header.type,0);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_CAMORGAN:
            tel_id = item_header.ident; // Telescope ID is in the header
            if ( (itel = find_tel_idx(tel_id)) < 0 )
            {
               char msg[256];
               snprintf(msg,sizeof(msg)-1,
                  "Camera organisation for unknown telescope ID %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_camorgan(iobuf,&hsdata->camera_org[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_camorgan(), rc = %d (tel. ID=%d, itel=%d)\n",rc,tel_id,itel);
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
                  "Pixel settings for unknown telescope ID %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_pixelset(iobuf,&hsdata->pixel_set[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_pixelset(), rc = %d (tel. ID=%d, itel=%d)\n",rc,tel_id,itel);
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
                  "Pixel disable block for unknown telescope ID %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_pixeldis(iobuf,&hsdata->pixel_disabled[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_pixeldis(), rc = %d (tel. ID=%d, itel=%d)\n",rc,tel_id,itel);
            /* No print function available */
            if ( showdata )
               print_hess_pixeldis(iobuf);
            set_disabled_pixels(hsdata,itel,broken_pixels_fraction);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_CAMSOFTSET:
            tel_id = item_header.ident; // Telescope ID is in the header
            if ( (itel = find_tel_idx(tel_id)) < 0 )
            {
               char msg[256];
               snprintf(msg,sizeof(msg)-1,
                  "Camera software settings for unknown telescope ID %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_camsoftset(iobuf,&hsdata->cam_soft_set[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_camsoftset(), rc = %d (tel. ID=%d, itel=%d)\n",rc,tel_id,itel);
            /* No print function available */
            if ( showdata )
               printf("\nCamera software settings for telescope ID %d ...\n", tel_id);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_POINTINGCOR:
            tel_id = item_header.ident; // Telescope ID is in the header
            if ( (itel = find_tel_idx(tel_id)) < 0 )
            {
               char msg[256];
               snprintf(msg,sizeof(msg)-1,
                  "Pointing correction for unknown telescope ID %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_pointingcor(iobuf,&hsdata->point_cor[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_pointingco(), rc = %d (tel. ID=%d, itel=%d)\n",rc,tel_id,itel);
            /* No print function available */
            if ( showdata )
               print_hess_pointingcor(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_TRACKSET:
            tel_id = item_header.ident; // Telescope ID is in the header
            if ( (itel = find_tel_idx(tel_id)) < 0 )
            {
               char msg[256];
               snprintf(msg,sizeof(msg)-1,
                  "Tracking settings for unknown telescope ID %d.", tel_id);
               Warning(msg);
               exit(1);
            }
            rc = read_hess_trackset(iobuf,&hsdata->tracking_set[itel]);
	    if ( verbose || rc != 0 )
               printf("read_hess_trackset(), rc = %d (tel. ID=%d, itel=%d)\n",rc,tel_id,itel);
            /* No print function available */
            if ( showdata )
               print_hess_trackset(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_EVENT:
            if ( skip_run )
               continue;
            if ( hsdata->mc_shower.energy < min_true_energy )
               continue;
            rc = read_hess_event(iobuf,&hsdata->event,-1);
	    if ( verbose || rc != 0 )
               printf("read_hess_event(), rc = %d\n",rc);
            events++;
#ifdef CHECK_MISSING_PE_LIST
            if ( check_missing_pe_list && ! got_pe_list )
               fprintf(stderr,"Event %d triggered without preceding p.e. list\n",
                  hsdata->mc_event.event);
#endif
            /* Not always all telescopes should get used in the analysis or the DST output */
            if ( num_only > 0 || num_not > 0 )
            {
               for (itel=0; itel<hsdata->run_header.ntel; itel++)
               {
                  if ( hsdata->event.teldata[itel].known )
                  {
                     size_t jtel, jimg;
                     for (jtel=0; jtel<num_not; jtel++)
                     {
                        if ( hsdata->event.teldata[itel].tel_id == not_telescope[jtel] )
			{
                           hsdata->event.teldata[itel].known = 0;
                           if ( hsdata->event.teldata[itel].img != NULL )
			      for ( jimg=0; jimg<(size_t)hsdata->event.teldata[itel].num_image_sets; jimg++  )
			         if ( hsdata->event.teldata[itel].img[jimg].known )
			            hsdata->event.teldata[itel].img[jimg].known = 0;
			   if ( hsdata->event.teldata[itel].raw != NULL )
			      if ( hsdata->event.teldata[itel].raw->known )
			         hsdata->event.teldata[itel].raw->known = 0;
			}
                     }
                     if ( num_only > 0 )
                     {
		     	int keep_known = 0;
                        for (jtel=0; jtel<num_only; jtel++)
                           if ( hsdata->event.teldata[itel].tel_id == only_telescope[jtel] )
                              keep_known = 1;
			if ( !keep_known )
			{
                           if ( hsdata->event.teldata[itel].known )
                              hsdata->event.teldata[itel].known = keep_known;
                           if ( hsdata->event.teldata[itel].img != NULL )
			      for ( jimg=0; jimg<(size_t)hsdata->event.teldata[itel].num_image_sets; jimg++  )
			         if ( hsdata->event.teldata[itel].img[jimg].known )
			            hsdata->event.teldata[itel].img[jimg].known = keep_known;
			   if ( hsdata->event.teldata[itel].raw != NULL )
			      if ( hsdata->event.teldata[itel].raw->known )
			         hsdata->event.teldata[itel].raw->known = keep_known;
                        }
                     }
                  }
               }
            }

            /* Is there a subset that requires hard stereo ? */
            if ( nhard_st > 0 ) 
            {
               int j, jimg;
               int have_st = 0;
               int itel_single = -1;
               for ( j=0; j<nhard_st; j++ )
               {
                  for (itel=0; itel<hsdata->run_header.ntel; itel++)
                  {
                     if ( hsdata->event.teldata[itel].tel_id == hard_stereo[j] &&
                          hsdata->event.teldata[itel].known )
                     {
                        have_st++;
                        itel_single = itel;
                     }
                  }
               }
               /* If there is only a single telescope in the hard stereo subset, we kill it now. */
               if ( have_st == 1 && itel_single >= 0 )
               {
                  itel = itel_single;
                  tel_id = hsdata->event.teldata[itel].tel_id;
                  if ( verbose )
                     printf("Identified single telescope ID %d in hardware stereo subset discarded\n", 
                        tel_id);
                  hsdata->event.teldata[itel].known = 0;
                  if ( hsdata->event.teldata[itel].img != NULL )
		     for ( jimg=0; jimg<(int)hsdata->event.teldata[itel].num_image_sets; jimg++  )
			if ( hsdata->event.teldata[itel].img[jimg].known )
			   hsdata->event.teldata[itel].img[jimg].known = 0;
		  if ( hsdata->event.teldata[itel].raw != NULL )
		     if ( hsdata->event.teldata[itel].raw->known )
			hsdata->event.teldata[itel].raw->known = 0;
                  hsdata->event.central.teltrg_type_mask[itel] = 0;
                  for ( j=0; j<hsdata->event.central.num_teltrg; j++ )
                  {
                     if ( hsdata->event.central.teltrg_list[j] == tel_id )
                     {
                        if ( j+1 < hsdata->event.central.num_teltrg )
                           hsdata->event.central.teltrg_list[j] = 
                              hsdata->event.central.teltrg_list[hsdata->event.central.num_teltrg-1];
                        hsdata->event.central.teltrg_list[hsdata->event.central.num_teltrg-1] = -1;
                        hsdata->event.central.num_teltrg--;
                     }
                  }
                  for ( j=0; j<hsdata->event.central.num_teldata; j++ )
                  {
                     if ( hsdata->event.central.teldata_list[j] == tel_id )
                     {
                        if ( j+1 < hsdata->event.central.num_teldata )
                           hsdata->event.central.teldata_list[j] = 
                              hsdata->event.central.teldata_list[hsdata->event.central.num_teldata-1];
                        hsdata->event.central.teldata_list[hsdata->event.central.num_teldata-1] = -1;
                        hsdata->event.central.num_teldata--;
                     }
                  }
               }
            }

            if ( showdata )
               print_hess_event(iobuf);

            /* First of all fix trigger type bit masks, if necessary */
            if ( ths != NULL && ths->run == hsdata->run_header.run )
            {
               struct trgmask_entry *tme;
               for (itel=0; itel<hsdata->event.central.num_teltrg; itel++ )
               {
                  tel_id = hsdata->event.central.teltrg_list[itel];
                  tme = find_trgmask(ths,item_header.ident,tel_id);
                  if ( tme == NULL )
                  { 
                     if ( verbose )
                        printf("No extra trigger type mask found for telescope %d, setting from %d to 0.\n",
                           tel_id, hsdata->event.central.teltrg_type_mask[itel]);
                     hsdata->event.central.teltrg_type_mask[itel] = 0;
                  }
                  else
                  {
                     if ( verbose )
                        printf("Setting trigger type mask for telescope %d from %d to %d.\n",
                           tel_id, hsdata->event.central.teltrg_type_mask[itel], tme->trg_mask);
                     hsdata->event.central.teltrg_type_mask[itel] = tme->trg_mask;
                  }
               }
            }
            if ( dead_time_fraction > 0. )
            {
               for (itel=0; itel<hsdata->run_header.ntel; itel++)
               {
                  if ( hsdata->event.teldata[itel].known )
                  {
#ifdef WITH_RANDFLAT
                     if ( RandFlat() < dead_time_fraction )
#else
                     if ( drand48() < dead_time_fraction )
#endif
                        hsdata->event.teldata[itel].known = 0;
                  }
               }
            }
            if ( trg_req > 0 )
            {
               int jtel; /* Used to iterate in the list of triggered telescopes */
               for (jtel=0; jtel<hsdata->event.central.num_teltrg; jtel++ )
               {
                  int trg_req_tel = trg_req;
                  tel_id = hsdata->event.central.teltrg_list[jtel];
                  itel = find_tel_idx(tel_id);
                  if ( itel < 0 || itel >= H_MAX_TEL )
                     continue;
                  if ( type_set )
                  {
                     int itype = user_get_type(itel);
                     struct user_parameters *up = user_get_parameters(itype);
                     trg_req_tel = up->i.trg_req;
                  }
                  if ( (trg_req_tel & hsdata->event.central.teltrg_type_mask[jtel]) == 0 )
                  {
                     printf("Telescope %d has trigger pattern %d but requirement is %d. Discarding now.\n",
                        tel_id, hsdata->event.central.teltrg_type_mask[jtel], trg_req_tel);
                     hsdata->event.teldata[itel].known = 0;
                  }
               }
            }
            /* Count number of telescopes (still) present in data and triggered */
            ntel_trg = 0;
            for (itel=0; itel<hsdata->run_header.ntel; itel++)
            {
               if ( hsdata->event.teldata[itel].known )
               {
                  /* If non-triggered telescopes record data (like HEGRA),
                     we may have to check the central trigger bit as well,
                     but ignore this for now. */
                  ntel_trg++;
               }
            }
	    if ( hsdata->event.shower.known )
	       hsdata->event.shower.num_trg = ntel_trg;
            if ( ntel_trg < min_tel_trg )
               continue;
            if ( max_tel_trg > 0 && ntel_trg > max_tel_trg )
               continue;
            if ( user_ana )
               do_user_ana(hsdata,item_header.type,0);
            wsum_trg += pow(hsdata->mc_shower.energy,
               plidx-hsdata->mc_run_header.spectral_index);
            ntrg++;
            if ( reco_flag >= 5 ) 
	       for ( itel=0; itel<hsdata->run_header.ntel; itel++ )
	          if ( hsdata->event.teldata[itel].known )
                  {
                     if ( show_true_pe && hsdata->mc_event.mc_pe_list[itel].npe > 0 )
                        hesscam_ps_plot(ps_fname, hsdata, itel, -1, 3, 0.);
	             hesscam_ps_plot(ps_fname, hsdata, itel, -1, flag_amp_tm, 0.);
                  }

            if ( reco_flag > 1 )
               reconstruct(hsdata, reco_flag, min_amp_tel, min_pix_tel,
	       		tailcut_low_tel, tailcut_high_tel, lref_tel, minfrac_tel, 0 /* -2 */, 
                        flag_amp_tm, cleaning);
            if ( user_ana )
               do_user_ana(hsdata,item_header.type,1);

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
               {
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
               }

               for (itel=0; itel<hsdata->run_header.ntel; itel++)
               {
                  if ( hsdata->event.teldata[itel].known &&
                       hsdata->event.teldata[itel].num_image_sets > 0 )
                  {
                     int jimg;
                     for ( jimg=0; jimg<hsdata->event.teldata[itel].num_image_sets; jimg++ )
                     {
                        if ( hsdata->event.teldata[itel].img[jimg].known )
                        {
                           /* Using only the first image set here. */
                           double direction = (180./M_PI) * hsdata->event.teldata[itel].img[jimg].phi;
                           double orientation = (180./M_PI) * atan2(hsdata->event.teldata[itel].img[jimg].y,
                                                               hsdata->event.teldata[itel].img[jimg].x);
                           double alpha = (180./M_PI) * (hsdata->event.teldata[itel].img[jimg].phi -
                              atan2(hsdata->event.teldata[itel].img[jimg].y,
                                    hsdata->event.teldata[itel].img[jimg].x));
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
                           if ( hsdata->event.teldata[itel].img[jimg].pixels == 0 &&
                                hsdata->event.teldata[itel].image_pixels.pixels > 0 )
                              hsdata->event.teldata[itel].img[jimg].pixels =
                                 hsdata->event.teldata[itel].image_pixels.pixels;
                           printf("@+ %d %d %6.4f %7.2f %7.2f %d %7.5f %7.5f %7.5f %7.5f %7.5f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f   %d  %4.2f %4.2f %4.2f %4.2f %4.2f    %d %d %f\n",
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
                              CALIB_SCALE*hsdata->event.teldata[itel].img[jimg].amplitude,
                              hsdata->event.teldata[itel].img[jimg].pixels,
                              hsdata->event.teldata[itel].img[jimg].w*(180./M_PI),
                              hsdata->event.teldata[itel].img[jimg].l*(180./M_PI),
                              sqrt(hsdata->event.teldata[itel].img[jimg].x*hsdata->event.teldata[itel].img[jimg].x +
                                hsdata->event.teldata[itel].img[jimg].y*hsdata->event.teldata[itel].img[jimg].y) *
                                (180./M_PI),
                              -1., alpha, orientation, direction,
                              hsdata->event.teldata[itel].img[jimg].x*(180./M_PI),
                              hsdata->event.teldata[itel].img[jimg].y*(180./M_PI),
                              hsdata->mc_shower.xmax, hsdata->mc_shower.hmax,
                              hsdata->mc_event.mc_pesum.num_pe[itel],
                              hsdata->event.teldata[itel].img[jimg].num_hot>0 ? CALIB_SCALE*hsdata->event.teldata[itel].img[jimg].hot_amp[0] : -1.,
                              hsdata->event.teldata[itel].img[jimg].num_hot>1 ? CALIB_SCALE*hsdata->event.teldata[itel].img[jimg].hot_amp[1] : -1.,
                              hsdata->event.teldata[itel].img[jimg].num_hot>2 ? CALIB_SCALE*hsdata->event.teldata[itel].img[jimg].hot_amp[2] : -1.,
                              hsdata->event.teldata[itel].img[jimg].num_hot>3 ? CALIB_SCALE*hsdata->event.teldata[itel].img[jimg].hot_amp[3] : -1.,
                              hsdata->event.teldata[itel].img[jimg].num_hot>4 ? CALIB_SCALE*hsdata->event.teldata[itel].img[jimg].hot_amp[4] : -1.,
                              jimg,
                              hsdata->event.teldata[itel].img[jimg].cut_id,
                              hsdata->event.teldata[itel].img[jimg].clip_amp);
                              
                           if ( hsdata->event.teldata[itel].img[jimg].num_hot>0 &&
                                hsdata->event.teldata[itel].raw != NULL && 
                                hsdata->event.teldata[itel].raw ->known )
                           {
                              int ihp = hsdata->event.teldata[itel].img[jimg].hot_pixel[0];
                              printf("Calibration cross check for hottest pixel of telescope %d: "
                                 "%4.2f versus %4.2f peak p.e. in pixel %d.\n",
                                 hsdata->camera_set[itel].tel_id,
                                 CALIB_SCALE*hsdata->event.teldata[itel].img[jimg].hot_amp[0],
                                 calibrate_pixel_amplitude(hsdata, itel, ihp, flag_amp_tm, -1, 0.), 
                                 ihp);
                           }
                        }
                     }
                  }
               }
            }

	    for (itel=0; itel<hsdata->run_header.ntel; itel++)
            {
	       if ( hsdata->event.teldata[itel].known )
               {
                  if ( show_true_pe && hsdata->mc_event.mc_pe_list[itel].npe > 0 )
                     hesscam_ps_plot(ps_fname, hsdata, itel, -1, 3, 0.);
	          hesscam_ps_plot(ps_fname, hsdata, itel, -1, flag_amp_tm, last_clip_amp);
               }
            }
            if ( show_type_sum )
            {
               int itype;
               if ( ntypes == 0 )
               {
                  for (itel=0; itel<hsdata->run_header.ntel; itel++)
                  {
                     itype = user_get_type(itel);
                     if ( itype >= ntypes )
                        ntypes = itype+1;
                  }
               }
               for ( itype=0; itype<ntypes; itype++ )
                  hesscam_type_sum_plot(ps_fname, hsdata, itype);
            }
            if ( (dst_level >= 0 && dst_level <= 3) ||
                 (dst_level >= 10 && dst_level <= 13 && user_ana && 
                    (last_event_selected = user_selected_event()) != 0 ) ||
                  dst_level >= 100 )
            {
               if ( dst_level != 100 && !no_mc_data )
               {
                  if ( mc_shower_stored == 0 && dst_level > 1 )
                  {
                     if ( write_hess_mc_shower(iobuf,&hsdata->mc_shower) < 0 )
                     {
                        fprintf(stderr,"Delayed writing of MC shower data failed.\n");
                        exit(1);
                     }
                     mc_shower_stored = 1;
                  }
                  if ( mc_event_stored == 0 && dst_level > 1  )
                  {
                     if ( write_hess_mc_event(iobuf,&hsdata->mc_event) < 0 )
                     {
                        fprintf(stderr,"Delayed writing of MC event data failed.\n");
                        exit(1);
                     }
                     mc_event_stored = 1;
                  }
                  if ( mc_pe_sum_stored == 0 && (dst_level%10) == 2 &&
                       last_mc_pe_sum_id == item_header.ident )
                  {
                     if ( write_hess_mc_pe_sum(iobuf,&hsdata->mc_event.mc_pesum) < 0 )
                     {
                        fprintf(stderr,"Delayed writing of MC p.e. sum data failed.\n");
                        exit(1);
                     }
                     mc_pe_sum_stored = 1;
                  }
               }

               if ( dst_level >= 0 )
               {
                  for (itel=0; itel<hsdata->run_header.ntel; itel++)
                  {
                     TelEvent *te = &hsdata->event.teldata[itel];
                     if ( te->known && te->raw != NULL && te->raw->known )
                     {
                        /* Apply "pure raw data" cleaning only where raw data is actually available */
                        if ( pure_raw )
                        {
                           int j;
                           if ( te->pixtm != NULL )
                              te->pixtm->known = 0;
                           if ( te->pixcal != NULL )
                              te->pixcal->known = 0;
                           if ( te->img != NULL )
                           {
                              for ( j=0; j<te->num_image_sets; j++ )
                                 te->img[j].known = 0;
                              te->num_image_sets = 0;
                           }
                           /* We keep the list of triggered pixels but discard image selected pixels */
                           te->image_pixels.pixels = 0;
                        }
                        /* Calibrated pixel intensities requested to be stored */
                        if ( do_calibrate && te->pixcal != NULL )
                        {
                           if ( ! te->pixcal->known ) /* Not filled already in reconstruction ? */
                           {
                              int tel_type = user_get_type(itel);
                              struct user_parameters *up = user_get_parameters(tel_type);
                              double clip_amp = up->d.clip_amp;
                              calibrate_amplitude(hsdata, itel, flag_amp_tm, clip_amp);
                           }
                           if ( only_calibrated )
                              te->raw->known = 0;
                        }
                        if ( dst_level != 0 && dst_level != 10 && dst_level < 100 )
                           te->raw->known = 0;
                     }
                  }
                  if ( pure_raw )
                     hsdata->event.shower.known = 0;
               }
               if ( (dst_level%10) >= 2 )
               {
                  double lgE = log10(hsdata->mc_shower.energy);
                  fill_histogram_by_ident(11001,lgE,0.,
                     pow(hsdata->mc_shower.energy, 
                         plidx-hsdata->mc_run_header.spectral_index));
                  fill_histogram_by_ident(11101,lgE,0.,1.);
               }
               if ( zero_suppression >= 0 )
               {
                  for (itel=0; itel<hsdata->run_header.ntel; itel++)
                     if ( hsdata->event.teldata[itel].known &&
                          hsdata->event.teldata[itel].raw != NULL )
                     {
                        hsdata->event.teldata[itel].raw->zero_sup_mode = zero_suppression & 0x1f;
                     }
               }
               if ( dst_level == 0 )
                  rc = write_hess_event(iobuf,&hsdata->event,-1 ^ RAWDATA_FLAG );
               else if ( dst_level == 10 || cleaning > 0 || dst_level >= 100 )
                  rc = write_hess_event(iobuf,&hsdata->event,-1);
               else
                  rc = write_hess_event(iobuf,&hsdata->event,
                        -1 ^ (RAWDATA_FLAG|RAWSUM_FLAG|TIME_FLAG|TRACKDATA_FLAG));
               if ( rc != 0 )
               {
                  fprintf(stderr,"Writing of event data failed (rc=%d).\n", rc);
                  exit(1);
               }
            }

            if ( ntuple_file != NULL && hsdata->event.shower.known && 
                 bnt.n_img >= 2 && bnt.lg_e > -3. )
            {
               list_ntuple(ntuple_file,&bnt,1);
            }
            break;

         /* =================================================== */
         case IO_TYPE_HESS_CALIBEVENT:
            if ( skip_run )
               continue;
         {
            int type = -1;
            rc = read_hess_calib_event(iobuf,&hsdata->event,-1,&type);
	    if ( verbose || rc != 0 )
               printf("read_hess_calib_event(), rc = %d, type=%d\n",rc,type);
            if ( showdata )
               print_hess_calib_event(iobuf);
	    for (itel=0; itel<hsdata->run_header.ntel; itel++)
	       if ( hsdata->event.teldata[itel].known )
	          hesscam_ps_plot(ps_fname, hsdata, itel, type, 0, 0.);
         }
            break;

         /* =================================================== */
         case IO_TYPE_HESS_MC_SHOWER:
            if ( skip_run )
               continue;
            rc = read_hess_mc_shower(iobuf,&hsdata->mc_shower);
	    if ( verbose || rc != 0 )
               printf("read_hess_mc_shower(), rc = %d (shower %d)\n",
                  rc, hsdata->mc_shower.shower_num);
            if ( showdata )
               print_hess_mc_shower(iobuf);
            if ( user_ana )
               do_user_ana(hsdata,item_header.type,0);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_MC_EVENT:
#ifdef CHECK_MISSING_PE_LIST
            got_pe_list = 0;
#endif
            if ( skip_run )
               continue;
            rc = read_hess_mc_event(iobuf,&hsdata->mc_event);
	    if ( verbose || rc != 0 )
               printf("read_hess_mc_event(), rc = %d (event %d, shower %d)\n",
                  rc, hsdata->mc_event.event, hsdata->mc_event.shower_num);
            if ( showdata )
               print_hess_mc_event(iobuf);
            if ( user_ana )
               do_user_ana(hsdata,item_header.type,0);
            if ( (dst_level%10) >= 2 )
               mc_event_fill(hsdata,plidx-hsdata->mc_run_header.spectral_index);
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
#ifdef CHECK_MISSING_PE_LIST
            got_pe_list = 1;
#endif
            if ( skip_run )
               continue;
            if ( hsdata && hsdata->run_header.ntel > 0 )
            {
               rc = read_hess_mc_phot(iobuf,&hsdata->mc_event);
	       if ( verbose || rc != 0 )
                  printf("read_hess_mc_phot(), rc = %d\n",rc);
            }
            if ( showdata )
               print_hess_mc_phot(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_MC_TELARRAY_HEAD:
            if ( showdata )
               printf("\nStart of split photon bunch blocks (id=%ld)\n",
                  item_header.ident);
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
         case IO_TYPE_MC_TELARRAY_END:
            if ( showdata )
               printf("\nEnd of split photon bunch blocks (id=%ld)\n",
                  item_header.ident);
            break;

         /* =================================================== */
         case IO_TYPE_MC_LONGI:
            if ( showdata )
               printf("\nShower longitudinal distributions (not shown in detail) ...\n");
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
         case IO_TYPE_MC_TELPOS:
            if ( showdata )
               print_tel_pos(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_MC_TELOFF:
            if ( showdata )
               print_tel_offset(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_MC_PE_SUM:
            if ( skip_run )
               continue;
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
               printf("read_hess_tel_monitor(), rc = %d (tel. ID=%d, itel=%d)\n",
                  rc, tel_id, itel);
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
               printf("read_hess_laser_calib(), rc = %d (tel. ID=%d, itel=%d)\n",
                  rc, tel_id, itel);
            if ( showdata )
               print_hess_laser_calib(iobuf);
            /* Optionally add random pixel calibration error */
            if ( calerror > 0.0 )
            {
               int igain, ipix;
               LasCalData *cal = &hsdata->tel_lascal[itel];
               for ( igain=0; igain<cal->num_gains; igain++ )
               {
                  for ( ipix=0; ipix<cal->num_pixels; ipix++ )
                  {
#ifdef WITH_RANDFLAT
                     double g = RandGauss(0.,calerror);
#else
                     double g = grand48(0.,calerror);
#endif
                     cal->calib[igain][ipix] *= exp(g);
                  }
               }
            }
            break;

         /* =================================================== */
         case IO_TYPE_HESS_RUNSTAT:
            if ( skip_run )
               continue;
            rc = read_hess_run_stat(iobuf,&hsdata->run_stat);
	    if ( verbose || rc != 0 )
               printf("read_hess_run_stat(), rc = %d\n",rc);
            if ( showdata )
               print_hess_run_stat(iobuf);
            if ( user_ana )
               do_user_ana(hsdata,item_header.type,0);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_MC_RUNSTAT:
            if ( skip_run )
               continue;
            rc = read_hess_mc_run_stat(iobuf,&hsdata->mc_run_stat);
	    if ( verbose || rc != 0 )
               printf("read_hess_mc_run_stat(), rc = %d\n",rc);
            if ( showdata )
               print_hess_mc_run_stat(iobuf);
            if ( user_ana )
               do_user_ana(hsdata,item_header.type,0);
            if ( (dst_level%10) >= 2 )
               write_dst_histos(iobuf2);
            break;

         /* =================================================== */
         /* In-data trigger pattern block, or in extra file (not auto-loaded) */
         case IO_TYPE_HESS_XTRGMASK:
            if ( tms == NULL )
               tms = (struct trgmask_set *) calloc(1,sizeof(struct trgmask_set));
            if ( ths == NULL )
               ths = (struct trgmask_hash_set *) calloc(1,sizeof(struct trgmask_hash_set));
            rc = read_trgmask(iobuf, tms);
	    if ( verbose || rc != 0 )
               printf("read_trgmask(), rc = %d\n",rc);
            if ( showdata )
               print_trgmask(iobuf);
            trgmask_fill_hashed(tms,ths);
            break;

         /* =================================================== */
         /* (End-of-job or DST) histograms */
         case 100:
            {
               int nhist;
#if 0
               HISTOGRAM *h_list[10];
               /* If the file contained output from multiple jobs, */
               /* existing histograms (e.g. ID #1) will be replaced. */
               if ( get_histogram_by_ident(1) != NULL )
                  printf("Existing histograms will be replaced!\n");
               nhist = read_histograms(h_list,10,iobuf);
#else
               HISTOGRAM *h_list[1] = { NULL };
               /* Histograms are added to existing ones of the same ID. */
               nhist = read_histograms(h_list,-1,iobuf);
#endif
               if ( !quiet )
                  printf("\nread_histograms() got %d histograms\n",nhist);
               if ( showdata )
               {
                  print_histograms(iobuf);
                  printf("\n");
               }

               /* Histogram 99 is actually a list of configuration parameters. We cannot add up those. */
               if ( nhist == 1 && h_list[0] != NULL && h_list[0]->ident == 99 )
               {
                  double plidx_new = plidx;
                  read_histograms(h_list,1,iobuf); /* Read again, replacing histogram */
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
                  if ( h_list[0]->nbins >= 5 &&
                       (fabs(h_list[0]->extension->ddata[3]-tailcut_low) > 0.01 ||
                        fabs(h_list[0]->extension->ddata[4]-tailcut_high) > 0.01) )
                  {
                     fprintf(stderr,"\nYou asked for images with tail-cuts %3.1f,%3.1f but "
                        "you are reading MC DST data produced with tail-cuts %3.1f,%3.1f\n\n",
                        tailcut_low, tailcut_high,
                        h_list[0]->extension->ddata[3], h_list[0]->extension->ddata[4]);
                  }
               }
            }
            break;

         default:
            /* Some data block types normally embedded in type 2010. Only showing contents. */
            if ( item_header.type == 2009 )
            {
               if ( showdata )
               {
                  printf("Stand-alone central trigger data block:\n");
                  print_hess_centralevent(iobuf);
               }
               rc = read_hess_centralevent(iobuf,&hsdata->event.central);
	       if ( verbose || rc != 0 )
                  printf("read_hess_centralevent(), rc = %d\n",rc);
               continue;
            }
            else if ( item_header.type > 2000 && item_header.type < 10000 )
            {
               if ( item_header.type%1000 >= 100 && item_header.type%1000 <= 199 )
               {
                  int tid = (item_header.type - IO_TYPE_HESS_TRACKEVENT)%100 +
                            100 * ((item_header.type - IO_TYPE_HESS_TRACKEVENT)/1000);
                  if ( showdata )
                  {
                     printf("Stand-alone tracking data block for telescope %d:\n", tid);
                     print_hess_trackevent(iobuf);
                  }
	          if ( (itel = find_tel_idx(tid)) < 0 || itel >= H_MAX_TEL )
	          {
	             fprintf(stderr,"Telescope number out of range for tracking data (telescope ID = %d)\n", tid);
	             continue;
	          }
                  rc = read_hess_trackevent(iobuf,&hsdata->event.trackdata[itel]);
	          if ( verbose || rc != 0 )
                     printf("read_hess_trackevent(), rc = %d\n",rc);
                  continue;
               }
               else if ( item_header.type%1000 >= 200 && item_header.type%1000 <= 299 )
               {
                  int tid = (item_header.type - IO_TYPE_HESS_TELEVENT)%100 +
                            100 * ((item_header.type - IO_TYPE_HESS_TELEVENT)/1000);
                  if ( showdata )
                  {
                     printf("Stand-alone telescope event data block for telescope %d:\n", tid);
                     print_hess_televent(iobuf);
                  }
	          if ( (itel = find_tel_idx(tid)) < 0 || itel >= H_MAX_TEL )
	          {
	             fprintf(stderr,"Telescope number out of range for telescope event data (telescope ID = %d)\n", tid);
	             continue;
	          }
                  rc = read_hess_televent(iobuf,&hsdata->event.teldata[itel],-1);
	          if ( verbose || rc != 0 )
                     printf("read_hess_televent(), rc = %d\n",rc);

                  if ( ps_fname != NULL ) 
	             if ( hsdata->event.teldata[itel].known )
                     {
                        if ( show_true_pe && hsdata->mc_event.mc_pe_list[itel].npe > 0 )
                           hesscam_ps_plot(ps_fname, hsdata, itel, -1, 3, 0.);
	                hesscam_ps_plot(ps_fname, hsdata, itel, -1, flag_amp_tm, 0.);
                     }
                  continue;
               }
            }
            if ( !ignore )
               fprintf(stderr,"Ignoring unknown data block type %ld\n",item_header.type);
      }
    } /* Loop over all data in a file */
    
    /* ================ Done with this input data file ============== */

    if ( iobuf->input_file != NULL && iobuf->input_file != stdin )
      fileclose(iobuf->input_file);
    iobuf->input_file = NULL;
    reset_io_block(iobuf);

    if ( hsdata != NULL && !quiet )
      show_run_summary(hsdata,nev,ntrg,plidx,wsum_all,wsum_trg,rmax_x,rmax_y,rmax_r);
    else if ( nev > 0 )
      printf("%d of %d events triggered\n", ntrg, nev);

    if ( user_ana && hsdata != NULL )
       do_user_ana(hsdata,0,0);

    /* If the RUNSTAT block was missing we should still record our stats. */
    if ( (dst_level%10) >= 2 )
       write_dst_histos(iobuf2);

    if ( hsdata != NULL )
       hsdata->run_header.run = 0;
   } /* Loop over all input files */

   /* ============== Done with all input data files ============== */

   if ( user_ana && hsdata != NULL )
      do_user_ana(hsdata,0,1);

   /* Close output files or pipes. */
   if ( iobuf2 != NULL && iobuf != NULL && 
        iobuf2->output_file != NULL && 
        iobuf2->output_file != iobuf->output_file )
      fileclose(iobuf2->output_file);
   if ( iobuf != NULL && iobuf->output_file != NULL )
      fileclose(iobuf->output_file);
   if ( ntuple_file != NULL )
      fileclose(ntuple_file);

   if ( iobuf2 != NULL )
   {
      free_io_buffer(iobuf2);
      iobuf2 = NULL;
   }
   if ( iobuf != NULL )
   {
      free_io_buffer(iobuf);
      iobuf = NULL;
   }

#ifdef TEST_MEMORY_LEAKS
   /* Optionally free most of the allocated memory to make sure */
   /* we have no leaks anywhere (using valgrind or such). */
   /* It is not actually needed though. */

   free_all_histograms();
   histogram_hashing(0);

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
      if ( hsdata->event.teldata[itel].pixcal != NULL )
      {
         free(hsdata->event.teldata[itel].pixcal);
         hsdata->event.teldata[itel].pixcal = NULL;
      }
      if ( hsdata->event.teldata[itel].img != NULL )
      {
         free(hsdata->event.teldata[itel].img);
         hsdata->event.teldata[itel].img = NULL;
      }
      
      deallocate_nb_list(itel);
   }
   if ( hsdata->run_header.target != NULL )
   {
      free(hsdata->run_header.target);
      hsdata->run_header.target = NULL;
   }
   if ( hsdata->run_header.observer != NULL )
   {
      free(hsdata->run_header.observer);
      hsdata->run_header.observer = NULL;
   }
   {
      int j;
      for ( j=0; j<H_MAX_PROFILE; j++)
      {
         if ( hsdata->mc_shower.profile[j].content != NULL )
         {
            free(hsdata->mc_shower.profile[j].content);
            hsdata->mc_shower.profile[j].content = NULL;
         }
      }
   }
   free(hsdata);
   hsdata = NULL;

   initpath("");
   initexepath(NULL);

   { /* Among the still allocated memory is a linked list of data files processed */
      NextFile *this_file = &first_file;
      for ( next_file = NULL; this_file != NULL; this_file = next_file )
      {
         if ( this_file->fname != NULL )
         {
            free(this_file->fname);
            this_file->fname = NULL;
         }
         next_file = this_file->next;
         this_file->next = NULL;
         if ( this_file != &first_file )
            free(this_file);
      }
   }
   
   /* At this point the only still allocated memory blocks should be the
      command line history and the user analysis parameters. */ 
#endif

   return 0;
}

/** @} */
