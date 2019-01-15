/* ============================================================================

Copyright (C) 2001, 2003, 2008, 2009  Konrad Bernloehr

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

/** @file camera_image.c
 *  @short Plot a camera image from H.E.S.S./CTA data.
 *
 *  This code is derived from sim_conv2hess.c but now getting
 *  the relevant data from the data structure filled after
 *  reading the eventio based data, rather than from the internal
 *  data structures of sim_hessarray. As a consequence not all
 *  information available in the sim_hessarray generated plots
 *  is available in the plots generated here. Also some flexibility
 *  is lost, concerning for example the pixel shape which is not
 *  included in the data.
 *
 *  @author  Konrad Bernloehr
 *  @date    @verbatim CVS $Date: 2018/09/18 15:10:17 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.42 $ @endverbatim
 */

#include "initial.h"      /* This file includes others as required. */
#include "io_basic.h"     /* This file includes others as required. */
#include "mc_tel.h"
#include "history.h"
#include "io_hess.h"
#include "fileopen.h"
#include "rec_tools.h"
#include "reconstruct.h"
#include "camera_image.h"
#include "user_analysis.h"

static char ps_head1a[] =
"%!PS-Adobe-2.0\n"
#ifdef CTA
"%%Title: CTA Telescope Simulation";
#else
"%%Title: H.E.S.S. Telescope Simulation";
#endif
static char ps_head1b[] =
"\n%%Creator:";

static char ps_head2[] =
#ifdef CTA
"%%For: Cherenkov Telescope Array\n"
#else
"%%For: High Energy Stereoscopic System\n"
#endif
"%%Orientation: Portrait\n"
"%%BoundingBox: 0 0 595 841\n"
"%%Pages: atend\n"
/* "%%BeginSetup\n" */
"%%IncludeFeature: *PageSize A4\n"
/* "%%EndSetup\n" */
"%%Magnification: 1.0000\n"
"%%EndComments\n"
"\n"
"%%BeginSetup\n"
"0.0 842.0 translate\n"
"1 -1 scale\n"
"\n"
"/black {0.0 0.0 0.0 srgb} bind def\n"
"/red {1.0 0.0 0.0 srgb} bind def\n"
"/blue {0.0 0.0 1.0 srgb} bind def\n"
"/yellow {1.0 1.0 0.0 srgb} bind def\n"
"/magenta {1.0 0.0 1.0 srgb} bind def\n"
"/green {0.0 1.0 0.0 srgb} bind def\n"
"/yxcr {gs -80 -80 rm 160 160 rl cp yellow s n gr \n"
" -80 80 rm 160 -160 rl cp gs yellow s gr} bind def\n"
"/mxcr {gs -80 -80 rm 160 160 rl cp magenta s n gr \n"
" -80 80 rm 160 -160 rl cp gs magenta s gr} bind def\n"
"/rxcr {gs -100 -100 rm 200 200 rl cp 30 slw red s n gr \n"
" -100 100 rm 200 -200 rl cp gs 30 slw red s gr} bind def\n"
"/gcr {gs -80 0 rm 160 0 rl cp green s n gr \n"
" 0 80 rm 0 -160 rl cp gs green s gr} bind def\n"
"/txt10 {/Helvetica-Bold-iso ff 1000.00 scf sf} bind def\n"
"/txt8 {/Helvetica-Bold-iso ff 800.00 scf sf} bind def\n"
"/txt6 {/Helvetica-Bold-iso ff 550.00 scf sf} bind def\n"
"/txt5 {/Helvetica-Bold-iso ff 480.00 scf sf} bind def\n"
"/txt4 {/Helvetica-Bold-iso ff 400.00 scf sf} bind def\n"
"/txt3 {/Helvetica-Bold-iso ff 320.00 scf sf} bind def\n"
"/txt2 {/Helvetica-iso ff 250.00 scf sf} bind def\n"
"/txt1 {/Helvetica-iso ff 100.00 scf sf} bind def\n"
"/txt125 {/Helvetica-iso ff 125.00 scf sf} bind def\n"
"/txt100 {/Helvetica-iso ff 100.00 scf sf} bind def\n"
"/txt80 {/Helvetica-iso ff 80.00 scf sf} bind def\n"
"/txt70 {/Helvetica-iso ff 70.00 scf sf} bind def\n"
"/txt60 {/Helvetica-iso ff 60.00 scf sf} bind def\n"
"/txt50 {/Helvetica-iso ff 50.00 scf sf} bind def\n"
"/mtxt {m gs 1 -1 sc} bind def\n"
"/tblack {dup sw pop 2 div neg 0 rm black sh gr} bind def\n"
"/tblue {dup sw pop 2 div neg 0 rm blue sh gr} bind def\n"
"/tred {dup sw pop 2 div neg 0 rm red sh gr} bind def\n"
"/rtblack {gs 1 -1 sc black sh gr} bind def\n"
"/cp {closepath} bind def\n"
"/ef {eofill} bind def\n"
"/gr {grestore} bind def\n"
"/gs {gsave} bind def\n"
"/sa {save} bind def\n"
"/rs {restore} bind def\n"
"/l {lineto} bind def\n"
"/rl {rlineto} bind def\n"
"/m {moveto} bind def\n"
"/rm {rmoveto} bind def\n"
"/n {newpath} bind def\n"
"/s {stroke} bind def\n"
"/sh {show} bind def\n"
"/slc {setlinecap} bind def\n"
"/slj {setlinejoin} bind def\n"
"/slw {setlinewidth} bind def\n"
"/srgb {setrgbcolor} bind def\n"
"/rot {rotate} bind def\n"
"/sc {scale} bind def\n"
"/sd {setdash} bind def\n"
"/ff {findfont} bind def\n"
"/sf {setfont} bind def\n"
"/scf {scalefont} bind def\n"
"/sw {stringwidth} bind def\n"
"/tr {translate} bind def\n"
"/tnt {dup dup currentrgbcolor\n"
"  4 -2 roll dup 1 exch sub 3 -1 roll mul add\n"
"  4 -2 roll dup 1 exch sub 3 -1 roll mul add\n"
"  4 -2 roll dup 1 exch sub 3 -1 roll mul add srgb}\n"
"  bind def\n"
"/shd {dup dup currentrgbcolor 4 -2 roll mul 4 -2 roll mul\n"
"  4 -2 roll mul srgb} bind def\n";

static char ps_head3[] = 
"/reencdict 12 dict def /ReEncode { reencdict begin\n"
"/newcodesandnames exch def /newfontname exch def /basefontname exch def\n"
"/basefontdict basefontname findfont def /newfont basefontdict maxlength dict def\n"
"basefontdict { exch dup /FID ne { dup /Encoding eq\n"
"{ exch dup length array copy newfont 3 1 roll put }\n"
"{ exch newfont 3 1 roll put } ifelse } { pop pop } ifelse } forall\n"
"newfont /FontName newfontname put newcodesandnames aload pop\n"
"128 1 255 { newfont /Encoding get exch /.notdef put } for\n"
"newcodesandnames length 2 idiv { newfont /Encoding get 3 1 roll put } repeat\n"
"newfontname newfont definefont pop end } def\n"
"/isovec [\n"
"8#055 /minus 8#200 /grave 8#201 /acute 8#202 /circumflex 8#203 /tilde\n"
"8#204 /macron 8#205 /breve 8#206 /dotaccent 8#207 /dieresis\n"
"8#210 /ring 8#211 /cedilla 8#212 /hungarumlaut 8#213 /ogonek 8#214 /caron\n"
"8#220 /dotlessi 8#230 /oe 8#231 /OE\n"
"8#240 /space 8#241 /exclamdown 8#242 /cent 8#243 /sterling\n"
"8#244 /currency 8#245 /yen 8#246 /brokenbar 8#247 /section 8#250 /dieresis\n"
"8#251 /copyright 8#252 /ordfeminine 8#253 /guillemotleft 8#254 /logicalnot\n"
"8#255 /hyphen 8#256 /registered 8#257 /macron 8#260 /degree 8#261 /plusminus\n"
"8#262 /twosuperior 8#263 /threesuperior 8#264 /acute 8#265 /mu 8#266 /paragraph\n"
"8#267 /periodcentered 8#270 /cedilla 8#271 /onesuperior 8#272 /ordmasculine\n"
"8#273 /guillemotright 8#274 /onequarter 8#275 /onehalf\n"
"8#276 /threequarters 8#277 /questiondown 8#300 /Agrave 8#301 /Aacute\n"
"8#302 /Acircumflex 8#303 /Atilde 8#304 /Adieresis 8#305 /Aring\n"
"8#306 /AE 8#307 /Ccedilla 8#310 /Egrave 8#311 /Eacute\n"
"8#312 /Ecircumflex 8#313 /Edieresis 8#314 /Igrave 8#315 /Iacute\n"
"8#316 /Icircumflex 8#317 /Idieresis 8#320 /Eth 8#321 /Ntilde 8#322 /Ograve\n"
"8#323 /Oacute 8#324 /Ocircumflex 8#325 /Otilde 8#326 /Odieresis 8#327 /multiply\n"
"8#330 /Oslash 8#331 /Ugrave 8#332 /Uacute 8#333 /Ucircumflex\n"
"8#334 /Udieresis 8#335 /Yacute 8#336 /Thorn 8#337 /germandbls 8#340 /agrave\n"
"8#341 /aacute 8#342 /acircumflex 8#343 /atilde 8#344 /adieresis 8#345 /aring\n"
"8#346 /ae 8#347 /ccedilla 8#350 /egrave 8#351 /eacute\n"
"8#352 /ecircumflex 8#353 /edieresis 8#354 /igrave 8#355 /iacute\n"
"8#356 /icircumflex 8#357 /idieresis 8#360 /eth 8#361 /ntilde 8#362 /ograve\n"
"8#363 /oacute 8#364 /ocircumflex 8#365 /otilde 8#366 /odieresis 8#367 /divide\n"
"8#370 /oslash 8#371 /ugrave 8#372 /uacute 8#373 /ucircumflex\n"
"8#374 /udieresis 8#375 /yacute 8#376 /thorn 8#377 /ydieresis] def\n"
"/Helvetica-Bold /Helvetica-Bold-iso isovec ReEncode\n"
"/Helvetica /Helvetica-iso isovec ReEncode\n"
"%%EndSetup\n"
"\n";

static char ps_begin_page1[] =
"%%Page: ";
static char ps_begin_page2[] =
"save\n"
"10 setmiterlimit\n"
"n -1000 31000 m -1000 -1000 l 22000 -1000 l 22000 31000 l cp clip\n"
" 0.02835 0.02835 sc\n"
"gs\n"
"7.500 slw\n"
"black\n";
static char ps_end_page[] =
"gr\n"
"showpage\n";

static char ps_trailer[] =
"rs\n";
/*
static char polyline_begin[] = "n %d %d m";
static char polyline_next[] = " %d %d l";
static char polyline_end[] = 
"\n cp gs %5.3f %5.3f %5.3f srgb 1.00 shd ef gr gs black s gr\n";
*/

static char alt_az_arrow[] =
      "n 18000 26000 m "
      "0 100 rl 200 -100 rl -200 -100 rl 0 100 rl -1000 0 rl "
      "cp gs 20 slw black s gr\n"
      "txt5 18700 26100 mtxt (Az) tblack\n"
      "n 17000 25000 m "
      "100 0 rl -100 -200 rl -100 200 rl 100 0 rl 0 1000 rl "
      "cp gs 20 slw black s gr\n"
      "txt5 17000 24600 mtxt (Alt) tblack\n"
      "gs 17800 25500 tr %f rot -17800 -25500 tr\n"
        "n 17800 25500 m "
        "0 100 rl 200 -100 rl -200 -100 rl 0 100 rl -300 0 rl "
        "cp gs 10 slw black s gr\n"
        "txt2 17950 25350 mtxt (y) tblack\n"
        "n 17500 25200 m "
        "100 0 rl -100 -200 rl -100 200 rl 100 0 rl 0 300 rl "
        "cp gs 10 slw black s gr\n"
        "txt2 17700 25200 mtxt (x) tblack\n"
      "gr\n";

static int guessed_pixel_shape_type(CameraSettings *camset, int itel);

static double dist2 (double x, double y) { return sqrt(x*x + y*y); }

static void print_pix_col(double n_o_r, FILE *psfile, double gamma_coeff);

/** Print a false-colour RGB value for a pixel intensity. */

static void print_pix_col(double n_o_r, FILE *psfile, double gamma_coeff)
{
   if ( gamma_coeff > 0. )
   {
      /* Colour image */
      if ( n_o_r <= 0. )
	 fprintf(psfile,"%5.3f %5.3f %5.3f pxe\n",0.,0.,0.);
      else if ( n_o_r >= 1. )
	 fprintf(psfile,"%5.3f %5.3f %5.3f pxe\n",
	    1.,pow(1./n_o_r,0.8*gamma_coeff),pow(1./n_o_r,1.2*gamma_coeff));
      else
	 fprintf(psfile,"%5.3f %5.3f %5.3f pxe\n",
	    pow(n_o_r,0.8*gamma_coeff),
	    pow(n_o_r,0.6*gamma_coeff),
	    pow(n_o_r,0.4*gamma_coeff) );
   }
   else
   {
      /* Gray-scale image */
      gamma_coeff *= -1.;
      if ( n_o_r <= 0. )
	 fprintf(psfile,"%5.3f %5.3f %5.3f pxe\n",0.,0.,0.);
      else if ( n_o_r >= 3.0 )
	 fprintf(psfile,"%5.3f %5.3f %5.3f pxe\n",1.,1.,1.);
      else
	 fprintf(psfile,"%5.3f %5.3f %5.3f pxe\n",
	    pow(n_o_r/3.,0.5*gamma_coeff),
	    pow(n_o_r/3.,0.5*gamma_coeff),
	    pow(n_o_r/3.,0.5*gamma_coeff) );
   }
}

static int ps_num_page = 0;

/** 
 *  @short Write PostScript of camera sum image or sample image to a dedicated file.
 *
 *  Also controlled via environment variables 
 *  GAMMA_COEFF, GRAY_IMAGE, IMAGE_RANGE, IMAGE_OFFSET for image colors, 
 *  PLOT_WITH_PIXEL_ID, PLOT_WITH_PIXEL_AMP, PLOT_WITH_PIXEL_PE for overlay text,
 *  SHOW_TRUE_PE for showing color for true p.e. number in place of calibrated amplitude.
 *
 *  @param image_fname The name of the postscript image file. Opened for appending new images.
 *  @param hsdata      Pointer to the structure containing all data.
 *  @param itel        The telescope index number.
 *  @param type        Event type (<0: MC events, >=0: various type of calib data).
 *  @param amp_tm      0: Use normal integrated amplitude.
 *                     1: Use integration around global peak position from
 *                        pulse shape analysis. May include all pixels or only selected.
 *                     2: Use integration around local peak position from
 *                        pulse shape analysis. Return 0 for pixels without
 *                        a fairly significant peak.
 *                     3: Show only true p.e. content as amplitude (no samples).
 *  @param clip_amp: if >0, any calibrated amplitude is clipped not to exceed this value [mean p.e.].
 */

void hesscam_ps_plot(const char *image_fname, AllHessData *hsdata, int itel, 
   int type, int amp_tm, double clip_amp);

void hesscam_ps_plot(const char *image_fname, AllHessData *hsdata, int itel, 
   int type, int amp_tm, double clip_amp)
{
   int run = hsdata->run_header.run;
   int event = (type<0 ?
                hsdata->event.central.glob_count :
                hsdata->event.teldata[itel].loc_count*100);
   int tel = hsdata->camera_set[itel].tel_id;
   int itel_pe = itel;
   int npix = hsdata->camera_set[itel].num_pixels;
   double flen = hsdata->camera_set[itel].flen;
   double *x = &hsdata->camera_set[itel].xpix[0];
   double *y = &hsdata->camera_set[itel].ypix[0];
   TelEvent *teldata;
   FILE *psfile;
   // int ngt;
   int i, j, n;
   double xc, yc;
   double scale = 150.;
   double body_diameter = 0.;
   double gamma_coeff = (getenv("GAMMA_COEFF")==NULL)?0.65:atof(getenv("GAMMA_COEFF"));
   double img_gamma = (getenv("GRAY_IMAGE")==NULL)?gamma_coeff:-0.85*gamma_coeff;
   double img_range = (getenv("IMAGE_RANGE")==NULL)?20.:atof(getenv("IMAGE_RANGE"));
   double img_off = (getenv("IMAGE_OFFSET")==NULL)?4.:atof(getenv("IMAGE_OFFSET"));
   int with_id = (getenv("PLOT_WITH_PIXEL_ID")==NULL)?0:atoi(getenv("PLOT_WITH_PIXEL_ID"));
   int with_amp = (getenv("PLOT_WITH_PIXEL_AMP")==NULL)?0:atoi(getenv("PLOT_WITH_PIXEL_AMP"));
   int with_npe = (getenv("PLOT_WITH_PIXEL_PE")==NULL)?0:atoi(getenv("PLOT_WITH_PIXEL_PE"));
   int with_sum_only = (getenv("PLOT_WITH_SUM_ONLY")==NULL)?0:atoi(getenv("PLOT_WITH_SUM_ONLY"));
   int without_reco = (getenv("PLOT_WITHOUT_RECO")==NULL)?0:atoi(getenv("PLOT_WITHOUT_RECO"));
   double range;
   double hex_dx[6] = { 1.155, 0.577, -0.577, -1.155, -0.577, 0.577 };
   double hex_dy[6] = { 0.0, 1.0, 1.0, 0., -1.0, -1.0 };
   double sqr_dx[4] = { 1.0, -1.0, -1.0, 1.0 };
   double sqr_dy[4] = { 1.0, 1.0, -1.0, -1.0 };
   double showval[] = { 0., 0.2, 0.5, 1.0, 2.0, 5.0, 10., 20., 50. };
   // double sum = 0.;
   double sumsel = 0.;
   int npsum = 0;
   char primary[40];
   int ns = 0;
   int has_triggered[H_MAX_PIX], is_in_image[H_MAX_PIX];
   double r, rs = 0.;
   int last_time = -1, itime;
   double rot_deg = hsdata->camera_set[itel].cam_rot * (180./M_PI);
   int has_data = 0;
   int npe_true_total = 0, npe_true_sel = 0;
   int show_true_pe = (getenv("SHOW_TRUE_PE")==NULL)?(amp_tm==3?2:0):atoi(getenv("SHOW_TRUE_PE"));
# ifdef CTA
   const char *default_title = "CTA Telescope Simulation";
# else
   const char *default_title = "H.E.S.S. Telescope Simulation";
# endif
   const char *title = (getenv("PLOT_WITH_TITLE")==NULL) ? default_title : getenv("PLOT_WITH_TITLE");
   int its = (npix < 1000) ? 125 : ( (npix < 3000) ? 100 : ((npix < 5000) ? 80 : (( npix < 10000) ? 60 : 50)));

   if ( strcmp(image_fname,"none") == 0 ||
      	strcmp(image_fname,"/dev/null") == 0 )
      return;

   if ( itel < 0 || itel >= H_MAX_TEL )
      return;
   if ( itel_pe < 0 || itel_pe >= H_MAX_TEL )
   {
      show_true_pe = with_npe = 0;
   }

   if ( show_true_pe )
      img_off = 0.;

   teldata = &hsdata->event.teldata[itel];
   if ( !teldata->known )
      return;
   if ( teldata->raw == NULL && teldata->pixcal == NULL )
      return;
   if ( teldata->raw != NULL && teldata->raw->known )
   {
      has_data = 1;
      if ( teldata->raw->num_samples > 1 && amp_tm == 0 && !with_sum_only )
         last_time = teldata->raw->num_samples-1;
   }
   else if ( teldata->pixcal != NULL && teldata->pixcal->known )
      has_data = 2;
   if ( !has_data )
      return;

   if ( itel_pe >= 0 && itel_pe < H_MAX_TEL )
      npe_true_total = hsdata->mc_event.mc_pe_list[itel_pe].npe;
   if ( npe_true_total < 0 )
      show_true_pe = 0;
   if ( show_true_pe == 2 )
      last_time = -1;

   for ( i=0; i<npix; i++ )
   {
      has_triggered[i] = is_in_image[i] = 0;
      r = sqrt(x[i]*x[i] + y[i]*y[i]);
      rs += r;
   }
   for ( i=0; i<teldata->trigger_pixels.pixels; i++ )
      has_triggered[teldata->trigger_pixels.pixel_list[i]] = 1;
   for ( i=0; i<teldata->image_pixels.pixels; i++ )
      is_in_image[teldata->image_pixels.pixel_list[i]] = 1;

   body_diameter = 4.4*rs/(npix+0.1);
   scale *= 170./body_diameter;

   range = (img_range >= 1. ? img_range : 1.);
   gamma_coeff = (img_gamma != 0. ? img_gamma : 0.65);

   if ( image_fname[0] == '+' )
      psfile = fileopen(image_fname+1,"a");
   else if ( ps_num_page == 0 )
      psfile = fileopen(image_fname,"w");
   else
      psfile = fileopen(image_fname,"a");
   if ( psfile == NULL )
   {
      perror(image_fname);
      return;
   }

   if ( ps_num_page == 0 )
   {
      pid_t p = getpid();
      char prog[1024], proc[128];
      FILE *fproc;
      char date[60];
      time_t t = time(NULL);
      struct tm *ttm = localtime(&t);
      const char *bsnm = strrchr(image_fname,'/');
      if ( bsnm == NULL )
         bsnm = image_fname;
      else
         bsnm++;
      
      fputs(ps_head1a,psfile);
      fprintf(psfile," (%s)",bsnm);
      fputs(ps_head1b,psfile);
      sprintf(proc,"/proc/%d/cmdline",p);
      strcpy(prog,"program");
      if ( (fproc = fopen(proc,"r")) != NULL )
      {
         fgets(prog,sizeof(prog)-1,fproc);
         fclose(fproc);
      }
      fprintf(psfile," %s", prog);
      if ( getlogin() != NULL )
      {
      	 char host[60];
      	 fprintf(psfile," run by %s",getlogin());
	 if ( gethostname(host,sizeof(host)-1) == 0 )
	    fprintf(psfile," on %s",host);
      }
      fprintf(psfile,"\n");
      
      strftime(date,sizeof(date)-1,"%Y-%m-%d %T %Z",ttm);
      fprintf(psfile,"%%%%CreationDate: %s\n",date);
      fputs(ps_head2,psfile);
      fputs(ps_head3,psfile);
   }

   for ( itime=-1; itime<=last_time; itime++ )
   {
      ps_num_page++;
      fputs(ps_begin_page1,psfile);
      fprintf(psfile,"%s%d.%d:%d %d\n",
         1?"+":"-",event/100,event%100,tel,ps_num_page);
      fputs(ps_begin_page2,psfile);

   //   for (i=0; i<el->num_trigger_groups; i++)
   //      if ( el->trg_group[i].min_pixels < mintrg )
   //      	 mintrg = el->trg_group[i].min_pixels;

      for (i=0; i<1; i++)
      {
         // Only works if all pixels are of the same type.
         double hs = hsdata->camera_set[itel].size[0] / 2.;
         int pixel_shape = guessed_pixel_shape_type(&hsdata->camera_set[itel],itel);
         if ( pixel_shape < 0 )
            pixel_shape = 0;

         switch ( pixel_shape )
         {
      	    int idx, idy, ndx, ndy;
	    case 0:
	       if ( hs <= 0. )
      	          continue;
	       fprintf(psfile,"%% Pixel type %d has shape %d\n",
	          i, pixel_shape);
               ndx = Nint(hs*scale);
               fprintf(psfile,
                  "/pt%d { %d 0 rm currentpoint exch %d add exch %d 0 360 arc cp gs} def\n",
                  i, ndx, -ndx, ndx);
	       break;
      	    case 1:
	       ndx = idx = Nint(hs*scale*hex_dx[0]);
	       ndy = idy = Nint(hs*scale*hex_dy[0]);
	       fprintf(psfile,"/pt%d { %d %d rm", i, idx, idy);
	       for ( j=1; j<6; j++ )
	       {
	          idx = Nint(hs*scale*hex_dx[j]) - ndx;
	          idy = Nint(hs*scale*hex_dy[j]) - ndy;
	          ndx += idx;
	          ndy += idy;
	          fprintf(psfile," %d %d rl",idx,idy);
	       }
	       fprintf(psfile," cp gs } def\n");
	       break;
      	    case 2: /* square */
	       ndx = idx = Nint(hs*scale*sqr_dx[0]);
	       ndy = idy = Nint(hs*scale*sqr_dy[0]);
	       fprintf(psfile,"/pt%d { %d %d rm", i, idx, idy);
	       for ( j=1; j<4; j++ )
	       {
	          idx = Nint(hs*scale*sqr_dx[j]) - ndx;
	          idy = Nint(hs*scale*sqr_dy[j]) - ndy;
	          ndx += idx;
	          ndy += idy;
	          fprintf(psfile," %d %d rl",idx,idy);
	       }
	       fprintf(psfile," cp gs } def\n");
	       break;
      	    case 3:
	       ndx = idx = Nint(hs*scale*hex_dy[0]);
	       ndy = idy = Nint(hs*scale*hex_dx[0]);
	       fprintf(psfile,"/pt%d { %d %d rm", i, idx, idy);
	       for ( j=1; j<6; j++ )
	       {
	          idx = Nint(hs*scale*hex_dy[j]) - ndx;
	          idy = Nint(hs*scale*hex_dx[j]) - ndy;
	          ndx += idx;
	          ndy += idy;
	          fprintf(psfile," %d %d rl",idx,idy);
	       }
	       fprintf(psfile," cp gs } def\n");
	       break;
	    default:
	       if ( hs <= 0. )
      	          continue;
	       fprintf(psfile,"%% Pixel type %d has shape %d\n",
	          i, pixel_shape);
         }
      }
      fprintf(psfile,"/pxe { srgb 1.00 shd ef gr gs black s gr } def\n");

      // sum = hsdata->mc_event.mc_pesum.num_pe[itel];
      if ( hsdata->event.teldata[itel].img != NULL && 
           hsdata->event.teldata[itel].num_image_sets > 0 &&
           hsdata->event.teldata[itel].img->known &&
      	   hsdata->event.teldata[itel].img->amplitude > 0. &&
           !without_reco )
      {
         sumsel = hsdata->event.teldata[itel].img->amplitude;
         npsum  = hsdata->event.teldata[itel].img->pixels;
      }
      else
      {
         sumsel = 0.;
         npsum = 0;
      }
      if ( hsdata->event.teldata[itel].raw != NULL )
      {
   #ifdef DEBUG_NPE
      printf("\nNumber of samples: %d (%d)\n",
          hsdata->event.teldata[itel].raw->num_samples,
          hsdata->tel_moni[itel].num_ped_slices);
      printf("Zero suppression mode: %d, data reduction mode: %d\n",
          hsdata->event.teldata[itel].raw->zero_sup_mode,
          hsdata->event.teldata[itel].raw->data_red_mode);
      printf("Amplitude method: %d\n", amp_tm);
   #endif

       /* Rotate all pixels and pixel-related symbols by given camera rotation */
       fprintf(psfile,"gs\n10800 16200 tr\n%f rot\n-10800 -16200 tr\n", rot_deg);

       ns = 0;
       for (i=0; i<npix; i++)
       {
         /* Calibration the standard way. */
         double npe = calibrate_pixel_amplitude(hsdata,itel,i,amp_tm,itime,clip_amp);
         if ( itime >= 0 )
            npe *= 5.; // FIXME: need proper scaling to p.e. amplitude.
         else if ( show_true_pe )
            npe = hsdata->mc_event.mc_pe_list[itel_pe].pe_count[i];
         // double hs;
         /* Exactly zero amplitude indicates (most likely) we didn't know any amplitude. */
         int significant = (npe != 0.);
         if ( significant )
            npe += img_off; /* Offset used in color scale */

         if ( significant )
            ns++;

   #ifdef DEBUG_NPE
         if ( lg_known && !hg_known )
            printf("Alarm 1 in hesscam_ps_plot!\n");
         if ( significant && !hg_known )
            printf("Alarm 2 in hesscam_ps_plot!\n");
         if ( !significant && hg_known )
            printf("Alarm 3 in hesscam_ps_plot!\n");

         if ( getenv("SHOW_NPE") != 0 )
         { 
            double npe_thresh = atof(getenv("SHOW_NPE"));

            double hg2lg = (sig_lg != 0) ? sig_hg/sig_lg : 999.;

            double raw_hg = hsdata->event.teldata[itel].raw->adc_sum[HI_GAIN][i];
            double ped_hg = hsdata->tel_moni[itel].pedestal[HI_GAIN][i];
            double raw_lg = hsdata->event.teldata[itel].raw->adc_sum[LO_GAIN][i];
            double ped_lg = hsdata->tel_moni[itel].pedestal[LO_GAIN][i];

            if ( npe_hg > npe_thresh )
               printf("%7.2f %7.3f %7.3f %7.3f %7.3f %d %d  %7.1f %7.2f  %7.1f %7.2f\n", 
                  sig_hg, sig_lg, hg2lg, npe_hg, npe_lg, hg_known, lg_known,
	          raw_hg, ped_hg, raw_lg, ped_lg);
         }
   #endif

         /* Assume all pixels to be of the same type (and thus size) */
         // hs = hsdata->camera_set[itel].size[0] / 2.;

         xc = 10800. + scale*y[i];
         yc = 16200. - scale*x[i];

         fprintf(psfile,"n %ld %ld m pt%d ",
      	    Nint(xc),Nint(yc),/*pixtype*/ 0);

         print_pix_col(npe/range,psfile,gamma_coeff);
         if ( has_triggered[i] )
      	    fprintf(psfile,"n %ld %ld m gcr\n", Nint(xc), Nint(yc));
         if ( is_in_image[i] )
         {
      	    fprintf(psfile,"n %ld %ld m yxcr\n", Nint(xc), Nint(yc));
            npe_true_sel += hsdata->mc_event.mc_pe_list[itel_pe].pe_count[i];
         }

         if ( with_id )
         {
            /* Print the pixel ID on top of each pixel in very small font size. */
            fprintf(psfile,"txt%d %ld %ld mtxt (%d) tblack\n", its, Nint(xc), Nint(yc), i);
         }
         else if ( with_amp && !show_true_pe )
         {
           /* Print the pixel amplitude on top of each pixel in very small font size. */
            fprintf(psfile,"txt%d %ld %ld mtxt (%3.1f) tblack\n", its, Nint(xc), Nint(yc), npe-img_off);
         }
         else if ( with_npe || (with_amp && show_true_pe) )
         {
           /* Print the pixel amplitude on top of each pixel in very small font size. */
            fprintf(psfile,"txt%d %ld %ld mtxt (%d) tblack\n", its, Nint(xc), Nint(yc), 
               hsdata->mc_event.mc_pe_list[itel_pe].pe_count[i]);
         }

       }
       /* End of camera rotation */
       fprintf(psfile,"gr\n");
      }
      // ngt = 0;
   //   for (i=0; i<el->num_trigger_groups; i++)
   //      if ( el->trg_group[i].group_triggered )
   //      	 ngt++;
#ifdef LANG_DE
      if ( itime < 0 )
       fprintf(psfile,
# ifdef CTA
         "txt10 10800 2200 mtxt (CTA-Teleskop-Simulation) tblue\n");
# else
         "txt10 10800 2200 mtxt (H.E.S.S. Teleskop-Simulation) tblue\n");
# endif
      else
       fprintf(psfile,
# ifdef CTA
         "txt10 10800 2200 mtxt (CTA-Teleskop-Simulation \\(Interval %d\\)) tblue\n", itime);
# else
         "txt10 10800 2200 mtxt (H.E.S.S. Teleskop-Simulation \\(Interval %d\\)) tblue\n", itime);
# endif
      fprintf(psfile,
         "txt4 10800 3400 mtxt (Nr. %d, Ereignis %d, System %d, Teleskop %d) tblack\n",
         run,event/100,event%100,tel);
      if ( teldata->trigger_pixels.pixels > 0 )
         fprintf(psfile,"txt4 10800 4600 mtxt\n (Anzahl angesprochener Pixel: %d of %d)\n tblack\n",
            teldata->trigger_pixels.pixels, npix);
      if ( teldata->image_pixels.pixels > 0 )
         fprintf(psfile,"txt4 10800 5200 mtxt\n (Anzahl Pixel nach S\\344uberung: %d)\n tblack\n",
            teldata->image_pixels.pixels);
      fprintf(psfile, "txt4 10800 5800 mtxt (Anzahl signifikanter Pixel: %d) tblack\n", ns);
      fprintf(psfile,
         "txt4 10800 6400 mtxt\n"
         " (Summe der Signale in %d selektierten Pixeln: %3.1f p.e.%s)\n"
         " tblack\n", npsum, sumsel, amp_tm==0?"":(amp_tm==1?" am globalen Peak":" an lokalen Peaks"));

      if ( type < 0 )
      {
         switch ( hsdata->mc_shower.primary_id )
         {
            case 1:
	       strcpy(primary,"Gamma"); break;
            case 2:
            case 3:
	       strcpy(primary,"Elektron"); break;
            case 5:
            case 6:
	       strcpy(primary,"Myon"); break;
            case 14:
	       strcpy(primary,"Proton"); break;
            default:
	       sprintf(primary,"Typ %d",hsdata->mc_shower.primary_id);
         }
         fprintf(psfile,
            "txt6 10800 27800 mtxt\n"
            " (Prim\\344rteilchen: %s von %5.3f TeV Energie in %1.0f m Abstand) tred\n",
            primary, hsdata->mc_shower.energy,
            dist2(hsdata->mc_event.xcore - hsdata->run_header.tel_pos[itel][0],
                  hsdata->mc_event.ycore - hsdata->run_header.tel_pos[itel][1])
            );
      }
      else if ( type == 0 )
         fprintf(psfile,
            "txt6 10800 27800 mtxt\n"
            " (Kalibrationsdaten: Deckel geschlossen) tred\n");
      else if ( type == 1 )
         fprintf(psfile,
            "txt6 10800 27800 mtxt\n"
            " (Kalibrationsdaten: Deckel offen) tred\n");
      else if ( type == 2 )
         fprintf(psfile,
            "txt6 10800 27800 mtxt\n"
            " (Kalibrationsdaten: Lid-LED \\(Deckel geschlossen\\)) tred\n");
      else if ( type >= 3 )
         fprintf(psfile,
            "txt6 10800 27800 mtxt\n"
            " (Kalibrationsdaten: Laser in Einstellung %d) tred\n",
            type-2);
#else
      if ( itime < 0 )
       fprintf(psfile,
         "txt10 10800 2200 mtxt (%s) tblue\n", title);
      else
       fprintf(psfile,
         "txt10 10800 2200 mtxt (%s \\(Sample %d\\)) tblue\n", title, itime);
      fprintf(psfile,
         "txt4 10800 3400 mtxt (Run %d, event %d, array %d, telescope %d) tblack\n",
         run,event/100,event%100,tel);
      if ( teldata->trigger_pixels.pixels > 0 )
         fprintf(psfile,"txt4 10800 4600 mtxt\n (Number of triggered pixels: %d of %d)\n tblack\n",
            teldata->trigger_pixels.pixels, npix);
      fprintf(psfile, "txt4 10800 5200 mtxt (Number of significant pixels: %d) tblack\n", ns);
      if ( teldata->image_pixels.pixels > 0 )
         fprintf(psfile, "txt4 10800 5800 mtxt\n (Number of pixels after cleaning: %d)\n tblack\n",
            teldata->image_pixels.pixels);
      if ( show_true_pe )
      {
         fprintf(psfile,
            "txt4 10800 6400 mtxt\n"
            " (Actual p.e. count in %d selected pixels: %d of %d)\n"
            " tblack\n", npsum, npe_true_sel, npe_true_total);
      }
      else if ( npe_true_total > 0 )
      {
         if ( itime < 0 )
         {
            fprintf(psfile,
               "txt4 10800 6400 mtxt\n"
               " (Sum of signals in %d selected pixels: %3.1f p.e.%s \\(from %d of %d true p.e.\\))\n"
               " tblack\n", npsum, sumsel, amp_tm==0?"":(amp_tm==1?" at global peak":" at local peaks"),
               npe_true_sel, npe_true_total);
         }
         else
         {
            fprintf(psfile,
               "txt4 10800 6400 mtxt\n"
               " (Sum of signals in %d selected pixels: %3.1f p.e.%s)\n"
               " tblack\n", npsum, sumsel, amp_tm==0?"":(amp_tm==1?" at global peak":" at local peaks"));
         }
      }
      else
      {
         fprintf(psfile,
            "txt4 10800 6400 mtxt\n"
            " (Sum of signals in %d selected pixels: %3.1f p.e.%s)\n"
            " tblack\n", npsum, sumsel, amp_tm==0?"":(amp_tm==1?" at global peak":" at local peaks"));
      }

      if ( type < 0 )
      {
         switch ( hsdata->mc_shower.primary_id )
         {
            case 0:
	       strcpy(primary,"gamma"); break;
            case -1:
	       strcpy(primary,"positron"); break;
            case 1:
	       strcpy(primary,"electron"); break;
            case -2:
            case 2:
	       strcpy(primary,"muon"); break;
            case 101:
	       strcpy(primary,"proton"); break;
            default:
	       sprintf(primary,"type %d",hsdata->mc_shower.primary_id);
         }
         fprintf(psfile,
            "txt6 10800 27800 mtxt\n"
            " (Primary: %s of %5.3f TeV energy at %1.0f m distance) tred\n",
            primary, hsdata->mc_shower.energy,
            dist2(hsdata->mc_event.xcore - hsdata->run_header.tel_pos[itel][0],
                  hsdata->mc_event.ycore - hsdata->run_header.tel_pos[itel][1])
            );
      }
      else if ( type == 0 )
         fprintf(psfile,
            "txt6 10800 27800 mtxt\n"
            " (Calibration data: lid closed) tred\n");
      else if ( type == 1 )
         fprintf(psfile,
            "txt6 10800 27800 mtxt\n"
            " (Calibration data: lid open) tred\n");
      else if ( type == 2 )
         fprintf(psfile,
            "txt6 10800 27800 mtxt\n"
            " (Calibration data: lid-LED \\(lid closed\\)) tred\n");
      else if ( type >= 3 )
         fprintf(psfile,
            "txt6 10800 27800 mtxt\n"
            " (Calibration data: Laser at amplitude level %d) tred\n",
            type-2);
#endif

      if ( itime < 0 )
      {
         n = sizeof(showval)/sizeof(showval[0]);
         for (i=0; i<n; i++)
         {
            double npe = range * showval[i] + img_off;
            // double hs;

            // hs = hsdata->camera_set[itel].size[0] / 2.;

            xc = 10800 + 1000.*(i-0.5*n);
            yc = 26000;

            fprintf(psfile,"n %ld %ld m pt%d ",
      	       Nint(xc),Nint(yc),0);

            print_pix_col(npe/range,psfile,gamma_coeff);

            fprintf(psfile, "txt5 %d %d mtxt (%1.0f) tblue\n",
      	       (int)(xc),(int)(yc+800.),range * showval[i]);
         }
         fprintf(psfile, "txt5 %d %d mtxt (p.e.) tblue\n",
            (int)(xc+1000.),(int)(yc+800.));
      }

      /* Plot Az/Alt coordinate arrows (with x/y arrows rotated, if necessary) */
      fprintf(psfile,alt_az_arrow, rot_deg);

      /* Plot symbol description */
#ifdef LANG_DE
      fprintf(psfile,
         "%% Symbole:\n"
         "n 1900 23800 m 5600 23800 l 5600 26600 l 1900 26600 l 1900 23800 l cp\n"
         "  gs 0.82 0.92 1.0 srgb 1.00 shd ef gr gs black s gr gr\n"
         "gs 20 slw\n"
         "  n 2300 24200 m gcr  txt2 2600 24250 m (Pixel \\374ber Schwelle) rtblack\n"
         "  n 2300 24600 m yxcr txt2 2600 24650 m (Pixel selektiert \\374ber Bild) rtblack\n"
         "  n 2300 25000 m rxcr txt2 2600 25050 m (simulierte Richtung) rtblack\n"
         "  n 2300 25400 150 0 360 arc cp gs 30 slw red s gr\n"
         "    txt2 2600 25450 m (rekonstrierte Richtung) rtblack\n"
         "  gs 2700 26000 translate 5.0 2.0 scale n -110 0 m 110 0 l cp gs yellow s gr\n"
         "    n 0 0 100 0 360 arc cp gs 10 slw yellow s gr gr\n"
         "    txt2 3500 25850 m (\\(1*/2*\\) Ellipse aus) rtblack\n"
         "    txt2 3500 26250 m (2. Momenten) rtblack\n"
         "gr\n");
#else
      fprintf(psfile,
         "%% Symbol description:\n"
         "n 1900 23800 m 5600 23800 l 5600 26600 l 1900 26600 l 1900 23800 l cp\n"
         "  gs 0.82 0.92 1.0 srgb 1.00 shd ef gr gs black s gr gr\n"
         "gs 20 slw\n"
         "  n 2300 24200 m gcr  txt2 2600 24250 m (pixel above threshold) rtblack\n"
         "  n 2300 24600 m yxcr txt2 2600 24650 m (pixel in selected image) rtblack\n"
         "  n 2300 25000 m rxcr txt2 2600 25050 m (simulated direction) rtblack\n"
         "  n 2300 25400 150 0 360 arc cp gs 30 slw red s gr\n"
         "    txt2 2600 25450 m (reconstructed direction) rtblack\n"
         "  gs 2700 26000 translate 5.0 2.0 scale n -110 0 m 110 0 l cp gs yellow s gr\n"
         "    n 0 0 100 0 360 arc cp gs 10 slw yellow s gr gr\n"
         "    txt2 3500 25850 m (second moments) rtblack\n"
         "    txt2 3500 26250 m (ellipse \\(*1/*2\\)) rtblack\n"
         "gr\n");
#endif

      if ( hsdata->event.teldata[itel].img != NULL && 
           hsdata->event.teldata[itel].num_image_sets > 0 &&
           hsdata->event.teldata[itel].img->known &&
      	   hsdata->event.teldata[itel].img->amplitude > 0. &&
           !without_reco )
      {
         double ang_to_len = flen;
         double xe = hsdata->event.teldata[itel].img->x * ang_to_len;
         double ye = hsdata->event.teldata[itel].img->y * ang_to_len;
         double angle = 90. + (180./M_PI)*hsdata->event.teldata[itel].img->phi;
         double xs = 2. * hsdata->event.teldata[itel].img->l * ang_to_len;
         double ys = 2. * hsdata->event.teldata[itel].img->w * ang_to_len;

         /* // Note: we do not rotate here: fprintf(psfile,"gs\n10800 16200 tr\n%f rot\n-10800 -16200 tr\n", rot_deg); */
         xc = 10800. + scale*ye;
         yc = 16200. - scale*xe;
         if ( angle < 0. )
      	    angle += 360.;
         else if ( angle > 360. )
      	    angle -= 360.;
         fprintf(psfile,
            "%% Second moments ellipse at one and twice the actual size:\n"
      	    "gs 2 slw %d %d translate %f rotate %f %f scale\n"
	    "  n -200 0 m 200 0 l cp gs yellow s gr\n"
      	    "  n 0 0 100 0 360 arc cp gs yellow s gr gr\n"
      	    "gs 4 slw %d %d translate %f rotate %f %f scale\n"
	    "  n -200 0 m 200 0 l cp gs yellow s gr\n"
      	    "  n 0 0 100 0 360 arc cp gs yellow s gr gr\n",
	    (int)xc, (int)yc, angle, 0.01*scale*xs, 0.01*scale*ys,
	    (int)xc, (int)yc, angle, 0.005*scale*xs, 0.005*scale*ys);
      }

      if ( hsdata->event.shower.known && !without_reco )
      {
         double axs=0., ays=0.;
         angles_to_offset(hsdata->event.shower.Az, 
	    hsdata->event.shower.Alt, 
	    hsdata->event.trackdata[itel].azimuth_raw,
	    hsdata->event.trackdata[itel].altitude_raw,
	    flen, &axs, &ays);

      	 xc = 10800. + scale*ays;
      	 yc = 16200. - scale*axs;

	 fprintf(psfile,
	    "%% Reconstructed shower direction:\n"
	    "n %d %d 150 0 360 arc cp gs 30 slw red s gr\n",
	    (int)xc, (int)yc);
      }

      {
         double axs=0., ays=0.;
         angles_to_offset(hsdata->mc_shower.azimuth, 
	    hsdata->mc_shower.altitude, 
	    hsdata->event.trackdata[itel].azimuth_raw,
	    hsdata->event.trackdata[itel].altitude_raw,
	    flen, &axs, &ays);

         xc = 10800. + scale*ays;
         yc = 16200. - scale*axs;

         fprintf(psfile,
	    "%% Simulated shower direction:\n"
	    "n %d %d m %d %d l cp gs 30 slw red s gr\n"
	    "n %d %d m %d %d l cp gs 30 slw red s gr\n",
	    (int)(xc-100.), (int)(yc-100.), (int)(xc+100.), (int)(yc+100.),
	    (int)(xc-100.), (int)(yc+100.), (int)(xc+100.), (int)(yc-100.));
      }

      fputs(ps_end_page,psfile);
      fputs(ps_trailer,psfile);
   }
   fileclose(psfile);
}

/* ------------------------- hesscam_type_sum_plot ----------------------- */

void hesscam_type_sum_plot(const char *image_fname, AllHessData *hsdata, int teltype);

void hesscam_type_sum_plot(const char *image_fname, AllHessData *hsdata, int teltype)
{
   int run = hsdata->run_header.run;
   int event = hsdata->event.central.glob_count;
   int itel, itel_pe=-1, ktel=-1, ntel = hsdata->run_header.ntel;
   int nteltype = 0, iteltype, itlist[H_MAX_TEL];
   double npe_reco[H_MAX_PIX], npe_true[H_MAX_PIX], sum_npe_reco=0., sum_npe_true=0.;
   int npix, ipix;

   for ( itel=0; itel<ntel; itel++ )
   {
      TelEvent *teldata;
      int itype = user_get_type(itel);
      if ( itype != teltype )
         continue;
      teldata = &hsdata->event.teldata[itel];
      if ( !teldata->known )
         continue;
      if ( teldata->raw == NULL && teldata->pixcal == NULL )
         continue;
      itlist[nteltype] = itel;
      nteltype++;
   }
   if ( nteltype < 2 )
      return;

   ktel=itlist[0];
   npix = hsdata->camera_set[ktel].num_pixels;
   for ( ipix=0; ipix<npix; ipix++ )
      npe_reco[ipix] = npe_true[ipix] = 0.;
   
   for ( iteltype=0; iteltype<nteltype; iteltype++ )
   {
      int i;
      itel = itlist[iteltype];
      itel_pe = itel;
      if ( itel < 0 || itel >= H_MAX_TEL )
         continue;
      TelEvent *teldata = &hsdata->event.teldata[itel];
      if ( hsdata->camera_set[itel].num_pixels != npix )
         continue;
      for (i=0; i<npix; i++)
      {
         double npe = hsdata->mc_event.mc_pe_list[itel_pe].pe_count[i];
         npe_true[i] += npe;
         sum_npe_true += npe;
      }
      if ( (teldata->raw != NULL && teldata->raw->known) ||
           (teldata->pixcal != NULL && teldata->pixcal->known) )
      {
         for (i=0; i<npix; i++)
         {
            /* Calibration the standard way. */
            double npe = calibrate_pixel_amplitude(hsdata,itel,i,0,-1,0.);
            npe_reco[i] += npe;
            sum_npe_reco += npe;
         }
      }
   }
   itel = ktel;
   // int show_true_pe = (getenv("SHOW_TRUE_PE")==NULL)?0:atoi(getenv("SHOW_TRUE_PE"));
   int show_true_pe = 1;
   if ( sum_npe_true == 0. )
      show_true_pe = 0;

   double flen = hsdata->camera_set[itel].flen;
   double *x = &hsdata->camera_set[itel].xpix[0];
   double *y = &hsdata->camera_set[itel].ypix[0];
   // TelEvent *teldata = &hsdata->event.teldata[itel];
   FILE *psfile;
   // int ngt;
   int i, j, n;
   double xc, yc;
   double scale = 150.;
   double body_diameter = 0.;
   double gamma_coeff = (getenv("GAMMA_COEFF")==NULL)?0.65:atof(getenv("GAMMA_COEFF"));
   double img_gamma = (getenv("GRAY_IMAGE")==NULL)?gamma_coeff:-0.85*gamma_coeff;
   double img_range = (getenv("IMAGE_RANGE")==NULL)?20.:atof(getenv("IMAGE_RANGE"));
   double img_off = (getenv("IMAGE_OFFSET")==NULL)?4.:atof(getenv("IMAGE_OFFSET"));
   int with_id = (getenv("PLOT_WITH_PIXEL_ID")==NULL)?0:atoi(getenv("PLOT_WITH_PIXEL_ID"));
   int with_amp = (getenv("PLOT_WITH_PIXEL_AMP")==NULL)?0:atoi(getenv("PLOT_WITH_PIXEL_AMP"));
   int with_npe = (getenv("PLOT_WITH_PIXEL_PE")==NULL)?0:atoi(getenv("PLOT_WITH_PIXEL_PE"));
   int without_reco = (getenv("PLOT_WITHOUT_RECO")==NULL)?0:atoi(getenv("PLOT_WITHOUT_RECO"));
   double range;
   double hex_dx[6] = { 1.155, 0.577, -0.577, -1.155, -0.577, 0.577 };
   double hex_dy[6] = { 0.0, 1.0, 1.0, 0., -1.0, -1.0 };
   double sqr_dx[4] = { 1.0, -1.0, -1.0, 1.0 };
   double sqr_dy[4] = { 1.0, 1.0, -1.0, -1.0 };
   double showval[] = { 0., 0.2, 0.5, 1.0, 2.0, 5.0, 10., 20., 50. };
   char primary[40];
   double r, rs = 0.;
   double rot_deg = hsdata->camera_set[itel].cam_rot * (180./M_PI);
# ifdef CTA
   const char *default_title = "CTA Telescope Simulation";
# else
   const char *default_title = "H.E.S.S. Telescope Simulation";
# endif
   const char *title = (getenv("PLOT_WITH_TITLE")==NULL) ? default_title : getenv("PLOT_WITH_TITLE");
   int its = (npix < 1000) ? 125 : ( (npix < 3000) ? 100 : ((npix < 5000) ? 80 : (( npix < 10000) ? 60 : 50)));

   if ( strcmp(image_fname,"none") == 0 ||
      	strcmp(image_fname,"/dev/null") == 0 )
      return;

   if ( itel < 0 || itel >= H_MAX_TEL )
      return;
   if ( itel_pe < 0 || itel_pe >= H_MAX_TEL )
   {
      show_true_pe = with_npe = 0;
   }

   for ( i=0; i<npix; i++ )
   {
      r = sqrt(x[i]*x[i] + y[i]*y[i]);
      rs += r;
   }

   body_diameter = 4.4*rs/(npix+0.1);
   scale *= 170./body_diameter;

   range = (img_range >= 1. ? img_range : 1.);
   gamma_coeff = (img_gamma != 0. ? img_gamma : 0.65);

   if ( image_fname[0] == '+' )
      psfile = fileopen(image_fname+1,"a");
   else if ( ps_num_page == 0 )
      psfile = fileopen(image_fname,"w");
   else
      psfile = fileopen(image_fname,"a");
   if ( psfile == NULL )
   {
      perror(image_fname);
      return;
   }

   if ( ps_num_page == 0 )
   {
      pid_t p = getpid();
      char prog[1024], proc[128];
      FILE *fproc;
      char date[60];
      time_t t = time(NULL);
      struct tm *ttm = localtime(&t);
      const char *bsnm = strrchr(image_fname,'/');

      if ( bsnm == NULL )
         bsnm = image_fname;
      else
         bsnm++;
      
      fputs(ps_head1a,psfile);
      fprintf(psfile," (%s)",bsnm);
      fputs(ps_head1b,psfile);
      sprintf(proc,"/proc/%d/cmdline",p);
      strcpy(prog,"program");
      if ( (fproc = fopen(proc,"r")) != NULL )
      {
         fgets(prog,sizeof(prog)-1,fproc);
         fclose(fproc);
      }
      fprintf(psfile," %s", prog);
      if ( getlogin() != NULL )
      {
      	 char host[60];
      	 fprintf(psfile," run by %s",getlogin());
	 if ( gethostname(host,sizeof(host)-1) == 0 )
	    fprintf(psfile," on %s",host);
      }
      fprintf(psfile,"\n");
      
      strftime(date,sizeof(date)-1,"%Y-%m-%d %T %Z",ttm);
      fprintf(psfile,"%%%%CreationDate: %s\n",date);
      fputs(ps_head2,psfile);
      fputs(ps_head3,psfile);
   }

   {
    int ixtrue;
    for ( ixtrue=(show_true_pe?1:0); ixtrue>=0; ixtrue-- )
    {
      ps_num_page++;
      fputs(ps_begin_page1,psfile);
      fprintf(psfile,"%s%d.%d:T%d %d\n",
         1?"+":"-",event/100,event%100,teltype,ps_num_page);
      fputs(ps_begin_page2,psfile);

      if ( ixtrue )
         img_off = 0.;

      for (i=0; i<1; i++)
      {
         // Only works if all pixels are of the same type.
         double hs = hsdata->camera_set[itel].size[0] / 2.;
         int pixel_shape = guessed_pixel_shape_type(&hsdata->camera_set[itel],itel);
         if ( pixel_shape < 0 )
            pixel_shape = 0;

         switch ( pixel_shape )
         {
      	    int idx, idy, ndx, ndy;
	    case 0:
	       if ( hs <= 0. )
      	          continue;
	       fprintf(psfile,"%% Pixel type %d has shape %d\n",
	          i, pixel_shape);
               ndx = Nint(hs*scale);
               fprintf(psfile,
                  "/pt%d { %d 0 rm currentpoint exch %d add exch %d 0 360 arc cp gs} def\n",
                  i, ndx, -ndx, ndx);
	       break;
      	    case 1:
	       ndx = idx = Nint(hs*scale*hex_dx[0]);
	       ndy = idy = Nint(hs*scale*hex_dy[0]);
	       fprintf(psfile,"/pt%d { %d %d rm", i, idx, idy);
	       for ( j=1; j<6; j++ )
	       {
	          idx = Nint(hs*scale*hex_dx[j]) - ndx;
	          idy = Nint(hs*scale*hex_dy[j]) - ndy;
	          ndx += idx;
	          ndy += idy;
	          fprintf(psfile," %d %d rl",idx,idy);
	       }
	       fprintf(psfile," cp gs } def\n");
	       break;
      	    case 2: /* square */
	       ndx = idx = Nint(hs*scale*sqr_dx[0]);
	       ndy = idy = Nint(hs*scale*sqr_dy[0]);
	       fprintf(psfile,"/pt%d { %d %d rm", i, idx, idy);
	       for ( j=1; j<4; j++ )
	       {
	          idx = Nint(hs*scale*sqr_dx[j]) - ndx;
	          idy = Nint(hs*scale*sqr_dy[j]) - ndy;
	          ndx += idx;
	          ndy += idy;
	          fprintf(psfile," %d %d rl",idx,idy);
	       }
	       fprintf(psfile," cp gs } def\n");
	       break;
      	    case 3:
	       ndx = idx = Nint(hs*scale*hex_dy[0]);
	       ndy = idy = Nint(hs*scale*hex_dx[0]);
	       fprintf(psfile,"/pt%d { %d %d rm", i, idx, idy);
	       for ( j=1; j<6; j++ )
	       {
	          idx = Nint(hs*scale*hex_dy[j]) - ndx;
	          idy = Nint(hs*scale*hex_dx[j]) - ndy;
	          ndx += idx;
	          ndy += idy;
	          fprintf(psfile," %d %d rl",idx,idy);
	       }
	       fprintf(psfile," cp gs } def\n");
	       break;
	    default:
	       if ( hs <= 0. )
      	          continue;
	       fprintf(psfile,"%% Pixel type %d has shape %d\n",
	          i, pixel_shape);
         }
      }
      fprintf(psfile,"/pxe { srgb 1.00 shd ef gr gs black s gr } def\n");

       /* Rotate all pixels and pixel-related symbols by given camera rotation */
       fprintf(psfile,"gs\n10800 16200 tr\n%f rot\n-10800 -16200 tr\n", rot_deg);

       for (i=0; i<npix; i++)
       {
         /* Calibration the standard way. */
         double npe = (ixtrue ? npe_true[i] : npe_reco[i] ) + img_off;

         xc = 10800. + scale*y[i];
         yc = 16200. - scale*x[i];

         fprintf(psfile,"n %ld %ld m pt%d ",
      	    Nint(xc),Nint(yc),/*pixtype*/ 0);

         print_pix_col(npe/range,psfile,gamma_coeff);


         if ( with_id )
         {
            /* Print the pixel ID on top of each pixel in very small font size. */
            fprintf(psfile,"txt%d %ld %ld mtxt (%d) tblack\n", its, Nint(xc), Nint(yc), i);
         }
         else if ( with_amp>1 && !ixtrue )
         {
           /* Print the pixel amplitude on top of each pixel in very small font size. */
            fprintf(psfile,"txt%d %ld %ld mtxt (%3.1f) tblack\n", its, Nint(xc), Nint(yc), 
               npe_reco[i]);
         }
         else if ( with_npe>1 || (with_amp>1 && show_true_pe) )
         {
           /* Print the pixel amplitude on top of each pixel in very small font size. */
            fprintf(psfile,"txt%d %ld %ld mtxt (%1.0f) tblack\n", its, Nint(xc), Nint(yc), 
               npe_true[i]);
         }

       }
       /* End of camera rotation */
       fprintf(psfile,"gr\n");

#ifdef LANG_DE
       fprintf(psfile,
# ifdef CTA
         "txt10 10800 2200 mtxt (CTA-Teleskop-Simulation) tblue\n");
# else
         "txt10 10800 2200 mtxt (H.E.S.S. Teleskop-Simulation) tblue\n");
# endif
      fprintf(psfile,
         "txt4 10800 3400 mtxt (Nr. %d, Ereignis %d, System %d, Summe von %d Teleskopen vom Typ %d) tblack\n",
         run,event/100,event%100, nteltype, teltype);

         switch ( hsdata->mc_shower.primary_id )
         {
            case 1:
	       strcpy(primary,"Gamma"); break;
            case 2:
            case 3:
	       strcpy(primary,"Elektron"); break;
            case 5:
            case 6:
	       strcpy(primary,"Myon"); break;
            case 14:
	       strcpy(primary,"Proton"); break;
            default:
	       sprintf(primary,"Typ %d",hsdata->mc_shower.primary_id);
         }
         fprintf(psfile,
            "txt6 10800 27800 mtxt\n"
            " (Prim\\344rteilchen: %s von %5.3f TeV Energie) tred\n",
            primary, hsdata->mc_shower.energy);
#else
       fprintf(psfile,
         "txt10 10800 2200 mtxt (%s) tblue\n", title);
      fprintf(psfile,
         "txt4 10800 3400 mtxt (Run %d, event %d, array %d, sum of %d telescopes of type %d) tblack\n",
         run,event/100,event%100,nteltype,teltype);

         switch ( hsdata->mc_shower.primary_id )
         {
            case 0:
	       strcpy(primary,"gamma"); break;
            case -1:
	       strcpy(primary,"positron"); break;
            case 1:
	       strcpy(primary,"electron"); break;
            case -2:
            case 2:
	       strcpy(primary,"muon"); break;
            case 101:
	       strcpy(primary,"proton"); break;
            default:
	       sprintf(primary,"type %d",hsdata->mc_shower.primary_id);
         }
         fprintf(psfile,
            "txt6 10800 27800 mtxt\n"
            " (Primary: %s of %5.3f TeV energy) tred\n",
            primary, hsdata->mc_shower.energy);
#endif

         n = sizeof(showval)/sizeof(showval[0]);
         for (i=0; i<n; i++)
         {
            double npe = range * showval[i] + img_off;
            // double hs;

            // hs = hsdata->camera_set[itel].size[0] / 2.;

            xc = 10800 + 1000.*(i-0.5*n);
            yc = 26000;

            fprintf(psfile,"n %ld %ld m pt%d ",
      	       Nint(xc),Nint(yc),0);

            print_pix_col(npe/range,psfile,gamma_coeff);

            fprintf(psfile, "txt5 %d %d mtxt (%1.0f) tblue\n",
      	       (int)(xc),(int)(yc+800.),range * showval[i]);
         }
         fprintf(psfile, "txt5 %d %d mtxt (p.e.) tblue\n",
            (int)(xc+1000.),(int)(yc+800.));

      /* Plot Az/Alt coordinate arrows (with x/y arrows rotated, if necessary) */
      fprintf(psfile,alt_az_arrow, rot_deg);

      if ( hsdata->event.shower.known && !without_reco )
      {
         double axs=0., ays=0.;
         angles_to_offset(hsdata->event.shower.Az, 
	    hsdata->event.shower.Alt, 
	    hsdata->event.trackdata[itel].azimuth_raw,
	    hsdata->event.trackdata[itel].altitude_raw,
	    flen, &axs, &ays);

      	 xc = 10800. + scale*ays;
      	 yc = 16200. - scale*axs;

	 fprintf(psfile,
	    "%% Reconstructed shower direction:\n"
	    "n %d %d 150 0 360 arc cp gs 30 slw red s gr\n",
	    (int)xc, (int)yc);
      }

      {
         double axs=0., ays=0.;
         angles_to_offset(hsdata->mc_shower.azimuth, 
	    hsdata->mc_shower.altitude, 
	    hsdata->event.trackdata[itel].azimuth_raw,
	    hsdata->event.trackdata[itel].altitude_raw,
	    flen, &axs, &ays);

         xc = 10800. + scale*ays;
         yc = 16200. - scale*axs;

         fprintf(psfile,
	    "%% Simulated shower direction:\n"
	    "n %d %d m %d %d l cp gs 30 slw red s gr\n"
	    "n %d %d m %d %d l cp gs 30 slw red s gr\n",
	    (int)(xc-100.), (int)(yc-100.), (int)(xc+100.), (int)(yc+100.),
	    (int)(xc-100.), (int)(yc+100.), (int)(xc+100.), (int)(yc-100.));
      }

      fputs(ps_end_page,psfile);
      fputs(ps_trailer,psfile);
    }
   }
   fileclose(psfile);
}

/* --------------------------- find_neighbours ---------------------------- */

#define H_MAX_NB1 8
#define H_MAX_NB2 24
static int neighbours1[H_MAX_TEL][H_MAX_PIX][H_MAX_NB1]; 
//static int neighbours2[H_MAX_TEL][H_MAX_PIX][H_MAX_NB2];
static int nnb1[H_MAX_TEL][H_MAX_PIX];
//static int nnb2[H_MAX_TEL][H_MAX_PIX];
static int has_nblist[H_MAX_TEL];
static int px_shape_type[H_MAX_TEL];

/** Find the list of neighbours for each pixel. */

static int find_neighbours(CameraSettings *camset, int itel);

static int find_neighbours(CameraSettings *camset, int itel)
{
   int npix = camset->num_pixels;
   int i, j;
   int stat_st[6] = {0, 0, 0, 0, 0, 0 };
   double asum = 0., dsum = 0., aod2 = 0.;

   for (i=0; i<npix; i++)
   {
#if 0
      if ( pixel_disabled[itel][i] )
         continue;
#endif
      asum += camset->area[i];
      dsum += camset->size[i];
      for (j=0; j<H_MAX_NB1; j++)
         neighbours1[itel][i][j] = -1;
      nnb1[itel][i] = 0;
      for (j=0; j<i; j++)
      {
         double ds, dx, dy, d2, a;
         int ia;
#if 0
         if ( pixel_disabled[itel][j] )
            continue;
#endif
         ds = camset->size[i] + camset->size[j];
         dx = camset->xpix[i] - camset->xpix[j];
         dy = camset->ypix[i] - camset->ypix[j];
         d2 = dx*dx + dy*dy;
         if ( d2 < 0.5*ds*ds ) 
         {
            if ( nnb1[itel][i] < H_MAX_NB1 )
               neighbours1[itel][i][nnb1[itel][i]++] = j;
            if ( nnb1[itel][j] < H_MAX_NB1 )
               neighbours1[itel][j][nnb1[itel][j]++] = i;
            a  = (180./M_PI) * atan2(dy,dx);
            if ( a < -1. )
               a += 180.;
            ia = ((int) ((a+0.5)/5.)) * 5;
            if ( ia == 0 )
               stat_st[0]++;
            else if ( ia == 60 )
               stat_st[1]++;
            else if ( ia == 90 )
               stat_st[2]++;
            else if ( ia == 120 )
               stat_st[3]++;
            else if ( ia == 30 )
               stat_st[4]++;
            else if ( ia == 150 )
               stat_st[5]++;
         }
      }
   }
   has_nblist[itel] = 1;

   asum /= (((double) npix)+1e-10);
   dsum /= (((double) npix)+1e-10);
   aod2 = asum/(dsum*dsum);

   if ( stat_st[0] > 0 && stat_st[2] > 0 && 
        stat_st[1] == 0 && stat_st[3] == 0 )
   {
      px_shape_type[itel] = 2; /* probably square pixels */
      if ( aod2 < 0.99 || aod2 > 1.01 )
      {
         fprintf(stderr,
            "Pixel positions in telescope %d indicate square pixels but area/size^2 does not match.\n",
            camset->tel_id);
      }
   }
   else /* most likely hexagonal, orientation? */
   {
      if ( 4*stat_st[2] < (stat_st[1] + stat_st[3]) )
         px_shape_type[itel] = 1;
      else if ( stat_st[2] > 0 && stat_st[0] == 0 )
         px_shape_type[itel] = 3;
      else
      {
         px_shape_type[itel] = 0; /* Fall back to circular */
         if ( aod2 < 0.99*M_PI/4. || aod2 > 1.01*M_PI/4. )
         {
            fprintf(stderr,
               "Pixel positions in telescope %d indicate round pixels but area/size^2 does not match.\n",
               camset->tel_id);
         }
      }
      if ( px_shape_type[itel] != 0 )
      {
         if ( aod2 < 0.99*sqrt(3.)/2. || aod2 > 1.01*sqrt(3.)/2. )
         {
            if ( aod2 >= 0.99*M_PI/4. && aod2 <= 1.01*M_PI/4. )
            {
               px_shape_type[itel] = 0; // Round pixels on hexagonal pattern
            }
            else
               fprintf(stderr,
                  "Pixel positions in telescope %d indicate hexagonal pixels but area/size^2 does not match.\n",
                  camset->tel_id);
         }
      }
   }
#if 0
   fprintf(stderr,"Pixel shape type of telescope #%d seems to be %d (stat = %d, %d, %d, %d)\n",
      itel, px_shape_type[itel], stat_st[0], stat_st[1], stat_st[2], stat_st[3]);
#endif
   return 0;
}

/* ------------------------ guessed_pixel_shape_type --------------------- */

static int guessed_pixel_shape_type(CameraSettings *camset, int itel)
{
   if ( itel >= 0 && itel < H_MAX_TEL )
   {
      if ( has_nblist[itel] == 0 )
         find_neighbours(camset, itel);
      if ( has_nblist[itel] == 0 ) /* did not work out ? */
      {
         fprintf(stderr, "Cannot guess the pixel type of telescope #%d yet.\n",
            itel);
         return 0;
      }
      if ( camset != NULL && camset->pixel_shape[0] >= 0 && 
           px_shape_type[itel] != camset->pixel_shape[0] )
      {
         fprintf(stderr,
            "Warning: configured pixel shape (%d) does not match neighbour geometry (%d)!\n",
            camset->pixel_shape[0], px_shape_type[itel]);
         return camset->pixel_shape[0];
      }
      return px_shape_type[itel];
   }
   else
      return 0;
}

