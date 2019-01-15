/* ============================================================================

Copyright (C) 2014, 2016  Konrad Bernloehr

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

/** @file split_hessio.c
 *  @short Rip out data for each telescope into individual files.
 *
 *
@verbatim
Syntax: split_hessio [ options ] [ - | input_fname ... ]
Options:
   -x              (Extract TelescopeEvent data from Event.)
   -X              (Extract TelescopeEvent raw data (samples or sum).)
   -i|--ignore     (Ignore unknown data block types.)
   -q|--quiet      (More quiet on standard output.)
   -v|--verbose    (More verbose on standard output.)
   --max-events n  (Skip remaining data after so many triggered events.)
   --pure-raw      (Discard any sub-items of TelescopeEvent which are not raw data.)
   --clean-history (Drop previous history data blocks)
   --output-path d (Create output files in given directory instead of current.)
   --only-telescope[s] (Only data for the given telescopes IDs is written.)
   --not-telescope[s]  (No data for the given telescopes IDs is written.)
@endverbatim
 *
 *  @author  Konrad Bernloehr
 *  @date    @verbatim CVS $Date: 2017/05/16 12:31:52 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.8 $ @endverbatim
 */

/** @defgroup split_hessio_c The split_hessio program */
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

/** Show program syntax */

static void syntax (char *program);

static void syntax (char *program)
{
   printf("Syntax: %s [ options ] [ - | input_fname ... ]\n",program);
   printf("Options:\n");
   printf("   -x              (Extract TelescopeEvent data from Event.)\n");
   printf("   -X              (Extract TelescopeEvent raw data (samples or sum).)\n");
   printf("   -i|--ignore     (Ignore unknown data block types.)\n");
   printf("   -q|--quiet      (More quiet on standard output.)\n");
   printf("   -v|--verbose    (More verbose on standard output.)\n");
   printf("   --max-events n  (Skip remaining data after so many triggered events.)\n");
   printf("   --pure-raw      (Discard any sub-items of TelescopeEvent which are not raw data.)\n");
   printf("   --clean-history (Drop previous history data blocks)\n");
   printf("   --keep-compression (Keep same compression scheme as input file, .gz etc.)\n");
   printf("   --output-path d (Create output files in given directory instead of current.)\n");
   printf("   --only-telescope[s] (Only data for the given telescopes IDs is written.)\n");
   printf("   --not-telescope[s]  (No data for the given telescopes IDs is written.)\n");

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
   int verbose = 0, ignore = 0, quiet = 0;
   int ntel_trg = 0, min_tel_trg = 0;
   int nev = 0, ntrg = 0;
   char *program = argv[0];
   size_t events = 0, max_events = 0;
   int iarg;
   static FILE *f_teldata[H_MAX_TEL], *f_telaux[H_MAX_TEL], *f_telaux2[H_MAX_TEL];
   static FILE *f_mc = NULL, *f_head = NULL, *f_tail = NULL, *f_central = NULL;
   static int tel_used[H_MAX_TEL];
   size_t num_only = 0, num_not = 0;
   int only_telescope[H_MAX_TEL], not_telescope[H_MAX_TEL];
   char input_fname_current[1024], tmp_fname[2600], tmp_fname_head[2200], tmp_fname_tail[300];
   char output_path[1024];
   int active_tel[H_MAX_TEL];
   int num_active;
   int extract_televent = 0;
   int pure_raw = 0, clean_history = 0;
   int with_compression = 0;

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

   /* In case we write a DST file, we want to record how we did it. */
   push_command_history(argc,argv);

   if ( argc < 2 )
      input_fname = "iact.out";

   if ( (iobuf = allocate_io_buffer(5000000L)) == NULL )
   {
      Error("Cannot allocate I/O buffer");
      exit(1);
   }
   iobuf->max_length = 100000000L;
   output_path[0] = '\0';

   /* Command line options */
   while ( argc > 1 )
   {
      if ( strcmp(argv[1],"-i") == 0 || strcmp(argv[1],"--ignore") == 0 )
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
      else if ( strcmp(argv[1],"-x") == 0 )
      {
       	 extract_televent = 1;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"-X") == 0 )
      {
       	 extract_televent = 2;
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
      else if ( strcmp(argv[1],"--output-path") == 0 && argc > 2 )
      {
         size_t l = strlen(argv[2]);
       	 strncpy(output_path,argv[2],sizeof(output_path)-2);
         if ( l > 0 && argv[2][l-1] != '/' )
            strcat(output_path,"/");
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
      else if ( strcmp(argv[1],"--pure-raw") == 0 )
      {
       	 pure_raw = 1;
	 argc--;
	 argv++;
      }
      else if ( strcmp(argv[1],"--clean-history") == 0 )
      {
       	 clean_history = 1;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--keep-compression") == 0 )
      {
       	 with_compression = 1;
	 argc--;
	 argv++;
	 continue;
      }
      else if ( strcmp(argv[1],"--help") == 0 )
      {
        printf("\nhessio_extract_tel: Extract data to telescope-specific files.\n\n");
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
    strncpy(input_fname_current, input_fname, sizeof(input_fname_current)-1);
    tmp_fname_tail[0] = '\0';
    {
       const char *s;
       char *t;
       if ( (s = strrchr(input_fname,'/')) != NULL )
          strncpy(tmp_fname_head, s+1, sizeof(tmp_fname_head));
       else
          strncpy(tmp_fname_head, input_fname, sizeof(tmp_fname_head));
       if ( (t = strstr(tmp_fname_head,".sim")) != NULL )
       {
          strncpy(tmp_fname_tail,t,sizeof(tmp_fname_tail)-1);
          *t = '\0';
       }
    }
    if ( !with_compression )
    {
       char *u = strchr(tmp_fname_tail+4,'.');
       if ( u != NULL )
       	  *u = '\0';
    }
    input_fname = NULL;
    if ( strlen(output_path)+strlen(tmp_fname_head)+30+strlen(tmp_fname_tail) >= sizeof(tmp_fname) )
    {
       fprintf(stderr,"%s: File name (including paths) gets too long.\n", input_fname);
       continue;
    }
    
    if ( f_head != NULL )
    {
       fileclose(f_head);
       f_head = NULL;
    }
    if ( f_mc != NULL )
    {
       fileclose(f_mc);
       f_mc = NULL;
    }
    if ( f_central != NULL )
    {
       fileclose(f_central);
       f_central = NULL;
    }
    if ( f_tail != NULL )
    {
       fileclose(f_tail);
       f_tail = NULL;
    }
    
    for ( itel=0; itel<H_MAX_TEL; itel++ )
    {
       if ( f_telaux[itel] != NULL )
       {
          fileclose(f_telaux[itel]);
          f_telaux[itel] = NULL;
       }
       if ( f_telaux2[itel] != NULL )
       {
          fileclose(f_telaux2[itel]);
          f_telaux2[itel] = NULL;
       }
       if ( f_teldata[itel] != NULL )
       {
          fileclose(f_teldata[itel]);
          f_teldata[itel] = NULL;
       }
    }

    if ( f_head == NULL )
    {
       strcpy(tmp_fname,output_path);
       strcat(tmp_fname,tmp_fname_head);
       strcat(tmp_fname,".headers");
       strcat(tmp_fname,tmp_fname_tail);
       if ( (f_head = fileopen(tmp_fname,"w")) == NULL )
       {
          perror(tmp_fname);
          exit(1);
       }
       iobuf->output_file = f_head;
       write_history(0,iobuf);  /* Save the command line history */
    }
    if ( f_tail == NULL )
    {
       strcpy(tmp_fname,output_path);
       strcat(tmp_fname,tmp_fname_head);
       strcat(tmp_fname,".tail");
       strcat(tmp_fname,tmp_fname_tail);
       if ( (f_tail = fileopen(tmp_fname,"w")) == NULL )
       {
          perror(tmp_fname);
          exit(1);
       }
    }
    if ( f_mc == NULL )
    {
       strcpy(tmp_fname,output_path);
       strcat(tmp_fname,tmp_fname_head);
       strcat(tmp_fname,".mc");
       strcat(tmp_fname,tmp_fname_tail);
       if ( (f_mc = fileopen(tmp_fname,"w")) == NULL )
       {
          perror(tmp_fname);
          exit(1);
       }
    }
    if ( f_central == NULL )
    {
       strcpy(tmp_fname,output_path);
       strcat(tmp_fname,tmp_fname_head);
       strcat(tmp_fname,".central");
       strcat(tmp_fname,tmp_fname_tail);
       if ( (f_central = fileopen(tmp_fname,"w")) == NULL )
       {
          perror(tmp_fname);
          exit(1);
       }
    }

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
            nev = ntrg = 0;
 
            if ( (iobuf->output_file = f_head) != NULL )
            {
               write_io_block(iobuf);   /* Re-write the run header as is */
            }

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

            for (itel=0; itel<hsdata->run_header.ntel; itel++)
            {
               size_t j;
               tel_id = hsdata->run_header.tel_id[itel];
               if ( num_only == 0 )
                  tel_used[itel] = 1;
               else
               {
                  tel_used[itel] = 0;
                  for ( j=0; j<num_only; j++ )
                  {
                     if ( only_telescope[j] == tel_id )
                        tel_used[itel] = 1;
                  }
               }
               for ( j=0; j<num_not; j++ )
                  if ( not_telescope[j] == tel_id )
                     tel_used[itel] = 0;
               
               if ( tel_used[itel] )
               {
                  char exn[100];
                  snprintf(exn,sizeof(exn),"%d",tel_id);
                  if ( f_teldata[itel] == NULL )
                  {
                     strcpy(tmp_fname,output_path);
                     strcat(tmp_fname,tmp_fname_head);
                     strcat(tmp_fname,".teldata-");
                     strcat(tmp_fname,exn);
                     strcat(tmp_fname,tmp_fname_tail);
                     if ( (f_teldata[itel] = fileopen(tmp_fname,"w")) == NULL )
                     {
                        perror(tmp_fname);
                        exit(1);
                     }
                  }
                  if ( f_telaux[itel] == NULL )
                  {
                     strcpy(tmp_fname,output_path);
                     strcat(tmp_fname,tmp_fname_head);
                     strcat(tmp_fname,".telaux-");
                     strcat(tmp_fname,exn);
                     strcat(tmp_fname,tmp_fname_tail);
                     if ( (f_telaux[itel] = fileopen(tmp_fname,"w")) == NULL )
                     {
                        perror(tmp_fname);
                        exit(1);
                     }
                  }
                  if ( f_telaux2[itel] == NULL )
                  {
                     strcpy(tmp_fname,output_path);
                     strcat(tmp_fname,tmp_fname_head);
                     strcat(tmp_fname,".telaux2-");
                     strcat(tmp_fname,exn);
                     strcat(tmp_fname,tmp_fname_tail);
                     if ( (f_telaux2[itel] = fileopen(tmp_fname,"w")) == NULL )
                     {
                        perror(tmp_fname);
                        exit(1);
                     }
                  }
               }

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
            if ( (iobuf->output_file = f_mc) != NULL )
               write_io_block(iobuf);
            break;

         /* =================================================== */
	 case IO_TYPE_MC_INPUTCFG:
            if ( (iobuf->output_file = f_mc) != NULL )
               write_io_block(iobuf);
            break;

         /* =================================================== */
         case 70: /* How sim_hessarray was run and how it was configured. */
            if ( !clean_history && (iobuf->output_file = f_head) != NULL )
               write_io_block(iobuf);
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
            if ( (iobuf->output_file = f_telaux[itel]) != NULL )
               write_io_block(iobuf);
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
            if ( (iobuf->output_file = f_telaux[itel]) != NULL )
               write_io_block(iobuf);
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
            if ( (iobuf->output_file = f_telaux[itel]) != NULL )
               write_io_block(iobuf);
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
            if ( (iobuf->output_file = f_telaux2[itel]) != NULL )
               write_io_block(iobuf);
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
            if ( (iobuf->output_file = f_telaux[itel]) != NULL )
               write_io_block(iobuf);
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
            if ( (iobuf->output_file = f_telaux[itel]) != NULL )
               write_io_block(iobuf);
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
            if ( (iobuf->output_file = f_telaux[itel]) != NULL )
               write_io_block(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_EVENT:
            /* This is the only block type which we need to read and process */
            rc = read_hess_event(iobuf,&hsdata->event,-1);
	    if ( verbose || rc != 0 )
               printf("read_hess_event(), rc = %d\n",rc);
            events++;
            /* Count number of telescopes (still) present in data and triggered */
            ntel_trg = 0;
            num_active = 0;
            for (itel=0; itel<hsdata->run_header.ntel; itel++)
            {
               if ( hsdata->event.teldata[itel].known )
               {
                  active_tel[num_active++] = itel;
                  hsdata->event.teldata[itel].known = 0;
                  hsdata->event.trackdata[itel].raw_known = 0;
                  if ( pure_raw )
                  {
                     int j;
                     if ( hsdata->event.teldata[itel].pixtm != NULL )
                        hsdata->event.teldata[itel].pixtm->known = 0;
                     if ( hsdata->event.teldata[itel].pixcal != NULL )
                        hsdata->event.teldata[itel].pixcal->known = 0;
                     if ( hsdata->event.teldata[itel].img != NULL )
                     {
                        for ( j=0; j<hsdata->event.teldata[itel].num_image_sets; j++ )
                           hsdata->event.teldata[itel].img[j].known = 0;
                        hsdata->event.teldata[itel].num_image_sets = 0;
                     }
                     /* We keep the list of triggered pixels but discard image selected pixels */
                     hsdata->event.teldata[itel].image_pixels.pixels = 0;
                  }
                  /* If non-triggered telescopes record data (like HEGRA),
                     we may have to check the central trigger bit as well,
                     but ignore this for now. */
                  ntel_trg++;
               }
            }
            if ( pure_raw )
               hsdata->event.shower.known = 0;
	    if ( hsdata->event.shower.known )
	       hsdata->event.shower.num_trg = ntel_trg;
            if ( ntel_trg < min_tel_trg )
               continue;
            ntrg++;
            if ( (iobuf->output_file = f_central) != NULL )
            {
               rc = write_hess_centralevent(iobuf,&hsdata->event.central);
            }

            { int j;
               for ( j=0; j<num_active; j++ )
               {
                  itel = active_tel[j];
                  if ( tel_used[itel] )
                  {
                     /* We temporarily mark just this one telescope as having participated and with data */
                     hsdata->event.num_teldata = 1;
                     hsdata->event.teldata_list[0] = hsdata->event.teldata[itel].tel_id;
                     hsdata->event.teldata[itel].known = 1;
                     hsdata->event.trackdata[itel].raw_known = 1;

                     if ( (iobuf->output_file = f_teldata[itel]) != NULL )
                     {
                        if ( extract_televent == 1 )
                        {
                           /* In this mode we store the tracking data first and then the telescope event */
                           if ( hsdata->event.trackdata[itel].raw_known ||
                                hsdata->event.trackdata[itel].cor_known )
                              (void) write_hess_trackevent(iobuf,&hsdata->event.trackdata[itel]);
	                   rc = write_hess_televent(iobuf,&hsdata->event.teldata[itel],-1);
                        }
                        else if ( extract_televent == 2 )
                        {
                           /* In this mode the tracking and telescope event header get lost */
                           TelEvent *te = &hsdata->event.teldata[itel];
                           if ( te->raw != NULL && te->raw->known )
                           {
                              if ( te->readout_mode > 0 ) /* Write only samples */
                                 rc = write_hess_teladc_samples(iobuf,te->raw);
                              else /* or write only sums, never both */
                                 rc = write_hess_teladc_sums(iobuf,te->raw);
                           }
                        }
                        else
                        {
                           /* In this mode we write the full event data block, although with just one telescope */
                           rc = write_hess_event(iobuf,&hsdata->event,-1);
                        }
                     }
                     hsdata->event.teldata[itel].known = 0;
                     hsdata->event.trackdata[itel].raw_known = 0;
                  }
               }
            }

            break;

         /* =================================================== */
         case IO_TYPE_HESS_CALIBEVENT:
            if ( (iobuf->output_file = f_head) != NULL )
               write_io_block(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_MC_SHOWER:
            if ( (iobuf->output_file = f_mc) != NULL )
               write_io_block(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_MC_EVENT:
            if ( (iobuf->output_file = f_mc) != NULL )
               write_io_block(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_MC_TELARRAY:
            if ( (iobuf->output_file = f_mc) != NULL )
               write_io_block(iobuf);
            break;

         /* =================================================== */
         /* With extended output option activated, the particles
            arriving at ground level would be stored as seemingly
            stray photon bunch block. */
         case IO_TYPE_MC_PHOTONS:
            if ( (iobuf->output_file = f_mc) != NULL )
               write_io_block(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_MC_RUNH:
         case IO_TYPE_MC_EVTH:
         case IO_TYPE_MC_EVTE:
         case IO_TYPE_MC_RUNE:
            if ( (iobuf->output_file = f_mc) != NULL )
               write_io_block(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_MC_PE_SUM:
            if ( (iobuf->output_file = f_mc) != NULL )
               write_io_block(iobuf);
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
            if ( (iobuf->output_file = f_telaux2[itel]) != NULL )
               write_io_block(iobuf);
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
            if ( (iobuf->output_file = f_telaux2[itel]) != NULL )
               write_io_block(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_RUNSTAT:
            if ( (iobuf->output_file = f_tail) != NULL )
               write_io_block(iobuf);
            break;

         /* =================================================== */
         case IO_TYPE_HESS_MC_RUNSTAT:
            if ( (iobuf->output_file = f_mc) != NULL )
               write_io_block(iobuf);
            break;

         /* (End-of-job or DST) histograms */
         case 100:
            if ( (iobuf->output_file = f_tail) != NULL )
               write_io_block(iobuf);
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

    if ( nev > 0 )
      printf("%d of %d events triggered\n", ntrg, nev);

    if ( hsdata != NULL )
       hsdata->run_header.run = 0;
   }

   return 0;
}

/** @} */
