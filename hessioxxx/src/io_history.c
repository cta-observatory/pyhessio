/* ============================================================================

   Copyright (C) 1997, 2001, 2009, 2010  Konrad Bernloehr

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

/** @file io_history.c
 *  @short Record history of configuration settings/commands.
 *
 *  This code has not been adapted for multi-threading.
 *
 *  @author  Konrad Bernloehr 
 *  @date 1997 to 2010
 *  @date    @verbatim CVS $Date: 2014/02/20 11:40:42 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.8 $ @endverbatim
 */
/* ================================================================ */

#include "initial.h"
#include "io_basic.h"
#include "history.h"
#include "current.h"
#include "hconfig.h"

/** Use to build a linked list of configuration history. */
struct history_struct
{
   char *text;  /** Configuration test */
   time_t time; /** Time when the configuration was entered */
   struct history_struct *next; /* Next element */
};
typedef struct history_struct HSTRUCT;

static char *cmdline = NULL; /**< A copy of the program's command line. */
static time_t cmdtime;       /**< The time when the program was started. */
static HSTRUCT *configs = NULL; /**< Start of configuration history. */

static void listtime (time_t t,FILE * f);

/* ------------------------ push_command_history --------------------- */
/*
@@@push_command_history()
@@    Save the command line for later output in a history block.
@*
@*    Arguments:  argc  --  Number of command line arguments (incl. command)
@*                argv  --  Pointers to argument text strings
@*
@*    Return value:  0 (o.k.), -1 (invalid argument or no memory)
@*
 */

int push_command_history (int argc, char **argv)
{
   int i;
   size_t len;

   for (i=0, len=0; i<argc && argv[i] != (char *) NULL; i++)
      len += strlen(argv[i])+1;
   if ( len < 1 )
      return -1;

   if ( cmdline != (char *) NULL )
      free((void *)cmdline);
   if ( (cmdline = (char *) malloc(len)) == (char *) NULL )
      return -1;
   strcpy(cmdline,argv[0]);
   for (i=1; i<argc && argv[i] != (char *) NULL; i++)
   {
      strcat(cmdline," ");
      strcat(cmdline,argv[i]);
   }
   cmdtime = current_time();

   return 0;
}

/* -------------------------- push_config_history ---------------------- */
/*
@@@push_config_history()
@@    Save a line of configuration text for later output in a history block.
@@    If any configuration text for the same keyword was present before,
@@    it may be replaced by the new text, even if the new setting does not
@@    replace all old settings.
@*
@*    Arguments:   line    --  Configuration text line.
@*                 replace --  Replace old text for same keyword (1) or not (0).
@*
@*    Return value:  0 (o.k.),  -1 (memory allocation failed)
@*
 */

int push_config_history (const char *line, int noreplace)
{
   char word[32], word2[32];
   char *text;
   int ipos;
   int replace;
   HSTRUCT *hp, *htmp;
   if ( noreplace )
      replace = 0;
   else
      replace = 1;

   if ( line == (char *) NULL )
      return -1;

   ipos = 0;
   if ( getword(line,&ipos,word,sizeof(word)-1,' ','%') <= 0 )
      return 0;

   if ( (text = (char*) malloc(strlen(line)+1)) == (char *) NULL )
      return -1;
   strcpy(text,line);

   if ( configs == (HSTRUCT *) NULL )
   {
      
      if ( (configs = (HSTRUCT *) calloc(1,sizeof(HSTRUCT))) == 
           (HSTRUCT *) NULL )
      {
         free((void *)text);
         return -1;
      }
      configs->text = text;
      configs->time = current_time();
      configs->next = (HSTRUCT *) NULL;
   }
   else
   {
      for (hp=configs; hp->next != (HSTRUCT *) NULL; hp=hp->next)
      {
         if ( replace )
         {
            ipos = 0;
            if ( getword(hp->text,&ipos,word2,sizeof(word2)-1,' ','%') <= 0 )
               continue;
            /* If a line for the same keyword exists, remove it. */
            if ( strcmp(word,word2) == 0 )
            {
               free((void *)hp->text);
               hp->text = text;
               hp->time = current_time();
               return 0;
            }
         }
      }
      if ( replace )
      {
         ipos = 0;
         if ( getword(hp->text,&ipos,word2,sizeof(word2)-1,' ','%') > 0 )
            /* If a line for the same keyword exists, remove it. */
            if ( strcmp(word,word2) == 0 )
            {
               free((void *)hp->text);
               hp->text = text;
               hp->time = current_time();
               return 0;
            }
      }

      if ( (htmp = (HSTRUCT *) calloc(1,sizeof(HSTRUCT))) == 
           (HSTRUCT *) NULL )
      {
         free((void *)text);
         return -1;
      }
      hp->next = htmp;
      htmp->text = text;
      htmp->time = current_time();
      htmp->next = (HSTRUCT *) NULL;
   }

   return 0;
}

/* -------------------------- write_history ---------------------- */
/*
@@@write_history()
@@    Write the block of accumulated history lines (command line and
@@    configuration lines) to an I/O buffer.
@*
@*    Arguments:   id   --  Identifier (detector number)
@*                 iobuf -- I/O buffer descriptor
@*
@*    Return value:  0 (o.k.), <0 (error with I/O buffer)
@*
 */

int write_history (long id, IO_BUFFER *iobuf)
{
   IO_ITEM_HEADER item_header, sub_item_header;
   HSTRUCT *hp;

   item_header.type = 70;           /* History block is type 70 */
   item_header.version = 1;         /* Version 1 */
   item_header.ident = id;

   if ( put_item_begin(iobuf,&item_header) < 0 )
      return -1;

   if ( cmdline != (char *) NULL )
   {
      sub_item_header.type = 71;         /* Command line history is type 71 */
      sub_item_header.version = 1;       /* Version 1 */
      sub_item_header.ident = id;
      if ( put_item_begin(iobuf,&sub_item_header) < 0 )
         return -1;
      put_long(cmdtime,iobuf);
      put_string(cmdline,iobuf);
      put_item_end(iobuf,&sub_item_header);
   }

   if ( configs != (HSTRUCT *) NULL )
      for (hp=configs; hp!=(HSTRUCT *) NULL; hp=hp->next)
         if ( write_config_history(hp->text,hp->time,id,iobuf) < 0 )
            break;

   return(put_item_end(iobuf,&item_header));
}

/* ------------------------ write_config_history -------------------- */
/*
@@@write_config_history()
@@    Write a configuration history line to an I/O buffer.
@*
@*    Arguments:   htext -- Text of configuration line
@*                 htime -- Time when the configuration was set.
@*                 id    -- Identifier (detector number)
@*                 iobuf -- I/O buffer descriptor
@*
@*    Return value:  0 (o.k.), <0 (error with I/O buffer)
@*
 */

int write_config_history (const char *htext, long htime, long id, IO_BUFFER *iobuf)
{
   IO_ITEM_HEADER item_header;

   item_header.type = 72;           /* Configuration line is type 72 */
   item_header.version = 1;         /* Version 1 */
   item_header.ident = id;

   if ( put_item_begin(iobuf,&item_header) < 0 )
      return -1;
   put_long(htime,iobuf);
   put_string(htext,iobuf);

   return(put_item_end(iobuf,&item_header));
}

/* ------------------------- listtime ----------------------- */

static void listtime (time_t t, FILE *f)
{
   struct tm *tmp;
   
   if ( (tmp = localtime(&t)) == (struct tm *) NULL )
      return;
   fprintf(f,"%04d-%02d-%02d %02d:%02d:%02d",
      tmp->tm_year+1900, tmp->tm_mon+1, tmp->tm_mday, 
      tmp->tm_hour, tmp->tm_min, tmp->tm_sec);
}

/* ------------------------ list_history ----_------------------------- */
/*
@@@list_history()
@@    List history block contents on standard output
@*
@*    Arguments:  iobuf -- I/O buffer descriptor
@*
@*    Return value:  0 (o.k.), <0 (error with I/O buffer)
@*
 */

int list_history (IO_BUFFER *iobuf, FILE *file)
{
   IO_ITEM_HEADER item_header, sub_item_header;
   int rc;
   char line[10240];
   char *s;
   time_t t;
   
   if ( file == NULL )
      file = stdout;
  
   if ( iobuf == (IO_BUFFER *) NULL )
      return -1;
      
   item_header.type = 70;
   if ( (rc = get_item_begin(iobuf,&item_header)) < 0 )
      return rc;
   if ( item_header.version != 1 )
   {
      Warning("Wrong version number of history item to be read.");
      return -1;
   }
   if ( item_header.ident > 0 )
      fprintf(file,"\nHistory block (ident %ld):\n",item_header.ident);
   else
      fprintf(file,"\nHistory block:\n");
   sub_item_header.type = 71;
   if ( search_sub_item(iobuf,&item_header,&sub_item_header) == 0 )
   {
      if ( (rc = get_item_begin(iobuf,&sub_item_header)) < 0 )
      {
         get_item_end(iobuf,&item_header);
         return rc;
      }
      if ( sub_item_header.version != 1 )
      {
         Warning("Wrong version number of command line history item to be read.");
         get_item_end(iobuf,&sub_item_header);
         get_item_end(iobuf,&item_header);
         return -1;
      }
      t = get_long(iobuf);
      get_string(line,sizeof(line)-1,iobuf);
      fprintf(file,"   Command line (dated ");
      listtime(t,file);
      if ( sub_item_header.ident > 0 )
         fprintf(file,", ident %ld):\n      %s\n",sub_item_header.ident,line);
      else
         fprintf(file,"):\n      %s\n",line);
      get_item_end(iobuf,&sub_item_header);
   }
   rewind_item(iobuf,&item_header);
   fprintf(file,"   Configuration data:\n");
   for (;;)
   {
      sub_item_header.type = 72;
      if ( search_sub_item(iobuf,&item_header,&sub_item_header) != 0 )
         break;
      if ( (rc = get_item_begin(iobuf,&sub_item_header)) < 0 )
      {
         get_item_end(iobuf,&item_header);
         return rc;
      }
      if ( sub_item_header.type != 72 )
      {
         Warning("Wrong item type instead of configuration history item");
         return -1;
      }
      if ( sub_item_header.version != 1 )
      {
         Warning("Wrong version number of configuration history item to be read.");
         get_item_end(iobuf,&sub_item_header);
         get_item_end(iobuf,&item_header);
         return -1;
      }
      *line = '\0';
      t = get_long(iobuf);
      get_string(line,sizeof(line)-1,iobuf);
      fprintf(file,"      ");
      for (s=line; *s!='\0'; s++)
      {
         fputc(*s,file);
         if ( *s=='\n' && *(s+1)!='\0' )
            fputs("      ",file);
      }
      if ( s!=line && *(s-1) != '\n' )
         fputc('\n',file);

      get_item_end(iobuf,&sub_item_header);
   }
   fputs("End of history block\n",file);

   return (get_item_end(iobuf,&item_header));
}
