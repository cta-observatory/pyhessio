/* ============================================================================

   Copyright (C) 1995, 2000, 2007, 2010  Konrad Bernloehr

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

/** @file current.c
 *  @short Code to insert current time string into warnings.
 *
 *  This code is meant for inserting time strings into warnings
 *  passed through the code of warning.c.
 *  It is not currently used in my code and is not yet 
 *  multi-threading safe. 
 *  It is here mainly for improved backward-compatibility with config.c.
 *
 *  @author  Konrad Bernloehr
 *  @date    1995, 2000, 2007
 *  @date    $Date: 2018/05/11 11:54:29 $
 *  @version $Revision: 1.8 $
 */

#include "initial.h"
#define __Current_Module__ 1
#include "current.h"

/* Note: the following parameters are not used in a multi-threading */
/* safe way and should be modified only at program startup. */
static long tcor_parm[3];
static long local_offset = DEFAULT_LOCAL_OFFSET;
static int local_set =0;

static long time_correction (time_t now);

/* ----------------------- current_time -------------------------- */
/**
 *  @short Get the current time in seconds since 1970.0 GMT.
 *
 *  The resulting time includes the last time correction
 *  with respect to the server. Therefore, as long as the
 *  clock on the local computer is not much slower or faster
 *  than the clock on the I/O server, it is the current
 *  Greenwich Mean Time on the I/O server.
 *
 *  @return  Time in seconds since 0h UT on January 1, 1970.
 */

time_t current_time ()
{
   time_t now;
   (void) time(&now);
   return (now+time_correction(now));
}

/* -------------------- current_localtime ------------------------ */
/**
 *  @short Like current_time() but should return time in the local time zone.
 *
 *  The offset of the time zone to GMT must be set by set_local_offset()
 *  or it is derived from the machine's internal time zone setup.
 */

time_t current_localtime ()
{
   if ( !local_set )
   {
      time_t now;
      struct tm *local;
      (void) time(&now);
      if ( (local = localtime(&now)) != (struct tm *) NULL )
         /* local_offset = local->tm_gmtoff; */ /* Works on UNIX */
         local_offset = (long) (mkgmtime(local) - now);
      local_set = 1;
   }
   return (current_time()+local_offset);
}

/* -------------------- set_current_offset --------------------- */
/**
 *  @short Set current time offset.
 *
 *  Set the offset between the time on the time server and the
 *  local time (in seconds in the sense 'remote-local').
 *
 *  Note: in a multi-threaded program this function should be
 *  called only at program startup.
 *
 *  @param  off   Time offset in seconds
 *
 *  @return (none)
 */

void set_current_offset (long off)
{
   tcor_parm[0] = off;
}

/* --------------------- set_local_offset --------------------- */
/**
 *  @short Set offset of local time zone. 
 *
 *  Set the offset between the local time zone and GMT
 *  (in seconds in the sense 'local zone - GMT').
 *
 *  Note: in a multi-threaded program this function should be
 *  called only at program startup.
 *
 *  @param  off  Time offset in seconds
 *
 *  @return  (none)
 *
 */

void set_local_offset (long off)
{
   local_offset = off;
   local_set = 1;
}

/* -------------------- reset_local_offset -------------------- */
/**
 *  @short Reset any previous local time offset.
 *
 *  Reset any previously set local time offset. The next call
 *  to current_localtime() will therefore set the offset to
 *  present system value.
 *
 *  Note: in a multi-threaded program this function should be
 *  called only at program startup.
 *
 *  @return (none)
 */

void reset_local_offset ()
{
   local_set = 0;
   local_offset = 0;
}

/* --------------------- time_correction ----------------------- */

static long time_correction (time_t now)
{
   return tcor_parm[0];
}

/* ---------------------------- time_string ---------------------- */
/**
 *  @short Return a pointer to a formatted time-and-date string. 
 *  This string is reused (changed) on the next call.
 *
 *  @return Time/date character string pointer.
 */

char *time_string ()
{
   struct tm *tmp;
   time_t now;
   static char text[80];

   now = current_localtime();
   *text = '\0';
   if ( (tmp = gmtime(&now)) == (struct tm *) NULL )
#ifdef OS_VAXVMS
     /* On VAX/VMS the timezone is unknown and gmtime() always fails. */
     if ( (tmp = localtime(&now)) == (struct tm *) NULL )
#endif
      return text;
   snprintf(text,sizeof(text)-1,"%02d.%02d.%02d %02d:%02d:%02d",
      tmp->tm_mday, tmp->tm_mon+1, tmp->tm_year,
      tmp->tm_hour, tmp->tm_min, tmp->tm_sec);
   return text;
}

/* ------------------------ mkgmtime ------------------------ */
/**
 *  @short Inverse to gmtime() library function.
 *
 *  Inverse to gmtime() library function without correction for
 *  timezone and daylight saving time.
 *
 *  @param  tms   Pointer to time structure as filled by gmtime().
 *
 *  @return  Time in seconds since 1970.0
 *
 */

time_t mkgmtime (struct tm *tms)
{
   int year, month, day;
   int thour, tmin, tsec;
   long ndays;
   int i;
   time_t result;

   static int mlen[] = {31,28,31,30,31,30,31,31,30,31,30,31};

   if ( tms == (struct tm *) NULL )
      return((time_t) 0);

   year  = tms->tm_year;     /* 70 ... 199 */
   month = tms->tm_mon;      /*  0 ... 11  */
   day   = tms->tm_mday;     /*  1 ... 31  */
   thour = tms->tm_hour;     /*  0 ... 23  */
   tmin  = tms->tm_min;      /*  0 ... 59  */
   tsec  = tms->tm_sec;      /*  0 ... 59  */

   if ( year < 1000 )
      year += 1900;

   if ( year < 1970 || year >= 2100 )
      return((time_t) 0);

   ndays = 0;
   for ( i=0; i<month && i<12; i++ )
      ndays += mlen[i];
   if ( month >= 2 && (year/4)*4 == year ) /* leap year in 1901...2099 */
      ndays ++;
   ndays += day - 1;
   tms->tm_yday = (int) ndays;

   for ( i=1970; i<year; i++ )
   {
      if ( (i/4)*4 == i )
        ndays += 366;
      else
        ndays += 365;
   }
   tms->tm_wday = (int) ((ndays+4)%7);

   result = (time_t)(tsec + 60*(tmin + 60*(thour + 24*ndays)));
   return(result);
}
