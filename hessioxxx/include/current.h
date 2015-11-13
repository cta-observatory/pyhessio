/* ============================================================================

   Copyright (C) 1993, 2007, 2010  Konrad Bernloehr

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

/** @file current.h
 *  @short Header file for optional current time add-on to warning.c.
 *
 *  @author  Konrad Bernloehr
 *  @date 1993 (original version)
 *  $Date: 2010/07/20 13:37:45 $
 *  $Revision: 1.4 $
 *  
 */  

#ifndef CURRENT_H__LOADED

#define CURRENT_H__LOADED 1

#ifdef __cplusplus
extern "C" {
#endif

/* Include file <time.h> before this file. */

#ifdef __Current_Module__
time_t last_data_time = 0;
#else
extern time_t last_data_time;
#endif

#define DEFAULT_LOCAL_OFFSET 3600

time_t current_time (void);
time_t current_localtime (void);
void set_current_offset (long _toffset);
void set_local_offset (long _local_offset);
void reset_local_offset (void);
char *time_string (void);
time_t mkgmtime (struct tm *tms);

#ifdef __cplusplus
}
#endif

#endif
