/* ============================================================================

   Copyright (C) 1997, 2001, 2007, 2010  Konrad Bernloehr

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

/** @file io_history.h
 *  @short Record history of configuration settings/commands.
 *
 *  @author  Konrad Bernloehr 
 *  @date 1997 to 2010
 *  @date    @verbatim CVS $Date: 2014/02/20 11:40:42 $ @endverbatim
 *  @version @verbatim CVS $Revision: 1.5 $ @endverbatim
 */

#ifndef IO_HISTORY_H__LOADED
#define IO_HISTORY_H__LOADED 1

#ifndef INITIAL_H__LOADED
#include "initial.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* io_history.c */
int push_command_history (int argc, char **argv);
int push_config_history (const char *line, int noreplace);
int write_history (long id, IO_BUFFER *iobuf);
int write_config_history (const char *htext, long htime, 
    long id, IO_BUFFER *iobuf);
int list_history (IO_BUFFER *iobuf, FILE *file);

#ifdef __cplusplus
}
#endif

#endif
