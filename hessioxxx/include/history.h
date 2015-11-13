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

/** @file history.h
    @short Keep blocks of history in the data (like command line of
       programs operating on the data, ...)
    
    @author  Konrad Bernloehr
    @date    1997 to 2010
    @date    @verbatim $Date: 2014/02/20 11:40:42 $ @endverbatim
    @version @verbatim $Revision: 1.5 $ @endverbatim
*/

#ifndef __HISTORY_LOADED

#define __HISTORY_LOADED 1

#ifndef INITIAL_H__LOADED
#include "initial.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* history.c */
int push_command_history (int argc, char **argv);
int push_config_history (const char *line, int replace);
int write_history (long id, IO_BUFFER *iobuf);
int write_config_history (const char *htext, long htime, 
    long id, IO_BUFFER *iobuf);
int list_history (IO_BUFFER *iobuf, FILE *file);

#ifdef __cplusplus
}
#endif

#endif
