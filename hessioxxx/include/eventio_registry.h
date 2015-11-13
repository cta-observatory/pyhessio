/* ============================================================================

   Copyright (C) 2014  Konrad Bernloehr

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

/** @file eventio_registry.h
 *  @short Register and enquire about well-known I/O block types.
 *
 *
 *  @author  Konrad Bernloehr 
 *  @date    2014
 *  @date    CVS $Date: 2014/06/01 11:33:04 $
 *  @version CVS $Revision: 1.2 $
 */

#ifndef EV_REG_H__LOADED            /* Ignore if included a second time */

#define EV_REG_H__LOADED 1

#ifndef INITIAL_H__LOADED
#include "initial.h"
#endif
#include "io_basic.h"

#ifdef __cplusplus
extern "C" {
#endif

int read_eventio_registry(const char *fname);
struct ev_reg_entry *find_ev_reg_std(unsigned long t);
void set_ev_reg_std(void);

#ifdef __cplusplus
}
#endif

#endif
