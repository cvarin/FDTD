/* AUTORIGHTS
Copyright (C) 2006-2012  Ilker R. Capoglu

    This file is part of the Angora package.

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
*/

#ifndef READ_NFFFT_PD_H
#define READ_NFFFT_PD_H

//use the libconfig library
#include "config/cfg_utils.h"

//only the declaration of Cnffft_pd needed: use forward declaration
class Cnffft_pd;
//only the declaration of Ctfsf needed: use forward declaration
class Ctfsf;
//only the declaration of Cpointsources needed: use forward declaration
class Cpointsources;

void read_nffft_pd(Cnffft_pd &NFFFT_pd, const Config& fdtdconfig, const Config& validsettings, const Cpointsources &PointSources, const Ctfsf &TFSF);		//modifies the NFFFT_pd object

#endif
