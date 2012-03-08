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

#ifndef READ_NFFFT_TD_H
#define READ_NFFFT_TD_H

//use the libconfig library
#include "config/cfg_utils.h"

//only the declaration of Cnffft_td needed: use forward declaration
class Cnffft_td;
//only the declaration of Cpointsources needed: use forward declaration
class Cpointsources;

void read_nffft_td(Cnffft_td &NFFFT_td, const Config& fdtdconfig, const Config& validsettings, const Cpointsources &PointSources);		//modifies the NFFFT_td object

#endif
