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

#ifndef READ_LOGGING_DATA_H
#define READ_LOGGING_DATA_H

//use the libconfig library
#include "config/cfg_utils.h"

//forward declaration
class Cestimator;


void read_logging_data(Cestimator &Estimator, const Config& fdtdconfig);		//modifies the Estimator object

#endif // READ_LOGGING_DATA_H
