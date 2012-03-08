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

#ifndef HEADERS_H
#define HEADERS_H

//common headers used in all files

//GNU autoconf config header
#if HAVE_CONFIG_H
#include <config.h>
#endif

// #pragma warning( disable : 4800 4805 4996)

//use the blitz++ array library
#include <blitz/array.h>
#define _USE_MATH_DEFINES
#ifdef _DEBUG
	#define BZ_DEBUG
#endif

using namespace blitz;
using namespace std;

//Include common physical and mathematical constants
#include "constants.h"

//the minimum representable double-precision number
#define LIBSTD_DBL_EPSILON numeric_limits<double>::epsilon()
//the maximum representable double-precision number
#define LIBSTD_DBL_MAX numeric_limits<double>::max()

#endif
