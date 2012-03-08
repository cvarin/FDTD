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

#ifndef CWF_SHARED_PTR_H
#define CWF_SHARED_PTR_H

//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>

class Cwf;
typedef boost::shared_ptr<Cwf> Cwf_shared_ptr;
typedef boost::shared_ptr<const Cwf> const_Cwf_shared_ptr;

#endif // CWF_SHARED_PTR_H
