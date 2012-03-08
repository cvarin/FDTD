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

//Variables related to the mapping between time steps and actual time

#include "headers.h"

#include "time_axis.h"

static double initial_time_value = 0;
//static bool initial_time_value_can_be_set = true;


double get_initial_time_value()
{
	return initial_time_value;
}

void set_initial_time_value(const double&time_value)
{
//	if (initial_time_value_can_be_set)
//	{
		initial_time_value = time_value;
//	}
//	else
//	{
//		throw AngoraTimeAxisException();
//	}
}

//void allow_setting_initial_time_value()
//{
//	initial_time_value_can_be_set = true;
//}
//
//void disallow_setting_initial_time_value()
//{
//	initial_time_value_can_be_set = false;
//}
