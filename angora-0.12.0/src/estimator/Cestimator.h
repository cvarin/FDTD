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

#ifndef CESTIMATOR_H
#define CESTIMATOR_H

//use the libconfig library
#include "config/cfg_utils.h"

//Declaration of the high-precision duration estimator object.

#include <fstream>

#ifdef _WIN32
#include <sys/types.h>
#include <sys/timeb.h>
#define FDTD_TIMEVAR timeb
#else
#include <sys/timex.h>
#define FDTD_TIMEVAR ntptimeval
#endif


class Cestimator
{
 public:
	 void Initialize(const string& my_logfilename, const int& my_number_of_steps);//initializes the estimator (needs to be called right before the time loop)

	 void MakeEstimate(const int& n);
	 void WriteActualFinishingTime(time_t& actual_finishing_time, time_t& actual_elapsed_time);

 private:
	 ofstream LogFile;
	 string logfilename;

	 //login name of user running Angora
	 string current_user;

	 int sample_simulation_size,total_simulation_size;
	 time_t simulation_starting_time,estimate;
	 FDTD_TIMEVAR starting_time,finishing_time,elapsed_time;

	 void write_log(const time_t& simulation_starting_time,
		 const time_t& estimated_finishing_time,
		 const time_t& estimated_duration);
	 int time_subtract(FDTD_TIMEVAR *finish, FDTD_TIMEVAR *start, FDTD_TIMEVAR *elapsed);
};

#endif
