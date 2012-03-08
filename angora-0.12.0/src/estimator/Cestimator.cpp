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

//Defines the high-precision duration estimator object.

#include "headers.h"

#include "Cestimator.h"

#ifdef _WIN32
#define FDTD_GETTIME ftime
#define FDTD_SEC time
#define FDTD_FRACSEC millitm
#define FDTD_TIMECONVERSION 1000
#else
#define FDTD_GETTIME ntp_gettime
#define FDTD_SEC time.tv_sec
#define FDTD_FRACSEC time.tv_usec
//FIXME: This is a quick, dirty (and incorrect) fix to the nanosecond-microsecond bug for this Linux kernel (tv_usec returns nanosec values)
//#ifdef STA_NANO  //should be removed once tv_usec returns microseconds again in the future (as it should!!!)
// #define FDTD_TIMECONVERSION 1000000000  //should be removed once tv_usec returns microseconds again in the future (as it should!!!)
//#else  //should be removed once tv_usec returns microseconds again in the future (as it should!!!)
 #define FDTD_TIMECONVERSION 1000000  //this is the only line needed if tv_usec returns microseconds again in the future (as it should!!!)
//#endif   //should be removed once tv_usec returns microseconds again in the future (as it should!!!)
//FIXME
#endif

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif


void Cestimator::Initialize(const string& my_logfilename, const int& my_number_of_steps)
{//initializes the time estimator (needs to be called right before the main time loop)
	logfilename = my_logfilename;

	// get user name running the process
//	current_user = getlogin();
	current_user = getenv("LOGNAME");

	sample_simulation_size = min(20,my_number_of_steps);
	total_simulation_size = my_number_of_steps;

	simulation_starting_time = time(NULL);

	//start high-precision timing
	FDTD_GETTIME(&starting_time);
}

void Cestimator::MakeEstimate(const int& n)
{
	if (n==sample_simulation_size-1)
	{
		FDTD_GETTIME(&finishing_time); //finish high-precision timing

		time_subtract(&finishing_time,&starting_time,&elapsed_time);   //CAREFUL: Modifies finishing_time and starting_time!!

		estimate = (time_t)ceil((FDTD_TIMECONVERSION*elapsed_time.FDTD_SEC+elapsed_time.FDTD_FRACSEC)
								/(double)FDTD_TIMECONVERSION*((double)total_simulation_size/(double)sample_simulation_size));
												//FDTD_TIMECONVERSION converts from fraction of a second to second

		write_log(simulation_starting_time, simulation_starting_time + estimate,estimate);	//write out timing data
	}
}

void Cestimator::WriteActualFinishingTime(time_t& actual_finishing_time, time_t& actual_elapsed_time)
{
	LogFile.open(logfilename.c_str(),ios::out|ios::app);
	if (!LogFile)
	{
		cout << "Error opening log file!" << endl << endl;
//		exit(-1);
	}

	actual_finishing_time = time(NULL);	//get the time in basic time_t type
	actual_elapsed_time = (int)difftime(actual_finishing_time,simulation_starting_time);  //get the running time of the simulation

	//Write out the finishing time and date
	tm *struct_actual_finishing_time = localtime(&actual_finishing_time);
	const int time_size=40;
	char *str_actual_finishing_time = new char[time_size];	//time & date in string format
	strftime(str_actual_finishing_time,time_size,"%x %I:%M:%S%p",struct_actual_finishing_time);
	LogFile << "    Simulation finished on " << str_actual_finishing_time << endl;
	//Write out the elapsed time
	// Obtain coordinated universal time:
	tm *struct_actual_elapsed_time = gmtime(&actual_elapsed_time);
	if (struct_actual_elapsed_time->tm_hour>0)
	{
		LogFile << "    Elapsed time : " << struct_actual_elapsed_time->tm_hour << " hours, "
			<< struct_actual_elapsed_time->tm_min << " minutes and "
			<< struct_actual_elapsed_time->tm_sec << " seconds." << endl;
	}
	else if (struct_actual_elapsed_time->tm_min>0)
	{
		LogFile << "    Elapsed time : " << struct_actual_elapsed_time->tm_min << " minutes and "
			<< struct_actual_elapsed_time->tm_sec << " seconds." << endl;
	}
	else
	{
		LogFile << "    Elapsed time : " << struct_actual_elapsed_time->tm_sec << " seconds." << endl;
	}

	LogFile.close();
}

void Cestimator::write_log(const time_t& simulation_starting_time,
						 const time_t& estimated_finishing_time,
						 const time_t& estimated_duration)
{
	LogFile.open(logfilename.c_str(),ios::out|ios::app);
	if (!LogFile)
	{
		cout << "Error opening log file!" << endl << endl;
		exit(-1);
	}

	tm *struct_simulation_starting_time = localtime(&simulation_starting_time);
	const int time_size=40;
	char *str_simulation_starting_time = new char[time_size];	//starting time & date in string format
	strftime(str_simulation_starting_time,time_size,"%x %I:%M:%S%p",struct_simulation_starting_time);
	LogFile << current_user << " started Angora on " << str_simulation_starting_time << endl;

	tm *struct_estimated_finishing_time = localtime(&estimated_finishing_time);
	char *str_finishing_time = new char[time_size];	//finishing time & date in string format
	strftime(str_finishing_time,time_size,"%x %I:%M:%S%p",struct_estimated_finishing_time);
	LogFile << "    Estimated to finish on " << str_finishing_time << endl;

	tm *struct_estimated_time = gmtime(&estimated_duration);
	if (struct_estimated_time->tm_hour>0)
	{
		LogFile << "    Estimated duration : " << struct_estimated_time->tm_hour << " hours, "
			<< struct_estimated_time->tm_min << " minutes and "
			<< struct_estimated_time->tm_sec << " seconds." << endl;
	}
	else if (struct_estimated_time->tm_min>0)
	{
		LogFile << "    Estimated duration : " << struct_estimated_time->tm_min << " minutes and "
			<< struct_estimated_time->tm_sec << " seconds." << endl;
	}
	else
	{
		LogFile << "    Estimated duration : " << struct_estimated_time->tm_sec << " seconds." << endl;
	}

	LogFile.close();
}

int Cestimator::time_subtract(FDTD_TIMEVAR *finish, FDTD_TIMEVAR *start, FDTD_TIMEVAR *elapsed)
{
/* Subtract the `FDTD_TIMEVAR' values 'finish' and 'start',
storing the result in 'elapsed'.
Return 1 if the difference is negative, otherwise 0.  */

/* Perform the carry for the later subtraction by updating start. */
if (finish->FDTD_FRACSEC < start->FDTD_FRACSEC) {
 int nsec = (start->FDTD_FRACSEC - finish->FDTD_FRACSEC) / FDTD_TIMECONVERSION + 1;
 start->FDTD_FRACSEC -= FDTD_TIMECONVERSION * nsec;
 start->FDTD_SEC += nsec;
}
if (finish->FDTD_FRACSEC - start->FDTD_FRACSEC > FDTD_TIMECONVERSION) {
 int nsec = (finish->FDTD_FRACSEC - start->FDTD_FRACSEC) / FDTD_TIMECONVERSION;
 start->FDTD_FRACSEC += FDTD_TIMECONVERSION * nsec;
 start->FDTD_SEC -= nsec;
}

/* Compute the time remaining to wait.
  FDTD_FRACSEC is certainly positive. */
elapsed->FDTD_SEC = finish->FDTD_SEC - start->FDTD_SEC;
elapsed->FDTD_FRACSEC = finish->FDTD_FRACSEC - start->FDTD_FRACSEC;

/* Return 1 if result is negative. */
return finish->FDTD_SEC < start->FDTD_SEC;
}
