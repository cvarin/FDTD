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

#ifndef TIME_AXIS_H
#define TIME_AXIS_H

//Variables related to the mapping between time steps and actual time

////base Angora exception class
//#include "angora_excp.h"


double get_initial_time_value(); //get the time value corresponding to the start of the simulation
void set_initial_time_value(const double&time_value); //set the time value corresponding to the start of the simulation

//class AngoraTimeAxisException: public AngoraException
//{// exception class for signaling that setting the time axis is not allowed
//public:
//  virtual ~AngoraTimeAxisException() throw() {};
//
//  virtual const string getError() const
//  {
//  	return "Error: Changing the initial time value is not allowed once the simulation has started.";
//  };
//
//  virtual const char* what() const throw()
//  {
//  	return getError().c_str();
//  }
//};


#endif
