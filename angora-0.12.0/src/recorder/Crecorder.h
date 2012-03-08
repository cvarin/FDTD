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

#ifndef CRECORDER_H
#define CRECORDER_H

//Declaration of the recorder class "Crecorder"

//Use MPI if not disabled
#ifndef MPI_DISABLE
#include <mpi.h>
#endif

#include <fstream>

//for the vector STL class
#include <vector>
//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>

//only the declaration of Crec needed: use forward declaration
class Crec;


class Crecorder
{//generic recorder class (may include movie or field value recorders)
 public:
	 int AddMovieRecorder(const string& mysection, const int& Pos,
		const string& mycomponent, const string& myscale, const string& myrecordingtype, const string& MovieFileName,const bool& OnlyRecordsGeometry);
	 int AddLineRecorder(const string& myorientation, const int& x1Pos, const int& x2Pos,
		const string& mycomponent, const string& myscale, const string& LineFileName);
	 int AddFieldValueRecorder(const int& xPos, const int& yPos, const int& zPos,
		const string& mycomponent, const string& myscale, const string& FieldValueFileName);

	 //record everything that needs to be recorded at the current time step
	 void Record();

 private:

	 string RecorderOutputDir;	//recorder output directory

	 void RecordMovieFrames();
	 void RecordLines();
	 void RecordFieldValues();

	 //movie recorders
	 vector<boost::shared_ptr<Crec> > MovieRecorders;
	 //line recorders
	 vector<boost::shared_ptr<Crec> > LineRecorders;
	 //field value recorders
	 vector<boost::shared_ptr<Crec> > FieldValueRecorders;
};

#endif
