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

//Defines the recorder class "Crecorder"

#include "headers.h"

#include "Crecorder.h"

#include "Crecorder_movie.h"
#include "Crecorder_line.h"
#include "Crecorder_fieldvalue.h"

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif


int Crecorder::AddMovieRecorder(const string& mysection, const int& Pos,
				const string& mycomponent, const string& myscale, const string& myrecordingtype, const string& MovieFileName,
				const bool& OnlyRecordsGeometry)
{//adds a movie recorder to the recorder container
	boost::shared_ptr<Crec> new_movie_recorder_ptr(new Crecorder_movie(mysection, Pos, mycomponent, myscale, myrecordingtype, MovieFileName, MovieRecorders.size(),OnlyRecordsGeometry));
	MovieRecorders.push_back(new_movie_recorder_ptr);
}

int Crecorder::AddLineRecorder(const string& myorientation, const int& x1Pos, const int& x2Pos,
				const string& mycomponent, const string& myscale, const string& LineFileName)
{//adds a line recorder to the recorder container
	boost::shared_ptr<Crec> new_line_recorder_ptr(new Crecorder_line(myorientation, x1Pos, x2Pos, mycomponent, myscale, LineFileName, LineRecorders.size()));
	LineRecorders.push_back(new_line_recorder_ptr);
}

int Crecorder::AddFieldValueRecorder(const int& xPos, const int& yPos, const int& zPos,
				const string& mycomponent, const string& myscale, const string& FieldValueFileName)
{//adds a field value recorder to the recorder container
	boost::shared_ptr<Crec> new_fieldvalue_recorder_ptr(new Crecorder_fieldvalue(xPos, yPos, zPos, mycomponent, myscale, FieldValueFileName, FieldValueRecorders.size()));
	FieldValueRecorders.push_back(new_fieldvalue_recorder_ptr);
}

void Crecorder::Record()
{//records everything that needs to be recorded at the current time step
	RecordMovieFrames();
	RecordLines();
	RecordFieldValues();
}

void Crecorder::RecordMovieFrames()
{//records all movie frames
	for (int i=0; i<MovieRecorders.size(); i++)
	{
		MovieRecorders[i]->Record();		//each movie recorder adds movie frame
	}
}

void Crecorder::RecordLines()
{//records all lines
	for (int i=0; i<LineRecorders.size(); i++)
	{
		LineRecorders[i]->Record();		//each line recorder adds line
	}
}

void Crecorder::RecordFieldValues()
{//records all field values
	for (int i=0; i<FieldValueRecorders.size(); i++)
	{
		FieldValueRecorders[i]->Record();		//each field value recorder adds field value
	}
}
