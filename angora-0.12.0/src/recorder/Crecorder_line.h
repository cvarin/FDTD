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

#ifndef CRECORDER_LINE_H
#define CRECORDER_LINE_H

//Declaration of the line recorder class "Crecorder_line"
//Derived from the abstract recorder base class "Crec".

#include "Crec.h"

//Use MPI if not disabled
#ifndef MPI_DISABLE
#include <mpi.h>
#endif

#include <fstream>


class Crecorder_line: public Crec
{//line recorder class
 public:
	 Crecorder_line(const string& myorientation, const int& x1Pos, const int& x2Pos,
		const string& mycomponent, const string& myscale, const string& LineFileName,
		const int& Index);		//constructor

	 //record the line
	 void Record() {
	 	RecordLine();
	 }

 private:
	 int LineRecorderIndex;	//index of the line recorder

	 void RecordLine();

	 double field_value;	//field value at a certain point (dB or linear)
	 double FieldMax,FieldMin;	//maximum and minimum values in the discretization of the fields
	 string orientation,component,scale;

	 void PlacePartialLine(Array<double,1> LineArray); //places partial lines in LineArray
	 ofstream LineFile;
	 Array<Array<double,1>,1> PartialLineBuf_recv;	//1-D array of partial line arrays
	 Array<double,1> PartialLineBuf_send;	//Buffer array for partial line transmitted via MPI
	 Array<double,1> LineArray;			//Array for the whole line
	 int relative_LinePos_x1,relative_LinePos_x2;	//relative position of the line in the perpendicular plane w.r.t. the origin cell
	 int LinePos_x1,LinePos_x2;		//absolute index of the line in the perpendicular plane
	 int TotalLineLength;			//total length of the line
	 int LineMin,LineMax;			//minimum and maximum value of the coordinate variable along the line
	 int X_extra;	//0 or 1, depending on whether there are N or N+1 field components along the orientation of the line, where N is the number of cells in that direction
	 bool PassesThroughNode;	//does the line pass through the node?
	 int line_nodes;		//total number of nodes along the line
	 Array<int,1> LineMaxArray,LineMinArray;	//upper and lower limits of the nodes along the line
	 int xPos,yPos,zPos;	//3-D position indices of the E-field on the main grid
	 int n1,s;	//iteration indices

	 //Line communicator variables
#ifndef MPI_DISABLE
	 MPI_Comm MPI_LineComm;
	 MPI_Status Status;
#endif

	 int line_size,line_rank;
};

#endif
