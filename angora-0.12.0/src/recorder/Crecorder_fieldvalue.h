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

#ifndef CRECORDER_FIELDVALUE_H
#define CRECORDER_FIELDVALUE_H

//Declaration of the field-valie recorder class "Crecorder_fieldvalue"
//Derived from the abstract recorder base class "Crec".

#include "Crec.h"

#include <fstream>


class Crecorder_fieldvalue: public Crec
{//field value recorder class
 public:
	 Crecorder_fieldvalue(const int& xPos, const int& yPos, const int& zPos,
		 const string& mycomponent, const string& myscale, const string& myfilename,
		 const int& Index);		//constructor
	 ~Crecorder_fieldvalue(); //dtor

	 //record the field value
	 void Record() {
	 	RecordFieldValue();
	 }

 private:
	 int FieldValueRecorderIndex;	//index of the field value recorder

	 void RecordFieldValue();

	 bool FieldComponentInNode;	//is the recorded field component in this node?

	 double field_value;	//field value at a certain point (dB or linear)
	 double FieldMax,FieldMin;	//maximum and minimum values in the discretization of the fields
	 string component,scale;

	 Array<double,1> field_value_array; //array that holds the field values to be recorded

	 int counter; //internal counter

//	 ofstream FieldValueFile;
//	 H5File field_value_file;
	 const string FieldValueFileName;

	 int relative_FieldPos_x,relative_FieldPos_y,relative_FieldPos_z;	//relative x,y,z positions of the recorded field point w.r.t. the origin cell
	 int FieldPos_x,FieldPos_y,FieldPos_z;	//absolute x,y,z indices of the recorded field point
};

#endif
