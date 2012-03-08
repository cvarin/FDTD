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

#ifndef CRECORDER_MOVIE_H
#define CRECORDER_MOVIE_H

//Declaration of the movie recorder class "Crecorder_movie"
//Derived from the abstract recorder base class "Crec".

#include "Crec.h"

//for the definition of MaterialId
#include "material_id.h"

//Use MPI if not disabled
#ifndef MPI_DISABLE
#include <mpi.h>
#endif

#include <fstream>

// MPI_USE_PARALLEL_IO is determined at build time by configure option


class Crecorder_movie: public Crec
{//movie recorder class
 public:
	 Crecorder_movie(const string& mysection, const int& Pos,
		 const string& mycomponent, const string& myscale, const string& myrecordingtype,
    	 const string& MovieFileName,
		 const int& Index,
		 const bool& OnlyRecordGeometry=false);		//constructor

	 //record the movie frame
	 void Record() {
	 	RecordMovieFrame();
	 }

#ifndef MPI_DISABLE
#ifdef MPI_USE_PARALLEL_IO
	 ~Crecorder_movie();
#endif
#endif


 private:
	 int MovieRecorderIndex;	//index of the movie recorder

	 const bool OnlyRecordGeometry;	//if set, the movie recorder only records the geometry (permittivity,conductivity, etc.), NOT the field values at each time step

	 int NMOVIESTEPS;	//number of movie steps (0 or NSTEPS)

	 void RecordMovieFrame();
	 //more general form of RecordMovieFrame(): records any array (field, material property, etc.)
	 void RecordMovieFrame_uchar1(const Array<double,3>* RecordedArray);	//for "unsigned char" (1 byte) recording
	 void RecordMovieFrame_dbl8(const Array<double,3>* RecordedArray);	//for "double" (8 byte) recording

 	 //Only the electric material properties are recorded
	 //for "double" (8 byte) recording of permittivity
	 void RecordPermittivity(const Array<ElectricMaterialIndexType_X,3>* PropertyIndexArray, const Array<double,1>* PermittivityArray); /** isotropy assumed! **/
	 //for "double" (8 byte) recording of conductivity
	 void RecordConductivity(const Array<ElectricMaterialIndexType_X,3>* PropertyIndexArray, const Array<double,1>* ConductivityArray); /** isotropy assumed! **/

	 void PlaceFrame(const Array<double,3>* RecordedArray, Array<unsigned char,2> FrameArray);  //places field values on the cross-section (yz, xz, or xy) into FrameArray [in "unsigned char" (1 byte) type]
	 void PlaceFrame(const Array<double,3>* RecordedArray, Array<double,2> FrameArray);  //places field values on the cross-section (yz, xz, or xy) into FrameArray [in "double" (8 byte) type]
	 void PlacePermittivityFrame(const Array<ElectricMaterialIndexType_X,3>* PropertyIndexArray, const Array<double,1>* PermittivityArray, Array<double,2> DestinationArray);  /** isotropy assumed! **/  //places the permittivity component (eps_x,eps_y or eps_z) defined by PermittivityArray into DestinationArray in "double" (8 byte) type
	 void PlaceConductivityFrame(const Array<ElectricMaterialIndexType_X,3>* PropertyIndexArray, const Array<double,1>* ConductivityArray, Array<double,2> DestinationArray);  /** isotropy assumed! **/  //places the conductivity component (cond_e_x,cond_e_y or cond_e_z) defined by MaterialArray into ConductivityArray in "double" (8 byte) type

	void WritePreamble(); //writes the preamble to the output file

 	 ofstream MovieFile;

	 double field_value;	//field value at a certain point (dB or linear)
	 double FieldMax,FieldMin;	//maximum and minimum values in the discretization of the fields
	 string section,component,scale,recordingtype;

#ifndef MPI_DISABLE
#ifndef MPI_USE_PARALLEL_IO
	 //buffer arrays for "unsigned char" (1 byte) recordings:
	 Array<Array<unsigned char,2>,2> uchar1_PartialFrameBuf_recv;	//2-D array of partial buffer arrays
	 Array<unsigned char,2> uchar1_PartialFrameBuf_send;	//Buffer arrays for partial frames transmitted via MPI
	 Array<unsigned char,2> uchar1_MovieFrameArray;			//Array for movie frame
	 //buffer arrays for "double" (8 byte) recordings:
	 Array<Array<double,2>,2> dbl8_PartialFrameBuf_recv;	//2-D array of partial buffer arrays
	 Array<double,2> dbl8_PartialFrameBuf_send;	//Buffer arrays for partial frames transmitted via MPI
	 Array<double,2> dbl8_MovieFrameArray;			//Array for movie frame

	 Array<int,1> FrameMaxArray_X1,FrameMinArray_X1,FrameMaxArray_X2,FrameMinArray_X2;	//upper and lower limits of the nodes on the cross-section
#endif

#ifdef MPI_USE_PARALLEL_IO
	 MPI_File MPIMovieFile;
	 MPI_Datatype MPIMovieFileType;
	 //Partial frame for "unsigned char" recording:
	 Array<unsigned char,2> uchar1_PartialFrame;
	 //Partial frame for "double" recording:
	 Array<double,2> dbl8_PartialFrame;
#endif
#else
	 // there is no MPI, so there is only one movie array (char or double)
	 Array<unsigned char,2> uchar1_MovieFrameArray;		//Array for movie frame ("unsigned char" recording)
	 Array<double,2> dbl8_MovieFrameArray;				//Array for movie frame ("double" recording)
#endif

	 const Array<double,3>* RecordedFieldComponent;	//pointer to the recorded component array (non-changeable by the pointer), e.g. Ex.
	 const Array<ElectricMaterialIndexType_X,3>* RecordedMaterialIndex;	//pointer to the recorded material-index array (non-changeable by the pointer) e.g. Media_Ex.
	 const Array<double,1>* RecordedPermittivityComponent;	//pointer to the recorded permittivity component
	 const Array<double,1>* RecordedConductivityComponent;	//pointer to the recorded conductivity component

	 int TotalFrameSize_X1,TotalFrameSize_X2;		//dimensions of the whole movie frame
	 int FrameMin_X1,FrameMax_X1,FrameMin_X2,FrameMax_X2;	//minimum and maximum values of the coordinate variables across the cross-section
	 int X1_extra,X2_extra;	//0 or 1, depending on whether there are N or N+1 field components in a certain direction, where N is the number of cells in that direction

	 int OriginX1,OriginX2;	//indices of the field components on the recorded cross section that correspond to the "origin"
	 						// (this means that the OriginX1'th field component has the same index along the x1 direction as the corresponding component in the "origin cell" of the grid -- for ex. OriginY)
	 Array<double,1> X1_coord_range,X2_coord_range;	//spatial coordinates (in grid cells) of the field components on the recorded cross section

	 bool PassesThroughNode;	//does the cross-section pass through the node?
	 int movie_nodes_X1,movie_nodes_X2;	//total number of nodes in the x1 and x2 direction on the cross-section
	 int xPos,yPos,zPos;	//3-D position indices of the E-field on the main grid
	 int n1,n2,x1,x2;	//iteration indices

	 int relative_FramePos;		//relative position of the movie frame w.r.t. the origin cell
	 int FramePos;				//position of the movie frame (absolute)

#ifndef MPI_DISABLE
	 //Movie communicator variables
	 MPI_Comm MPI_MovieComm,MPI_CartMovieComm;
	 MPI_Status Status;

#ifndef MPI_USE_PARALLEL_IO
	 MPI_Comm MPI_CartMovieComm_X1,MPI_CartMovieComm_X2;
#endif
#endif

	 int movie_size,movie_rank;
	 int movie_rank_X1,movie_rank_X2;
	 int recv_rank;
	 int recv_coord[2];

};

#endif
