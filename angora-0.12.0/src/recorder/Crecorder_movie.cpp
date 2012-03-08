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

//Defines the movie recorder class "Crecorder_movie"

#include "headers.h"

#include "Crecorder_movie.h"

//routines for getting/setting the initial time value in the simulation
#include "time_axis.h"

extern int angora_version_major,angora_version_minor,angora_version_revision;

extern double dx,dt;

extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;
extern int NSTEPS;
extern int OriginX,OriginY,OriginZ;

extern Array<double,3> Ex,Ey,Ez;
extern Array<double,3> Hx,Hy,Hz;

#ifndef MPI_DISABLE
extern MPI_Comm MPI_CartSubComm;
#endif
extern int GridIndex;
extern int rank;
extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern int nodes_x, nodes_y, nodes_z;

extern Array<ElectricMaterialIndexType_X,3> Media_Ex;
extern Array<ElectricMaterialIndexType_Y,3> Media_Ey;
extern Array<ElectricMaterialIndexType_Z,3> Media_Ez;
extern Array<MagneticMaterialIndexType_X,3> Media_Hx;
extern Array<MagneticMaterialIndexType_Y,3> Media_Hy;
extern Array<MagneticMaterialIndexType_Z,3> Media_Hz;

extern Array<double,1> eps_x,mu_y,mu_z;
extern Array<double,1> cond_e_x,cond_h_y,cond_h_z;
extern Array<double,1> mu_x,eps_y,eps_z;
extern Array<double,1> cond_h_x,cond_e_y,cond_e_z;

extern double max_field_value;
extern double dB_accuracy;


Crecorder_movie::Crecorder_movie(const string& mysection, const int& Pos,
				const string& mycomponent, const string& myscale, const string& myrecordingtype,
				const string& MovieFileName,
				const int& Index,
		 		const bool& OnlyRecordsGeometry)
					: section(mysection), relative_FramePos(Pos), component(mycomponent), scale(myscale), recordingtype(myrecordingtype), MovieRecorderIndex(Index), OnlyRecordGeometry(OnlyRecordsGeometry)
//NOTE : if the magnitude is recorded, "Pos" must denote the CELL INDEX (therefore, must not exceed the number of
// cells in that direction), since magnitude is evaluated at the center of a cell by averaging four surrounding values
// for each field component. For the other field values, "Pos" must denote the actual index of the cross-section in the
// 3D component array.
{//constructor for the movie recorder
	if (OnlyRecordGeometry)
	{
		NMOVIESTEPS = 0;	//don't record any field values
	}
	else
	{
		NMOVIESTEPS = NSTEPS;
	}

	//Extract the maximum and minimum values in the discretization of the fields
	if (scale=="dB")
	{
		FieldMax = 20.0*log10(abs(max_field_value)+1e-20);
		FieldMin = FieldMax+dB_accuracy;
	}
	else if (scale=="linear")
	{
		FieldMax = max_field_value;
		FieldMin = -max_field_value;
	}
	else if (scale=="absolute")
	{
		FieldMax = max_field_value;
		FieldMin = 0;
	}
	else
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<string>(func_name,scale,
			"(valid arguments are \"dB\", \"linear\", or \"absolute\")");
	}

	// check if the recording type is valid
	if ((recordingtype!="uchar1")&&(recordingtype!="dbl8"))
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<string>(func_name,recordingtype,
			"(valid arguments are \"uchar1\" or \"dbl8\")");
	}

	//Initialize movie frame variables (every node does these, including node 0)

	//first,determine whether there are N or N+1 field components in a certain direction, where N is the number of cells in that direction.
	//for example, X1_extra=1 if there are N+1 field components in the 1st dimension in the 2D frame, and =0 if there are N.
	X1_extra = 0;
	X2_extra = 0;
	if (component=="E")
	{
		X1_extra = 0;
		X2_extra = 0;
	}
	else if (component=="Ex")
	{
		RecordedFieldComponent = &Ex;
		RecordedMaterialIndex = &Media_Ex;
		RecordedPermittivityComponent = &eps_x;
		RecordedConductivityComponent = &cond_e_x;
		if (section=="yz")
		{
			X1_extra = 1;
			X2_extra = 1;
		}
		else if (section=="xz")
		{
			X1_extra = 0;
			X2_extra = 1;
		}
		else if (section=="xy")
		{
			X1_extra = 0;
			X2_extra = 1;
		}
	}
	else if (component=="Ey")
	{
		RecordedFieldComponent = &Ey;
		RecordedMaterialIndex = &Media_Ey;
		RecordedPermittivityComponent = &eps_y;
		RecordedConductivityComponent = &cond_e_y;
		if (section=="yz")
		{
			X1_extra = 0;
			X2_extra = 1;
		}
		else if (section=="xz")
		{
			X1_extra = 1;
			X2_extra = 1;
		}
		else if (section=="xy")
		{
			X1_extra = 1;
			X2_extra = 0;
		}
	}
	else if (component=="Ez")
	{
		RecordedFieldComponent = &Ez;
		RecordedMaterialIndex = &Media_Ez;
		RecordedPermittivityComponent = &eps_z;
		RecordedConductivityComponent = &cond_e_z;
		if (section=="yz")
		{
			X1_extra = 1;
			X2_extra = 0;
		}
		else if (section=="xz")
		{
			X1_extra = 1;
			X2_extra = 0;
		}
		else if (section=="xy")
		{
			X1_extra = 1;
			X2_extra = 1;
		}
	}
	else
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<string>(func_name,component,
			"(valid arguments are \"E\", \"Ex\", \"Ey\", or \"Ez\")");
	}

	//determine the total frame sizes, partial frame limits for each node, and whether the frame passes through the node
	PassesThroughNode = false;
	if (section=="yz")
	{
		OriginX1 = OriginY;
		OriginX2 = OriginZ;
		movie_nodes_X1 = nodes_y;
		movie_nodes_X2 = nodes_z;
		TotalFrameSize_X1 = NCELLS_Y + 2*NPML + X1_extra;
		TotalFrameSize_X2 = NCELLS_Z + 2*NPML + X2_extra;
		FrameMin_X1 = jleft;
		FrameMax_X1 = jright + X1_extra;
		FrameMin_X2 = klower;
		FrameMax_X2 = kupper + X2_extra;
		FramePos = relative_FramePos + OriginX;
		if ((iback<=FramePos)&&(ifront>=FramePos))	//does the frame pass through the node?
		{
			PassesThroughNode = true;
		}
	}
	else if (section=="xz")
	{
		OriginX1 = OriginX;
		OriginX2 = OriginZ;
		movie_nodes_X1 = nodes_x;
		movie_nodes_X2 = nodes_z;
		TotalFrameSize_X1 = NCELLS_X + 2*NPML + X1_extra;
		TotalFrameSize_X2 = NCELLS_Z + 2*NPML + X2_extra;
		FrameMin_X1 = iback;
		FrameMax_X1 = ifront + X1_extra;
		FrameMin_X2 = klower;
		FrameMax_X2 = kupper + X2_extra;
		FramePos = relative_FramePos + OriginY;
		if ((jleft<=FramePos)&&(jright>=FramePos))	//does the frame pass through the node?
		{
			PassesThroughNode = true;
		}
	}
	else if (section=="xy")
	{
		OriginX1 = OriginX;
		OriginX2 = OriginY;
		movie_nodes_X1 = nodes_x;
		movie_nodes_X2 = nodes_y;
		TotalFrameSize_X1 = NCELLS_X + 2*NPML + X1_extra;
		TotalFrameSize_X2 = NCELLS_Y + 2*NPML + X2_extra;
		FrameMin_X1 = iback;
		FrameMax_X1 = ifront + X1_extra;
		FrameMin_X2 = jleft;
		FrameMax_X2 = jright + X2_extra;
		FramePos = relative_FramePos + OriginZ;
		if ((klower<=FramePos)&&(kupper>=FramePos))	//does the frame pass through the node?
		{
			PassesThroughNode = true;
		}
	}
	else
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<string>(func_name,section,
			"(valid arguments are \"yz\", \"xz\", or \"xy\")");
	}

	//determine the spatial coordinates (in grid cells) of the field components on the recorded cross section
	X1_coord_range.resize(Range(1,TotalFrameSize_X1));
	X2_coord_range.resize(Range(1,TotalFrameSize_X2));
	for (int x1=1; x1<=TotalFrameSize_X1; x1++)
	{
			X1_coord_range(x1) = (x1-OriginX1+0.5-(X1_extra/2.0))*dx;
	}
	for (int x2=1; x2<=TotalFrameSize_X2; x2++)
	{
			X2_coord_range(x2) = (x2-OriginX2+0.5-(X2_extra/2.0))*dx;
	}

// cout << "I am rank " << rank << " and PassesThroughNode is " << PassesThroughNode << endl;
#ifndef MPI_DISABLE
	//Create movie communicator
	int MovieMembershipKey = MPI_UNDEFINED;
	if (PassesThroughNode)
	{
		MovieMembershipKey = 1;
	}
	int NewRankKey = rank;	//the movie communicator is ranked (in both directions) according to rank in original communicator
	MPI_Comm_split(MPI_CartSubComm,MovieMembershipKey,NewRankKey,&MPI_MovieComm);
	//MPI_MovieComm is :
	//	- a valid communicator if node is in movie frame
	//	- an invalid communicator if node is not in movie frame
	//Therefore, we use the flag "PassesThroughNode" to invoke communication processes, otherwise an error occurs.

	if (PassesThroughNode)
	{
		//Get rank and size in the movie communicator
		MPI_Comm_size(MPI_MovieComm,&movie_size);
		MPI_Comm_rank(MPI_MovieComm,&movie_rank);

#ifndef MPI_USE_PARALLEL_IO
		//Create movie communicator with cartesian topology
		int ndims = 2;
		int dims[2] = {movie_nodes_X1,movie_nodes_X2};
		int periods[2] = {0,0};	//non-periodic
		int reorder = 1;	//permit reorder
		MPI_Cart_create(MPI_MovieComm, ndims, dims, periods, reorder, &MPI_CartMovieComm);
		//Get the coordinates of the partial frame on the frame
		int movie_coordinate[2];
		MPI_Cart_coords(MPI_CartMovieComm,movie_rank,2,movie_coordinate);
		movie_rank_X1 = movie_coordinate[0];
		movie_rank_X2 = movie_coordinate[1];

		//Create two other communicators, for the two directions in the movie frame
		int remain_dims[2];
		remain_dims[0] = true; remain_dims[1] = false;
		MPI_Cart_sub(MPI_CartMovieComm,remain_dims,&MPI_CartMovieComm_X1);
		remain_dims[0] = false; remain_dims[1] = true;
		MPI_Cart_sub(MPI_CartMovieComm,remain_dims,&MPI_CartMovieComm_X2);

		//receive node limits for each node that sends partial frames.
		FrameMinArray_X1.resize(movie_nodes_X1);
		FrameMaxArray_X1.resize(movie_nodes_X1);
		FrameMinArray_X2.resize(movie_nodes_X2);
		FrameMaxArray_X2.resize(movie_nodes_X2);
		int recv_node = 0;	//only the master movie node receives the node limits
		if (movie_rank_X2==0)	//send-receive only within the first X1 column
		{
			MPI_Gather((void*)&FrameMin_X1,1,MPI_INT,FrameMinArray_X1.data(),1,MPI_INT,recv_node,MPI_CartMovieComm_X1);
			MPI_Gather((void*)&FrameMax_X1,1,MPI_INT,FrameMaxArray_X1.data(),1,MPI_INT,recv_node,MPI_CartMovieComm_X1);
		}
		if (movie_rank_X1==0)	//send-receive only within the first X2 column
		{
			MPI_Gather((void*)&FrameMin_X2,1,MPI_INT,FrameMinArray_X2.data(),1,MPI_INT,recv_node,MPI_CartMovieComm_X2);
			MPI_Gather((void*)&FrameMax_X2,1,MPI_INT,FrameMaxArray_X2.data(),1,MPI_INT,recv_node,MPI_CartMovieComm_X2);
		}

		if (movie_rank==0)
		{
			//Only the master node with movie_rank=0 has to do the following tasks, since it builds the movie frame.

			//Initialize the whole movie frame array
			// (both uchar1 and dbl8 arrays are initialized because dbl8 is always used for material recording)
			//unsigned char array:
			uchar1_MovieFrameArray.resize(Range(1,TotalFrameSize_X1),Range(1,TotalFrameSize_X2));	//whole frame
			uchar1_MovieFrameArray = 0;
			//double array:
			dbl8_MovieFrameArray.resize(Range(1,TotalFrameSize_X1),Range(1,TotalFrameSize_X2));	//whole frame
			dbl8_MovieFrameArray = 0;

			//construct the 2-D array of partial buffer arrays
			// (both uchar1 and dbl8 arrays are initialized because dbl8 is always used for material recording)
			//unsigned char arrays:
			uchar1_PartialFrameBuf_recv.resize(movie_nodes_X1,movie_nodes_X2);
			//double arrays:
			dbl8_PartialFrameBuf_recv.resize(movie_nodes_X1,movie_nodes_X2);

			for (n1=0; n1<movie_nodes_X1; n1++)
			{
				for (n2=0; n2<movie_nodes_X2; n2++)
				{
					// (both uchar1 and dbl8 arrays are initialized because dbl8 is always used for material recording)
					//unsigned char array:
					uchar1_PartialFrameBuf_recv(n1,n2).resize
						(Range(FrameMinArray_X1(n1),FrameMaxArray_X1(n1)),Range(FrameMinArray_X2(n2),FrameMaxArray_X2(n2)));
					//double array:
					dbl8_PartialFrameBuf_recv(n1,n2).resize
						(Range(FrameMinArray_X1(n1),FrameMaxArray_X1(n1)),Range(FrameMinArray_X2(n2),FrameMaxArray_X2(n2)));
				}
			}

			//open file
			MovieFile.open(MovieFileName.c_str(),ios::binary);
			if (!MovieFile)
			{
				/** throw exception **/
				if (rank==0) cout << "Error opening movie file " << MovieFileName << endl;
			}
/*			//record the package version
			MovieFile.write((char*)&package_version,sizeof(package_version));
			//record the number of bytes used to represent each field value (ver. >=0.101)
			int num_of_bytes;
			if (recordingtype=="uchar1")
			{// "unsigned char" recording
				num_of_bytes = 1;
			}
			else if (recordingtype=="dbl8")
			{// "double" recording
				num_of_bytes = 8;
			}
			MovieFile.write((char*)&num_of_bytes,sizeof(num_of_bytes));
			//record the grid spacing (in m) (ver. >=0.101)
			MovieFile.write((char*)&dx,sizeof(dx));
			//record the time step (in sec) (ver. >=0.102)
			MovieFile.write((char*)&dt,sizeof(dt));
			//record frame dimensions
			MovieFile.write((char*)&FieldMax,sizeof(FieldMax));
			MovieFile.write((char*)&FieldMin,sizeof(FieldMin));
			MovieFile.write((char*)&TotalFrameSize_X1,sizeof(TotalFrameSize_X1));
			MovieFile.write((char*)&TotalFrameSize_X2,sizeof(TotalFrameSize_X2));
			MovieFile.write((char*)&NMOVIESTEPS,sizeof(NMOVIESTEPS));
			MovieFile.write((char*)&NPML,sizeof(NPML));
			//record the X1 and X2 coordinate ranges of the recorded cross section (ver >=0.105)
			MovieFile.write((char*)X1_coord_range.data(),X1_coord_range.size()*sizeof(X1_coord_range(0)));		//"double" type, of size TotalFrameSize_X1
			MovieFile.write((char*)X2_coord_range.data(),X2_coord_range.size()*sizeof(X2_coord_range(0)));		//"double" type, of size TotalFrameSize_X2*/
			WritePreamble();
		}
		else
		{
			//Initialize the partial movie frame array used for sending data to node 0
			// (both uchar1 and dbl8 arrays are initialized because dbl8 is always used for material recording)
			//unsigned char array:
			uchar1_PartialFrameBuf_send.resize(Range(FrameMin_X1,FrameMax_X1),Range(FrameMin_X2,FrameMax_X2));	//part of the frame
			//double array:
			dbl8_PartialFrameBuf_send.resize(Range(FrameMin_X1,FrameMax_X1),Range(FrameMin_X2,FrameMax_X2));	//part of the frame
		}

		//record the material properties
		//property arrays have the same size as the electric field array that is recorded
		//note that nothing is recorded if electric field is not recorded, in which case RecordedMaterialIndex is NULL
		RecordPermittivity(RecordedMaterialIndex,RecordedPermittivityComponent);
		RecordConductivity(RecordedMaterialIndex,RecordedConductivityComponent);
#endif

#ifdef MPI_USE_PARALLEL_IO
		//Initialize the partial movie frame array used for writing data
		// (both uchar1 and dbl8 arrays are initialized because dbl8 is always used for material recording)
		//unsigned char array:
		uchar1_PartialFrame.resize(Range(FrameMin_X1,FrameMax_X1),Range(FrameMin_X2,FrameMax_X2));
		//double array:
		dbl8_PartialFrame.resize(Range(FrameMin_X1,FrameMax_X1),Range(FrameMin_X2,FrameMax_X2));

		if (movie_rank==0)
		{
			//Only the master node with movie_rank=0 has to do the following tasks, since it builds the movie frame.
			//open file
			MovieFile.open(MovieFileName.c_str(),ios::binary);
			if (!MovieFile)
			{
				/** throw exception **/
				if (rank==0) cout << "Error opening movie file " << MovieFileName << endl;
			}
/*			//record the package version
			MovieFile.write((char*)&package_version,sizeof(package_version));
			//record the number of bytes used to represent each field value (ver. >=0.101)
			int num_of_bytes;
			if (recordingtype=="uchar1")
			{// "unsigned char" recording
				num_of_bytes = 1;
			}
			else if (recordingtype=="dbl8")
			{// "double" recording
				num_of_bytes = 8;
			}
			MovieFile.write((char*)&num_of_bytes,sizeof(num_of_bytes));
			//record the grid spacing (in m) (ver. >=0.101)
			MovieFile.write((char*)&dx,sizeof(dx));
			//record the time step (in sec) (ver. >=0.102)
			MovieFile.write((char*)&dt,sizeof(dt));
			//record frame dimensions
			MovieFile.write((char*)&FieldMax,sizeof(FieldMax));
			MovieFile.write((char*)&FieldMin,sizeof(FieldMin));
			MovieFile.write((char*)&TotalFrameSize_X1,sizeof(TotalFrameSize_X1));
			MovieFile.write((char*)&TotalFrameSize_X2,sizeof(TotalFrameSize_X2));
			MovieFile.write((char*)&NMOVIESTEPS,sizeof(NMOVIESTEPS));
			MovieFile.write((char*)&NPML,sizeof(NPML));
			//record the X1 and X2 coordinate ranges of the recorded cross section (ver >=0.105)
			MovieFile.write((char*)X1_coord_range.data(),X1_coord_range.size()*sizeof(X1_coord_range(0)));		//"double" type, of size TotalFrameSize_X1
			MovieFile.write((char*)X2_coord_range.data(),X2_coord_range.size()*sizeof(X2_coord_range(0)));		//"double" type, of size TotalFrameSize_X2*/
			WritePreamble();

			MovieFile.close();
		}
		//wait until the movie node is finished writing the preamble
//cout << "rank " << rank << " came here in Crecorder_movie" << endl;

		MPI_Barrier(MPI_MovieComm);

		//Open the file for parallel I/O
		int OpenMode = (MPI_MODE_WRONLY|MPI_MODE_APPEND);	//open for writing, and append at the end of file
		MPI_Info info=MPI_INFO_NULL;
		char *MPIMovieFileName = new char[MovieFileName.size()+1];	//one extra place for the null-termination
		strcpy(MPIMovieFileName, MovieFileName.c_str());
		//open the file for writing
		int opencode = MPI_File_open(MPI_MovieComm,MPIMovieFileName,OpenMode,info,&MPIMovieFile);
		delete [] MPIMovieFileName;
		//File opened for parallel I/O.

		//Set view for the current process
		int array_of_sizes[2] = {TotalFrameSize_X1,TotalFrameSize_X2};
		int array_of_subsizes[2] = {FrameMax_X1-FrameMin_X1+1,FrameMax_X2-FrameMin_X2+1};
		int array_of_starts[2] = {FrameMin_X1-1,FrameMin_X2-1};

		//First, set file view for "double" recording of material properties
		// Create filetype for parallel I/O
		MPI_Type_create_subarray(2,array_of_sizes,array_of_subsizes,array_of_starts,MPI_ORDER_C,
									MPI_DOUBLE,&MPIMovieFileType);
		MPI_Type_commit(&MPIMovieFileType);
		//Get the current size of the file in bytes, use it as offset for the file view
		MPI_Offset beginning;
		MPI_File_get_size(MPIMovieFile,&beginning);
		string datarepstr = "native";
		char *datarep = new char[datarepstr.size()+1];	//one extra place for the null-termination
		strcpy(datarep, datarepstr.c_str());
		MPI_File_set_view(MPIMovieFile,beginning,MPI_DOUBLE,MPIMovieFileType,datarep,info);
		//record the material properties using parallel I/O
		RecordPermittivity(RecordedMaterialIndex,RecordedPermittivityComponent);
		RecordConductivity(RecordedMaterialIndex,RecordedConductivityComponent);


		//When material property recording is finished, set file type for recording of field values
		if (recordingtype=="uchar1")
		{// "unsigned char" recording
			// Create filetype for parallel I/O
			MPI_Type_create_subarray(2,array_of_sizes,array_of_subsizes,array_of_starts,MPI_ORDER_C,
										MPI_UNSIGNED_CHAR,&MPIMovieFileType);
			MPI_Type_commit(&MPIMovieFileType);
			//Get the current size of the file in bytes, use it as offset for the file view
			MPI_File_get_size(MPIMovieFile,&beginning);
			MPI_File_set_view(MPIMovieFile,beginning,MPI_UNSIGNED_CHAR,MPIMovieFileType,datarep,info);
		}
		else if (recordingtype=="dbl8")
		{// "double" recording
			// Create filetype for parallel I/O
			MPI_Type_create_subarray(2,array_of_sizes,array_of_subsizes,array_of_starts,MPI_ORDER_C,
										MPI_DOUBLE,&MPIMovieFileType);
			MPI_Type_commit(&MPIMovieFileType);
			//Get the current size of the file in bytes, use it as offset for the file view
			MPI_File_get_size(MPIMovieFile,&beginning);
			MPI_File_set_view(MPIMovieFile,beginning,MPI_DOUBLE,MPIMovieFileType,datarep,info);
		}
#endif
	}
#else //MPI_DISABLE
	//no MPI, only one node
	movie_size = 1;
	movie_rank = 0;
	movie_rank_X1 = 0;
	movie_rank_X2 = 0;

	//Initialize the whole movie frame array
	// (both uchar1 and dbl8 arrays are initialized because dbl8 is always used for material recording)
	//unsigned char array:
	uchar1_MovieFrameArray.resize(Range(1,TotalFrameSize_X1),Range(1,TotalFrameSize_X2));	//whole frame
	uchar1_MovieFrameArray = 0;
	//double array:
	dbl8_MovieFrameArray.resize(Range(1,TotalFrameSize_X1),Range(1,TotalFrameSize_X2));	//whole frame
	dbl8_MovieFrameArray = 0;

	//open file
	MovieFile.open(MovieFileName.c_str(),ios::binary);
	if (!MovieFile)
	{
		/** throw exception **/
		if (rank==0) cout << "Error opening movie file " << MovieFileName << endl;
	}
/*	//record the package version
	MovieFile.write((char*)&package_version,sizeof(package_version));
	//record the number of bytes used to represent each field value (ver. >=0.101)
	int num_of_bytes;
	if (recordingtype=="uchar1")
	{// "unsigned char" recording
		num_of_bytes = 1;
	}
	else if (recordingtype=="dbl8")
	{// "double" recording
		num_of_bytes = 8;
	}
	MovieFile.write((char*)&num_of_bytes,sizeof(num_of_bytes));
	//record the grid spacing (in m) (ver. >=0.101)
	MovieFile.write((char*)&dx,sizeof(dx));
	//record the time step (in sec) (ver. >=0.102)
	MovieFile.write((char*)&dt,sizeof(dt));
	//record frame dimensions
	MovieFile.write((char*)&FieldMax,sizeof(FieldMax));
	MovieFile.write((char*)&FieldMin,sizeof(FieldMin));
	MovieFile.write((char*)&TotalFrameSize_X1,sizeof(TotalFrameSize_X1));
	MovieFile.write((char*)&TotalFrameSize_X2,sizeof(TotalFrameSize_X2));
	MovieFile.write((char*)&NMOVIESTEPS,sizeof(NMOVIESTEPS));
	MovieFile.write((char*)&NPML,sizeof(NPML));
	//record the X1 and X2 coordinate ranges of the recorded cross section (ver >=0.105)
	MovieFile.write((char*)X1_coord_range.data(),X1_coord_range.size()*sizeof(X1_coord_range(0)));		//"double" type, of size TotalFrameSize_X1
	MovieFile.write((char*)X2_coord_range.data(),X2_coord_range.size()*sizeof(X2_coord_range(0)));		//"double" type, of size TotalFrameSize_X2*/
	WritePreamble();

	//record the material properties
	//property arrays have the same size as the electric field array that is recorded
	//note that nothing is recorded if electric field is not recorded, in which case RecordedMaterialIndex is NULL
	RecordPermittivity(RecordedMaterialIndex,RecordedPermittivityComponent);
	RecordConductivity(RecordedMaterialIndex,RecordedConductivityComponent);
#endif //MPI_DISABLE
}

void Crecorder_movie::RecordMovieFrame()
{
	if (!OnlyRecordGeometry)//if only the geometry is recorded, do nothing
	{
		if (recordingtype=="uchar1")
		{// "unsigned char" recording
			RecordMovieFrame_uchar1(RecordedFieldComponent);
		}
		else if (recordingtype=="dbl8")
		{// "double" recording
			RecordMovieFrame_dbl8(RecordedFieldComponent);
		}
	}
}

void Crecorder_movie::RecordMovieFrame_uchar1(const Array<double,3>* RecordedArray)
{//records a single frame from the specified array in "unsigned char" (1 byte) type (used for recording field components, yields faster access but discretization is necessary)
#ifndef MPI_DISABLE
	if (PassesThroughNode)
	{
#ifndef MPI_USE_PARALLEL_IO
		if (movie_rank==0)
		{
			//receive data from the nodes
			for (n1=0; n1<movie_nodes_X1; n1++)
			{
				for (n2=0; n2<movie_nodes_X2; n2++)
				{
					if ((n1==0)&&(n2==0))
					//do NOT try to receive data from node 0 (yourself), place available data
					{
						PlaceFrame(RecordedArray,uchar1_MovieFrameArray);
					}
					else
					//receive data from other nodes
					{
						recv_coord[0] = n1; recv_coord[1] = n2;
						MPI_Cart_rank(MPI_CartMovieComm,recv_coord,&recv_rank);
						MPI_Recv(uchar1_PartialFrameBuf_recv(n1,n2).data(),uchar1_PartialFrameBuf_recv(n1,n2).size(),MPI_UNSIGNED_CHAR,
							recv_rank,0,MPI_CartMovieComm,&Status);
						//place buffer into place
						uchar1_MovieFrameArray(Range(FrameMinArray_X1(n1),FrameMaxArray_X1(n1)),
							Range(FrameMinArray_X2(n2),FrameMaxArray_X2(n2)))
							=uchar1_PartialFrameBuf_recv(n1,n2);
					}
				}
			}
			//write out the data
			//blitz++ does row-major ordering, therefore the 2nd dimension gets written first
			MovieFile.write((char*)uchar1_MovieFrameArray.data(),uchar1_MovieFrameArray.size()*sizeof(uchar1_MovieFrameArray(0)));
		}
		else
		{
			//send data to node 0
			PlaceFrame(RecordedArray,uchar1_PartialFrameBuf_send);
			MPI_Send(uchar1_PartialFrameBuf_send.data(),uchar1_PartialFrameBuf_send.size(),MPI_UNSIGNED_CHAR,0,0,MPI_CartMovieComm);
		}
#endif

#ifdef MPI_USE_PARALLEL_IO
		PlaceFrame(RecordedArray,uchar1_PartialFrame);
		//write out the data
		//blitz++ does row-major ordering, therefore the 2nd dimension gets written first
		MPI_File_write(MPIMovieFile,uchar1_PartialFrame.data(),uchar1_PartialFrame.size(),MPI_UNSIGNED_CHAR,&Status);
#endif
	}
#else //MPI_DISABLE
	//no MPI, only one node
	PlaceFrame(RecordedArray,uchar1_MovieFrameArray);
	MovieFile.write((char*)uchar1_MovieFrameArray.data(),uchar1_MovieFrameArray.size()*sizeof(uchar1_MovieFrameArray(0)));
#endif //MPI_DISABLE
}

void Crecorder_movie::RecordMovieFrame_dbl8(const Array<double,3>* RecordedArray)
{//records a single frame from the specified array in "double" (8 byte) type
#ifndef MPI_DISABLE
	if (PassesThroughNode)
	{
#ifndef MPI_USE_PARALLEL_IO
		if (movie_rank==0)
		{
			//receive data from the nodes
			for (n1=0; n1<movie_nodes_X1; n1++)
			{
				for (n2=0; n2<movie_nodes_X2; n2++)
				{
					if ((n1==0)&&(n2==0))
					//do NOT try to receive data from node 0 (yourself), place available data
					{
						PlaceFrame(RecordedArray,dbl8_MovieFrameArray);
					}
					else
					//receive data from other nodes
					{
						recv_coord[0] = n1; recv_coord[1] = n2;
						MPI_Cart_rank(MPI_CartMovieComm,recv_coord,&recv_rank);
						MPI_Recv(dbl8_PartialFrameBuf_recv(n1,n2).data(),dbl8_PartialFrameBuf_recv(n1,n2).size(),MPI_DOUBLE,
							recv_rank,0,MPI_CartMovieComm,&Status);
						//place buffer into place
						dbl8_MovieFrameArray(Range(FrameMinArray_X1(n1),FrameMaxArray_X1(n1)),
							Range(FrameMinArray_X2(n2),FrameMaxArray_X2(n2)))
							=dbl8_PartialFrameBuf_recv(n1,n2);
					}
				}
			}
			//write out the data
			//blitz++ does row-major ordering, therefore the 2nd dimension gets written first
			MovieFile.write((char*)dbl8_MovieFrameArray.data(),dbl8_MovieFrameArray.size()*sizeof(dbl8_MovieFrameArray(0)));
		}
		else
		{
			//send data to node 0
			PlaceFrame(RecordedArray,dbl8_PartialFrameBuf_send);
			MPI_Send(dbl8_PartialFrameBuf_send.data(),dbl8_PartialFrameBuf_send.size(),MPI_DOUBLE,0,0,MPI_CartMovieComm);
		}
#endif

#ifdef MPI_USE_PARALLEL_IO
		PlaceFrame(RecordedArray,dbl8_PartialFrame);
		//write out the data
		//blitz++ does row-major ordering, therefore the 2nd dimension gets written first
		MPI_File_write(MPIMovieFile,dbl8_PartialFrame.data(),dbl8_PartialFrame.size(),MPI_DOUBLE,&Status);
#endif
	}
#else //MPI_DISABLE
	//no MPI, only one node
	PlaceFrame(RecordedArray,dbl8_MovieFrameArray);
	MovieFile.write((char*)dbl8_MovieFrameArray.data(),dbl8_MovieFrameArray.size()*sizeof(dbl8_MovieFrameArray(0)));
#endif //MPI_DISABLE
}

void Crecorder_movie::RecordPermittivity(const Array<ElectricMaterialIndexType_X,3>* PropertyIndexArray, const Array<double,1>* PermittivityArray)
{
#ifndef MPI_DISABLE
	if (PassesThroughNode)
	{
#ifndef MPI_USE_PARALLEL_IO
		if (movie_rank==0)
		{
			//receive data from the nodes
			for (n1=0; n1<movie_nodes_X1; n1++)
			{
				for (n2=0; n2<movie_nodes_X2; n2++)
				{
					if ((n1==0)&&(n2==0))
					//do NOT try to receive data from node 0 (yourself), place available data
					{
						PlacePermittivityFrame(PropertyIndexArray,PermittivityArray,dbl8_MovieFrameArray);
					}
					else
					//receive data from other nodes
					{
						recv_coord[0] = n1; recv_coord[1] = n2;
						MPI_Cart_rank(MPI_CartMovieComm,recv_coord,&recv_rank);
						MPI_Recv(dbl8_PartialFrameBuf_recv(n1,n2).data(),dbl8_PartialFrameBuf_recv(n1,n2).size(),MPI_DOUBLE,
							recv_rank,0,MPI_CartMovieComm,&Status);
						//place buffer into place
						dbl8_MovieFrameArray(Range(FrameMinArray_X1(n1),FrameMaxArray_X1(n1)),
							Range(FrameMinArray_X2(n2),FrameMaxArray_X2(n2)))
							=dbl8_PartialFrameBuf_recv(n1,n2);
					}
				}
			}
			//write out the data
			//blitz++ does row-major ordering, therefore the 2nd dimension gets written first
			MovieFile.write((char*)dbl8_MovieFrameArray.data(),dbl8_MovieFrameArray.size()*sizeof(dbl8_MovieFrameArray(0)));
		}
		else
		{
			//send data to node 0
			PlacePermittivityFrame(PropertyIndexArray,PermittivityArray,dbl8_PartialFrameBuf_send);
			MPI_Send(dbl8_PartialFrameBuf_send.data(),dbl8_PartialFrameBuf_send.size(),MPI_DOUBLE,0,0,MPI_CartMovieComm);
		}
#endif

#ifdef MPI_USE_PARALLEL_IO
		PlacePermittivityFrame(PropertyIndexArray,PermittivityArray,dbl8_PartialFrame);
		//write out the data
		//blitz++ does row-major ordering, therefore the 2nd dimension gets written first
		//use collective write to ensure that all nodes write their material properties
		MPI_File_write_all(MPIMovieFile,dbl8_PartialFrame.data(),dbl8_PartialFrame.size(),MPI_DOUBLE,&Status);
#endif
	}
#else //MPI_DISABLE
	//no MPI, only one node
	PlacePermittivityFrame(PropertyIndexArray,PermittivityArray,dbl8_MovieFrameArray);
	MovieFile.write((char*)dbl8_MovieFrameArray.data(),dbl8_MovieFrameArray.size()*sizeof(dbl8_MovieFrameArray(0)));
#endif //MPI_DISABLE
}

void Crecorder_movie::RecordConductivity(const Array<ElectricMaterialIndexType_X,3>* PropertyIndexArray, const Array<double,1>* ConductivityArray)
{
#ifndef MPI_DISABLE
	if (PassesThroughNode)
	{
#ifndef MPI_USE_PARALLEL_IO
		if (movie_rank==0)
		{
			//receive data from the nodes
			for (n1=0; n1<movie_nodes_X1; n1++)
			{
				for (n2=0; n2<movie_nodes_X2; n2++)
				{
					if ((n1==0)&&(n2==0))
					//do NOT try to receive data from node 0 (yourself), place available data
					{
						PlaceConductivityFrame(PropertyIndexArray,ConductivityArray,dbl8_MovieFrameArray);
					}
					else
					//receive data from other nodes
					{
						recv_coord[0] = n1; recv_coord[1] = n2;
						MPI_Cart_rank(MPI_CartMovieComm,recv_coord,&recv_rank);
						MPI_Recv(dbl8_PartialFrameBuf_recv(n1,n2).data(),dbl8_PartialFrameBuf_recv(n1,n2).size(),MPI_DOUBLE,
							recv_rank,0,MPI_CartMovieComm,&Status);
						//place buffer into place
						dbl8_MovieFrameArray(Range(FrameMinArray_X1(n1),FrameMaxArray_X1(n1)),
							Range(FrameMinArray_X2(n2),FrameMaxArray_X2(n2)))
							=dbl8_PartialFrameBuf_recv(n1,n2);
					}
				}
			}
			//write out the data
			//blitz++ does row-major ordering, therefore the 2nd dimension gets written first
			MovieFile.write((char*)dbl8_MovieFrameArray.data(),dbl8_MovieFrameArray.size()*sizeof(dbl8_MovieFrameArray(0)));
		}
		else
		{
			//send data to node 0
			PlaceConductivityFrame(PropertyIndexArray,ConductivityArray,dbl8_PartialFrameBuf_send);
			MPI_Send(dbl8_PartialFrameBuf_send.data(),dbl8_PartialFrameBuf_send.size(),MPI_DOUBLE,0,0,MPI_CartMovieComm);
		}
#endif

#ifdef MPI_USE_PARALLEL_IO
		PlaceConductivityFrame(PropertyIndexArray,ConductivityArray,dbl8_PartialFrame);
		//write out the data
		//blitz++ does row-major ordering, therefore the 2nd dimension gets written first
		//use collective write to ensure that all nodes write their material properties
		MPI_File_write_all(MPIMovieFile,dbl8_PartialFrame.data(),dbl8_PartialFrame.size(),MPI_DOUBLE,&Status);
#endif
	}
#else //MPI_DISABLE
	//no MPI, only one node
	PlaceConductivityFrame(PropertyIndexArray,ConductivityArray,dbl8_MovieFrameArray);
	MovieFile.write((char*)dbl8_MovieFrameArray.data(),dbl8_MovieFrameArray.size()*sizeof(dbl8_MovieFrameArray(0)));
#endif //MPI_DISABLE
}

void Crecorder_movie::PlaceFrame(const Array<double,3>* RecordedArray, Array<unsigned char,2> DestinationArray)
{//for "unsigned char" (1 byte) recording
//places the partial frame in DestinationArray, which may be the same size as the partial frame (uchar1_PartialFrameBuf_send for rank!=0)
// or the same size as the whole cross section (uchar1_MovieFrameArray for rank=0).
	for (x1=FrameMin_X1; x1<=FrameMax_X1; x1++)
	{
		for (x2=FrameMin_X2; x2<=FrameMax_X2; x2++)
		{
			//find the position on the main grid, depending on the cross-section
			if (section=="yz")
			{
				xPos = FramePos;
				yPos = x1;
				zPos = x2;
			}
			else if (section=="xz")
			{
				xPos = x1;
				yPos = FramePos;
				zPos = x2;
			}
			else if (section=="xy")
			{
				xPos = x1;
				yPos = x2;
				zPos = FramePos;
			}

			//take the E-field component at the position found above
			if (component=="E")
			{
				field_value = sqrt(pow((Ex(xPos,yPos,zPos)+Ex(xPos,yPos+1,zPos)
										+Ex(xPos,yPos,zPos+1)+Ex(xPos,yPos+1,zPos+1))/4.0,2)
								+pow((Ey(xPos,yPos,zPos)+Ey(xPos+1,yPos,zPos)
										+Ey(xPos,yPos,zPos+1)+Ey(xPos+1,yPos,zPos+1))/4.0,2)
								+pow((Ez(xPos,yPos,zPos)+Ez(xPos,yPos+1,zPos)
										+Ez(xPos+1,yPos,zPos)+Ez(xPos+1,yPos+1,zPos))/4.0,2));
			}
			else
			{
				field_value = (*RecordedArray)(xPos,yPos,zPos);
			}

			//if necessary, convert to dB
			if (scale=="dB")
			{
				field_value = 20.0*log10(abs(field_value)+1e-20);
			}
			else if (scale=="absolute")
			{
				field_value = abs(field_value);
			}
			//clip, discretize and record the field value
			if (field_value>=FieldMax)
			{
				DestinationArray(x1,x2) = 255;
			}
			else if (field_value<=FieldMin)
			{
				DestinationArray(x1,x2) = 0;
			}
			else
			{
				DestinationArray(x1,x2) = (unsigned char)
					(255/(FieldMax-FieldMin)*(field_value-FieldMin)+.5);
			}
		}
	}
}

void Crecorder_movie::PlaceFrame(const Array<double,3>* RecordedArray, Array<double,2> DestinationArray)
{//for "double" (8 byte) recording
//places the partial frame in DestinationArray, which may be the same size as the partial frame (dbl8_PartialFrameBuf_send for rank!=0)
// or the same size as the whole cross section (dbl8_MovieFrameArray for rank=0).
	for (x1=FrameMin_X1; x1<=FrameMax_X1; x1++)
	{
		for (x2=FrameMin_X2; x2<=FrameMax_X2; x2++)
		{
			//find the position on the main grid, depending on the cross-section
			if (section=="yz")
			{
				xPos = FramePos;
				yPos = x1;
				zPos = x2;
			}
			else if (section=="xz")
			{
				xPos = x1;
				yPos = FramePos;
				zPos = x2;
			}
			else if (section=="xy")
			{
				xPos = x1;
				yPos = x2;
				zPos = FramePos;
			}

			//take the E-field component at the position found above
			if (component=="E")
			{
				field_value = sqrt(pow((Ex(xPos,yPos,zPos)+Ex(xPos,yPos+1,zPos)
										+Ex(xPos,yPos,zPos+1)+Ex(xPos,yPos+1,zPos+1))/4.0,2)
								+pow((Ey(xPos,yPos,zPos)+Ey(xPos+1,yPos,zPos)
										+Ey(xPos,yPos,zPos+1)+Ey(xPos+1,yPos,zPos+1))/4.0,2)
								+pow((Ez(xPos,yPos,zPos)+Ez(xPos,yPos+1,zPos)
										+Ez(xPos+1,yPos,zPos)+Ez(xPos+1,yPos+1,zPos))/4.0,2));
			}
			else
			{
				field_value = (*RecordedArray)(xPos,yPos,zPos);
			}

			//if necessary, convert to dB
			if (scale=="dB")
			{
				field_value = 20.0*log10(abs(field_value)+1e-20);
			}
			else if (scale=="absolute")
			{
				field_value = abs(field_value);
			}
			//no discretization necessary, since field value is recorded in "double" type
			DestinationArray(x1,x2) = field_value;
		}
	}
}

void Crecorder_movie::PlacePermittivityFrame(const Array<ElectricMaterialIndexType_X,3>* PropertyIndexArray, const Array<double,1>* PermittivityArray, Array<double,2> DestinationArray)
{
	for (x1=FrameMin_X1; x1<=FrameMax_X1; x1++)
	{
		for (x2=FrameMin_X2; x2<=FrameMax_X2; x2++)
		{
			//find the position on the main grid, depending on the cross-section
			if (section=="yz")
			{
				xPos = FramePos;
				yPos = x1;
				zPos = x2;
			}
			else if (section=="xz")
			{
				xPos = x1;
				yPos = FramePos;
				zPos = x2;
			}
			else if (section=="xy")
			{
				xPos = x1;
				yPos = x2;
				zPos = FramePos;
			}

			//take the E-field component at the position found above
			if (component=="E")
			{
				DestinationArray(x1,x2) = (eps_x(Media_Ex(xPos,yPos,zPos))+eps_x(Media_Ex(xPos,yPos+1,zPos))
										+eps_x(Media_Ex(xPos,yPos,zPos+1))+eps_x(Media_Ex(xPos,yPos+1,zPos+1))
								+eps_y(Media_Ey(xPos,yPos,zPos))+eps_y(Media_Ey(xPos+1,yPos,zPos))
										+eps_y(Media_Ey(xPos,yPos,zPos+1))+eps_y(Media_Ey(xPos+1,yPos,zPos+1))
								+eps_z(Media_Ez(xPos,yPos,zPos))+eps_z(Media_Ez(xPos,yPos+1,zPos))
										+eps_z(Media_Ez(xPos+1,yPos,zPos))+eps_z(Media_Ez(xPos+1,yPos+1,zPos)))/12.0;
								//doesn't really make much sense if the material is anisotropic, but that's for later
			}
			else
			{
				DestinationArray(x1,x2) = (*PermittivityArray)((*PropertyIndexArray)(xPos,yPos,zPos));
			}
		}
	}
}

void Crecorder_movie::PlaceConductivityFrame(const Array<ElectricMaterialIndexType_X,3>* PropertyIndexArray, const Array<double,1>* ConductivityArray, Array<double,2> DestinationArray)
{
	for (x1=FrameMin_X1; x1<=FrameMax_X1; x1++)
	{
		for (x2=FrameMin_X2; x2<=FrameMax_X2; x2++)
		{
			//find the position on the main grid, depending on the cross-section
			if (section=="yz")
			{
				xPos = FramePos;
				yPos = x1;
				zPos = x2;
			}
			else if (section=="xz")
			{
				xPos = x1;
				yPos = FramePos;
				zPos = x2;
			}
			else if (section=="xy")
			{
				xPos = x1;
				yPos = x2;
				zPos = FramePos;
			}

			//take the E-field component at the position found above
			if (component=="E")
			{
				DestinationArray(x1,x2) = (cond_e_x(Media_Ex(xPos,yPos,zPos))+cond_e_x(Media_Ex(xPos,yPos+1,zPos))
										+cond_e_x(Media_Ex(xPos,yPos,zPos+1))+cond_e_x(Media_Ex(xPos,yPos+1,zPos+1))
								+cond_e_y(Media_Ey(xPos,yPos,zPos))+cond_e_y(Media_Ey(xPos+1,yPos,zPos))
										+cond_e_y(Media_Ey(xPos,yPos,zPos+1))+cond_e_y(Media_Ey(xPos+1,yPos,zPos+1))
								+cond_e_z(Media_Ez(xPos,yPos,zPos))+cond_e_z(Media_Ez(xPos,yPos+1,zPos))
										+cond_e_z(Media_Ez(xPos+1,yPos,zPos))+cond_e_z(Media_Ez(xPos+1,yPos+1,zPos)))/12.0;
								//doesn't really make much sense if the material is anisotropic, but that's for later
			}
			else
			{
				DestinationArray(x1,x2) = (*ConductivityArray)((*PropertyIndexArray)(xPos,yPos,zPos));
			}
		}
	}
}

void Crecorder_movie::WritePreamble()
{//writes the preamble to the output file
	//record the package version
//	MovieFile.write((char*)&package_version,sizeof(package_version));
	MovieFile.write((char*)&angora_version_major,sizeof(angora_version_major));
	MovieFile.write((char*)&angora_version_minor,sizeof(angora_version_minor));
	MovieFile.write((char*)&angora_version_revision,sizeof(angora_version_revision));
	//record the number of bytes used to represent each field value (ver. >=0.101)
	int num_of_bytes;
	if (recordingtype=="uchar1")
	{// "unsigned char" recording
		num_of_bytes = 1;
	}
	else if (recordingtype=="dbl8")
	{// "double" recording
		num_of_bytes = 8;
	}
	MovieFile.write((char*)&num_of_bytes,sizeof(num_of_bytes));
	//record the grid spacing (in m) (ver. >=0.101)
	MovieFile.write((char*)&dx,sizeof(dx));
	//record the time step (in sec) (ver. >=0.102)
	MovieFile.write((char*)&dt,sizeof(dt));
	//record the initial time value (ver. >=0.12.0)
	double initial_time_value = get_initial_time_value();
	MovieFile.write((char*)&initial_time_value,sizeof(initial_time_value));
	//record frame dimensions
	MovieFile.write((char*)&FieldMax,sizeof(FieldMax));
	MovieFile.write((char*)&FieldMin,sizeof(FieldMin));
	MovieFile.write((char*)&TotalFrameSize_X1,sizeof(TotalFrameSize_X1));
	MovieFile.write((char*)&TotalFrameSize_X2,sizeof(TotalFrameSize_X2));
	MovieFile.write((char*)&NMOVIESTEPS,sizeof(NMOVIESTEPS));
	MovieFile.write((char*)&NPML,sizeof(NPML));
	//record the X1 and X2 coordinate ranges of the recorded cross section (ver >=0.105)
	MovieFile.write((char*)X1_coord_range.data(),X1_coord_range.size()*sizeof(X1_coord_range(0)));		//"double" type, of size TotalFrameSize_X1
	MovieFile.write((char*)X2_coord_range.data(),X2_coord_range.size()*sizeof(X2_coord_range(0)));		//"double" type, of size TotalFrameSize_X2
}

#ifndef MPI_DISABLE
#ifdef MPI_USE_PARALLEL_IO
Crecorder_movie::~Crecorder_movie()
{
	if (PassesThroughNode)
	{
		MPI_File_close(&MPIMovieFile);	//close the MPI I/O file, since it is not done automatically
	}
}
#endif
#endif
