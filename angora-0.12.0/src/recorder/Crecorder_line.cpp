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

//Defines the line recorder class "Crecorder_line"

#include "headers.h"

#include "Crecorder_line.h"

//routines for getting/setting the initial time value in the simulation
#include "time_axis.h"

extern int angora_version_major,angora_version_minor,angora_version_revision;

extern double dt;

extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;
extern int OriginX,OriginY,OriginZ;
extern int NSTEPS;

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

extern double max_field_value;
extern double dB_accuracy;


Crecorder_line::Crecorder_line(const string& myorientation, const int& x1Pos, const int& x2Pos,
			const string& mycomponent, const string& myscale, const string& LineFileName,
			const int& Index)
		 : orientation(myorientation), relative_LinePos_x1(x1Pos), relative_LinePos_x2(x2Pos), component(mycomponent), scale(myscale),
		 LineRecorderIndex(Index)
{//constructor for the line recorder
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

	//Initialize line recorder variables (every node does these, including node 0)

	//first,determine whether there are N or N+1 field components along the orientation of the line, where N is the number of cells in that direction.
	//for example, X_extra=1 if there are N+1 field components along the orientation of the line, and =0 if there are N.
	X_extra = 0;
	if ((component!="E")&&(component!="Ex")&&(component!="Ey")&&(component!="Ez"))
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//			InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<string>(func_name,component,
			"(valid arguments are \"E\", \"Ex\", \"Ey\", or \"Ez\")");
	}
	if ((orientation!="x_directed")&&(orientation!="y_directed")&&(orientation!="z_directed"))
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<string>(func_name,orientation,
			"(valid arguments are \"x_directed\", \"y_directed\", or \"z_directed\")");
	}
	if ((component=="E")||
		((component=="Ex")&&(orientation=="x_directed"))||((component=="Ey")&&(orientation=="y_directed"))||((component=="Ez")&&(orientation=="z_directed")))
	{
		X_extra = 0;
	}
	else
	{
		X_extra = 1;
	}

	//determine the total line length, partial line lengths for each node, and whether the line passes through the node
	PassesThroughNode = false;
	if (orientation=="x_directed")
	{
		line_nodes = nodes_x;
		TotalLineLength = NCELLS_X + 2*NPML + X_extra;
		LineMin = iback;
		LineMax = ifront + X_extra;
		LinePos_x1 = relative_LinePos_x1 + OriginY;
		LinePos_x2 = relative_LinePos_x2 + OriginZ;
		if ((LinePos_x1<=jright)&&(LinePos_x1>=jleft)
		&&(LinePos_x2<=kupper)&&(LinePos_x2>=klower))	//does the line pass through the node?
		{
			PassesThroughNode = true;
		}
	}
	else if (orientation=="y_directed")
	{
		line_nodes = nodes_y;
		TotalLineLength = NCELLS_Y + 2*NPML + X_extra;
		LineMin = jleft;
		LineMax = jright + X_extra;
		LinePos_x1 = relative_LinePos_x1 + OriginX;
		LinePos_x2 = relative_LinePos_x2 + OriginZ;
		if ((LinePos_x1<=ifront)&&(LinePos_x1>=iback)
		&&(LinePos_x2<=kupper)&&(LinePos_x2>=klower))	//does the line pass through the node?
		{
			PassesThroughNode = true;
		}
	}
	else if (orientation=="z_directed")
	{
		line_nodes = nodes_z;
		TotalLineLength = NCELLS_Z + 2*NPML + X_extra;
		LineMin = klower;
		LineMax = kupper + X_extra;
		LinePos_x1 = relative_LinePos_x1 + OriginX;
		LinePos_x2 = relative_LinePos_x2 + OriginY;
		if ((LinePos_x1<=ifront)&&(LinePos_x1>=iback)
		&&(LinePos_x2<=jright)&&(LinePos_x2>=jleft))	//does the line pass through the node?
		{
			PassesThroughNode = true;
		}
	}

	//Create line recorder communicator
#ifndef MPI_DISABLE
	int LineMembershipKey = MPI_UNDEFINED;
	if (PassesThroughNode)
	{
		LineMembershipKey = 1;
	}
	int NewRankKey = rank;	//the line communicator is ranked according to rank in original communicator
	MPI_Comm_split(MPI_CartSubComm,LineMembershipKey,NewRankKey,&MPI_LineComm);
	//MPI_LineComm is :
	//	- a valid communicator if the line passes through node
	//	- an invalid communicator if the line does not pass through node
	//Therefore, we use the flag "PassesThroughNode" to invoke communication processes, otherwise an error occurs.
#endif

	if (PassesThroughNode)
	{
		//Get rank and size in the line communicator
#ifndef MPI_DISABLE
		MPI_Comm_size(MPI_LineComm,&line_size);
		MPI_Comm_rank(MPI_LineComm,&line_rank);
#else
		//no MPI, only one node
		line_size = 1;
		line_rank = 0;
#endif

		//receive node limits for each node that sends partial lines.
		LineMinArray.resize(line_nodes);
		LineMaxArray.resize(line_nodes);
#ifndef MPI_DISABLE
		int recv_node = 0;	//only the master node receives the node limits
		MPI_Gather((void*)&LineMin,1,MPI_INT,LineMinArray.data(),1,MPI_INT,recv_node,MPI_LineComm);
		MPI_Gather((void*)&LineMax,1,MPI_INT,LineMaxArray.data(),1,MPI_INT,recv_node,MPI_LineComm);
/*		if (line_rank==recv_node)
		{
			cout << "My rank is " << rank << endl;
			cout << LineMinArray << endl;
			cout << LineMaxArray << endl;
		}*/
// 		cout << "I am rank " << line_rank << " with min: " << LineMin << ", max: " << LineMax << endl;
#else
		//no MPI: node limits are trivial
		LineMinArray(0) = LineMin;
		LineMaxArray(0) = LineMax;
#endif

		if (line_rank==0)
		{
			//Only the master node with line_rank=0 has to do the following tasks, since it builds the whole line.
			//open file
			LineFile.open(LineFileName.c_str(),ios::binary);
			if (!LineFile)
			{
				/** Throw exception **/
			}
			//record the package version
//			LineFile.write((char*)&package_version,sizeof(package_version));
			LineFile.write((char*)&angora_version_major,sizeof(angora_version_major));
			LineFile.write((char*)&angora_version_minor,sizeof(angora_version_minor));
			LineFile.write((char*)&angora_version_revision,sizeof(angora_version_revision));
			//record the time step (in sec) (ver. >=0.12.0)
			LineFile.write((char*)&dt,sizeof(dt));
			//record the initial time value (ver. >=0.12.0)
			double initial_time_value = get_initial_time_value();
			LineFile.write((char*)&initial_time_value,sizeof(initial_time_value));
			//record line dimensions
//			LineFile.write((char*)&FieldMax,sizeof(FieldMax));
//			LineFile.write((char*)&FieldMin,sizeof(FieldMin));
			LineFile.write((char*)&TotalLineLength,sizeof(TotalLineLength));
			LineFile.write((char*)&NSTEPS,sizeof(NSTEPS));
			LineFile.write((char*)&NPML,sizeof(NPML));

			//Initialize the whole line array
			LineArray.resize(Range(1,TotalLineLength));	//whole line
			LineArray = 0;

			//construct the 1-D array of partial line buffer arrays
			PartialLineBuf_recv.resize(line_nodes);
			for (n1=0; n1<line_nodes; n1++)
			{
				PartialLineBuf_recv(n1).resize(Range(LineMinArray(n1),LineMaxArray(n1)));
			}
		}
		else
		{
			//Initialize the partial line array used for sending data to node 0
			PartialLineBuf_send.resize(Range(LineMin,LineMax));	//send the partial line that belongs to the current node
		}
	}
}

void Crecorder_line::RecordLine()
{//records a single line
	if (PassesThroughNode)
	{
		if (line_rank==0)
		{
			//do NOT try to receive data from node 0 (yourself), place available data
			PlacePartialLine(LineArray);
			for (n1=1; n1<line_nodes; n1++)
			{
#ifndef MPI_DISABLE
				//receive data from other nodes
				MPI_Recv(PartialLineBuf_recv(n1).data(),PartialLineBuf_recv(n1).size(),MPI_DOUBLE,
					n1,0,MPI_LineComm,&Status);
				// place buffer into place
				LineArray(Range(LineMinArray(n1),LineMaxArray(n1))) = PartialLineBuf_recv(n1);
#else
				// no MPI, just place field values into array
				PlacePartialLine(PartialLineBuf_send);
				LineArray(Range(LineMinArray(n1),LineMaxArray(n1))) = PartialLineBuf_send;
#endif
			}
			//write out the data
			LineFile.write((char*)LineArray.data(),LineArray.size()*sizeof(LineArray(0)));
		}
		else
		{
#ifndef MPI_DISABLE
			//send data to node 0
			PlacePartialLine(PartialLineBuf_send);
			MPI_Send(PartialLineBuf_send.data(),PartialLineBuf_send.size(),MPI_DOUBLE,0,0,MPI_LineComm);
#endif
		}
	}
}

void Crecorder_line::PlacePartialLine(Array<double,1> DestinationArray)
{//places the partial line in DestinationArray, which may be the same size as the partial line (PartialLineBuf_send for rank!=0)
// or the same size as the whole line (LineArray for rank=0).
	for (s=LineMin; s<=LineMax; s++)
	{
		//find the position on the main grid, depending on the orientation and position
		if (orientation=="x_directed")
		{
			xPos = s;
			yPos = LinePos_x1;
			zPos = LinePos_x2;
		}
		else if (orientation=="y_directed")
		{
			xPos = LinePos_x1;
			yPos = s;
			zPos = LinePos_x2;
		}
		else if (orientation=="z_directed")
		{
			xPos = LinePos_x1;
			yPos = LinePos_x2;
			zPos = s;
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
		else if (component=="Ex")
		{
			field_value = Ex(xPos,yPos,zPos);
		}
		else if (component=="Ey")
		{
			field_value = Ey(xPos,yPos,zPos);
		}
		else if (component=="Ez")
		{
			field_value = Ez(xPos,yPos,zPos);
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
		//record the field value (no need for discretization, since only a line is recorded)
		DestinationArray(s) = field_value;
	}
}
