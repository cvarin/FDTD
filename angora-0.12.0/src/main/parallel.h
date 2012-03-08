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

#ifndef PARALLEL_H
#define PARALLEL_H

//Use MPI if not disabled
#ifndef MPI_DISABLE
#include <mpi.h>
#endif

//Parallel information exchange routine
void exchangeH();
void init_parallel();
void MPI_exit(const int& exitcode);


//class CMPI
//{
// public:
//	 CMPI(int argc, char *argv[])
//	 {
//#ifndef MPI_DISABLE
//		 MPI_Init(&argc,&argv);
//#endif
//	 }
//	 ~CMPI()
//	 {
//#ifndef MPI_DISABLE
//		 MPI_Barrier(MPI_COMM_WORLD);	//wait until all groups are finished
//		 MPI_Finalize();
//#endif
//	 }
//};

#endif
