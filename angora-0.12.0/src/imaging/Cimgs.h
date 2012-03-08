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

#ifndef CIMGS_H
#define CIMGS_H

//Declaration of the collector class "Cimgs" for optical image objects Cimg

//only the declaration of Cimg needed: use forward declaration
class Cimg;
//only the declaration of ImgDataType needed: use forward declaration
class ImgDataType;

//for the vector STL class
#include <vector>
//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>

class Cimgs
{
 public:
	 int AddOpticalImage(const ImgDataType& MyData,
		 const string& InputFarFieldFileName = "");	//adds an optical image

	 void UpdateFarField(const int& n);		//Update far-field arrays

	 void ConstructImages();	//Construct optical image from far-field arrays

	 int NumberOfOpticalImages()
	 {//returns the number of images in the collector object
		 return Images.size();
	 }

 private:
	 vector<boost::shared_ptr<Cimg> > Images;		//array of smart pointers to image objects
};

#endif
