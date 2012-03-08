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

#ifndef CWFS_H
#define CWFS_H

//Declaration of the class "Cwfs" representing a collection of time waveforms
//A good reason to define this as a class is the automatic deletion of dynamically-allocated "Cwf" objects upon destruction.

#include "Cwfs_excp.h"

//declaration of shared-pointers to Cwf objects
#include "waveforms/Cwf_shared_ptr.h"

//for the C++ STL map class
#include <map>
//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>


class Cwfs
{
 public:
	 void AddGaussianWaveform(const double& myamplitude, const double& mytau, const double& mydelay, const string& waveform_tag);	//adds a Gaussian waveform
	 void AddDiffGaussianWaveform(const double& myamplitude, const double& mytau, const double& mydelay, const int& ndiff, const string& waveform_tag);	//adds a n'th order-differentiated Gaussian waveform
	 void AddSineModulatedGaussianWaveform(const double& myamplitude, const double& mytau, const double& myf_0, const double& myphase, const double& mydelay, const string& waveform_tag);	//adds a sine-modulated Gaussian waveform
	 void AddCosineModulatedGaussianWaveform(const double& myamplitude, const double& mytau, const double& myf_0, const double& myphase, const double& mydelay, const string& waveform_tag);	//adds a cosine-modulated Gaussian waveform
	 void AddDiffSineModulatedGaussianWaveform(const double& myamplitude, const double& mytau, const double& myf_0, const double& myphase, const double& mydelay, const string& waveform_tag);	//adds a sine-modulated, then differentiated Gaussian waveform
	 void AddDiffCosineModulatedGaussianWaveform(const double& myamplitude, const double& mytau, const double& myf_0, const double& myphase, const double& mydelay, const string& waveform_tag);	//adds a cosine-modulated, then differentiated Gaussian waveform

//	 void AddWaveformTag(const string& Tag, const int& Index);	//Add the string "Tag" for the waveform with index "Index" into the tag array

//	 bool TagExists(const string& tag);	//returns true if the tag already exists

//	 bool lookupWaveform(const string& tag, int& waveformindex) const;	//returns in "waveformindex" the index to the Cwf object with the specified tag

	 const_Cwf_shared_ptr operator[] (const string& waveform_tag) const; //returns the smart pointer to the waveform corresponding to the waveform tag, throws exception if not found

//	 const_Cwf_shared_ptr Waveform(const int& index) const
//	 {//the pointer to the i'th waveform object
//	 	if ((index<0)||(index>=NumberOfWaveforms()))
//	 	{
//	 		return const_Cwf_shared_ptr();
//	 	}
//	 	else
//		{
////			return WaveformPtrs[index].get();
//			return WaveformPtrs[index];
//		}
//	 }

//	 int NumberOfWaveforms() const
//	 {
//		 return WaveformPtrs.size();
//	 }

private:
	 //C++ STL map object that holds the named-waveform information
	 map<string,Cwf_shared_ptr> NamedWaveforms;

	 bool WaveformTagExists(const string& waveform_tag);

//	 vector<const_Cwf_shared_ptr> WaveformPtrs;	//array of pointers to waveform objects
};

#endif
