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

//Definiton of the class "Cwfs" representing a collection of time waveforms
//A good reason to define this as a class is that the automatic deletion of dynamically-allocated "Cwf" objects upon destruction.

#include "headers.h"

#include "Cwfs.h"

#include "Cwf_gauss.h"


bool Cwfs::WaveformTagExists(const string& waveform_tag)
{//returns true if a waveform with the given tag exists.
	return (NamedWaveforms.find(waveform_tag)!=NamedWaveforms.end());
}

const_Cwf_shared_ptr Cwfs::operator[] (const string& waveform_tag) const
{//returns the smart pointer to the waveform corresponding to the waveform tag, throws exception if not found
	map<string,Cwf_shared_ptr>::const_iterator map_it = NamedWaveforms.find(waveform_tag);
	if (map_it==NamedWaveforms.end())
	{//string tag does not correspond to any waveform, throw exception
		throw NamedWaveformNotFoundException(waveform_tag);
	}
	else
	{
		return map_it->second;
	}
}

void Cwfs::AddGaussianWaveform(const double& myamplitude, const double& mytau, const double& mydelay, const string& waveform_tag)
{//adds a Gaussian waveform
	if (!WaveformTagExists(waveform_tag))
	{
		//create waveform
		Cwf_shared_ptr new_wf_ptr(new Cwf_gaussian(myamplitude,mytau,mydelay));
		//add the shape with the given tag
		NamedWaveforms.insert(pair<string,Cwf_shared_ptr>(waveform_tag,new_wf_ptr));
	}
	else
	{//string tag already exists, throw exception
		throw NamedWaveformExistsException(waveform_tag);
	}
}

void Cwfs::AddDiffGaussianWaveform(const double& myamplitude, const double& mytau, const double& mydelay, const int& ndiff, const string& waveform_tag)
{//adds a n'th order-differentiated Gaussian waveform
	if (!WaveformTagExists(waveform_tag))
	{
		//create waveform
		Cwf_shared_ptr new_wf_ptr(new Cwf_diffgaussian(myamplitude,mytau,mydelay,ndiff));
		//add the shape with the given tag
		NamedWaveforms.insert(pair<string,Cwf_shared_ptr>(waveform_tag,new_wf_ptr));
	}
	else
	{//string tag already exists, throw exception
		throw NamedWaveformExistsException(waveform_tag);
	}
}

void Cwfs::AddSineModulatedGaussianWaveform(const double& myamplitude, const double& mytau, const double& myf_0, const double& myphase, const double& mydelay, const string& waveform_tag)
{//adds a sine-modulated Gaussian waveform
	if (!WaveformTagExists(waveform_tag))
	{
		//create waveform
		Cwf_shared_ptr new_wf_ptr(new Cwf_sinemodulatedgaussian(myamplitude,mytau,myf_0,myphase,mydelay));
		//add the shape with the given tag
		NamedWaveforms.insert(pair<string,Cwf_shared_ptr>(waveform_tag,new_wf_ptr));
	}
	else
	{//string tag already exists, throw exception
		throw NamedWaveformExistsException(waveform_tag);
	}
}

void Cwfs::AddCosineModulatedGaussianWaveform(const double& myamplitude, const double& mytau, const double& myf_0, const double& myphase, const double& mydelay, const string& waveform_tag)
{//adds a cosine-modulated Gaussian waveform
	if (!WaveformTagExists(waveform_tag))
	{
		//create waveform
		Cwf_shared_ptr new_wf_ptr(new Cwf_cosinemodulatedgaussian(myamplitude,mytau,myf_0,myphase,mydelay));
		//add the shape with the given tag
		NamedWaveforms.insert(pair<string,Cwf_shared_ptr>(waveform_tag,new_wf_ptr));
	}
	else
	{//string tag already exists, throw exception
		throw NamedWaveformExistsException(waveform_tag);
	}
}

void Cwfs::AddDiffSineModulatedGaussianWaveform(const double& myamplitude, const double& mytau, const double& myf_0, const double& myphase, const double& mydelay, const string& waveform_tag)
{//adds a sine-modulated, then differentiated Gaussian waveform
	if (!WaveformTagExists(waveform_tag))
	{
		//create waveform
		Cwf_shared_ptr new_wf_ptr(new Cwf_diffsinemodulatedgaussian(myamplitude,mytau,myf_0,myphase,mydelay));
		//add the shape with the given tag
		NamedWaveforms.insert(pair<string,Cwf_shared_ptr>(waveform_tag,new_wf_ptr));
	}
	else
	{//string tag already exists, throw exception
		throw NamedWaveformExistsException(waveform_tag);
	}
}

void Cwfs::AddDiffCosineModulatedGaussianWaveform(const double& myamplitude, const double& mytau, const double& myf_0, const double& myphase, const double& mydelay, const string& waveform_tag)
{//adds a cosine-modulated, then differentiated Gaussian waveform
	if (!WaveformTagExists(waveform_tag))
	{
		//create waveform
		Cwf_shared_ptr new_wf_ptr(new Cwf_diffcosinemodulatedgaussian(myamplitude,mytau,myf_0,myphase,mydelay));
		//add the shape with the given tag
		NamedWaveforms.insert(pair<string,Cwf_shared_ptr>(waveform_tag,new_wf_ptr));
	}
	else
	{//string tag already exists, throw exception
		throw NamedWaveformExistsException(waveform_tag);
	}
}

//void Cwfs::AddWaveformTag(const string& Tag, const int& Index)
//{//Add the string "Tag" for the waveform with index "Index" into the tag array
//	//increase size by 1
//	WaveformTags.resizeAndPreserve(WaveformTags.size()+1);
//	//insert the string tag
//	WaveformTags(WaveformTags.size()-1).Tag = Tag;
//	//insert the index
//	WaveformTags(WaveformTags.size()-1).Index = Index;
//}

//bool Cwfs::TagExists(const string& tag)
//{//returns true if the tag already exists
//	bool tagfound=false;
//	for (int i=0; (tagfound==false)&&(i<WaveformTags.size()); i++)
//	{
//		if (WaveformTags(i).Tag==tag)
//		{
//			tagfound = true;
//		}
//	}
//	return tagfound;
//}

//bool Cwfs::lookupWaveform(const string& tag, int& waveformindex) const
//{//returns in "waveformindex" the index to the Cwf object with the specified tag
//	bool tagfound=false;
//	for (int i=0; (tagfound==false)&&(i<WaveformTags.size()); i++)
//	{
//		if (WaveformTags(i).Tag==tag)
//		{
//			tagfound = true;
//			waveformindex = i;	//return the index in "waveformindex"
//		}
//	}
//	return tagfound;	//return false if tag is not found
//}
