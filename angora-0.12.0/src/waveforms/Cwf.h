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

#ifndef CWF_H
#define CWF_H

//Declaration of the abstract base class "Cwf" representing a time waveform

//for the vector STL class
#include <vector>
//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>

//declaration of shared-pointers to Cwf objects
#include "Cwf_shared_ptr.h"


class Cwf
{
 public:
//	Cwf(const double& amp_factor = 1, const double& extra_delay = 0 ): _amp_factor(amp_factor), _extra_delay(extra_delay) {};
 	virtual ~Cwf() {};		//virtual destructor is needed when deleting derived objects using a base class pointer

	virtual Cwf_shared_ptr clone() const=0;  //"virtual" copy constructor (http://www.parashift.com/c++-faq-lite/virtual-functions.html#faq-20.8)
								   //also called the "Virtual Constructor Idiom" (http://www.parashift.com/c++-faq-lite/abcs.html#faq-22.5)

	virtual double Value(const double& t) const=0;	//value of the waveform at time t
	virtual double Derivative(const double& t) const=0;		//time derivative of the waveform at time t
	virtual double Integral(const double& t) const=0;		//time integral of the waveform at time t (should ->0 as t->-inf)
	complex<double> FourierComponent(const double& w) const
	{//non-virtual member that returns 1/(2*pi) times the temporal Fourier transform at angular frequency w
	 //The definition is similar to the inverse Fourier transform, except the 1/(2pi) factor is included in the result.
	 //The original waveform can be written as f(t) = Int{-inf}->{+inf}{FC(w)exp(iwt)dw}, which is a sum of monochromatic waveforms with weights FC(w). This member can be used to find the "weights" of the frequency components in a time waveform.
		return 1.0/(2*M_PI)*FourierTransform(w);
	};
	virtual complex<double> FourierTransform(const double& w) const=0;	//temporal Fourier transform of the waveform at angular frequency w
																		//defined as F[f(t)](w) = Int{-inf}->{+inf}{f(t)exp(-iwt)dt}
	virtual double HilbertTransform(const double& t) const=0;	//Hilbert transform of the waveform (applies a -90deg phase shift to all the frequencies in the waveform: converts sines into (-cosine)s, cosines into sines -- equivalent to multiplication by (-i)sgn(w) in the freq. domain
	virtual double PhaseShift(const double& extra_phase,const double& t) const=0;	//applies an arbitrary phase shift "extra_phase" (in rad) to all the frequencies in the waveform: converts sin(wt) into sin(wt+extra_phase), cos(wt) into cos(wt+extra_phase)
	//this operation reduces to the Hilbert transform for extra_phase=pi/2 (90deg)

	virtual double A_max() const=0;		//maximum amplitude in the waveform
	virtual double w_min_20() const=0;		//minimum angular frequency beyond which the spectral amplitude is above -20dB of the peak
	virtual double w_max_20() const=0;		//maximum angular frequency beyond which the spectral amplitude is above -20dB of the peak
	virtual double w_min_40() const=0;		//minimum angular frequency beyond which the spectral amplitude is above -40dB of the peak
	virtual double w_max_40() const=0;		//maximum angular frequency beyond which the spectral amplitude is above -40dB of the peak

	virtual double starting_time() const=0;

	///These should be made more clear
	virtual double F_0() const=0;
	virtual double W_0() const=0;
	///These should be made more clear

	// Virtual methods that generate smart pointers to new Cwf (or child) objects
	virtual Cwf_shared_ptr Derivative() const=0;
	virtual Cwf_shared_ptr Integral() const=0;
	virtual Cwf_shared_ptr HilbertTransform() const=0;
	virtual Cwf_shared_ptr PhaseShift(const double& extra_phase) const=0;

	//multiply the amplitude factor by the given value
	virtual void multiply_amplitude_by(const double& factor)=0;

	//add given value (in seconds) to the extra delay
	virtual void add_extra_delay(const double& additional_delay)=0;

// private:
// protected:
};

#endif
