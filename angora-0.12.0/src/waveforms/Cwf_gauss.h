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

#ifndef CWF_GAUSS_H
#define CWF_GAUSS_H

//Declaration of classes that represent certain time waveforms with a Gaussian envelope, and their time derivatives

#include "Cwf.h"

//base Angora exception class
#include "angora_excp.h"

const double n_tau = 5;


class Cwf_gaussian : public Cwf
{//base class for non-modulated Gaussian waveforms
 public:
	Cwf_gaussian(const double& myamplitude, const double& mytau, const double& mydelay)
		: amplitude(myamplitude), tau(mytau), delay(mydelay)
	{};
	double starting_time() const{return (-n_tau+delay)*tau;};
	///These should be made more clear
	double F_0() const {return 0;};
	double W_0() const {return 0;};
	///These should be made more clear
	//these are redefined in other types of Gaussian waveforms
	double Value(const double& t) const;	//value of the waveform at time t
	double Derivative(const double& t) const;		//time derivative of the waveform at time t
	double Integral(const double& t) const;		//time integral of the waveform at time t (should ->0 as t->-inf)
	complex<double> FourierTransform(const double& w) const;	//temporal Fourier transform of the waveform at angular frequency w
	double HilbertTransform(const double& t) const; //Hilbert transform of the waveform (applies a -90deg phase shift to all the frequencies in the waveform: converts sines into (-cosine)s, cosines into sines -- equivalent to multiplication by (-i)sgn(w) in the freq. domain
	double PhaseShift(const double& extra_phase,const double& t) const;	//applies an arbitrary phase shift "extra_phase" (in rad) to all the frequencies in the waveform: converts sin(wt) into sin(wt+extra_phase), cos(wt) into cos(wt+extra_phase)
	//this operation reduces to the Hilbert transform for extra_phase=pi/2 (90deg)

	double w_min_20() const {return 0;};	//-20 dB frequency in the waveform (minimum) (NOT DEFINED, left as 0)
	double w_max_20() const {return 2.146/tau;};	//-20 dB frequency in the waveform  (maximum)
	double w_min_40() const {return 0;};	//-40 dB frequency in the waveform (minimum) (NOT DEFINED, left as 0)
	double w_max_40() const {return 3.035/tau;};	//-40 dB frequency in the waveform  (maximum)
	double A_max() const {return amplitude;	/*this is exact.*/};	//maximum amplitude in the waveform

	virtual void multiply_amplitude_by(const double& factor)
	{//multiply the amplitude factor by the given value
		amplitude *= factor;
	}

	virtual void add_extra_delay(const double& additional_delay)
	{//add given value (in seconds) to the extra delay
		delay += additional_delay/tau;  //delay is in terms of tau
	}

	virtual Cwf_shared_ptr clone() const;
	virtual Cwf_shared_ptr Derivative() const;
	virtual Cwf_shared_ptr Integral() const;
	virtual Cwf_shared_ptr HilbertTransform() const;
	virtual Cwf_shared_ptr PhaseShift(const double& extra_phase) const;

 protected:
	double amplitude;	//amplitude of the waveform
	double tau;	//time constant for the waveform
	double delay;	//delay of the waveform (in tau)
};

class Cwf_diffgaussian : public Cwf_gaussian
{
 public:
	Cwf_diffgaussian(const double& myamplitude, const double& mytau, const double& mydelay, const int& mydifforder)
		: Cwf_gaussian(myamplitude,mytau,mydelay), ndiff(mydifforder)
	{
		if (ndiff<0)
		{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//			InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
			string func_name = __FUNCTION__;
#else
			string func_name = "";
#endif
			throw AngoraInvalidArgumentExceptionWithType<int>(func_name,ndiff,"(order of differentiation should be nonnegative)");
		}
	};
	double Value(const double& t) const;	//value of the waveform at time t
	double Derivative(const double& t) const;		//time derivative of the waveform at time t
	double Integral(const double& t) const;		//time integral of the waveform at time t (should ->0 as t->-inf)
	complex<double> FourierTransform(const double& w) const;	//temporal Fourier transform of the waveform at angular frequency w
	double HilbertTransform(const double& t) const; //Hilbert transform of the waveform (applies a -90deg phase shift to all the frequencies in the waveform: converts sines into (-cosine)s, cosines into sines -- equivalent to multiplication by (-i)sgn(w) in the freq. domain
	double PhaseShift(const double& extra_phase,const double& t) const;	//applies an arbitrary phase shift "extra_phase" (in rad) to all the frequencies in the waveform: converts sin(wt) into sin(wt+extra_phase), cos(wt) into cos(wt+extra_phase)
	//this operation reduces to the Hilbert transform for extra_phase=pi/2 (90deg)

	double w_min_20() const {return 0.0608/tau;};	//-20 dB frequency in the waveform (minimum) (defined, but extremely small)
	double w_max_20() const {return 2.764/tau;};	//-20 dB frequency in the waveform  (maximum)
	double w_min_40() const {return 0.0061/tau;};	//-40 dB frequency in the waveform (minimum) (defined, but extremely small)
	double w_max_40() const {return 3.57/tau;	/*this is only exact for ndiff=1.*/};	//-40 dB frequency in the waveform (maximum)
	double A_max() const {return (1/pow(tau,ndiff))*amplitude;	/*this is only correct up to an order of magnitude.*/};	//maximum amplitude in the waveform

	virtual Cwf_shared_ptr clone() const;
	virtual Cwf_shared_ptr Derivative() const;
	virtual Cwf_shared_ptr Integral() const;
	virtual Cwf_shared_ptr HilbertTransform() const;
	virtual Cwf_shared_ptr PhaseShift(const double& extra_phase) const;

 private:
	int ndiff;		//order of differentiation (ndiff>=0)
};

class Cwf_modulatedgaussian : public Cwf_gaussian
{//base class for modulated Gaussian waveforms
 public:
	Cwf_modulatedgaussian(const double& myamplitude, const double& mytau, const double& myf_0, const double& myphase, const double& mydelay)
		: Cwf_gaussian(myamplitude,mytau,mydelay), phase(myphase), f_0(myf_0), w_0(2*M_PI*f_0)
	{};
	///These should be made more clear
	double F_0() const {return f_0;};
	double W_0() const {return w_0;};
	double Phase() const {return phase;};
	///These should be made more clear

	double w_min_20() const {return max(0.0, 2*M_PI*f_0 - 2.146/tau);};	//-20 dB frequency in the waveform (minimum)
	double w_max_20() const {return 2*M_PI*f_0 + 2.146/tau;};	//-20 dB frequency in the waveform  (maximum)
	double w_min_40() const {return max(0.0, 2*M_PI*f_0 - 3.035/tau);};	//-40 dB frequency in the waveform (minimum)
	double w_max_40() const {return 2*M_PI*f_0 + 3.035/tau;};	//-40 dB frequency in the waveform (minimum)

 protected:
	double f_0;	//modulation frequency (Hz)
	double w_0;	//angular modulation frequency = 2*pi*f_0
	double phase;	//phase of the modulation sinusoid (0<phase<2pi)
};

class Cwf_sinemodulatedgaussian : virtual public Cwf_modulatedgaussian
{
 public:
	Cwf_sinemodulatedgaussian(const double& myamplitude, const double& mytau, const double& myf_0, const double& myphase, const double& mydelay)
		: Cwf_modulatedgaussian(myamplitude,mytau,myf_0,myphase,mydelay)
	{};
	double Value(const double& t) const;	//value of the waveform at time t
	double Derivative(const double& t) const;		//time derivative of the waveform at time t
	double Integral(const double& t) const;		//time integral of the waveform at time t (should ->0 as t->-inf)
	complex<double> FourierTransform(const double& w) const;	//temporal Fourier transform of the waveform at angular frequency w
	double HilbertTransform(const double& t) const; //Hilbert transform of the waveform (applies a -90deg phase shift to all the frequencies in the waveform: converts sines into (-cosine)s, cosines into sines -- equivalent to multiplication by (-i)sgn(w) in the freq. domain
	double PhaseShift(const double& extra_phase,const double& t) const;	//applies an arbitrary phase shift "extra_phase" (in rad) to all the frequencies in the waveform: converts sin(wt) into sin(wt+extra_phase), cos(wt) into cos(wt+extra_phase)
	//this operation reduces to the Hilbert transform for extra_phase=pi/2 (90deg)

	double A_max() const {return amplitude;	/*this is approximate.*/};	//maximum amplitude in the waveform

	virtual Cwf_shared_ptr clone() const;
	virtual Cwf_shared_ptr Derivative() const;
	virtual Cwf_shared_ptr Integral() const;
	virtual Cwf_shared_ptr HilbertTransform() const;
	virtual Cwf_shared_ptr PhaseShift(const double& extra_phase) const;

};

class Cwf_cosinemodulatedgaussian : virtual public Cwf_modulatedgaussian
{
 public:
	Cwf_cosinemodulatedgaussian(const double& myamplitude, const double& mytau, const double& myf_0, const double& myphase, const double& mydelay)
		: Cwf_modulatedgaussian(myamplitude,mytau,myf_0,myphase,mydelay)
	{};
	double Value(const double& t) const;	//value of the waveform at time t
	double Derivative(const double& t) const;		//time derivative of the waveform at time t
	double Integral(const double& t) const;		//time integral of the waveform at time t (should ->0 as t->-inf)
	complex<double> FourierTransform(const double& w) const;	//temporal Fourier transform of the waveform at angular frequency w
	double HilbertTransform(const double& t) const; //Hilbert transform of the waveform (applies a -90deg phase shift to all the frequencies in the waveform: converts sines into (-cosine)s, cosines into sines -- equivalent to multiplication by (-i)sgn(w) in the freq. domain
	double PhaseShift(const double& extra_phase,const double& t) const;	//applies an arbitrary phase shift "extra_phase" (in rad) to all the frequencies in the waveform: converts sin(wt) into sin(wt+extra_phase), cos(wt) into cos(wt+extra_phase)
	//this operation reduces to the Hilbert transform for extra_phase=pi/2 (90deg)

	double A_max() const {return amplitude;	/*this is exact.*/};	//maximum amplitude in the waveform

	virtual Cwf_shared_ptr clone() const;
	virtual Cwf_shared_ptr Derivative() const;
	virtual Cwf_shared_ptr Integral() const;
	virtual Cwf_shared_ptr HilbertTransform() const;
	virtual Cwf_shared_ptr PhaseShift(const double& extra_phase) const;

};

class Cwf_diffsinemodulatedgaussian : public Cwf_sinemodulatedgaussian, public Cwf_cosinemodulatedgaussian //also inherit from cos. modulated, for Hilbert-transform purposes
{
 public:
	Cwf_diffsinemodulatedgaussian(const double& myamplitude, const double& mytau, const double& myf_0, const double& myphase, const double& mydelay)
		: Cwf_sinemodulatedgaussian(myamplitude,mytau,myf_0,myphase,mydelay),
		Cwf_cosinemodulatedgaussian(myamplitude,mytau,myf_0,myphase,mydelay),
		Cwf_modulatedgaussian(myamplitude,mytau,myf_0,myphase,mydelay)	//because Cwf_modulatedgaussian is virtually inherited by Cwf_sinemodulatedgaussian and Cwf_cosinemodulatedgaussian, this is where it has to be initialized. Virtual inheritance means that objects of the intermediate classes Cwf_sinemodulatedgaussian and Cwf_cosinemodulatedgaussian contain only a vtable rather than the actual base-class object Cwf_modulatedgaussian, leaving the actual instantiation to the object of the most-derived class (Cwf_diffsinemodulatedgaussian in this case).
	{};

	double Value(const double& t) const;	//value of the waveform at time t
	double Derivative(const double& t) const;		//time derivative of the waveform at time t
	double Integral(const double& t) const;		//time integral of the waveform at time t (should ->0 as t->-inf)
	complex<double> FourierTransform(const double& w) const;	//temporal Fourier transform of the waveform at angular frequency w
	double HilbertTransform(const double& t) const; //Hilbert transform of the waveform (applies a -90deg phase shift to all the frequencies in the waveform: converts sines into (-cosine)s, cosines into sines -- equivalent to multiplication by (-i)sgn(w) in the freq. domain
	double PhaseShift(const double& extra_phase,const double& t) const;	//applies an arbitrary phase shift "extra_phase" (in rad) to all the frequencies in the waveform: converts sin(wt) into sin(wt+extra_phase), cos(wt) into cos(wt+extra_phase)
	//this operation reduces to the Hilbert transform for extra_phase=pi/2 (90deg)

	// w_max_40,w_min_40,w_max_20,w_min_20 is approximately the same as in Cwf_modulatedgaussian, if w does not change appreciably within the Gaussian spectral envelope.
	double A_max() const {return 2*M_PI*f_0*amplitude;	/*this is approximate.*/};	//maximum amplitude in the waveform

	virtual Cwf_shared_ptr clone() const;
	virtual Cwf_shared_ptr Derivative() const;
	virtual Cwf_shared_ptr Integral() const;
	virtual Cwf_shared_ptr HilbertTransform() const;
	virtual Cwf_shared_ptr PhaseShift(const double& extra_phase) const;

};

class Cwf_diffcosinemodulatedgaussian : public Cwf_cosinemodulatedgaussian, public Cwf_sinemodulatedgaussian //also inherit from sin. modulated, for Hilbert-transform purposes
{
 public:
	Cwf_diffcosinemodulatedgaussian(const double& myamplitude, const double& mytau, const double& myf_0, const double& myphase, const double& mydelay)
		: Cwf_cosinemodulatedgaussian(myamplitude,mytau,myf_0,myphase,mydelay),
		Cwf_sinemodulatedgaussian(myamplitude,mytau,myf_0,myphase,mydelay),
		Cwf_modulatedgaussian(myamplitude,mytau,myf_0,myphase,mydelay)	//because Cwf_modulatedgaussian is virtually inherited by Cwf_sinemodulatedgaussian and Cwf_cosinemodulatedgaussian, this is where it has to be initialized. Virtual inheritance means that objects of the intermediate classes Cwf_sinemodulatedgaussian and Cwf_cosinemodulatedgaussian contain only a vtable rather than the actual base-class object Cwf_modulatedgaussian, leaving the actual instantiation to the object of the most-derived class (Cwf_diffsinemodulatedgaussian in this case).
	{};

	double Value(const double& t) const;	//value of the waveform at time t
	double Derivative(const double& t) const;		//time derivative of the waveform at time t
	double Integral(const double& t) const;		//time integral of the waveform at time t (should ->0 as t->-inf)
	complex<double> FourierTransform(const double& w) const;	//temporal Fourier transform of the waveform at angular frequency w
	double HilbertTransform(const double& t) const; //Hilbert transform of the waveform (applies a -90deg phase shift to all the frequencies in the waveform: converts sines into (-cosine)s, cosines into sines -- equivalent to multiplication by (-i)sgn(w) in the freq. domain
	double PhaseShift(const double& extra_phase,const double& t) const;	//applies an arbitrary phase shift "extra_phase" (in rad) to all the frequencies in the waveform: converts sin(wt) into sin(wt+extra_phase), cos(wt) into cos(wt+extra_phase)
	//this operation reduces to the Hilbert transform for extra_phase=pi/2 (90deg)

	// w_max_40,w_min_40,w_max_20,w_min_20 is approximately the same as in Cwf_modulatedgaussian, if w does not change appreciably within the Gaussian spectral envelope.
	double A_max() const {return 2*M_PI*f_0*amplitude;	/*this is approximate.*/};	//maximum amplitude in the waveform

	virtual Cwf_shared_ptr clone() const;
	virtual Cwf_shared_ptr Derivative() const;
	virtual Cwf_shared_ptr Integral() const;
	virtual Cwf_shared_ptr HilbertTransform() const;
	virtual Cwf_shared_ptr PhaseShift(const double& extra_phase) const;

};

#endif
