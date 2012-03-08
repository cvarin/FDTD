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

//Definiton of classes that represent certain time waveforms with a Gaussian envelope, and their time derivatives

#include "headers.h"

#include "Cwf_gauss.h"

extern double Hermite(const int& N, const double& x);


//********************************************************************

double Cwf_gaussian::Value(const double& t) const
{
	return amplitude*exp(-pow2((t-delay*tau)/(sqrt(2)*tau)));
};

double Cwf_gaussian::Derivative(const double& t) const
{
	return amplitude*(-(t-delay*tau)/pow2(tau))*exp(-pow2((t-delay*tau)/(sqrt(2)*tau)));
};

double Cwf_gaussian::Integral(const double& t) const
{
	return 0;		///may change later to erf(t)
};

complex<double> Cwf_gaussian::FourierTransform(const double& w) const
{
	return amplitude*sqrt(2*M_PI)*tau*exp(-pow2(tau)*pow2(w)/2)*exp(-ii*(delay*tau*w));
};

double Cwf_gaussian::HilbertTransform(const double& t) const
{
	return 0;		///doubt that this will ever be used
}

double Cwf_gaussian::PhaseShift(const double& extra_phase,const double& t) const
{
	return 0;		///doubt that this will ever be used
}

Cwf_shared_ptr Cwf_gaussian::clone() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_gaussian(*this));
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_gaussian::Derivative() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_diffgaussian(amplitude,tau,delay,1)); //order of differentiation = 1
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_gaussian::Integral() const
{
	return Cwf_shared_ptr(); ///may be implemented later with the error function (erf)
}

Cwf_shared_ptr Cwf_gaussian::HilbertTransform() const
{
	return Cwf_shared_ptr();		///doubt that this will ever be used
}

Cwf_shared_ptr Cwf_gaussian::PhaseShift(const double& extra_phase) const
{
	return Cwf_shared_ptr();		///doubt that this will ever be used
}


//********************************************************************

double Cwf_diffgaussian::Value(const double& t) const
{
	return amplitude*pow(-1.0,ndiff)*1/pow(tau*sqrt(2),ndiff)*Hermite(ndiff,(t-delay*tau)/(sqrt(2)*tau))*exp(-pow2((t-delay*tau)/(sqrt(2)*tau)));
};

double Cwf_diffgaussian::Derivative(const double& t) const
{
	return amplitude*pow(-1.0,ndiff+1)*1/pow(tau*sqrt(2),ndiff+1)*Hermite(ndiff+1,(t-delay*tau)/(sqrt(2)*tau))*exp(-pow2((t-delay*tau)/(sqrt(2)*tau)));
};

double Cwf_diffgaussian::Integral(const double& t) const
{
	if (ndiff>0)
	{
		return amplitude*pow(-1.0,ndiff-1)*1/pow(tau*sqrt(2),ndiff-1)*Hermite(ndiff-1,(t-delay*tau)/(sqrt(2)*tau))*exp(-pow2((t-delay*tau)/(sqrt(2)*tau)));
	}
	else
	{
		return 0;		///may change later to erf(t)
	}
};

complex<double> Cwf_diffgaussian::FourierTransform(const double& w) const
{
	return pow(ii*w,ndiff)*  //n'th derivative -> factor of (iw)^n in Fourier space [remember that the kernel in the inverse Fourier integral is exp(iwt)]
		Cwf_gaussian::FourierTransform(w);
};

double Cwf_diffgaussian::HilbertTransform(const double& t) const
{
	return 0;		///doubt that this will ever be used
}

double Cwf_diffgaussian::PhaseShift(const double& extra_phase,const double& t) const
{
	return 0;		///doubt that this will ever be used
}

Cwf_shared_ptr Cwf_diffgaussian::clone() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_diffgaussian(*this));
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_diffgaussian::Derivative() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_diffgaussian(amplitude,tau,delay,ndiff+1)); //increment the order of differentiation by 1
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_diffgaussian::Integral() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_diffgaussian(amplitude,tau,delay,ndiff-1)); //decrement the order of differentiation by 1
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_diffgaussian::HilbertTransform() const
{
	return Cwf_shared_ptr();		///doubt that this will ever be used
}

Cwf_shared_ptr Cwf_diffgaussian::PhaseShift(const double& extra_phase) const
{
	return Cwf_shared_ptr();		///doubt that this will ever be used
}

//********************************************************************

double Cwf_sinemodulatedgaussian::Value(const double& t) const
{
	return amplitude*sin(2*M_PI*f_0*(t-delay*tau)+phase)*exp(-pow2((t-delay*tau)/tau)/2.0);
};

double Cwf_sinemodulatedgaussian::Derivative(const double& t) const
{
	return amplitude*(-(t-delay*tau)/pow2(tau)*sin(2*M_PI*f_0*(t-delay*tau)+phase)+2*M_PI*f_0*cos(2*M_PI*f_0*(t-delay*tau)+phase))*exp(-pow2((t-delay*tau)/tau)/2.0);
};

double Cwf_sinemodulatedgaussian::Integral(const double& t) const
{
	return 0;		///doubt that this will ever be used
};

complex<double> Cwf_sinemodulatedgaussian::FourierTransform(const double& w) const
{
	return amplitude*( (-0.5)*ii*sqrt(2*M_PI)*tau*(exp(-pow2(tau)*pow2(w-w_0)/2)-exp(-pow2(tau)*pow2(w+w_0)/2))*exp(-ii*(delay*tau*w))*cos(phase)
					+ 0.5*sqrt(2*M_PI)*tau*(exp(-pow2(tau)*pow2(w-w_0)/2)+exp(-pow2(tau)*pow2(w+w_0)/2))*exp(-ii*(delay*tau*w))*sin(phase));
};

double Cwf_sinemodulatedgaussian::HilbertTransform(const double& t) const
{
	return -amplitude*cos(2*M_PI*f_0*(t-delay*tau)+phase)*exp(-pow2((t-delay*tau)/tau)/2.0);  //Note the inverted amplitude, since Hilb. transform converts sine to (-cosine)
}

double Cwf_sinemodulatedgaussian::PhaseShift(const double& extra_phase,const double& t) const
{
	return amplitude*sin(2*M_PI*f_0*(t-delay*tau)+phase+extra_phase)*exp(-pow2((t-delay*tau)/tau)/2.0);  //converts sin(wt) to sin(wt+extra_phase)
}

Cwf_shared_ptr Cwf_sinemodulatedgaussian::clone() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_sinemodulatedgaussian(*this));
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_sinemodulatedgaussian::Derivative() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_diffsinemodulatedgaussian(amplitude,tau,f_0,phase,delay));
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_sinemodulatedgaussian::Integral() const
{
	return Cwf_shared_ptr();		///doubt that this will ever be used
}

Cwf_shared_ptr Cwf_sinemodulatedgaussian::HilbertTransform() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_cosinemodulatedgaussian(-amplitude,tau,f_0,phase,delay));//Note the inverted amplitude, since Hilb. transform converts sine to (-cosine)
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_sinemodulatedgaussian::PhaseShift(const double& extra_phase) const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_sinemodulatedgaussian(amplitude,tau,f_0,phase+extra_phase,delay)); //converts sin(wt) to sin(wt+extra_phase)
	return new_wf_ptr;
}

//********************************************************************

double Cwf_cosinemodulatedgaussian::Value(const double& t) const
{
	return amplitude*cos(2*M_PI*f_0*(t-delay*tau)+phase)*exp(-pow2((t-delay*tau)/tau)/2.0);
};

double Cwf_cosinemodulatedgaussian::Derivative(const double& t) const
{
	return amplitude*(-(t-delay*tau)/pow2(tau)*cos(2*M_PI*f_0*(t-delay*tau)+phase)-2*M_PI*f_0*sin(2*M_PI*f_0*(t-delay*tau)+phase))*exp(-pow2((t-delay*tau)/tau)/2.0);
};

double Cwf_cosinemodulatedgaussian::Integral(const double& t) const
{
	return 0;		///doubt that this will ever be used
};

complex<double> Cwf_cosinemodulatedgaussian::FourierTransform(const double& w) const
{
	return amplitude*( 0.5*sqrt(2*M_PI)*tau*(exp(-pow2(tau)*pow2(w-w_0)/2)+exp(-pow2(tau)*pow2(w+w_0)/2))*exp(-ii*(delay*tau*w))*cos(phase)
					- (-0.5)*ii*sqrt(2*M_PI)*tau*(exp(-pow2(tau)*pow2(w-w_0)/2)-exp(-pow2(tau)*pow2(w+w_0)/2))*exp(-ii*(delay*tau*w))*sin(phase));
};

double Cwf_cosinemodulatedgaussian::HilbertTransform(const double& t) const
{
	return amplitude*sin(2*M_PI*f_0*(t-delay*tau)+phase)*exp(-pow2((t-delay*tau)/tau)/2.0);  //converts cosine to sine
}

double Cwf_cosinemodulatedgaussian::PhaseShift(const double& extra_phase,const double& t) const
{
	return amplitude*cos(2*M_PI*f_0*(t-delay*tau)+phase+extra_phase)*exp(-pow2((t-delay*tau)/tau)/2.0);  //converts cos(wt) to cos(wt+extra_phase)
}

Cwf_shared_ptr Cwf_cosinemodulatedgaussian::clone() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_cosinemodulatedgaussian(*this));
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_cosinemodulatedgaussian::Derivative() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_diffcosinemodulatedgaussian(amplitude,tau,f_0,phase,delay));
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_cosinemodulatedgaussian::Integral() const
{
	return Cwf_shared_ptr();		///doubt that this will ever be used
}

Cwf_shared_ptr Cwf_cosinemodulatedgaussian::HilbertTransform() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_sinemodulatedgaussian(amplitude,tau,f_0,phase,delay)); //converts cosine to sine
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_cosinemodulatedgaussian::PhaseShift(const double& extra_phase) const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_cosinemodulatedgaussian(amplitude,tau,f_0,phase+extra_phase,delay));  //converts cos(wt) to cos(wt+extra_phase)
	return new_wf_ptr;
}

//********************************************************************

double Cwf_diffsinemodulatedgaussian::Value(const double& t) const
{
	return Cwf_sinemodulatedgaussian::Derivative(t);
};

double Cwf_diffsinemodulatedgaussian::Derivative(const double& t) const
{
	return 0;		///doubt that this will ever be used
};

double Cwf_diffsinemodulatedgaussian::Integral(const double& t) const
{
	return Cwf_sinemodulatedgaussian::Value(t);
};

complex<double> Cwf_diffsinemodulatedgaussian::FourierTransform(const double& w) const
{
	return (ii*w)*  //1st derivative -> factor of (iw) in Fourier space [remember that the kernel in the inverse Fourier integral is exp(iwt)]
		Cwf_sinemodulatedgaussian::FourierTransform(w);
};

double Cwf_diffsinemodulatedgaussian::HilbertTransform(const double& t) const
{
	return Cwf_sinemodulatedgaussian::Derivative()->HilbertTransform(t);  //takes derivative, applies Hilbert transform
}

double Cwf_diffsinemodulatedgaussian::PhaseShift(const double& extra_phase,const double& t) const
{
	return Cwf_sinemodulatedgaussian::Derivative()->PhaseShift(extra_phase,t);  //takes derivative, applies phase shift
}

Cwf_shared_ptr Cwf_diffsinemodulatedgaussian::clone() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_diffsinemodulatedgaussian(*this));
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_diffsinemodulatedgaussian::Derivative() const
{
	return Cwf_shared_ptr();		///doubt that this will ever be used
}

Cwf_shared_ptr Cwf_diffsinemodulatedgaussian::Integral() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_sinemodulatedgaussian(amplitude,tau,f_0,phase,delay));
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_diffsinemodulatedgaussian::HilbertTransform() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_diffcosinemodulatedgaussian(-amplitude,tau,f_0,phase,delay)); //Note the inverted amplitude, since Hilb. transform converts sine to (-cosine)
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_diffsinemodulatedgaussian::PhaseShift(const double& extra_phase) const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_diffsinemodulatedgaussian(amplitude,tau,f_0,phase+extra_phase,delay));  //converts sin(wt) to sin(wt+extra_phase)
	return new_wf_ptr;
}

//********************************************************************

double Cwf_diffcosinemodulatedgaussian::Value(const double& t) const
{
	return Cwf_cosinemodulatedgaussian::Derivative(t);
};

double Cwf_diffcosinemodulatedgaussian::Derivative(const double& t) const
{
	return 0;		///doubt that this will ever be used
};

double Cwf_diffcosinemodulatedgaussian::Integral(const double& t) const
{
	return Cwf_cosinemodulatedgaussian::Value(t);
};

complex<double> Cwf_diffcosinemodulatedgaussian::FourierTransform(const double& w) const
{
	return (ii*w)*  //1st derivative -> factor of (iw) in Fourier space [remember that the kernel in the inverse Fourier integral is exp(iwt)]
		Cwf_cosinemodulatedgaussian::FourierTransform(w);
};

double Cwf_diffcosinemodulatedgaussian::HilbertTransform(const double& t) const
{
	return Cwf_cosinemodulatedgaussian::Derivative()->HilbertTransform(t);  //takes derivative, applies Hilbert transform
}

double Cwf_diffcosinemodulatedgaussian::PhaseShift(const double& extra_phase,const double& t) const
{
	return Cwf_cosinemodulatedgaussian::Derivative()->PhaseShift(extra_phase,t);  //takes derivative, applies phase shift
}

Cwf_shared_ptr Cwf_diffcosinemodulatedgaussian::clone() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_diffcosinemodulatedgaussian(*this));
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_diffcosinemodulatedgaussian::Derivative() const
{
	return Cwf_shared_ptr();		///doubt that this will ever be used
}

Cwf_shared_ptr Cwf_diffcosinemodulatedgaussian::Integral() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_cosinemodulatedgaussian(amplitude,tau,f_0,phase,delay));
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_diffcosinemodulatedgaussian::HilbertTransform() const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_diffsinemodulatedgaussian(amplitude,tau,f_0,phase,delay)); //converts cosine to sine
	return new_wf_ptr;
}

Cwf_shared_ptr Cwf_diffcosinemodulatedgaussian::PhaseShift(const double& extra_phase) const
{
	Cwf_shared_ptr new_wf_ptr(new Cwf_diffcosinemodulatedgaussian(amplitude,tau,f_0,phase+extra_phase,delay));  //converts cos(wt) to cos(wt+extra_phase)
	return new_wf_ptr;
}
