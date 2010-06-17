/*
  Copyright (C) 2006 - 2008  Jeffrey M. McMahon

  This file is part of JSCIENCE.

  JSCIENCE is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  JSCIENCE is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with JSCIENCE.  If not, see <http://www.gnu.org/licenses/>.
*/


//************************************************************************
//------------------------------------------------------------------------
//
// NAME: JSCIENCE SPECIAL FUNCTION SUBROUTINES
//
// VERSION: 1.00
//		Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jmath_special.cpp
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
//	DESC:	Special subroutines are for calculating Bessel functions
//		of the first and second kinds and Hankel functions.
//
//
//	LIST OF SUBROUTINES:
//
//		BESSEL FUNCTIONS:
//
//			1) double bessel_1(int n, double arg) :: returns the  
//				Bessel function of the first kind of order n
//			2) complex<double> bessel_1_complex(int n, 
//				complex<double> arg) :: returns the Bessel 
//				function of the first kind of order n for 
//				complex arguments
//			3) double bessel_2(int n, double arg) :: returns the  
//				Bessel function of the second kind of order n
//			4) complex<double> bessel_2_complex(int n, 
//				complex<double> arg) :: returns the Bessel 
//				function of the second kind of order n for 
//				complex arguments
//			5) double h_s(int s) :: returns (1/1 + 1/2 + ... + 1/s)
//
//		BESSEL FUNCTION DERIVATIVES:
//
//			6) complex<double> bessel_1_deriv_complex(int n, 
//			 	complex<double> arg) :: returns the derivative 
//				of the first kind of order n for a complex
//				argument
//
//		HANKEL FUNCTIONS:
//
//			7) complex<double> hankel_1(int n, complex<double> arg) 
//				:: returns the Hankel function of the first kind 
//				of order n for complex arguments
//
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************

//************************************************************************
//			INCLUDE FILES
//************************************************************************
#include "jmath_special.h"


//************************************************************************
//			SUBROUTINES
//************************************************************************



//========================================================================
//========================================================================
//
//	NAME:	double bessel_1(int n, double arg)
//
//	DESC:	Calculates the Bessel function of the first kind (Jn) of order n.
//
//	INPUT:
//		int n:: Order of Bessel function
//		double arg: Bessel function argument
//	
//	OUTPUT:
//		Jn:: Bessel function of the first kind of order n
//		
//	NOTES:	1) The Bessel function of the first kind is defined as:
//			Jn(x) = ((x/2.0)^n)*sum{m=0, m=inf}[(((-x^2)/4.0)^m)/(m!*(n+m)!)
//
//		2) J-n = ((-1)^n)*Jn
//
//		3) The infinite series representing the Bessel function
//		converges for all arguments given that enough terms are taken:
//		k >> |arg| 
//
//========================================================================
//========================================================================
double bessel_1(int n, double arg)
{

  double tolerance = 100.0*depsilon();

  double Jn, Jn1;
  Jn = Jn1 = 0.0;

  // make sure we calculate the positive order Bessel function
  int m = abs(n);

  // !!! my concern with this do loop is that if a certain term contributes 0 then the 
  // loop may inappropriately exit, and I don't know if this will be large enough
  int k_max2 = int(2.0*arg + 1.0);

  int k = 0;

  do
  {
    Jn1 = Jn;

    Jn += pow(-arg*arg/4.0, k)/(dfactorial(k)*dfactorial(k+m));

    k++;
  } while ((fabs(Jn - Jn1) > tolerance) || (k < k_max2));

  Jn *= pow((arg/2.0), m);

  // NOW WE WILL MAKE USE OF THE BESSEL FUNCTION RELATION FOR NEGATIVE n
  // J(n-1) = (-1)^n*Jn
  if(n < 0)
  {
    Jn *= pow(-1.0, m);
  }


  return Jn;

}



//========================================================================
//========================================================================
//
//	NAME:	complex<double> bessel_1_complex(int n, complex<double> arg)
//
//	DESC:	Calculates the Bessel function of the first kind (Jn) for
//		complex arguments.
//
//	INPUT:
//		int n:: Order of Bessel function
//		complex<double> arg: Bessel function argument
//	
//	OUTPUT:
//		Jn:: Bessel function of the first kind of order n
//		
//	NOTES:	1) The Bessel function of the first kind is defined as:
//			Jn(x) = ((x/2.0)^n)*sum{m=0, m=inf}[(((-x^2)/4.0)^m)/(m!*(n+m)!)
//
//		2) J-n = ((-1)^n)*Jn
//
//		3) The infinite series representing the Bessel function
//		converges for all arguments given that enough terms are taken:
//		k >> |arg| 
//
//========================================================================
//========================================================================
complex<double> bessel_1_complex(int n, complex<double> arg)
{

  complex<double> Jn(0.0, 0.0);
  complex<double> Jn_series(0.0, 0.0);
  complex<double> Jn_seriesm1(0.0, 0.0);

  double tolerance = 100.0*depsilon();

  // Make sure that we calculate the positive order Bessel function
  int m = abs(n);

  // !!! I don't know if this will be large enough
  int k_max_r, k_max_im; 
  k_max_r = static_cast<int>(real(arg));
  k_max_im = static_cast<int>(imag(arg));

/*
  int k_max = 10*imax(k_max_r, k_max_im);

  for(int k = 0; k <= k_max; k++)
  {
    Jn_series = Jn_series + pow((0.0 - arg*arg/4.0), k)/(dfactorial(k)*dfactorial(k+m));
  }
*/

  // !!! my concern with this do loop is that if a certain term contributes 0 then the 
  // loop may inappropriately exit
  int k_max2 = static_cast<int>(2.0*static_cast<double>(imax(k_max_r, k_max_im))) + 1;

  int k = 0;

  do
  {
    Jn_seriesm1 = Jn_series;

    Jn_series = Jn_series + pow((0.0 - arg*arg/4.0), k)/(dfactorial(k)*dfactorial(k+m));

    k = k + 1;

  } while ((sqrt(real(Jn_series*conj(Jn_series) - Jn_seriesm1*conj(Jn_seriesm1))) > tolerance)   || (k < k_max2));


  Jn = pow((arg/2.0), m)*Jn_series;


  // NOW WE WILL MAKE USE OF THE BESSEL FUNCTION RELATION FOR NEGATIVE n:
  // J(n-1) = (-1)^n*Jn
  if(n < 0)
  {
    Jn = pow((0.0 - 1.0), m)*Jn;
  }


  return Jn;

}



//========================================================================
//========================================================================
//
//	NAME:	double bessel_2(int n, double arg)
//
//	DESC:	Calculates the Bessel function of the second kind (Yn).
//
//	INPUT:
//		int n:: Order of Bessel function
//		double arg: Bessel function argument
//	
//	OUTPUT:
//		Yn:: Bessel function of the second kind of order n
//		
//	NOTES:	1) The Bessel function of the second kind is defined as:
//			Yn(x) = (2.0*Jn(x)/pi)*(ln(x/2) + gamma_e) 
//				+ sum{m=0, m=inf}[(((-1)^(m-1))*(hs(m) + 
//					hs(m+n)) + x^(2*m))/(((2.0^(2*m+n))*m!*(m+n)!)
//				- sum{m=0, m=(n-1)}[((n-m-1)!*x^(2*m))/((2.0^(2*m-n))*m!)
//
//		2) Y-n = ((-1)^n)*Yn
//
//========================================================================
//========================================================================
double bessel_2(int n, double arg)
{

  // get a tolerance for the "infinite" sum
  double tol = 100.0*depsilon();

  int m, mmax;

  double Jn, Yn, sum1, sum1_prev, sum2;
  sum1 = sum1_prev = sum2 = 0.0;

  // Make sure that we calculate the positive order Bessel function
  int k = abs(n);

  // GET THE BESSEL FUNCTION OF THE FIRST KIND
  Jn = bessel_1(k, arg);


  // !!! my concern with this do loop is that if a certain term contributes 0 then the 
  // loop may inappropriately exit
  mmax = static_cast<int>(2.0*arg) + 1;
  m = 0;

  // GET TERM SUM 1
  do
  {
    sum1_prev = sum1;

    sum1 += pow((0.0-1.0), (m-1))*(h_s(m) + h_s(m+k))*pow(arg, (2*m))/(pow(2.0, (2*m+k))*dfactorial(m)*dfactorial(m+k));

    m++;

  } while((fabs(sum1 - sum1_prev) > tol) || (m < mmax));

  sum1 = pow(arg, k)*sum1/(1.0*PI);


  // GET TERM SUM 2
  for(m = 0; m <= (k-1); ++m)
  {
    sum2 = sum2 + dfactorial(k-m-1)*pow(arg, (2*m))/(pow(2.0, (2*m-k))*dfactorial(m));
  }
  
  sum2 = pow(arg, (0-k))*sum2/(1.0*PI);


  // NOW GET Yn
  Yn = (2.0*Jn/(1.0*PI))*(log(arg/2.0) + EULERC) + sum1 - sum2;


  // NOW WE WILL MAKE USE OF THE BESSEL FUNCTION RELATION FOR NEGATIVE n:
  // Y(n-1) = (-1)^n*Yn
  if(n < 0)
  {
    Yn = pow((0.0-1.0), k)*Yn;
  }

  return Yn;
}



//========================================================================
//========================================================================
//
//	NAME:	complex<double> bessel_2_complex(int n, complex<double> arg)
//
//	DESC:	Calculates the Bessel function of the second kind (Yn) for
//		complex arguments.
//
//	INPUT:
//		int n:: Order of Bessel function
//		complex<double> arg: Bessel function argument
//	
//	OUTPUT:
//		complex<double> Yn:: Bessel function of the second kind of
//			order n
//		
//	NOTES:	1) The Bessel function of the second kind is defined as:
//			Yn(x) = (2.0*Jn(x)/pi)*(ln(x/2) + gamma_e) 
//				+ sum{m=0, m=inf}[(((-1)^(m-1))*(hs(m) + 
//					hs(m+n)) + x^(2*m))/(((2.0^(2*m+n))*m!*(m+n)!)
//				- sum{m=0, m=(n-1)}[((n-m-1)!*x^(2*m))/((2.0^(2*m-n))*m!)
//
//		2) Y-n = ((-1)^n)*Yn
//
//========================================================================
//========================================================================
complex<double> bessel_2_complex(int n, complex<double> arg)
{

  complex<double> Yn(0.0, 0.0);

  // get a tolerance for the "infinite" sum
  double tolerance = 100.0*depsilon();

 
  // make sure we calculate the positive order Bessel function
  int k = abs(n);


  // GET BESSEL FUNCTIONS OF THE FIRST KIND
  complex<double> Jn = bessel_1_complex(k, arg);


  // GET TERM SUM 1
  complex<double> sum1(0.0,0.0), sum11(0.0,0.0);


/*
  for(int i = 0; i <= 20; i++)
  {
    sum1 = sum1 + pow((0.0-1.0), (i-1))*(h_s(i) + h_s(i+k))*pow(arg, (2*i))/(pow(2.0, (2*i+k))*dfactorial(i)*dfactorial(i+k));
  }
*/

  // !!! my concern with this do loop is that if a certain term contributes 0 then the 
  // loop may inappropriately exit
  int mm_max = static_cast<int>(2.0*real(arg)) + 1;

  int mm = 0;

  do
  {
    sum11 = sum1;

    sum1 = sum1 + pow((0.0-1.0), (mm-1))*(h_s(mm) + h_s(mm+k))*pow(arg, (2*mm))/(pow(2.0, (2*mm+k))*dfactorial(mm)*dfactorial(mm+k));

    mm = mm + 1;

  } while((fabs(real(sqrt(sum1*conj(sum1) - sum11*conj(sum11)))) > tolerance) || (mm < mm_max));

  sum1 = pow(arg, k)*sum1/(1.0*PI);


  // GET TERM SUM 2
  complex<double> sum2(0.0,0.0);

  for(int m = 0; m <= (k-1); m++)
  {
    sum2 = sum2 + dfactorial((k-m-1))*pow(arg, (2*m))/(pow(2.0, (2*m-k))*dfactorial(m));
  }
  
  sum2 = pow(arg, (0-k))*sum2/(1.0*PI);


  // NOW GET Yn
  Yn = (2.0*Jn/(1.0*PI))*(log(arg/2.0) + EULERC) + sum1 - sum2;


  // NOW WE WILL MAKE USE OF THE BESSEL FUNCTION RELATION FOR NEGATIVE n:
  // Y(n-1) = (-1)^n*Yn
  if(n < 0)
  {
    Yn = pow((0.0-1.0), k)*Yn;
  }


  return Yn;

}


//========================================================================
//========================================================================
//
//	NAME:	double h_s(int s)
//
//	DESC:	Calculates the hs function for use in the calculation
//		of the Bessel functions of the second kind.
//
//	INPUT:
//		int s:: 
//	
//	OUTPUT:
//		double hs:: 
//		
//	NOTES:	1) The hs function is calculated as:
//			hs = 1/1 + 1/2 + 1/3 + ... + 1/s
//		2) h0 = 0.0
//
//========================================================================
//========================================================================
double h_s(int s)
{

  double hs = 0.0;

  if(s == 0)
  {
    hs = 0.0;
  }
  else
  {

    for(int m = 1; m <= s; m++)
    {
      hs = hs + 1.0/static_cast<double>(m);
    }

  }


  return hs;

}




//========================================================================
//========================================================================
//
//	NAME:	complex<double> bessel_1_deriv_complex(int n, 
//			complex<double> arg)
//	DESC:	Calculates the derivative of the Bessel function of the 
//		first kind.
//
//	INPUT:
//		int n:: Function order
//		double arg:: Function argument
//	
//	OUTPUT:
//		complex<double> Jnd:: Derivative of Bessel function of
//			the first kind
//		
//	NOTES:	1) The derivative of the Bessel function is defined as:
//			Jnd(x) = Jn-1(x) - n*Jn(x)/x
//
//========================================================================
//========================================================================
complex<double> bessel_1_deriv_complex(int n, complex<double> arg)
{

  complex<double> Jnd(0.0, 0.0);

  Jnd = bessel_1_complex((n-1), arg) - static_cast<double>(n)*bessel_1_complex(n, arg)/arg;

  return Jnd;

}



//========================================================================
//========================================================================
//
//	NAME:	complex<double> hankel_1(int n, double arg)
//	DESC:	Calculates the Hankel function of the first kind of order n
//
//	INPUT:
//		int n:: Function order
//		double arg:: Function argument
//	
//	OUTPUT:
//		complex<double> Hn:: Hankel function of order n
//		
//	NOTES:	1) The Hankel function of the first kind is defined as:
//			Hn(x) = Jn(x) + j*Yn(x)
//
//========================================================================
//========================================================================
complex<double> hankel_1(int n, double arg)
{

  complex<double> Hn(0.0, 0.0);

  Hn = bessel_1(n, arg) + COMPLEXJ*bessel_2(n, arg);

  return Hn;

}


//========================================================================
//========================================================================
//
//	NAME:	static void get_g_2d(double k0, double xpos, double xposp, COMPLEX_DOUBLE& g)
//	DESC:	Gets the 2D vacuum Green's function.
//
//	INPUT:
//	
//	OUTPUT:
//		g:: the minimum value of two real numbers
//		
//	NOTES:	i. !!! this is a 2d Green's function, but it seems to depend only
//		on x
//		ii. this is the Green's function in vacuum so we use k0
//
//========================================================================
//========================================================================
void get_g(double k0, double xpos, double xposp, complex<double>& g)
{

  // get the free space wavelength
  double lambda = 2.0*PI/k0;

  // get the distance between x and x'
  double x = fabs(xpos - xposp);

  // get the n=0 Hankel function of the first kind with k0*x as the arguments
  complex<double> h0;
  get_hn(0, 1, k0*x, h0);


  g = (PI/lambda)*h0;

  return;
}


//========================================================================
//========================================================================
//
//	NAME:	static void get_hm(int n, int m, double z, COMPLEX_DOUBLE& hm)
//	DESC:	Gets the n order Hankel function of the mth kind with argument z
//
//	INPUT:
//	
//	OUTPUT:
//		g:: the minimum value of two real numbers
//		
//	NOTES:	i. the Hankel functions of the first and second kind are defined
//		as:
//
//			Hn_1(z) = Jn(z) + jYn(z)
//			Hn_2(z) = Jn(z) - jYn(z)
//
//========================================================================
//========================================================================
void get_hn(int n, int m, double z, complex<double>& hn)
{


  // GET THE BESSEL FUNCTIONS
  double jn, yn;
  jn = bessel_1(n, z);
  yn = bessel_2(n, z);


  // ASSIGN THE HANKEL FUNCTION

  // Hankel function of the first kind ...
  if(m == 1)
  {
    hn = jn + COMPLEXJ*yn;
  }
  // Hankel function of the second kind ...
  else if(m == 2)
  {
    hn = jn - COMPLEXJ*yn;
  }
 

  return;
}




double digamma(int m)
{

  double digamma;

  int k;


  digamma = -EULERC;


  for(k = 1; k <= m-1; ++k)
  {
    digamma += 1.0/static_cast<double>(k);
  }
 

  return digamma;
}


/*

double get_Pl(int l, double x);
double dfactorial2(int x);

  double Rinv;


cout << "here." << endl;
  if(R < 3.0)
  {

    double sr = get_length(r);
    double srp = get_length(r0);
    
    double rless, rgreat;
    if(sr > srp)
    {
      rless = srp;
      rgreat = sr;
    }
    else
    {
      rless = sr;
      rgreat = srp;
    }

    double s = sqrt(r.x*r.x + r.y*r.y);
    double sp = sqrt(r0.x*r0.x + r0.y*r0.y);

    double theta = acos(r.z/sr);
    double thetap = acos(r0.z/srp);

    double phi, phip;

    if(r.x < 0.0)
    {
      phi = PI - asin(r.y/s);
    }
    else
    {
      phi = asin(r.y/s);
    }

    if(r0.x < 0.0)
    {
      phip = PI - asin(r0.y/sp);
    }
    else
    {
      phip = asin(r0.y/sp);
    }

    int l, lmax = 10;
    Rinv = 0.0;
cout << "here2." << endl;
    double cosgamma = cos(theta)*cos(thetap) + sin(theta)*sin(thetap)*cos(phi - phip);
    for(l = 0; l < lmax; ++l)
    {
      cout << "get_Pl(l, cosgamma): " << get_Pl(l, cosgamma) << endl;
      Rinv += ( pow(rless, l)/pow(rgreat, l+1) )*get_Pl(l, cosgamma);
    }

  cout << "Rinv: " << Rinv << endl;
   cout << "1/r: " << 1.0/R << endl; 
  int gh;
  cin >> gh;
  }
  else
  {
    Rinv = 1.0/R;
  }
    //Rinv = 1.0/R;



//========================================================================
//========================================================================
//
//	NAME:	complex<double> get_Gmu(VECTOR_DOUBLE r, VECTOR_DOUBLE r0)
//	DESC:	Read in initial parameters
//
//	NOTES:
//
//
//========================================================================
//========================================================================
double get_Pl(int l, double x)
{
  double Pl = 0.0;

  int k;
  double dl = static_cast<double>(l);
  int krange = static_cast<int>(floor(dl/2.0));
 
  for(k = 0; k <= krange; ++k)
  {
    Pl += pow(-1.0, k)*(dfactorial2(2*(l - k)))*pow(x, l-2*k)/(dfactorial2(k)*dfactorial2(l-k)*dfactorial2(l-2*k));
  }

  Pl /= pow(2.0, l);

  return Pl;
}

//========================================================================
//========================================================================
//
//	NAME:	complex<double> get_Gmu(VECTOR_DOUBLE r, VECTOR_DOUBLE r0)
//	DESC:	Read in initial parameters
//
//	NOTES:
//
//
//========================================================================
//========================================================================
double dfactorial2(int x)
{

  if(x <= 1)
  {
    return 1.0;
  }

  return dfactorial2(x-1)*static_cast<double>(x);
}
*/
