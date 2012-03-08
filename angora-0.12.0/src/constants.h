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

#ifndef CONSTANTS_H
#define CONSTANTS_H

//general physical and mathematical constants

const complex<double> ii (0,1.0);		//imaginary number 0+1.0i
const double c = 299792458.;			//speed of light in free space
const double mu_0 = 4*M_PI*1e-7;		//permeability of free space
const double epsilon_0 = 1/(c*c*mu_0);	//permittivity of free space
//const double mu_0 = 1/c;		//permeability of free space
//const double epsilon_0 = 1/c;	//permittivity of free space
const double eta_0 = mu_0*c;			//wave impedance of free space

#endif
