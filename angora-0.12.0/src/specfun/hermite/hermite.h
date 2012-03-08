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

#ifndef HERMITE_H
#define HERMITE_H

//function for calculating (physicist's) Hermite polynomials : http://en.wikipedia.org/wiki/Hermite_polynomials


double Hermite(const int& N, const double& x)

void HermiteCoeff(const int& N, Array<int,1>& coeffs);

#endif // HERMITE_H
