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

#ifndef FLOAT_COMP_H
#define FLOAT_COMP_H


/** comparison functions for double type **/
/** from: http://stackoverflow.com/questions/17333/most-effective-way-for-float-and-double-comparison **/
inline bool approximatelyEqual(double a, double b, double epsilon)
{
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

inline bool essentiallyEqual(double a, double b, double epsilon)
{
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

inline bool definitelyGreaterThan(double a, double b, double epsilon)
{
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

inline bool definitelyLessThan(double a, double b, double epsilon)
{
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

/** by I. R. Capoglu **/
inline bool lessThanOrEqual(double a, double b, double epsilon)
{
	return (!definitelyGreaterThan(a,b,epsilon));
}

inline bool greaterThanOrEqual(double a, double b, double epsilon)
{
	return (!definitelyLessThan(a,b,epsilon));
}
/** by I. R. Capoglu **/

#endif // FLOAT_COMP_H
