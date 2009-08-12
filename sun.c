
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

/* Find the coordinates of the Sun, from Astronomical Formulae for Calculators,
 * by Jean Meeus, 4th edition, chapter 18.
 *
 * Meeus also gives formulae for computing the Sun's position relative
 * to the standard ecliptic of 1950, but I don't feel like doing that
 * right now.
 *
 * These functions return the mean position of the Sun.  I have not
 * written any code for the apparent position.
 *
 *
 * void
 * SunEcliptic(double jdate, double *lat, double *lon, double *rad)
 *	Return lat,lon,radius from the date.  lat,lon are in degrees,
 *	radius is in AU.  This is the Sun's coordinates relative to
 *	the Earth in ecliptic coordinates.  Take the complement of
 *	lat,lon to get the Earth's coordinates relative to the Sun.
 *	Note: by definition, lat is always 0 for the mean ecliptic,
 *	but I may change this function later to include perturbations.
 *
 * void
 * SunEquatorial(double jdate, double *decl, double *RA, double *rad)
 *	Return declination, right ascension and distance of the
 *	Sun for a given date.
 */


void
SunEcliptic(jdate, lat,lon,rad)
	double	jdate ;
	double	*lat, *lon, *rad ;
{
	double	T,T2,T3 ;	/* time, in centuries */
	double	L ;		/* geometric mean longitude of Sun */
	double	M ;		/* Sun's mean anomoly */
	double	e ;		/* eccentricity of Earth's orbit */
	double	v ;		/* true anomoly */
	double	C ;		/* Sun's equation of center */
	double	R ;		/* Sun's radius vector */

	T = (jdate-2415020.0)/36525 ; T2 = T*T ; T3 = T*T*T ;
	L = 279.69668 + 36000.76892*T + .0003025*T2 ;
	M = 358.47583 + 35999.04975*T - .000150*T2 - .0000033*T3 ;
	e = .01675104 - .0000418*T - .000000126*T2 ;

	L *= RAD ;
	M *= RAD ;

	C = + (1.919460 - .004789*T - .000014*T2) * sin(M)
	    + (0.020094 - .000100*T) * sin(2*M)
	    + 0.000293 * sin(3*M) ;
	C *= RAD ;

	L += C ;
	v = M + C ;

	R = ( 1.000002 * (1 - e*e) ) / ( 1 + e*cos(v) ) ;

	L *= DEG ;
	if( L > 360. ) L -= (int)L/360*360 ;

	*lat = 0. ;
	*lon = L ;
	*rad = R ;
}


void
SunEquatorial(jdate, decl,RA,rad)
	double	jdate ;
	double	*decl, *RA, *rad ;
{
	double	T,T2,T3 ;	/* time, in centuries */
	double	lat,lon ;
	double	ecl ;		/* obliquity of ecliptic */
	double	d,r ;		/* declination, right ascension */

	T = (jdate-2415020.0)/36525 ; T2 = T*T ; T3 = T*T*T ;

	SunEcliptic(jdate, &lat,&lon,rad) ;
	lat *= RAD ; lon *= RAD ;

	ecl = obliquity(jdate) ;
	ecl *= RAD ;

	r = atan2( cos(ecl)*sin(lon), cos(lon) ) * DEG ;
	d = asin(sin(ecl)*sin(lon)) * DEG ;
	if( r < 0. ) r += 360. ;

	*decl = d ;
	*RA = r*24/360 ;
}
