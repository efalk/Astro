
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

/* Find the coordinates of the Sun, from Astronomical Formulae for Calculators,
 * by Jean Meeus, 4th edition, chapter 18.
 *
 * References:
 *
 * [1] Astronomical Formulae for Calculators, * by Jean Meeus, 4th edition
 * [2] Atronomical Algorithms by Jean Meeus, 2nd edition.
 *
 * Meeus also gives formulae for computing the Sun's position relative
 * to the standard ecliptic of 1950, but I don't feel like doing that
 * right now.
 *
 * These functions return the mean position of the Sun.  I have not
 * written any code for the apparent position.
 *
 *
 *
 */

static inline double sind(double a) { return sin(a*RAD); }
static inline double cosd(double a) { return cos(a*RAD); }
static inline double asind(double x) { return asin(x) * DEG; }
static inline double acosd(double x) { return acos(x) * DEG; }
static inline double atan2d(double y, double x) { return atan2(y,x) * DEG; }


/**
 * Return lat,lon,radius from the date.  lat,lon are in degrees,
 * radius is in AU.  This is the Sun's coordinates relative to
 * the Earth in ecliptic coordinates.  Take the complement of
 * lat,lon to get the Earth's coordinates relative to the Sun.
 * Note: by definition, lat is always 0 for the mean ecliptic,
 * but I may change this function later to include perturbations.
 */
void
SunEcliptic(double jdate, double *lat, double *lon, double *rad)
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

/**
 * Return declination, right ascension and distance of the
 * Sun for a given date.  See [2], ch 25.  Accurate to about 0.01
 * degree.  If you need more accuracy, see [2], ch 26 or use the Vosp87
 * values.  RA given in hours, declination given in degrees.
 */
void
SunEquatorial(double jdate, double *decl, double *RA, double *rad)
{
	double	T,T2,T3 ;	/* time, in centuries, squared, cubed */
	double L0, M;		/* Mean longitude, Mean anomoly */
	double e;		/* Eccentricity */
	double C;		/* Equation of Center */
	double lon, v;		/* True longitude, true anomoly */
	double obl;		/* Obliqiuty of the ecliptic, see [2] 22.2 */
	double psi, eps;

	T = (jdate-JD2000)/36525 ; T2 = T*T ; T3 = T*T*T ;
	L0 = limitAngle(280.46646 + 36000.76983 * T + 0.0003032 * T2);
	M = limitAngle(357.52911 + 35999.05029 * T - 0.0001537 * T2);
	e = 0.016708634 - 0.000042037 * T - 0.0000001267 * T2;
	C = (1.914602 - 0.004817*T - 0.000014 * T2) * sind(M) +
	    (0.019993 - 0.000101 * T) * sind(2*M) +
	    0.000289 * sind(3*M);
	lon = L0 + C;
	v = M + C;
	*rad = 1.000001018 * (1 - e*e) / (1 + e * cosd(v));
	obl = 23.439291111 - (46.8150/3600) * T - (.00059/3600) * T2 + 
		(.001813/3600) * T3;
	nutation(&psi, &eps, jdate);	/* TODO: eps is not right */
	obl += eps/3600;
	*RA = limitAngle(atan2d(cosd(obl)*sind(lon), cosd(lon))) * (24./360.);
	*decl = asind(sind(obl) * sind(lon));
}
