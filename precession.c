
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

/* Precession is the slow rotation of the Earth's axis of about 3 seconds
 * of right ascension per year.  One full rotation of the equinoxes takes
 * about 26,000 years.  Source: Meeus, chap. 14
 *
 * Nutation is the small elliptical wobble in the Earth's axis caused
 * by the Moon and other influences.  Nutation has a period of about 18.6
 * years, and an amplitude of about 9.2 seconds of arc.
 * Source: Meeus, chap 15.
 *
 * void
 * precessionRate(double decl,RA, double *dd,*dr, double jdate)
 *	Accepts declination and right ascension in degrees and hours
 *	respectively, and returns the annual differences, also in
 *	degrees and hours
 *
 * precession(double decl0,RA0,jdate0, double *decl1,*RA1,jdate1)
 *	Accepts declination and right ascension for one date, and
 *	returns declination and right ascension for another date.
 *
 * nutation(double *psi, *eps, double jdate)
 *	Given a julian date, return the nutation in longitude and
 *	nutation in obliquity.  I have no idea what these really
 *	mean, but they're used as inputs to other functions.
 *	Return values are in seconds of arc.
 */



void
precessionRate(decl,RA, dd,dr, jdate)
	double	decl,RA, *dd,*dr, jdate ;
{
	double	T = (jdate-2415020)/36525 ;
	double	m,n ;
	double	d,r ;

	m = 3.07234 + .00186*T ;	/* seconds RA/year */
	n = 20.0468 - .0085*T ;		/* seconds of arc */

	/* convert to radians */
	decl *= RAD ;
	RA *= RAD*360/24 ;
	n *= (1./3600.)*RAD ;
	m *= (1./3600.)*RAD*360/24 ;

	r = m + n*sin(RA)*tan(decl) ;
	d = n*cos(RA) ;

	/* convert back to degrees and hours */
	*dd = d*DEG ;
	*dr = r*DEG*24/360 ;
}


void
precessionRad(decl0,RA0, jdate0, decl1,RA1, jdate1)
	double	decl0, RA0, jdate0 ;
	double	*decl1, *RA1, jdate1 ;
{
	double	T0, T ;		/* note: tropical centuries here */
	double	tau, z, theta ;
	double	T2,T3 ;
	double	A,B,C ;

	/* convert dates to tropical centuries relative to 1900 */
	T0 = (jdate0 - 2415020.313) / 36524.2199 ;
	T = (jdate1 - jdate0) / 36524.2199 ;
	T2 = T*T ;
	T3 = T2*T ;

	/* these units are seconds of arc. */
	tau = (2304.250 + 1.396*T0)*T + .302*T2 + .018*T3 ;
	z = tau + .791*T2 + .001*T3 ;
	theta = (2004.682 - .83*T0)*T - .426*T2 - .042*T3 ;

	tau *= (1./3600.)*RAD ;
	z *= (1./3600.)*RAD ;
	theta *= (1./3600.)*RAD ;

	A = cos(decl0)*sin(RA0 + tau) ;
	B = cos(theta)*cos(decl0)*cos(RA0+tau) - sin(theta)*sin(decl0) ;
	C = sin(theta)*cos(decl0)*cos(RA0+tau) + cos(theta)*sin(decl0) ;

	RA0 = atan2(A,B) ; RA0 += z ;
	decl0 = asin(C) ;

	*decl1 = decl0 ;
	*RA1 = RA0 ;
}


void
precession(decl0,RA0, jdate0, decl1,RA1, jdate1)
	double	decl0, RA0, jdate0 ;
	double	*decl1, *RA1, jdate1 ;
{
	/* convert to radians */
	decl0 *= RAD ;
	RA0 *= RAD*360/24 ;

	precessionRad(decl0,RA0, jdate0, decl1,RA1, jdate1) ;

	/* convert back to degrees and hours */
	*decl1 *= DEG ;
	*RA1 *= DEG*24/360 ;
}


void
nutation(psi, eps, jdate)
	double	*psi, *eps, jdate ;
{
	double	T, T2 ;		/* centuries since 1900 */
	double	L,Lm ;		/* mean longitude of Sun and Moon */
	double	M,Mm ;		/* mean anomaly of Sun and Moon */
	double	ohm ;		/* longitude of Moon's ascending node */

	T = (jdate - 2415020.0)/36525. ;
	T2 = T*T ;

	/* Meeus gave these equations (constants in degrees) */
	L =   (279.6967 +  36000.7689*T + .000303*T2) * RAD ;
	Lm =  (270.4342 + 481267.8831*T - .001133*T2) * RAD ;
	M =   (358.4758 +  35999.0498*T - .000150*T2) * RAD ;
	Mm =  (296.1046 + 477198.8491*T + .009192*T2) * RAD ;
	ohm = (259.1833 -   1934.1420*T + .002078*T2) * RAD ;

	/* Meeus gave these equations (constants in seconds of arc) */

	*psi =	- (17.2327 + .01737*T)	* sin(ohm)	/* 6798 days */
		- (1.2729 + 0.00013*T)	* sin(2*L)	/* 182.62 days */
		+ 0.2088		* sin(2*ohm)
		- 0.2037		* sin(2*Lm)
		+ (0.1261 - 0.00031*T)	* sin(M)
		+ 0.0675		* sin(Mm)
		- (0.0497 - 0.00012*T)	* sin(2*L + M)
		- 0.0342		* sin(2*Lm - ohm)
		- 0.0261		* sin(2*Lm + Mm)
		+ 0.0214		* sin(2*L - M)
		- 0.0149		* sin(2*L - 2*Lm + Mm)
		+ 0.0124		* sin(2*L - ohm)
		+ 0.0114		* sin(2*Lm - Mm) ;

	*eps =	+ (9.2100 + 0.00091*T)	* cos(ohm)	/* 6798 days */
		+ (0.5522 - 0.00029*T)	* cos(2*L)	/* 182.62 days */
		- 0.0904		* cos(2*ohm)
		+ 0.0884		* cos(2*Lm)
		+ 0.0216		* cos(2*L + M)
		+ 0.0183		* cos(2*Lm - ohm)
		+ 0.0113		* cos(2*Lm + Mm)
		- 0.0093		* cos(2*L - M)
		- 0.0066		* cos(2*L - ohm) ;
}


	/* compatibility routines from xephem */

void
precess(mjd1, mjd2, ra, dec)
	double	mjd1, mjd2 ;
	double	*ra, *dec ;
{
	mjd1 += JD1900 ;
	mjd2 += JD1900 ;
	precessionRad(*dec,*ra,mjd1, dec,ra,mjd2) ;
}


#ifdef	STANDALONE

main()
{
	double	date,time ;
	double	dd,dr, decl,RA ;
	double	decl1,RA1, date1 ;
	double	psi, epsilon ;

	date = date2julian(1978,1,1) ;
	RA = hms2h(10,5,42.7) ;
	decl = hms2h(12,12,45.) ;
	precessionRate(decl,RA, &dd,&dr, date) ;
	printf("%.6f, %.6f => %f,%f @ %f\n", RA,decl, dr*3600,dd*3600, date) ;

	RA = hms2h(2,40,46.276) ;
	decl = hms2h(49,01,06.45) ;
#ifdef	COMMENT
	date = date2julian(1950,1,1) ;
#endif	/* COMMENT */
	date = JD1950 ;
	date1 = date2julian(1978,11,13) + .19 ;
	RA += 28.8665*.0342/3600. ;
	decl -= 28.8665*.083/3600. ;
	precession(decl,RA,date, &decl1,&RA1,date1) ;
	printf("%f,%f @ %f => %s", RA,decl,date, convertHms(RA1)) ;
	printf(",%s @ %f\n", convertHms(decl1),date1) ;

	date = date2julian(1978,11,13) + hms2h(4,35,0)/24. ;
	nutation(&psi, &epsilon, date) ;
	printf("%f, %f\n", psi, epsilon) ;

	exit(0) ;
}

#endif

