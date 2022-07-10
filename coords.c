
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

/* Coordinate conversion routines, from Astronomical Formulae for Calculators,
 * by Jean Meeus, 4th edition, chapter 8.
 *
 *
 * double
 * obliquity(double date)
 *	Find obliquity of eccliptic for a specified date.
 *
 * void
 * equat2ecliptic(double decl,RA, double *lat,*lon, double jdate)
 *	convert equatorial coordinates to ecliptic coordinates.
 *
 * void
 * ecliptic2equat(double lat,lon, double *decl,*RA, double jdate)
 *	convert ecliptic coordinates to equatorial
 *
 * void
 * equat2bearings(double decl,RA, double lat,lon, double A,h, double jdate)
 *	convert equatorial coordinates to observer's local horizontal coords
 *
 *
 * TODO: rise/set times from chapter 8.
 */

/**
 * Find obliquity of eccliptic for a specified date.  Accurate to within
 * 1" for +/- 2000 years.  See Meeus, ch. 22 for a more accurate formula
 * if you need it.  Does not account for nutation.
 *
 * @param jdate  Julian day
 * @return  Obliquity, degrees
 */
double
obliquity(double jdate)
{
	double T, T2, T3;

	/* Meeus, formula 22.2 */

	T = (jdate-JD2000)/36525 ; T2 = T*T ; T3 = T*T*T ;
	return EarthTilt - (46.8150/3600) * T - (.00059/3600) * T2 + 
		(.001813/3600) * T3;
}

/**
 * convert equatorial coordinates to ecliptic coordinates.
 */
void
equat2ecliptic(double decl, double RA, double *lat, double *lon, double jdate)
{
	double	e ;	/* obliquity of ecliptic */
	double	la, lo ;

	e = obliquity(jdate) ;

	/* convert to radians */
	e *= RAD ;
	decl *= RAD ;
	RA *= RAD*360/24 ;


	lo = atan2( sin(RA)*cos(e)+tan(decl)*sin(e), cos(RA) ) ;
	la = asin( sin(decl)*cos(e) - cos(decl)*sin(e)*sin(RA) ) ;

	/* convert to degrees */
	lo *= DEG ;
	la *= DEG ;

	*lon = lo ;
	*lat = la ;
}


/**
 * Convert ecliptic coordinates to equatorial
 */
void
ecliptic2equat(double lat, double lon, double *decl, double *RA, double jdate)
{
	double	e ;	/* obliquity of ecliptic */
	double	d,r ;

	e = obliquity(jdate) ;

	/* convert to radians */
	e *= RAD ;
	lat *= RAD ;
	lon *= RAD ;

	/* rotate everything 23 degrees about X axis.  In 3-d coords,
	 * we can do this as:
	 *
	 *	x = cos(lon) * cos(lat) ;
	 *	y = sin(lon) * cos(lat) ;
	 *	z =	       sin(lat) ;
	 *
	 *	xe = x ;
	 *	ye = y*cos(e) - z*sin(e) ;
	 *	ze = y*sin(e) + z*cos(e) ;
	 *
	 *	r = atan2(ye,xe) ;
	 *	d = atan2(ze, sqrt(xe*xe+ye*ye)) ;
	 *
	 * but that can be simplified.  The following two lines come
	 * from Meeus.
	 */

	r = atan2( sin(lon)*cos(e) - tan(lat)*sin(e), cos(lon) ) ;
	d = asin( sin(lat)*cos(e) + cos(lat)*sin(e)*sin(lon) ) ;

	/* convert to hours, degrees */
	r *= DEG*24/360 ;
	d *= DEG ;

	*decl = d ;
	*RA = r ;
}

/**
 * Convert equatorial coordinates to observer's local horizontal coords
 * @param decl, RA  object's equatorial coords
 * @param lat, lon  observer's lat,lon
 * @param A, H      return: observer's bearings to object, azimuth, elevation
 * @param jdate
 * @param time      time, in hours GMT
 */
void
equat2bearings(double decl, double RA, double lat, double lon,
	double *A, double *H, double jdate, double time)
{
	double	a, h ;
	double	ha ;		/* hour angle */
	double	st ;		/* sidereal time at Greenwich */

	/* convert to radians */
	decl *= RAD ;
	RA *= RAD*360/24 ;
	lat *= RAD ;
	lon *= RAD ;

	/* get sidereal time at Greenwich, convert to radians */
	st = julianTime2sidereal(jdate,time) ;
	st *= RAD*360/24 ;

	/* convert to local hour angle */
	ha = st - lon - RA ;

	a = atan2( sin(ha), cos(ha)*sin(lat) - tan(decl)*cos(lat) ) ;
	h = asin( sin(lat)*sin(decl) + cos(lat)*cos(decl)*cos(ha) ) ;

	/* convert to degrees */
	a *= DEG ;
	h *= DEG ;

	*A = a ;
	*H = h ;
}



#ifdef	STANDALONE

main()
{
	double	date,time ;
	double	lat,lon, decl,RA ;
	double	A,h ;

	date = date2julian(1950,1,1) ;
	decl = hms2h(28,8,55.11) ;
	RA = hms2h(7,42,15.525) ;
	printf("%.6f, %.6f\n", decl,RA) ;

	equat2ecliptic(decl,RA, &lat,&lon, date) ;
	printf("%.6f, %.6f\n", lat,lon) ;

	ecliptic2equat(lat,lon, &decl,&RA, date) ;
	printf("%.6f, %.6f\n", decl,RA) ;

	RA = hms2h(10,57,35.681) ;
	decl = hms2h(8,25,58.1) ;
	lon = -hms2h(0,17,25.94)*360/24 ;
	lat = hms2h(50,47,55.0) ;
	date = date2julian(1978,11,13) ;
	time = hms2h(4,34,0.) ;

	equat2bearings(decl,RA, lat,lon, &A,&h, date,time) ;
	printf("%.6f, %.6f\n", A,h) ;

	exit(0) ;
}

#endif
