
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

/* Precession is the slow rotation of the Earth's axis of about 3 seconds
 * of right ascension per year.  One full rotation of the equinoxes takes
 * about 26,000 years.  Source: Meeus, chap. 14
 *
 * References:
 *
 * [1] Astronomical Formulae for Calculators, * by Jean Meeus, 4th edition
 * [2] Astronomical Algorithms by Jean Meeus, 2nd edition.
 *
 * Nutation is the small elliptical wobble in the Earth's axis caused
 * by the Moon and other influences.  Nutation has a period of about 18.6
 * years, and an amplitude of about 9.2 seconds of arc.
 * Source: [1], chap 15, [2] chap 22
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

#define	NA(a)	(sizeof(a)/sizeof(a[0]))


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

#if HIGH_PRECISION
/* Table 22.A from [2]; units are .0001".  Coefficients smaller than
 * 0.0003" have been omitted.
 */
typedef struct {
  int D, M, Mm, F, omega;
  int s0;
  double s1;
  int c0;
  double c1;
} NutationCoeffs;

static const NutationCoeffs nutationCoeffs[] = {
  {  0,  0,  0,  0,  1, -171996,  -174.2, 92025,  8.9 },
  { -2,  0,  0,  2,  2,  -13187,    -1.6,  5736, -3.1 },
  {  0,  0,  0,  2,  2,   -2274,    -0.2,   977, -0.5 },
  {  0,  0,  0,  0,  2,    2062,     0.2,  -895,  0.5 },
  {  0,  1,  0,  0,  0,    1426,    -3.4,    54, -0.1 },
  {  0,  0,  1,  0,  0,     712,     0.1,    -7,    0 },
  { -2,  1,  0,  2,  2,    -517,     1.2,   224, -0.6 },
  {  0,  0,  0,  2,  1,    -386,    -0.4,   200,    0 },
  {  0,  0,  1,  2,  2,    -301,       0,   129, -0.1 },
  { -2, -1,  0,  2,  2,     217,    -0.5,   -95,  0.3 },
  { -2,  0,  1,  0,  0,    -158,       0,     0,    0 },
  { -2,  0,  0,  2,  1,     129,     0.1,   -70,    0 },
  {  0,  0, -1,  2,  2,     123,       0,   -53,    0 },
  {  2,  0,  0,  0,  0,      63,       0,     0,    0 },
  {  0,  0,  1,  0,  1,      63,     0.1,   -33,    0 },
  {  2,  0, -1,  2,  2,     -59,       0,    26,    0 },
  {  0,  0, -1,  0,  1,     -58,    -0.1,    32,    0 },
  {  0,  0,  1,  2,  1,     -51,       0,    27,    0 },
  { -2,  0,  2,  0,  0,      48,       0,     0,    0 },
  {  0,  0, -2,  2,  1,      46,       0,   -24,    0 },
  {  2,  0,  0,  2,  2,     -38,       0,    16,    0 },
  {  0,  0,  2,  2,  2,     -31,       0,    13,    0 },
  {  0,  0,  2,  0,  0,      29,       0,     0,    0 },
  { -2,  0,  1,  2,  2,      29,       0,   -12,    0 },
  {  0,  0,  0,  2,  0,      26,       0,     0,    0 },
  { -2,  0,  0,  2,  0,     -22,       0,     0,    0 },
  {  0,  0, -1,  2,  1,      21,       0,   -10,    0 },
  {  0,  2,  0,  0,  0,      17,    -0.1,     0,    0 },
  {  2,  0, -1,  0,  1,      16,       0,    -8,    0 },
  { -2,  2,  0,  2,  2,     -16,     0.1,     7,    0 },
  {  0,  1,  0,  0,  1,     -15,       0,     9,    0 },
  { -2,  0,  1,  0,  1,     -13,       0,     7,    0 },
  {  0, -1,  0,  0,  1,     -12,       0,     6,    0 },
  {  0,  0,  2, -2,  0,      11,       0,     0,    0 },
  {  2,  0, -1,  2,  1,     -10,       0,     5,    0 },
  {  2,  0,  1,  2,  2,     -8,        0,     3,    0 },
  {  0,  1,  0,  2,  2,      7,        0,    -3,    0 },
  { -2,  1,  1,  0,  0,     -7,        0,     0,    0 },
  {  0, -1,  0,  2,  2,     -7,        0,     3,    0 },
  {  2,  0,  0,  2,  1,     -7,        0,     3,    0 },
  {  2,  0,  1,  0,  0,      6,        0,     0,    0 },
  { -2,  0,  2,  2,  2,      6,        0,    -3,    0 },
  { -2,  0,  1,  2,  1,      6,        0,    -3,    0 },
  {  2,  0, -2,  0,  1,     -6,        0,     3,    0 },
  {  2,  0,  0,  0,  1,     -6,        0,     3,    0 },
  {  0, -1,  1,  0,  0,      5,        0,     0,    0 },
  { -2, -1,  0,  2,  1,     -5,        0,     3,    0 },
  { -2,  0,  0,  0,  1,     -5,        0,     3,    0 },
  {  0,  0,  2,  2,  1,     -5,        0,     3,    0 },
  { -2,  0,  2,  0,  1,      4,        0,     0,    0 },
  { -2,  1,  0,  2,  1,      4,        0,     0,    0 },
  {  0,  0,  1, -2,  0,      4,        0,     0,    0 },
  { -1,  0,  1,  0,  0,     -4,        0,     0,    0 },
  { -2,  1,  0,  0,  0,     -4,        0,     0,    0 },
  {  1,  0,  0,  0,  0,     -4,        0,     0,    0 },
  {  0,  0,  1,  2,  0,      3,        0,     0,    0 },
  {  0,  0, -2,  2,  2,     -3,        0,     0,    0 },
  { -1, -1,  1,  0,  0,     -3,        0,     0,    0 },
  {  0,  1,  1,  0,  0,     -3,        0,     0,    0 },
  {  0, -1,  1,  2,  2,     -3,        0,     0,    0 },
  {  2, -1, -1,  2,  2,     -3,        0,     0,    0 },
  {  0,  0,  3,  2,  2,     -3,        0,     0,    0 },
  {  2, -1,  0,  2,  2,     -3,        0,     0,    0 },
};
#endif

/**
 * Given a julian date, return the nutation in longitude and
 * nutation in obliquity.  The first is is along the ecliptic
 * and the second is perpindicular to it.
 * Return values are in seconds of arc.
 * Source: [2], ch. 22
 *
 * @param psi    Returned nutation in longitude, arcseconds
 * @param eps    Returned nutation of obliquity of the eliptic, arcseconds
 * @param jdate  Julian day
 */
void
nutation(double *psi, double *eps, double jdate)
{
	double	T, T2, T3;	/* centuries since 2000 */
	double	D;		/* Mean elongation of the Moon from the Sun */
	double	M;		/* Mean anomoly of the Sun */
	double	Mm;		/* Mean anomoly of the Moon */
	double	F;		/* Moon's argument of longitude */
	double	om;		/* Longitude of ascending node of the moon */

	/* TODO: convert JD to JDE */

	T = (jdate - JD2000)/36525.; T2 = T*T; T3 = T2*T;
	D = limitAngle(297.85036 + 445267.111480*T - 0.0019142*T2 + T3/189474);
	M = limitAngle(357.52772 + 35999.050340*T - 0.0001603*T2 - T3/300000);
	Mm = limitAngle(134.96298 + 477198.867398*T + 0.0086972*T2 + T3/56250);
	F = limitAngle(93.27191 + 483202.017538*T - 0.0036825*T2 + T3/327270);
	om = limitAngle(125.04452 - 1934.136261*T + 0.0020708*T2 + T3/450000);

#if HIGH_PRECISION
	{
	    double p=0, e=0;
	    int i;
	    const NutationCoeffs *nc = nutationCoeffs;
	    for( i=0; i < NA(nutationCoeffs); ++i, ++nc) {
		double arg = nc->D * D + nc->M * M + nc->Mm * Mm +
		    nc->F * F + nc->omega * om;
		p += (nc->s0 + nc->s1*T) * sind(arg);
		e += (nc->c0 + nc->c1*T) * cosd(arg);
	    }
	    *psi = p * 0.0001;
	    *eps = e * 0.0001;
	}
#else
	{
	    double L, LL;		/* Mean longitudes of Sun and Moon */
	    L = limitAngle(280.4665 + 36000.7698 * T);
	    LL = limitAngle(218.3165 + 481267.8813 * T);
	    *psi = -17.2*sind(om) - 1.32*sind(2*L) -
		    0.23*sind(2*LL) + 0.21*sind(2*om);

	    *eps = 9.2*cosd(om) + 0.57*cosd(2*L) +
		    0.10*cosd(2*LL) + 0.09*cosd(2*om);
	}
#endif
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

