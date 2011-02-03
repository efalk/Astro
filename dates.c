
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

/* Date conversion routines, from Astronomical Formulae for Calculators,
 * by Jean Meeus, 4th edition [1] and Astronimical Algorithms by Jean Meeus,
 * 2nd edition [2].
 *
 * See below for more detailed explanations.
 *
 * double
 * date2julian(int y, int m, int d)
 *	Convert (year, month, day) to a julian date
 *
 * double
 * time2julian(int y, int m, int d, int h, int m, double s)
 *	Convert (year, month, day, h, m, s) to a julian date
 *
 * void
 * julian2date(double jdate, int *y, int *m, int *d)
 *	Convert a julian date to y/m/d
 *
 * void
 * julian2time(double jdate, int *y, int *m, int *d, int *h, int *m, double *s)
 *	Convert a julian date to y/m/d h:m:s GMT
 *
 * double
 * unix2julian(time_t)
 *	Convert unix time to julian
 *
 * double
 * jnow()
 *	Return the current julian date.
 *
 * double
 * julian2sidereal(double jdate)
 *	Return sidereal time in hours for a specific julian date
 *	(Note: this is the sidereal time of midnight GMT.  I.e. draw a
 *	line from the Sun through the center of the Earth and out
 *	to the stars.  Wherever the line passes through the far side of
 *	the Earth, that's where it's midnight.  Finding the sidereal
 *	time at other longitudes will require further conversion.)
 *	TODO: double-check this.  Where did I get it?
 *
 * double
 * julianTime2sidereal(double jdate, double time)
 *	Return sidereal time in hours for a specific julian date
 *	and time other than 0h.  Time is specified in hours since midnight
 *	TODO: double-check this.  Where did I get it?
 *
 * double
 * time2sidereal(double jdate)
 *	Return sidereal time at Grenwich
 *
 * double
 * gmst2gast(double gmst, double JD)
 *	Convert mean sidereal time to apparent sidereal time
 *
 * double
 * siderealMean2Apparent(double jdate)
 *	Return a correction factor to convert mean sidereal time in
 *	hours to apparent sidereal time (correcting for nutation).
 *	meanTime + siderealMean2Apparent() = apparentTime
 *
 * int
 * date2yday(int y, int m, int d)
 *	Obtain the year of the day from y/m/d
 *
 * void
 * yday2date(int y, int yday,  int *m, int *d)
 *	Obtain month and day from year and year of the day
 *
 *
 *
 * A few notes:
 *  GMT is Greenwich Mean Time.
 *  UT is Universal Time.
 *   GMT and UT are both based on the rotation of the Earth, which is not
 *   constant.  These are not uniform times.
 *
 *  Ephemeris time is based on the motion of the planets instead of the
 *   rotation of the Earth.  It has been replaced by Dynamical Time.
 *
 *  Dynamical Time (TD) is based on atomic clocks.  There are actually
 *   two Dynamical Times: Barycentric Dynamical Time, (TDB) based
 *   on the center of mass of the solar system, and Terrestrial
 *   Dynamical Time (TDT).  The difference is caused by relativistic effects
 *   of the Earth's orbit.  Most of the time, you can ignore the difference.
 *   The difference between TD and UT is determined by astronomical
 *   observation.  As of 2004, TD - UT = 64.6 seconds.
 *
 *  Julian Days (JD) are days since the start of the year -4712.
 *
 *  Sidereal time (ST) measures the rotation of the Earth relative to
 *   the vernal equinox rather than the Sun.  One sidereal day is 23 hours,
 *   56 minutes, 4.091 seconds.  This is about .008 seconds shorter than
 *   a day relative to the stars, due to the precession of the equinoxes.
 *
 *
 * Overview:
 *
 * In astronomy, dates are recorded as Julian dates, which are simply
 * days since noon, 1 Jan 4712 BC.  In the Julian system, the new day
 * starts at noon (GMT) instead of midnight.
 *
 * Converting between Julian dates and standard (Gregorian) dates is
 * non-trivial, because of the various changes that have been made in
 * the civil calendar over the years.  For example, in September 1752,
 * 11 days were removed from the calendar to make up for lack of leap year
 * adjustments (see cal(1) manpage).  More importantly, the modern
 * Gregorian calendar was created by Pope Gregory in 1582, which included
 * leap years.
 *
 * Leap years themselves are fairly complex.  The length of the year
 * is not a whole number of days, but closer to 365.2424 days in
 * length (this varies, as the rate of the Earth's rotation is not
 * constant.)  To compensate, an extra day in inserted into the year
 * every fourth year, except for century years when it is not added,
 * except for fourth centuries when it is.  Thus, 1900 was not a leap
 * year, but 2000 will be.
 *
 * Theory of Operations:
 *  Meeus's book doesn't explain much, but this is how I think it
 *  works.  The averge length of a month, ignoring February, is 30.60
 *  days.  To get around the February problem, Meeus "rotates" the
 *  year to put February at the end.  This is done by making January
 *  and February the 13th and 14th months of the year and then subtracting
 *  1 from the year.
 *
 *  Next, if the year is in the Gregorian calendar, the century #, a, is
 *  computed from y/100.  The leap century # is a/4.  The number of skipped
 *  leap years, b, is the number of centuries minus the number of leap
 *  centuries, minus a fudge factor.
 *
 *  Finally, the Julian day is the year*365.25 + month*30.6001 + day +
 *	julian(0,0,0) + b
 *
 * BUG?  This code does not seem to allow for the September 1752 adjustment,
 * unless it's built into the Gregorian conversion.
 */


/**
 * Return Julian day for given date.  Algorithm from [2], ch 7
 */
double
date2julian(int yy, int mm, double d)
{
	int	y = yy ;
	int	m = mm ;
	int	a,b ;

	if( m <= 2 ) {
	  --y ;
	  m += 12 ;
	}

	/* Gregorian calendar started 15 Oct 1582 */
	if( yy*10000 + mm*100 + d > 15821015 ) {
	  a = y/100 ;
	  b = 2 - a + a/4;
	}
	else
	  b = 0 ;

	return (int)(365.25*(y+4716)) + (int)(30.6001*(m+1)) + d + b - 1524.5;
}


double
time2julian(int yy, int mm, int dd,  int hr, int mn, double s)
{
	return date2julian(yy,mm,dd) + hms2h(hr,mn,s)/24. ;
}


void
julian2date(double jdate, int *yy, int *mm, int *dd)
{
	int	z = jdate+.5 ;
	int	a,b,c,d,e ;

	if( z <= 2299161 )
	  a = z ;
	else {
	  a = ((float)z-1867216.25)/36524.25 ;
	  a = z + 1 + a - a/4 ;
	}

	b = a + 1524 ;
	c = (b-122.1)/365.25 ;
	d = 365.25*c ;
	e = (float)(b-d)/30.6001 ;
	*dd = b - d - (int)(30.6001*e) ;
	*mm = e<=13 ? e-1 : e-13 ;
	*yy = *mm > 2 ? c-4716 : c-4715 ;
}


void
julian2time(double jdate,
	int *yy, int *mm, int *dd, int *hr, int *mn, double *s)
{
	double	j ;
	julian2date(jdate, yy,mm,dd) ;
	j = jdate - date2julian(*yy,*mm,*dd) ;
	j *= 24. ;
	*hr = j ; j -= *hr ; j *= 60. ;
	*mn = j ; j -= *mn ; j *= 60. ;
	*s = j ;
}

/**
 * Return the GMT hour for this jdate.
 */
double
julian2hour(double jdate)
{
	int yy, mm, dd;
	julian2date(jdate, &yy,&mm,&dd);
	return (jdate - date2julian(yy,mm,dd)) * 24;
}

/**
 * Print a julian date as yyy-mm-dd hh:mm:ss
 */
void
printDate(double jdate)
{
	int	y3,m3,d3 ;
	int	h,m ;
	double	s ;

	julian2time(jdate, &y3,&m3,&d3, &h,&m,&s) ;
	printf("%4d-%2.2d-%2.2d %d:%2.2d:%.1f\n", y3,m3,d3, h,m,s);
}



double
unix2julian(time_t t)
{
	return JDUnix + (float)t/(24.*60.*60.) ;
}



double
jnow()
{
	time_t	tnow = time(NULL) ;	/* seconds since 1/1/1970 */
	return unix2julian(tnow) ;
}


/**
 * Given a Julian date (not including time),
 * compute the Sidereal time at midnight UT.  Jdate should end in .5.
 * Based on [2], ch. 12
 * @return  Sidereal time in hours
 */
double
julian2sidereal(double jdate)
{
	double	T = (jdate-JD2000)/36525 ;
	double	st ;

	/* Degrees */
	st = 100.46061837 + 36000.770053608 * T +
		0.000387993 * T*T + T*T*T / 38710000;
	st = limitAngle(st);
	return st*24/360;
}


/**
 * Given a Julian date (not including time), and GMT in float hours,
 * compute the Sidereal time.
 */
double
julianTime2sidereal(double jdate, double time)
{
	return limitHour(julian2sidereal(jdate) + time*(366.2422/365.2422));
}


double
time2sidereal(double jdate)
{
	double	rval ;
	rval = julian2sidereal(jdate) ;
	jdate -= (int)jdate ;
	return limitHour( rval + jdate*24 - 12. ) ;
}


double
siderealMean2Apparent(double jdate)
{
	double	dpsi, deps ;
#define	cos_eps	0.9175		/* cos(23d26'30") */

	nutation(&dpsi, &deps, jdate) ;
	return dpsi*cos_eps/15./3600. ;
}

/**
 * @brief convert Greenwich Mean Sidereal Time to Greenwich Apparent
 * Sidereal Time.
 */
double
gmst2gast(double gmst, double JD)
{
	double omega;	/* ascending node of Moon */
	double L;	/* mean longitude of Sun */
	double psi;	/* nutation in longitude */
	double epsilon;	/* obliquity */
	double eqeq;	/* equation of equinoxes */
	double D = JD - 2451545.0;

	omega = 125.04 - 0.052954 * D;
	L = 280.47 + 0.98565 * D;
	epsilon = 23.4393 - 0.0000004 * D;
	psi = -0.000319 * sin(omega*RAD) - 0.000024 * sin(2*L*RAD);
	eqeq = psi * cos(epsilon*RAD);
	return gmst + eqeq;
}


int
date2yday(int y, int m, int d)
{
	/* if it's not a leap year */
	if( y%4 != 0  ||  (y%100 == 0 && y%400 != 0) )
	  return (275*m/9) - 2*((m+9)/12) + d - 30 ;
	else
	  return (275*m/9) - ((m+9)/12) + d - 30 ;
}

void
yday2date(int y, int yday, int *m, int *d)
{
	int	a,b,c,e ;

	/* if it's not a leap year */
	if( y%4 != 0  ||  (y%100 == 0 && y%400 != 0) )
	  a = 1889 ;
	else
	  a = 1523 ;

	b = ((float)yday + a - 122.1)/365.25 ;
	c = yday + a - floor(365.25*b) ;
	e = c/30.6001 ;
	if( e <= 13 )
	  *m = e-1 ;
	else
	  *m = e-13 ;
	*d = c - floor(30.6001*e) ;
}



	/* compatibility routines for xephem */

void
cal_mjd(int m, double d, int y, double *mjd)
{
	int	di ;
	double	df ;

	di = d ;
	df = d-di ;
	*mjd = date2julian(y, m, di) + df ;
}


#ifdef	STANDALONE

#include <stdio.h>
#include <stdlib.h>

static void
test2(double jdate)
{
	int	y3,m3,d3 ;
	int	h,m ;
	double	s ;

	julian2time(jdate, &y3,&m3,&d3, &h,&m,&s) ;
	printf( "%.1f = %4d-%2.2d-%2.2d %d:%2.2d:%.1f\n",
	  jdate, y3,m3,d3, h,m,s) ;
}

static void
test(int y, int m, int d)
{
	double	rval ;
	int	yday ;
	int	y2,m2,d2 ;
	int	y3,m3,d3 ;

	rval = time2julian(y,m,d,0,0,0.) ;
	yday = date2yday(y,m,d) ;
	yday2date(y2=y,yday, &m2,&d2) ;
	julian2date(rval, &y3,&m3,&d3) ;
	printf(
"%4d-%2.2d-%2.2d = %.1f = %4d-%2.2d-%2.2d = %4d-%3.3d = %4d-%2.2d-%2.2d %s",
	  y,m,d, rval, y3,m3,d3, y,yday, y2,m2,d2,
	  y2!=y || m2!=m || d2!=d || y3!=y || m3!=m || d3!=d ? "!" : "") ;
	test2(rval) ;
}


static void
testSt(int y, int m, int d)
{
	double	jd,T, st ;

	jd = date2julian(y,m,d) ;
	T = (jd-JD1900)/36525 ;
	st = julian2sidereal(jd) ;
	printf("jd = %f, T=%f, st=", jd, T) ;
	printHms(st) ;
	putchar('\n') ;
}


int
main()
{
	double	jd,t,st ;

	test2(jnow()) ;
	test(1900,1,1) ;
	test(-4712,1,1) ;
	test(1957,10,4) ;
	test(222,8,10) ;
	test(1978,11,14) ;
	test(1980,4,22) ;
	test(333,1,27) ;
	test(1835,11,16) ;
	test(1910,4,20) ;
	test(-584,5,28) ;
	test(1999,5,17) ;

	testSt(1978,11,13) ;
	testSt(1900,11,13) ;
	testSt(1978,3,21) ;
	testSt(1900,3,21) ;

	jd = date2julian(1978,11,13) ;
	t = hms2h(4,34,0.) ;
	st = julianTime2sidereal(jd,t) ;
	st += siderealMean2Apparent(jd) ;
	printf("jd = %f, st=", jd) ;
	printHms(st) ;
	putchar('\n') ;

	exit(0) ;
}
#endif	/* STANDALONE */
