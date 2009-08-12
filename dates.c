
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

/* date conversion routines, from Astronomical Formulae for Calculators,
 * by Jean Meeus, 4th edition.
 *
 *
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
 *	Convert a julian date to y/m/d h:m:s
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
 *	(Note: this is the sidereal time of midnight.  I.e. draw a
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


double
date2julian(yy,mm,dd)
	int	yy,mm,dd ;
{
	int	y = yy ;
	int	m = mm ;
	int	d = dd ;
	int	a,b ;

	if( m <= 2 ) {
	  --y ;
	  m += 12 ;
	}

	/* Gregorian calendar started 15 Oct 1582 */
	if( yy*10000 + mm*100 + dd > 15821015 ) {
	  a = y/100 ;
	  b = a - a/4 - 2 ;
	}
	else
	  b = 0 ;

	return floor(365.25*y) + floor(30.6001*(m+1)) + d + 1720994.5 - b ;
}


double
time2julian(yy,mm,dd, hr,mn,s)
	int	yy,mm,dd ;
	int	hr,mn ;
	double	s ;
{
	return date2julian(yy,mm,dd) + hms2h(hr,mn,s)/24. ;
}


void
julian2date(jdate, yy,mm,dd)
	double	jdate ;
	int	*yy ;
	int	*mm ;
	int	*dd ;
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
julian2time(jdate, yy,mm,dd, hr,mn,s)
	double	jdate ;
	int	*yy ;
	int	*mm ;
	int	*dd ;
	int	*hr, *mn ;
	double	*s ;
{
	double	j ;
	julian2date(jdate, yy,mm,dd) ;
	j = jdate - date2julian(*yy,*mm,*dd) ;
	j *= 24. ;
	*hr = j ; j -= *hr ; j *= 60. ;
	*mn = j ; j -= *mn ; j *= 60. ;
	*s = j ;
}



double
unix2julian(t)
	time_t	t ;
{
	return JDUnix + (float)t/(24.*60.*60.) ;
}



double
jnow()
{
	time_t	tnow = time(NULL) ;	/* seconds since 1/1/1970 */
	return unix2julian(tnow) ;
}


	/* sidereal time at midnight (far side of the Earth) */
double
julian2sidereal(jdate)
	double	jdate ;
{
	double	T = (jdate-2415020)/36525 ;
	double	st ;
	int	i ;

	st = 6.6460656 + 2400.051262*T + .00002581*T*T ;
	i = st/24 ;
	st -= i*24 ;
	return st ;
}


double
julianTime2sidereal(jdate, time)
	double	jdate, time ;
{
	/* TODO: double-check this.  Where did I get it? */

	return julian2sidereal(jdate) + time*1.002737908 ;
}


double
time2sidereal(jdate)
	double	jdate ;
{
	double	rval ;
	rval = julian2sidereal(jdate) ;
	jdate -= (int)jdate ;
	return limitHour( rval + jdate*24 - 12. ) ;
}


double
siderealMean2Apparent(jdate)
	double	jdate ;
{
	double	dpsi, deps ;
#define	cos_eps	0.9175		/* cos(23d26'30") */

	nutation(&dpsi, &deps, jdate) ;
	return dpsi*cos_eps/15./3600. ;
}


int
date2yday(y,m,d)
	int	y,m,d ;
{
	/* if it's not a leap year */
	if( y%4 != 0  ||  y%100 == 0 && y%400 != 0 )
	  return (275*m/9) - 2*((m+9)/12) + d - 30 ;
	else
	  return (275*m/9) - ((m+9)/12) + d - 30 ;
}

void
yday2date(y,yday, m,d)
	int	y,yday ;
	int	*m,*d ;
{
	int	a,b,c,e ;

	/* if it's not a leap year */
	if( y%4 != 0  ||  y%100 == 0 && y%400 != 0 )
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
cal_mjd(m,d,y,mjd)
	int	m,y ;
	double	d ;
	double	*mjd ;
{
	int	di ;
	double	df ;

	di = d ;
	df = d-di ;
	*mjd = date2julian(y, m, di) + df ;
}


#ifdef	STANDALONE

static
test2(jdate)
	double	jdate ;
{
	int	yday ;
	int	y2,m2,d2 ;
	int	y3,m3,d3 ;
	int	h,m ;
	double	s ;

	julian2time(jdate, &y3,&m3,&d3, &h,&m,&s) ;
	printf( "%.1f = %4d-%2.2d-%2.2d %d:%2.2d:%.1f\n",
	  jdate, y3,m3,d3, h,m,s) ;
}

static
test(y,m,d)
	int	y,m,d ;
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


static
testSt(y,m,d)
	int	y,m,d ;
{
	double	jd,T, st ;

	jd = date2julian(y,m,d) ;
	T = (jd-2415020)/36525 ;
	st = julian2sidereal(jd) ;
	printf("jd = %f, T=%f, st=", jd, T) ;
	printHms(st) ;
	putchar('\n') ;
}


main()
{
	double	jd,t,st ;
	int	yy,mm,dd, h,m ;
	double	s ;

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
