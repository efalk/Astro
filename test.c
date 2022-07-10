#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

static void checkDate(int y, int m, double d, double expected);
static const char * match(double a, double b, double epsilon);

int
main()
{
	double	date, time, st;
	double	lat,lon,rad ;
	double	decl, RA ;
	double	E ;
	double	psi, eps, obl;
	PlanetState p ;
	int	y,m,d ;
	int	H, M;
	double	S;
	double n;
	double elev, az;

	checkDate(1978,11,12, 2443824.5);
	checkDate(2000,1,1.5, JD2000);
	checkDate(1900,1,1, 2415020.5);
	checkDate(1600,1,1, 2305447.5);
	checkDate(1957,10,4.81, 2436116.31); /* Example [2]7.a */
	checkDate(333,1,27.5, 1842713.0); /* Example [2]7.a */
	checkDate(-4712,1,1.5, 0);
	putchar('\n');

	julian2date(JD2000, &y, &m, &d);
	printf("JD2000 = %d-%d-%d\n", y, m, d);
	printDate(JD2000);
	julian2date(2436116.81, &y, &m, &d);
	printf("2436116.81 = %d-%d-%d\n", y, m, d);	/* Example [2]7.c */
	printDate(2436116.81);
	putchar('\n');

	date = time2julian(1992, 10, 13, 0,0,0.);
	printf("date: %lf = ", date);
	printDate(date);
	putchar('\n');

	/* Sidereal time */
	date = date2julian(1987,4,10);	/* Example [2]12.a */
	st = julian2sidereal(date);
	h2hms(st, &H, &M, &S);
	printf("1987-4-10 = %lf = %lf ST = %d:%2.2d:%f\n",
	  date, st, H, M, S);
	date = date2julian(1987,4,10);	/* Example [2]12.b */
	time = hms2h(19,21,0);
	st = julianTime2sidereal(date, time);
	h2hms(st, &H, &M, &S);
	printf("1987-4-10 19:21:00 = %lf = %lf ST = %d:%2.2d:%f\n",
	  date, st, H, M, S);
	date = time2julian(1987,4,10, 19,21,0);
	st = time2sidereal(date);
	printf("1987-4-10 19:21:00 = %lf = %lf ST = %d:%2.2d:%f\n",
	  date, st, H, M, S);
	putchar('\n');

	/* Coordinate conversion */
	RA = hms2h(7,45,18.946);	/* Example [2]13.a */
	decl = hms2h(28,1,34.26);
	equat2ecliptic(decl, RA, &lat, &lon, JD2000);
	printf("RA=%lf, decl=%lf => lat = %lf (%s), lon = %lf (%s)\n", RA, decl,
	  lat, match(lat, 6.684170, .00001),
	  lon, match(lon, 113.215630, .00001));
	ecliptic2equat(lat, lon, &decl, &RA, JD2000);
	printf("lat=%lf, lon=%lf => RA = %lf (%s), decl = %lf (%s)\n", lat, lon,
	  RA, match(RA, 7.755263, .00001),
	  decl, match(decl, 28.026183, .00001));
	RA = hms2h(23,9,16.641);	/* Example [2]13.b */
	decl = -hms2h(6,43,11.61);
	date = time2julian(1987,4,10,19,21,0);
	st = time2sidereal(date);
	printf("Mean sidereal time: %s (%s)\n",
	  convertHms(st), match(st, hms2h(8,34,57.0896), .00000001));
	obl = obliquity(date);
	nutation(&psi, &eps, date);
	obl += eps/3600;
	printf("Nutation: %lf (%s), obliquity: %s (%s)\n",
	  psi, match(psi, -3.868, .00001),
	  convertHms(obl), match(obl, hms2h(8,34,56.853), .00000001));
	st += (eps/15) * cosd(obl) / 3600;
	printf("Apparent sidereal time: %s (%s)\n",
	  convertHms(st), match(st, hms2h(8,34,56.853), .000001));
	printf("Venus: RA=%lf, decl=%lf, ST=%lf\n",
	  RA, decl, st);
	putchar('\n');

	date = date2julian(1987, 4, 10);	/* Example [2] 22.a */
	nutation(&psi, &eps, date);
	obl = obliquity(date) + eps/3600;
	printf("Nutation on %lf = %lf (%s), %lf (%s)\n",
		date, psi, match(psi, -3.788, .0005),
		eps, match(eps, 9.443, .0005));
	printf("Obliquity on %lf = %s (%s)\n", date,
		convertHms(obl), match(obl, hms2h(23,26,36.850), .0001));
	putchar('\n');




	SunEquatorial(date, &decl, &RA, &rad);
	printf("decl=%lf, RA=%lf = %lf, rad=%lf\n", decl, RA, RA * 15, rad);
	lat = hm2h(37,23.1);
	lon = hm2h(122,4.9);

	SunPosition(jnow(), lat, lon, &az, &elev);
	printf("elev=%lf, az=%lf\n", elev, az);
	SunPosition(time2julian(2011,1,30,1,0,0.), lat, lon, &az, &elev);
	printf("@1700: elev=%lf, az=%lf\n", elev, az);
	SunPosition(time2julian(2011,1,30,1,25,0.), lat, lon, &az, &elev);
	printf("@1725: elev=%lf, az=%lf\n", elev, az);
	SunPosition(time2julian(2011,1,30,1,30,0.), lat, lon, &az, &elev);
	printf("@1730: elev=%lf, az=%lf\n", elev, az);
	SunPosition(time2julian(2011,1,30,1,35,0.), lat, lon, &az, &elev);
	printf("@1735: elev=%lf, az=%lf\n", elev, az);
	SunPosition(time2julian(2011,1,30,1,59,0.), lat, lon, &az, &elev);
	printf("@1759: elev=%lf, az=%lf\n", elev, az);

	printf("jnow = %lf\n", jnow());
	printf("Sun GHA = %lf\n", SunGHA(jnow()));
	n = jnow() - 2451545 - 0.0009 - lon/360.;
	printf("n = %lf\n", n);
	n = round(n);
	printf("n = %lf\n", n);
	date = 2451545 + 0.0009 + lon/360. + n;
	printf("noon = %lf, = ", date);
	printDate(date);
	date = SunNoon(jnow(), lat, lon);
	printf("noon = %lf, = ", date);
	printDate(date);
	printf("Sunset = %lf = ", SunSet(jnow(), lat, lon));
	printDate(SunSet(jnow(), lat, lon)-8./24.);

	julian2date(2415020.0, &y,&m,&d) ;
	printf("2415020.0 = %d - %d - %d\n", y,m,d) ;

	date = jnow() ;
	printf("Current julian date: %f\n", date) ;

	/* example 22a */
	E = keplerE(5.*RAD, 0.1) ;
	printf("kepler(5.,0.1) = %f\n", E*DEG) ;

	E = keplerE(2.*RAD, 0.99) ;
	printf("kepler(2.,0.99) = %f\n", E*DEG) ;

	/* example 18a */
	date = date2julian(1978,11,12) ;
	SunEcliptic(date, &lat,&lon,&rad) ;
	printf("sun @ %f = %f,%f,%f\n", date, lat,lon,rad) ;

	SunEquatorial(date, &decl,&RA,&rad) ;
	printf("sun @ %f = %f,%f,%f\n", date, decl,RA,rad) ;

	/* example 25a */
	date = date2julian(1978,11,12) ;
	Mercury(date, &p) ;
	printf("Mercury @ %f = %f,%f,%f\n", date, p.lat,p.lon,p.R) ;

	date = date2julian(1979,12,7) ;
	MoonPrecise(date, &p) ;
	printf("Moon @ %f = %f,%f,%f, par=%f\n",
		date, p.lat,p.lon,p.R, p.ad) ;

	Moon(date, &p) ;
	printf("Moon @ %f = %f,%f,%f, par=%f\n",
		date, p.lat,p.lon,p.R, p.ad) ;

	exit(0) ;
}

static const char *
match(double a, double b, double epsilon)
{
	static char rval[80];
	a -= b;
	if( a>= -epsilon && a < epsilon ) return "ok";
	snprintf(rval, sizeof(rval), "wrong, %lf should be %lf", a+b, b);
	return rval;
}

static void
checkDate(int y, int m, double d, double expected)
{
	double date = date2julian(y,m,d);
	printf("date(%d,%d,%lg) = %lf (%lf) %s\n",
		y,m,d, date, expected, date == expected ? "ok" : "wrong!");
}
