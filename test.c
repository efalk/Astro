#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"


main()
{
	double	date, time ;
	double	lat,lon,rad ;
	double	decl, RA ;
	double	E ;
	PlanetState p ;
	int	y,m,d ;

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
