static	char	usage[] =
"ephem - show positions of all planets\n"
"\n"
"  usage:  ephem [options] [yymmdd [hhmmss]]\n"
"	-a	include hour angles\n"
"	-l	local time\n"
;

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>

#include "astro.h"

static	void	showPlanet(), showSat(), showEarth() ;

static	PlanetState	earthState ;
static	double	jdate ;

static	int	hourAngles = 0 ;
static	double	aries ;		/* hour angle of aries */

main(argc, argv)
	int	argc ;
	char	**argv ;
{
	double	lat,lon,rad ;
	double	decl, RA ;
	double	E ;
	PlanetState p ;
	int	y,m,d ;
	int	hh, mm ;
	double	ss ;
	long	s ;
	int	haveDay = 0 ;
	int	local = 0 ;
	double	stime ;

	jdate = jnow() ;

	while( --argc > 0 ) {
	  ++argv ;
	  if( strcmp(*argv, "-a") == 0 )
	    hourAngles = 1 ;
	  else if( strcmp(*argv, "-l") == 0 )
	    local = 1 ;
	  else if( isdigit(**argv) )
	  {
	    if( !haveDay )
	    {
	      y = atoi(*argv) ;
	      d = y%100 ;
	      m = y/100 % 100 ;
	      y /= 10000 ;
	      if( y < 50 )
		y += 2000 ;
	      else if( y < 100 )
		y += 1900 ;
	      jdate = date2julian(y,m,d) ;
	      haveDay = 1 ;
	    }
	    else
	    {
	      s = atoi(*argv) ;
	      hh = s/10000 ;
	      mm = (s/100) % 100 ;
	      s %= 100 ;
	      jdate += hms2h(hh,mm,(double)s) / 24. ;
	    }
	  }
	  else {
	    fprintf(stderr, usage) ;
	    exit(2) ;
	  }
	}

	/* TODO: local to GMT conversion */

	printf("julian date: %f = %s\n", jdate, julian2str(jdate)) ;

	printf("%.4f = %s; sidereal = %s\n",
	  floor(jdate), julian2str(floor(jdate)),
	  hours2hmsStr(julian2sidereal(floor(jdate))) ) ;

	printf("%.4f = %s; sidereal = %s\n",
	  jdate, julian2str(jdate),
	  hours2hmsStr(julian2sidereal(jdate)) ) ;

	printf("%.4f = %s; sidereal = %s\n",
	  jdate+.5, julian2str(jdate+.5),
	  hours2hmsStr(julian2sidereal(jdate+.5)) ) ;

	printf("%.4f = %s; sidereal = %s\n",
	  ceil(jdate), julian2str(ceil(jdate)),
	  hours2hmsStr(julian2sidereal(ceil(jdate))) ) ;

	stime = julian2sidereal(jdate) ;
	printf("sidereal time (midnight) = %.4f = %s\n",
		stime, hours2hmsStr(stime) ) ;

	stime = time2sidereal(jdate) ;
	printf("sidereal time = %.4f = %s = %s\n\n",
		stime, hours2hmsStr(stime), deg2dmStr(stime*15) ) ;

	aries = 24. - stime ;
	printf("aries RA = %s", hours2hmsStr(aries)) ;
	printf(" = %s\n", hours2hmStr(aries)) ;
	printf("aries GHA = %s", deg2dmsStr(ra2sha(aries))) ;
	printf(" = %s\n", deg2dmStr(ra2sha(aries))) ;

	SunEquatorial(jdate, &decl,&RA,&rad) ;
	printf("sun GHA = %s", hours2hmsStr(RA)) ;
	printf(" = %s\n", hours2hmStr(RA)) ;

	printf("Object     lat         lon          r          decl        RA        dist") ;
	if( hourAngles )
	  printf("	SHA	GHA") ;
	putchar('\n') ;

	Earth(jdate, &earthState) ;

	SunEquatorial(jdate, &decl, &RA, &rad) ; showSat("Sun",decl,RA,rad) ;

	Moon(jdate, &p) ;
	ecliptic2equat(p.lat, p.lon, &decl, &RA, jdate) ;
	showSat("Moon", decl,RA,p.R) ;

	Mercury(jdate, &p) ; showPlanet("Mercury", &p) ;
	Venus(jdate, &p) ; showPlanet("Venus", &p) ;
	showEarth("Earth", &earthState) ;
	Mars(jdate, &p) ; showPlanet("Mars", &p) ;
	Jupiter(jdate, &p) ; showPlanet("Jupiter", &p) ;
	Saturn(jdate, &p) ; showPlanet("Saturn", &p) ;
	Uranus(jdate, &p) ; showPlanet("Uranus", &p) ;
	Neptune(jdate, &p) ; showPlanet("Neptune", &p) ;
#ifdef	COMMENT
	Pluto(jdate, &p) ; showPlanet("Pluto", &p) ;
#endif	/* COMMENT */



	printf("\n") ;
	printf("Notes: for planets, lat,lon,r are relative to the Sun, in\n") ;
	printf("eccliptic coordinates.  RA, decl are relative to the Earth,\n");
	printf("in celestial coordinates.  Distances in AU\n") ;

#ifdef	COMMENT
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
#endif	/* COMMENT */

	exit(0) ;
}


static	void
showPlanet(name, state)
	char	*name ;
	PlanetState *state ;
{
	int	lad,lam,las ;
	int	lod,lom,los ;
	int	rad, rah,ram,ras ;
	int	dd,dm,ds ;
	double	tmp ;
	char	sign = 'N' ;
	double	lat, lon, dist ;
	double	RA, decl ;
	double	sha, gha ;

	lat = state->lat ;
	if( lat < 0. ) {
	  sign = 'S' ;
	  lat = -lat ;
	}

	h2hms(lat, &lad, &lam, &tmp) ; las = tmp ;
	h2hms(state->lon, &lod, &lom, &tmp) ; los = tmp ;

	printf("%s	%2d°%2.2d'%2.2d%c   %3d°%2.2d'%2.2d   %8.4f",
		name, lad,lam,las,sign, lod,lom,los, state->R) ;

	deltaPolar(earthState.lat, earthState.lon, earthState.R,
		state->lat, state->lon, state->R,
		&lat, &lon, &dist) ;
	ecliptic2equat(lat, lon, &decl, &RA, jdate) ;
	if( decl < 0. ) {
	  sign = 'S' ;
	  decl = -decl ;
	}
	else
	  sign = 'N' ;
	if( RA < 0. ) RA += 24. ;

	h2hms(decl, &dd, &dm, &tmp) ; ds = tmp ;
	h2hms(RA, &rah, &ram, &tmp) ; ras = tmp ;

	printf("   %2d°%2.2d'%2.2d%c   %3d:%2.2d:%2.2d  %8.4f",
		dd,dm,ds,sign, rah,ram,ras, dist) ;

	if( hourAngles ) {
	  sha = ra2sha(RA) ;
	  gha = sha2gha(sha, jdate) ;
	  printf("	%s", deg2dmStr(sha) ) ;
	  printf(" %s", deg2dmStr(gha) ) ;
	}

	printf("\n") ;
}


static	void
showSat(name, decl,RA,dist)
	char	*name ;
	double	decl, RA, dist ;
{
	int	rad,rah,ram,ras ;
	int	dd,dm,ds ;
	double	tmp ;
	double	sha, gha ;
	char	sign = 'N' ;

	printf("%s					", name) ;

	if( decl < 0. ) {
	  sign = 'S' ;
	  decl = -decl ;
	}
	if( RA < 0. ) RA += 24. ;

	h2hms(decl, &dd, &dm, &tmp) ; ds = tmp ;
	h2hms(RA, &rah, &ram, &tmp) ; ras = tmp ;

	printf("   %2d°%2.2d'%2.2d%c   %3d:%2.2d:%2.2d  %8.4f",
		dd,dm,ds,sign, rah,ram,ras, dist) ;

	if( hourAngles ) {
	  sha = ra2sha(RA) ;
	  gha = sha2gha(sha, jdate) ;
	  printf("	%s", deg2dmStr(sha) ) ;
	  printf(" %s", deg2dmStr(gha) ) ;
	}

	printf("\n") ;
}


static	void
showEarth(name, state)
	char	*name ;
	PlanetState *state ;
{
	int	lad,lam,las ;
	int	lod,lom,los ;
	double	lat, tmp ;
	char	sign = 'N' ;

	lat = state->lat ;
	if( lat < 0. ) {
	  sign = 'S' ;
	  lat = -lat ;
	}

	h2hms(lat, &lad, &lam, &tmp) ; las = tmp ;
	h2hms(state->lon, &lod, &lom, &tmp) ; los = tmp ;

	printf("%s	%2d°%2.2d'%2.2d%c   %3d°%2.2d'%2.2d   %8.4f",
		name, lad,lam,las,sign, lod,lom,los, state->R) ;

	printf("\n") ;
}
