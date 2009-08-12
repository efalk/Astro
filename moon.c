
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

/* Find the coordinates of the Moon, from Astronomical Formulae for
 * Calculators, by Jean Meeus, 4th edition, chapter 30.
 */

#define	degrees	* RAD

static
range(x)
double *x ;
{
	while( *x < 0. ) *x += 360. ;
	if( *x >= 360. ) *x -= (int)*x/360*360 ;
}


	/* Return info about the moon.  Not all fields in planet state
	 * are filled in
	 */

void
MoonPrecise(date, m)
	double	date ;
	PlanetState *m ;
{
	double	T,T2,T3 ;
	double	Lm ;		/* moon's mean longitude */
	double	Mm ;		/* moon's mean anomoly */
	double	D ;		/* moon's mean elongation */
	double	F ;		/* mean distance of moon from ascending node */
	double	Ohm ;		/* longitude of moon's ascending node */
	double	M ;		/* sun's mean longitude */
	double	e,e2 ;		/* Meeus doesn't explain this */
	double	s ;
	double	w1,w2 ;
	double	par ;		/* parralax */

	m->date = date ;

	T = (date - 2415020.0) / 36525. ; T2 = T*T ; T3 = T*T*T ;

	Lm = 270.434164 + 481267.8831*T - .001133*T2 + .0000019*T3 ;
	Mm = 296.104608 + 477198.8491*T + .009192*T2 + .0000144*T3 ;
	D  = 350.737486 + 445267.1142*T - .001436*T2 + .0000019*T3 ;
	F  =  11.250889 + 483202.0251*T - .003211*T2 - .0000003*T3 ;
	Ohm= 259.183275 -   1934.1420*T + .002078*T2 + .0000022*T3 ;
	M  = 358.475833 +  35999.0498*T - .000150*T2 - .0000033*T3 ;

	range(&Ohm) ;

	/* additive terms */

	s = dsin(51.2 + 20.2*T) ;
	Lm += .000233 * s ;
	M  -= .001778 * s ;
	Mm += .000817 * s ;
	D  += .002011 * s ;

	/* "Great Venus Term" */
	s = dsin(346.560 + 132.870*T - .0091731*T2) ;
	Lm += .003964 * s ;
	Mm += .003964 * s ;
	D  += .003964 * s ;
	F  += .003964 * s ;

	s = dsin(Ohm) ;
	Lm += .001964 * s ;
	Mm += .002541 * s ;
	D  += .001964 * s ;
	F  -= .024691 * s ;

	F -= .004328 * dsin(Ohm + 275.05 - 2.3*T) ;

	m->M = Mm ;

	e = 1 - .002495*T - .00000752*T2 ; e2 = e*e ;


	range(&Lm) ;
	range(&Mm) ;
	range(&D) ;
	range(&F) ;
	range(&M) ;

	/* convert M,Mm,D,F to radians to make things quicker */
	M *= RAD ;
	Mm *= RAD ;
	D *= RAD ;
	F *= RAD ;

	m->lon = Lm
		+ 6.288750 * sin ( Mm )
		+ 1.274018 * sin ( 2*D - Mm )
		+ .658309 * sin ( 2*D )
		+ .213616 * sin ( 2*Mm )
		- .185596 * sin ( M ) * e
		- .114336 * sin ( 2*F )
		+ .058793 * sin ( 2*D - 2*Mm )
		+ .057212 * sin ( 2*D - M - Mm ) * e
		+ .053320 * sin ( 2*D + Mm )
		+ .045874 * sin ( 2*D - M ) * e
		+ .041024 * sin ( Mm - M ) * e
		- .034718 * sin ( D )
		- .030465 * sin ( M + Mm ) * e
		+ .015326 * sin ( 2*D - 2*F )
		- .012528 * sin ( 2*F + Mm )
		- .010980 * sin ( 2*F - Mm )
		+ .010674 * sin ( 4*D - Mm )
		+ .010034 * sin ( 3*Mm )
		+ .008548 * sin ( 4*D - 2*Mm )
		- .007910 * sin ( M - Mm + 2*D ) * e
		- .006783 * sin ( 2*D + M ) * e
		+ .005162 * sin ( Mm - D )
		+ .005000 * sin ( M + D ) * e
		+ .004049 * sin ( Mm - M + 2*D ) * e
		+ .003996 * sin ( 2*Mm + 2*D )
		+ .003862 * sin ( 4*D )
		+ .003665 * sin ( 2*D - 3*Mm )
		+ .002695 * sin ( 2*Mm - M ) * e
		+ .002602 * sin ( Mm - 2*F - 2*D )
		+ .002396 * sin ( 2*D - M - 2*Mm ) * e
		- .002349 * sin ( Mm + D )
		+ .002249 * sin ( 2*D - 2*M ) * e2
		- .002125 * sin ( 2*Mm + M ) * e
		- .002079 * sin ( 2*M ) * e2
		+ .002059 * sin ( 2*D - Mm - 2*M ) * e2
		- .001773 * sin ( Mm + 2*D - 2*F )
		- .001595 * sin ( 2*F + 2*D )
		+ .001220 * sin ( 4*D - M - Mm ) * e
		- .001110 * sin ( 2*Mm + 2*F )
		+ .000892 * sin ( Mm - 3*D )
		- .000811 * sin ( M + Mm + 2*D ) * e
		+ .000761 * sin ( 4*D - M - 2*Mm ) * e
		+ .000717 * sin ( Mm - 2*M )* e2
		+ .000704 * sin ( Mm - 2*M - 2*D ) * e2
		+ .000693 * sin ( M - 2*Mm + 2*D ) * e
		+ .000598 * sin ( 2*D - M - 2*F ) * e
		+ .000550 * sin ( Mm + 4*D )
		+ .000538 * sin ( 4*Mm )
		+ .000521 * sin ( 4*D - M ) * e
		+ .000486 * sin ( 2*Mm - D ) ;

	m->lat = 5.128189 * sin ( F )
		+ .280606 * sin ( Mm + F )
		+ .277693 * sin ( Mm - F )
		+ .173238 * sin ( 2*D - F )
		+ .055413 * sin ( 2*D + F - Mm )
		+ .046272 * sin ( 2*D - F - Mm )
		+ .032573 * sin ( 2*D + F )
		+ .017198 * sin ( 2*Mm + F )
		+ .009267 * sin ( 2*D + Mm - F )
		+ .008823 * sin ( 2*Mm - F )
		+ .008247 * sin ( 2*D - M - F )  * e
		+ .004323 * sin ( 2*D - F - 2*Mm )
		+ .004200 * sin ( 2*D + F + Mm )
		+ .003372 * sin ( F - M - 2*D )  * e
		+ .002472 * sin ( 2*D + F - M - Mm )  * e
		+ .002222 * sin ( 2*D + F - M )  * e
		+ .002072 * sin ( 2*D - F - M - Mm )  * e
		+ .001877 * sin ( F - M + Mm )  * e
		+ .001828 * sin ( 4*D - F - Mm )
		- .001803 * sin ( F + M )  * e
		- .001750 * sin ( 3*F )
		+ .001570 * sin ( Mm - M - F )  * e
		- .001487 * sin ( F + D )
		- .001481 * sin ( F + M + Mm )  * e
		+ .001417 * sin ( F - M - Mm )  * e
		+ .001350 * sin ( F - M )  * e
		+ .001330 * sin ( F - D )
		+ .001106 * sin ( F + 3*M )
		+ .001020 * sin ( 4*D - F )
		+ .000833 * sin ( F + 4*D - Mm )
		+ .000781 * sin ( Mm - 3*F )
		+ .000670 * sin ( F + 4*D - 2*Mm )
		+ .000606 * sin ( 2*D - 3*F )
		+ .000597 * sin ( 2*D + 2*Mm - F )
		+ .000492 * sin ( 2*D + Mm - M - F )  * e
		+ .000450 * sin ( 2*Mm - F - 2*D )
		+ .000439 * sin ( 3*Mm - F )
		+ .000423 * sin ( F + 2*D + 2*Mm )
		+ .000422 * sin ( 2*D - F - 3*Mm )
		- .000367 * sin ( M + F + 2*D - Mm )  * e
		- .000353 * sin ( M + F + 2*D )  * e
		+ .000331 * sin ( F + 4*D )
		+ .000317 * sin ( 2*D + F - M + Mm )  * e
		+ .000306 * sin ( 2*D - 2*M - F )  * e2
		- .000283 * sin ( Mm + 3*F ) ;

	w1 = .0004664 * dcos(Ohm) ;
	w2 = .0000754 * dcos(Ohm + 275.05 - 2.3*T) ;
	m->lat *= (1. - w1 - w2) ;

	par = 0.950724
		+ .051818 * cos ( Mm )
		+ .009531 * cos ( 2*D - Mm )
		+ .007843 * cos ( 2*D )
		+ .002824 * cos ( 2*Mm )
		+ .000857 * cos ( 2*D + Mm )
		+ .000533 * cos ( 2*D - M )  * e
		+ .000401 * cos ( 2*D - M - Mm )  * e
		+ .000320 * cos ( Mm - M )  * e
		- .000271 * cos ( D )
		- .000264 * cos ( M + Mm )  * e
		- .000198 * cos ( 2*F - Mm )
		+ .000173 * cos ( 3*Mm )
		+ .000167 * cos ( 4*D - Mm )
		- .000111 * cos ( M )  * e
		+ .000103 * cos ( 4*D - 2*Mm )
		- .000084 * cos ( 2*Mm - 2*D )
		- .000083 * cos ( 2*D + M )  * e
		+ .000079 * cos ( 2*D + Mm )
		+ .000072 * cos ( 4*D )
		+ .000064 * cos ( 2*D - M + Mm )
		- .000063 * cos ( 2*D + M - Mm )
		+ .000041 * cos ( M + D )
		+ .000035 * cos ( 2*Mm - M )
		- .000033 * cos ( 3*Mm - 2*D )
		- .000030 * cos ( Mm + D )
		- .000029 * cos ( 2*F - 2*D )
		- .000029 * cos ( 2*Mm + M )
		+ .000026 * cos ( 2*D - 2*M )
		- .000023 * cos ( 2*F - 2*D + Mm )
		+ .000019 * cos ( 4*D - M - Mm ) ;


	m->L = Lm ;
	m->dL = 0. ;	/* TODO */
	m->a = 0. ;	/* TODO */
	m->e = 0. ;	/* TODO */
	m->i = 0. ;	/* TODO */
	m->w = 0. ;	/* TODO */
	m->om = Ohm ;
	m->pi = 0. ;	/* TODO */
	m->ad = par ;	/* TODO */
	m->mag = 0. ;	/* TODO */
	m->v = 0. ;	/* TODO */
	m->R = 6378.14 / dsin(par) ;	/* kilometers, TODO: AU */
	range( &m->L ) ;
	range( &m->i ) ;
	range( &m->w ) ;
	range( &m->om ) ;
	range( &m->pi ) ;
	range( &m->M ) ;
	range( &m->lat ) ; if( m->lat > 180. ) m->lat -= 360. ;
	range( &m->lon ) ;
	m->year = 0. ;	/* TODO */
}




void
Moon(date, m)
	double	date ;
	PlanetState *m ;
{
	double	T,T2,T3 ;
	double	Lm ;		/* moon's mean longitude */
	double	Mm ;		/* moon's mean anomoly */
	double	D ;		/* moon's mean elongation */
	double	F ;		/* mean distance of moon from ascending node */
	double	M ;		/* sun's mean longitude */
	double	e,e2 ;		/* Meeus doesn't explain this */
	double	s ;
	double	w1,w2 ;
	double	par ;		/* parralax */

	m->date = date ;

	T = (date - 2415020.0) / 36525. ; T2 = T*T ; T3 = T*T*T ;

	Lm = 270.434164 + 481267.8831*T - .001133*T2 + .0000019*T3 ;
	Mm = 296.104608 + 477198.8491*T + .009192*T2 + .0000144*T3 ;
	D  = 350.737486 + 445267.1142*T - .001436*T2 + .0000019*T3 ;
	F  =  11.250889 + 483202.0251*T - .003211*T2 - .0000003*T3 ;
	M  = 358.475833 +  35999.0498*T - .000150*T2 - .0000033*T3 ;

	m->M = Mm ;

	e = 1 - .002495*T - .00000752*T2 ; e2 = e*e ;


	range(&Lm) ;
	range(&Mm) ;
	range(&D) ;
	range(&F) ;
	range(&M) ;

	/* convert M,Mm,D,F to radians to make things quicker */
	M *= RAD ;
	Mm *= RAD ;
	D *= RAD ;
	F *= RAD ;

	m->lon = Lm
		+ 6.288750 * sin ( Mm )
		+ 1.274018 * sin ( 2*D - Mm )
		+ .658309 * sin ( 2*D ) ;

	m->lat = 5.128189 * sin ( F )
		+ .280606 * sin ( Mm + F )
		+ .277693 * sin ( Mm - F ) ;

	w1 = .0004664 ;
	w2 = .0000754 ;
	m->lat *= (1. - w1 - w2) ;

	par = 0.950724
		+ .051818 * cos ( Mm )
		+ .009531 * cos ( 2*D - Mm )
		+ .007843 * cos ( 2*D ) ;


	m->L = Lm ;
	m->dL = 0. ;	/* TODO */
	m->a = 0. ;	/* TODO */
	m->e = 0. ;	/* TODO */
	m->i = 0. ;	/* TODO */
	m->w = 0. ;	/* TODO */
	m->om = 0. ;	/* TODO */
	m->pi = 0. ;	/* TODO */
	m->ad = par ;	/* TODO */
	m->mag = 0. ;	/* TODO */
	m->v = 0. ;	/* TODO */
	m->R = 6378.14 / dsin(par) ;	/* kilometers, TODO: AU */
	range( &m->L ) ;
	range( &m->i ) ;
	range( &m->w ) ;
	range( &m->om ) ;
	range( &m->pi ) ;
	range( &m->M ) ;
	range( &m->lat ) ; if( m->lat > 180. ) m->lat -= 360. ;
	range( &m->lon ) ;
	m->year = 0. ;	/* TODO */
}
