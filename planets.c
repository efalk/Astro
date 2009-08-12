
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

/* Find the coordinates of the Planets, from Astronomical Formulae for
 * Calculators, by Jean Meeus, 4th edition, chapter 23.
 *
 * Meeus also gives formulae for computing the positions relative
 * to the standard ecliptics of 1950 and 2000, but I don't feel
 * like doing that right now.
 *
 * Planetary orbits are described by these terms:
 *
 *  L = mean longitude of planet
 *  a = semimajor axis of the orbit (a constant for each planet)
 *  e = eccentricity of orbit
 *  i = inclination on the plane of the ecliptic
 *  w = argument of perihelion
 *  omega = longitude of ascending node
 *
 *  pi = longitude of perihelion = w + omega
 *  M = mean anomaly = L - pi
 *  q = perihelion distance = a*(1-e)
 *  Q = aphelion distance = a*(1+e)
 *
 * Note: elements for longitude of perihelion taken from "ephem"
 * by Elwood Charles Downey.  These are probably just the sum of the
 * longitude of ascending node and argument of perihelion, but I don't
 * feel like double-checking right now.  Also, I'm not sure these will
 * work; I think he uses centuries since JD0 while Meeus uses
 * centuries since 1900

 */


	/* describe a planet's motion.  All angles in degrees */

typedef	struct	{
	  double L0,L1,L2,L3 ;		/* coeffs of longitude */
	  double p0,p1,p2,p3 ;		/* longitude of perihelion */
	  double w0,w1,w2,w3 ;		/* arg. of perihelion */
	  double e0,e1,e2,e3 ;		/* eccentricity */
	  double i0,i1,i2,i3 ;		/* inclination */
	  double o0,o1,o2,o3 ;		/* longitude of ascending node */
	  double M0,M1,M2 ;		/* mean anomoly */
	  double a0 ;			/* semimajor axis */
	  double ad ;			/* angular diameter @ 1AU, seconds */
	  double mag ;			/* magnitude @ 1Au */
	} Elements ;


static
range(x)
double *x ;
{
	while( *x < 0. ) *x += 360. ;
	if( *x >= 360. ) *x -= (int)*x/360*360 ;
}


	/* utility: return orbital elements, units in degrees */

static
getElements(el, date, p)
	Elements *el ;
	double	date ;
	PlanetState *p ;
{
	double	T,T2,T3 ;

	p->date = date ;

	T = (date - 2415020.0) / 36525. ; T2 = T*T ; T3 = T*T*T ;

	p->L = el->L0 + el->L1*T + el->L2*T2 + el->L3*T3 ;
	p->dL = 0. ;	/* TODO */
	p->a = el->a0 ;
	p->e = el->e0 + el->e1*T + el->e2*T2 + el->e3*T3 ;
	p->i = el->i0 + el->i1*T + el->i2*T2 + el->i3*T3 ;
	p->w = el->w0 + el->w1*T + el->w2*T2 + el->w3*T3 ;
	p->om = el->o0 + el->o1*T + el->o2*T2 + el->o3*T3 ;
	p->pi = el->p0 + el->p1*T + el->p2*T2 + el->p3*T3 ;
	p->ad = el->ad ;
	p->mag = el->mag ;
	p->M = el->M0 + el->M1*T + el->M2*T2 ;
	range( &p->L ) ;
	range( &p->i ) ;
	range( &p->w ) ;
	range( &p->om ) ;
	range( &p->pi ) ;
	range( &p->M ) ;
	p->year = 360.*36525./el->L1 ;
}


	/* utility: compute mean anomaly, degrees */

static	double
anomaly(date, m0,m1,m2)
	double	date,m0,m1,m2 ;
{
	double	T = (date - 2415020.0) / 36525. ;
	double	M ;

	M = m0 + m1*T + m2*T*T ;
	range( &M ) ;
	return M ;
}

	/* utility: given planetary elements, compute lat,lon,R */

static
planet(p)
	PlanetState *p ;
{
	double	E,v ;		/* eccentric,true anomaly, radians */
	double	u ;		/* argument of latitude */
	double	e = p->e ;	/* eccentricity */
	double	i = p->i*RAD ;	/* inclination, radians */

	E = keplerE(p->M*RAD, e) ;
	v = atan( sqrt((1.+e)/(1.-e))*tan(E/2.) ) * 2 ;

	p->R = p->a*(1. - e*cos(E)) ;

	u = (p->L - p->M - p->om)*RAD + v ;

	p->lon = datan2(cos(i)*sin(u), cos(u)) + p->om ;
	range( &p->lon ) ;

	p->lat = dasin( sin(u) * sin(i) ) ;

	p->v = v*DEG ;
	range( &p->v ) ;
}



static	Elements mercury = {
	  178.179078,	149474.07078,	0.0003011,	0.,	/* lon */
	  75.899697,	1.5554889,	2.947e-4,	0.0,	/* lon of p */
	  28.753753,	0.3702806,	+0.0001208,	0.,	/* arg of p */
	  0.20561421,	0.00002046,	-0.000000030,	0.,	/* e */
	  7.002881,	0.0018608,	-0.0000183,	0.,	/* i */
	  47.145944,	1.1852083,	0.0001739,	0.,	/* lon of asc */
	  102.27938,	149472.51529,	0.000007,		/* M */
	  0.3870986,						/* a */
	  6.74,							/* diam @ 1au */
	  -0.42							/* mag @ 1au */
	} ;

void
Mercury(date, p)
	double	date ;
	PlanetState *p ;
{
	getElements(&mercury, date, p) ;

	/* TODO: compute perturbations */

	p->M = anomaly(date, 102.27938, 149472.51529, 0.000007) ;
	planet(p) ;
}



static	Elements venus = {
	  342.767053,	58519.21191,	0.0003097,	0.,
	  130.163833,	1.4080361,	-9.764e-4,	0.0,
	  54.384186,	0.5081861,	-0.0013864,	0.,
	  0.00682069,	-0.00004774,	0.000000091,	0.,
	  3.393631,	0.001058,	-0.0000010,	0.,
	  75.779647,	0.8998500,	0.0004100,	0.,
	  212.60322,	5817.80387,	0.001286,
	  0.7233316,
	  16.92,
	  -4.4,
	} ;

void
Venus(date, p)
	double	date ;
	PlanetState *p ;
{
	getElements(&venus, date, p) ;

	/* TODO: compute perturbations */

	p->M = anomaly(date, 212.60322, 5817.80387, 0.001286) ;
	planet(p) ;
}


void
Earth(date, p)
	double	date ;
	PlanetState *p ;
{
	double	T,T2,T3 ;
	T = (date - 2415020.0) / 36525. ; T2 = T*T ; T3 = T2*T ;

	p->date = date ;
	p->L = 279.69668 + 36000.76892*T + .0003025*T2 ;
	p->e = .01675104 - .0000418*T - .000000126*T2 ;
	p->i = 0. ;
	p->om = 0. ;
	p->w = 0. ;	/* TODO */
	p->pi = 0. ;	/* TODO */
	p->a = 1. ;
	p->M = 358.47583 + 35999.04975*T - .000150*T2 - .0000033*T3 ;
	SunEcliptic(date, &p->lat, &p->lon, &p->R) ;
	p->lat = -p->lat ;
	p->lon += 180. ;
	range(&p->lon) ;
	p->year = 365.2424 ;
}


static	Elements mars = {
	  293.737334,	19141.69551,	0.0003107,	0.,
	  334.218203,	1.8407584,	1.299e-4,	-1.19e-6,
	  285.431761,	1.0697667,	0.0001313,	0.00000414,
	  0.09331290,	0.000092064,	-0.000000077,	0.,
	  1.850333,	-0.0006750,	0.0000126,	0.,
	  48.786442,	0.7709917,	-0.0000014,	-0.00000533,
	  319.51913,	19139.85475,	0.000181,
	  1.5236883,
	  9.36,
	  -1.52
	} ;

void
Mars(date, p)
	double	date ;
	PlanetState *p ;
{
	getElements(&mars, date, p) ;

	/* TODO: compute perturbations */

	p->M = anomaly(date, 319.51913, 19139.85475, 0.000181) ;
	planet(p) ;
}



static	Elements jupiter = {
	  238.049257,	3036.301986,	0.0003347,	-0.00000165,
	  12.720972,	1.6099617,	1.05627e-3,	-3.43e-6,
	  273.277558,	0.5994317,	0.00070405,	0.00000508,
	  0.04833475,	0.000164180,	-0.0000004676,	-0.0000000017,
	  1.308736,	-0.0056961,	0.0000039,	0.,
	  99.443414,	1.0105300,	0.00035222,	-0.00000851,
	  225.32833,	3034.69202,	0.000722,
	  5.202561,
	  196.74,
	  -9.4
	} ;

void
Jupiter(date, p)
	double	date ;
	PlanetState *p ;
{
	getElements(&jupiter, date, p) ;

	/* TODO: compute perturbations */

	p->M = anomaly(date, 225.32833, 3034.69202, 0.000722) ;
	planet(p) ;
}



static	Elements saturn = {
	  266.564337,	1223.509884,	0.0003245,	-0.0000058,
	  91.098214,	1.9584158,	8.2636e-4,	4.61e-6,
	  338.307800,	1.0852207,	0.00097854,	0.00000992,
	  0.05589232,	-0.00034550,	-0.000000728,	0.00000000074,
	  2.492519,	-0.0039189,	-0.00001549,	0.00000004,
	  112.790414,	0.8731951,	-0.00015218,	-0.00000531,
	  175.46622,	1221.55147,	0.000502,
	  9.554747,
	  165.6,
	  -8.88
	} ;

void
Saturn(date, p)
	double	date ;
	PlanetState *p ;
{
	getElements(&saturn, date, p) ;

	/* TODO: compute perturbations */

	p->M = anomaly(date, 175.46622, 1221.55147, 0.000502) ;
	planet(p) ;
}



static	Elements uranus = {
	  244.197470,	429.863546,	0.0003160,	-0.00000060,
	  171.548692,	1.4844328,	2.372e-4,	-6.1e-7,
	  98.071581,	0.9857650,	-0.0010745,	-0.00000061,
	  0.0463444,	-0.00002658,	0.000000077,	0.,
	  0.772464,	0.0006253,	0.0000395,	0.,
	  73.477111,	0.4986678,	0.0013117,	0.,
	  72.64878,	428.37911,	.000079,
	  19.21814,
	  65.8,
	  -7.19
	} ;

void
Uranus(date, p)
	double	date ;
	PlanetState *p ;
{
	double	T = (date - 2415020.0) / 36525. ;
	double	u,P,Q,S,W ;
	double	G,H ;
	double	M ;
	double	tau, mu, theta ;
	double	A,B ;
	double	ec ;
	double	la,lo,r ;

	getElements(&uranus, date, p) ;

	/* Meeus doesn't even bother approximating the
	 * mean anomaly, due to large perturbations, so
	 * we compute them all out.
	 */

	u = T/5. + 0.1 ;
	P = ( 237.47555 + 3034.9061*T ) * RAD ;
	Q = ( 265.91650 + 1222.1139*T ) * RAD ;
	S = ( 243.51721 +  428.4677*T ) * RAD ;
	W = 2.*P - 6.*Q + 3.*S ;

	G = ( 83.76922 + 218.4901*T ) * RAD ;
	H = 2.*G - S ;

	getElements(&uranus, date, p) ;

	p->M = anomaly(date, 72.64878, 428.37911, 0.000079) ;

	tau = S - P ;
	mu = S - Q ;
	theta = G - S ;

	A = + (0.864319 - 0.001583*u) * sin(H)
	    + (0.082222 - 0.006833*u) * cos(H)
	    + 0.036017 * sin(2*H)
	    - 0.003019 * cos(2*H)
	    + 0.008122 * sin(W) ;

	B = + 0.120303 * sin(H)
	    + (0.019472 - 0.000947*u) * cos(H)
	    + 0.006197 * sin(2*H) ;

	p->M += A - B/p->e ;

	ec = + (-3349. + 163.*u) * sin(H)
	     + 20981. * cos(H)
	     + 1311. * cos(2*H) ;

	p->e += ec * 1e-7 ;

	p->a -= 0.003825 * cos(H) ;

	planet(p) ;

	p->lon += + (0.012122 - 0.000988*u) * sin(S+mu)
		  + (-0.038581 + 0.002031*u - 0.001910*u*u) * cos(S+mu)
		  + (0.034964 - 0.001038*u + 0.000868*u*u) * cos(2*S+mu)
		  + 0.005594 * sin(S+3*theta)
		  - 0.014808 * sin(tau)
		  - 0.005794 * sin(mu)
		  + 0.002347 * cos(mu)
		  + 0.009872 * sin(theta)
		  + 0.008803 * sin(2*theta)
		  - 0.004308 * sin(3*theta) ;

	p->lat += + (  0.000458 * sin(mu)
		     - 0.000642 * cos(mu)
		     - 0.000517 * cos(4*theta) ) * sin(S)
		  - (  0.000347 * sin(mu)
		     + 0.000853 * cos(mu)
		     + 0.000517 * sin(4*mu) ) * cos(S)
		  + 0.000403 * ( cos(2*theta)*sin(2*S) + sin(2*theta)*cos(2*S));

	r =	-25948
		+ (5795*cos(S) - 1165*sin(S) + 1388*cos(2*S)) * sin(mu)
		+ 4985*cos(tau)
		+ (1351*cos(S) + 5702*sin(S) + 1388*sin(2*S)) * cos(mu)
		- 1230*cos(S)
		+ 904*cos(2*theta)
		+ 3354*cos(mu)
		+ 894*(cos(theta) - cos(3*theta)) ;

	p->R += r * 1e-6 ;
}


	/* TODO: */

static	Elements neptune = {
	  84.457994,	219.885914,	.0003205,	-0.00000060,
	  46.727364,	1.4245744,	3.9082e-4,	-6.05e-7,
	  276.045975,	.3256394,	.00014095,	.000004113,
	  .00899704,	.000006330,	.000000002,	0.,
	  1.779242,	-9.5436e-3,	-9.1e-6,	0.0,
	  130.681389,	1.098935,	2.4987e-4,	-4.718e-6,
	  37.73063,	218.46134,	.000070,
	  30.10957,
	  62.2,
	  -6.87
	} ;


void
Neptune(date, p)
	double	date ;
	PlanetState *p ;
{
	getElements(&neptune, date, p) ;

	/* TODO: compute perturbations */

	p->M = anomaly(date, 37.73063, 218.46134, .000070) ;
	planet(p) ;
}

	/* TODO: */
static	Elements pluto = {
	  95.3113544,	.3980332167,	0.0,	0.0,
	  224.017,	0.0,		0.0,	0.0,
	  0,	0,	0,	0,
	  .25515,	0.0,		0.0,	0.0,
	  17.1329,	0.0,		0.0,	0.0,
	  110.191,	0.0,		0.0,	0.0,
	  0.,		0.,		0.,
	  39.8151,
	  8.2,
	  1.0
	} ;

