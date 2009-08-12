
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

/* Celestial navigation utilities.  Based on the Nautical Almanac,
 * 1999 Commercial Edition.
 *
 * double
 * ra2sha(double ra)
 *	convert Right Ascension (hours) to Sidereal Hour Angle (degrees)
 *
 * double
 * sha2gha(double sha, double jdate)
 *	convert Sidereal Hour Angle to Grenwich Hour Angle
 *
 * double
 * gha2lha(double sha, double longitude)
 *	convert Grenwich Hour Angle to Local Hour Angle
 *
 * double
 * interpolate(double a, double b, double hours)
 *	interpolate between a,b based on the fraction part of hours
 *
 * altaz(double lha, double decl, double lat, double *alt, double *az)
 *	Compute altitude (Hc) and Azimuth (Z) from lha, declination, latitude
 *
 * double
 * sext2obs(Hs, ie, h, T, P, HP)
 *	Compute sextant altitude.
 */


  /* convert right ascension (hours) to sidereal hour angle (degrees) */

double
ra2sha(ra)
  double ra ;
{
  return limitAngle( -15. * ra ) ;
}


  /* convert sidereal hour angle and time to grenwhich hour angle */

double
sha2gha(sha, jdate)
  double sha, jdate ;
{
  return limitAngle( sha + 15.*time2sidereal(jdate) ) ;
}


  /* convert grenwhich hour angle and longitude to local hour angle */

double
gha2lha(sha, lon)
  double sha, lon ;
{
  return sha + lon ;
}


  /* convert two hourly values and time in hours to intermediate value */

double
interpolate(a,b, time)
  double a,b, time ;
{
  time -= (int) time ;
  return a + (b-a)*time ;
}


/* given a local hour angle, declination, and observer's lattitude,
 * compute altitude and azimuth.  Do this with integers, and you
 * can re-create the sight reduction tables.
 */

void
altaz(lha, decl, lat, alt, az)
  double lha, decl, lat ;
  double *alt, *az ;
{
  double s,c,Hc ;
  double cHc, x, a, z ;

  lha *= RAD ;
  decl *= RAD ;
  lat *= RAD ;
  s = sin(decl) ;
  c = cos(decl) * cos(lha) ;
  Hc = asin( s*sin(lat) + c*cos(lat) ) ;
  *alt = Hc * DEG ;

  cHc = cos(Hc) ;
  if( cHc == 0. )
    a = 0. ;
  else {
    x = (s * cos(lat) - c * sin(lat) ) / cHc ;
    if( x > 1. ) x = 1. ;
    else if( x < -1. ) x = -1. ;
    a = acos(x) * DEG ;
  }
  if( lha > M_PI )
    z = a ;
  else
    z = 360. - a ;

  *az = z ;
}

#define	sind(a)	(sin((a)*RAD))
#define	cosd(a)	(cos((a)*RAD))
#define	tand(a)	(tan((a)*RAD))


  /* Given a raw sextant reading, index error, observer's height
   * above sea level, temperature, pressure and horizontal parallax,
   * compute altitude.
   */

double
sext2obs(Hs, ie, h, T, P, HP)
  double Hs ;	/* raw sextant reading */
  double ie ;	/* index error */
  double h ;	/* height above sea level, meters */
  double T ;	/* temperature, °C, if known, else 0. */
  double P ;	/* pressure, mb, if known, else 0. */
  double HP ;	/* horizontal parallax; see almanac.  0 for stars. */
{
  double D ;	/* compute dip from height */
  double H ;	/* compute apparent altitude */
  double Ro ;	/* refraction */
  double PA ;	/* parallax in altitude */

  D = 0.0293 * sqrt(h) ;	/* compute dip from height */
  H = (Hs + ie - D) ;		/* compute apparent altitude */
  Ro = 0.0167 / tand( H + 7.31/(H + 4.4) ) ;

  if( P > 0. )
    Ro *= 0.28 * P / (T + 273) ;

  PA = HP * cosd(H) ;

  return H - Ro + PA ;
}



/* Quick overview of spherical trig:
 *
 * Given a right triangle with sides a,b,c and angles A,B,C
 * (opposite of sides) and with angle C == 90 degrees, then Napier's
 * rules apply:  (Note that the lengths of sides are measured as
 * angles along the sphere, not linear lengths)
 *
 * sin(a) = tan(b) * cot(B) = sin(c) * sin(A)
 * sin(b) = tan(a) * cot(A) = sin(c) * sin(B)
 * cos(c) = cot(A) * cot(B) = cos(a) * cos(b)
 * cos(A) = tan(b) * cot(c) = cos(a) * sin(B)
 * cos(B) = tan(a) * cot(c) = cos(b) * sin(A)
 *
 * A quadrantal spherical triangle has one 90° side
 *
 * A biquadrantal spherical triangle has two 90° sides, is
 * isosceles and has two right angles opposite the 90° sides.
 *
 * A triquadrantal spherical triangle has three 90° sides,
 * is equilateral, has three right angles, and bounds
 * 1/8 of the sphere.
 */





#ifdef	STANDALONE
main()
{
  double jdate, time, ra, sha, gha, lha ;
  double decl, dist ;
  double lat, lon ;
  double alt, az ;
  double Hs, Ho ;

  jdate = time2julian(1999,11,16, 20,13,25.) ;
  time = hms2h(20,13,25.) ;
  ra = 0. ;
  sha = ra2sha(ra) ;
  gha = sha2gha(sha, jdate) ;
  printf("date = %.4f, ra=%.4f, sha = %s, gha = %s\n",
    jdate, ra, deg2dmStr(sha), deg2dmStr(gha)) ;

  SunEquatorial(jdate, &decl, &ra, &dist) ;
  sha = ra2sha(ra) ;
  gha = sha2gha(sha, jdate) ;
  printf("date = %.4f, ra=%.4f, dec = %.4f, sha = %.4f, gha = %.4f\n",
    jdate, ra, decl, sha, gha) ;

  gha = interpolate(123.8050, 138.8033, time) ;
  decl = interpolate(-18.7750, -18.7850, time) ;
  printf("interpolated: gha = %.4f, decl = %.4f\n", gha, decl) ;

  gha = interpolate(355.4450, 370.4867, time) ;
  gha += 80. + 46.4/60. ;
  gha = limitAngle(gha) ;
  printf("interpolated: gha = %.4f\n", gha) ;

  gha = 53. ;
  decl = -15. ;
  lat = 32. ;
  lon = -16. ;
  lha = gha2lha(gha,lon) ;
  altaz(lha, decl, lat, &alt, &az) ;
  printf("lha = %.4f, alt=%.4f, az=%.4f\n", lha, alt, az) ;

  /* exampe, p.281 */
  Ho = sext2obs(21.3283, 0., 5.4, -3., 982., .0024) + .27 ;
  printf("Ho = %.4f\n", Ho) ;

  Ho = sext2obs(3.3367, 0., 5.4, -3., 982., .0024) - .27 ;
  printf("Ho = %.4f\n", Ho) ;

  Ho = sext2obs(33.4600, 0., 5.4, -3., 982., .9317) + .2538 ;
  printf("Ho = %.4f\n", Ho) ;

  Ho = sext2obs(26.1117, 0., 5.4, -3., 982., .9317) - .2538 ;
  printf("Ho = %.4f\n", Ho) ;

  Ho = sext2obs(4.5433, 0., 5.4, -3., 982., .0033) ;
  printf("Ho = %.4f\n", Ho) ;

  Ho = sext2obs(49.6083, 0., 5.4, -3., 982., 0.) ;
  printf("Ho = %.4f\n", Ho) ;


  exit(0) ;
}
#endif	/* STANDALONE */
