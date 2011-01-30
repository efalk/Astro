
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

/* conversion and print utilities */


	/* convert hours to hh:mm:ss.ssss */
	/* works for degrees too */

char	*
convertHms(double hours)
{
static	char	obuf[80] ;

	int	h,m ;
	double	s ;

	h2hms(hours, &h,&m,&s) ;
	snprintf(obuf, sizeof(obuf), "%d:%2.2d:%f", h,m,s) ;
	return obuf ;
}



	/* convert hours:minutes:seconds (or dms) to hours (or degrees) */
double
hms2h(h,m,s)
	int	h,m ;
	double	s ;
{
	double	rval = h ;
	rval += (1./60.)*m ;
	rval += (1./3600.)*s ;
	return rval ;
}

	/* convert hours (or degrees) to hours:minutes:seconds */

void
h2hms(hours, h,m,s)
	double	hours ;
	int	*h,*m ;
	double	*s ;
{
	*h = hours ;
	*m = hours*60. - *h*60 ;
	*s = hours*3600. - *h*3600 - *m*60 ;
}



void
printHms(double hours)
{
	printf("%s", convertHms(hours)) ;
}


void
polar2rect(lat, lon, R, X,Y,Z)
	double	lat,lon,R ;
	double	*X,*Y,*Z ;
{
	lat *= RAD ;
	lon *= RAD ;

	*X = R*cos(lat)*cos(lon) ;
	*Y = R*cos(lat)*sin(lon) ;
	*Z = R*sin(lat) ;
}


void
rect2polar(X,Y,Z, lat,lon,R)
	double	X,Y,Z ;
	double	*lat,*lon,*R ;
{
	*R = sqrt(X*X+Y*Y+Z*Z) ;
	if( *R > 0. ) {
	  *lat = asin(Z / *R) * DEG ;
	  *lon = atan2(Y,X) * DEG ;
	}
	else
	  *lat = *lon = 0. ;
}



	/* given polar coordinates of two objects, find the bearing
	 * and distance of object2 relative to object1
	 */

void
deltaPolar(lat1,lon1,r1, lat2,lon2,r2, lat3,lon3,r3)
	double	lat1,lon1,r1 ;
	double	lat2,lon2,r2 ;
	double	*lat3, *lon3, *r3 ;
{
	double	x1,y1,z1, x2,y2,z2, x3,y3,z3 ;

	polar2rect(lat1, lon1, r1, &x1,&y1,&z1) ;
	polar2rect(lat2, lon2, r2, &x2,&y2,&z2) ;
	x3 = x2 - x1 ;
	y3 = y2 - y1 ;
	z3 = z2 - z1 ;
	rect2polar(x3,y3,z3, lat3,lon3,r3) ;
}


double
limitAngle(a)
  double a ;
{
  while( a < 0. ) a += 360. ;
  if( a >= 360. ) a -= (int)a/360*360 ;
  return a ;
}


double
limitHour(a)
  double a ;
{
  while( a < 0. ) a += 24. ;
  if( a >= 24. ) a -= (int)a/24*24 ;
  return a ;
}
