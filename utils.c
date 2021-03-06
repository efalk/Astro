
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

/* conversion and print utilities */


/**
 * convert hours to "hh:mm:ss.ssss"
 * Returns a buffer that will become invalid the next time this
 * function is called. Not thread-safe.
 */
const char *
convertHms(double hours)
{
static	char	obuf[80] ;

	int	h,m ;
	double	s ;

	h2hms(hours, &h,&m,&s) ;
	snprintf(obuf, sizeof(obuf), "%d:%2.2d:%f", h,m,s) ;
	return obuf ;
}


/**
 * convert hours (or degrees) to hours:minutes:seconds
 */
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
	printf("%s\n", convertHms(hours)) ;
}

/**
 * Convert spherical coordinates to rectangular
 */
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

/**
 * Convert rectangular coordinates to spherical
 */
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



/**
 * given polar coordinates of two objects, find the bearing
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


/**
 * Return a % 360
 */
double
limitAngle(double a)
{
  if( a >= 360 ) {
    int ia = a;
    a -= ia;
    return ia % 360 + a;
  } else if( a < 0 ) {
    return 360 - limitAngle(-a);
  }
  return a;
}


/**
 * Return h % 24
 */
double
limitHour(double h)
{
  if( h >= 24 ) {
    int ih = h;
    h -= ih;
    return ih % 24 + h;
  } else if( h < 0 ) {
    return 24 - limitHour(-h);
  }
  return h;
}

/**
 * Extract a float field from the given buffer
 * @param buffer - buffer to extract the field from
 * @param start  - index of first character, 1-based
 * @param end    - index of last character, inclusive
 */
double
recFloat(const char *buffer, int start, int end)
{
    double rval = 0;
    char tmp[80];
    recString(buffer, start, end, tmp);
    sscanf(tmp, "%lf", &rval);
    return rval;
}

/**
 * Extract a long field from the given buffer
 * @param buffer - buffer to extract the field from
 * @param start  - index of first character, 1-based
 * @param end    - index of last character, inclusive
 */
long
recLong(const char *buffer, int start, int end)
{
    long rval = 0;
    char tmp[80];
    recString(buffer, start, end, tmp);
    sscanf(tmp, "%ld", &rval);
    return rval;
}

/**
 * Extract a string from the given buffer
 * @param buffer - buffer to extract the field from
 * @param start  - index of first character, 1-based
 * @param end    - index of last character, inclusive
 * @param tmp    - caller-supplied buffer to receive string
 *
 * Caller is responsible for making sure tmp[] is at least
 * (end-start+1) characters.
 */
char *
recString(const char *buffer, int start, int end, char *tmp)
{
    memcpy(tmp, buffer+start-1, end-start+1);
    tmp[end-start+1] = '\0';
    return tmp;
}
