
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

/* I/O routines.
 *
 *
 * char *
 * julian2ymdStr(double jdate)
 *
 * char *
 * julian2hmsStr(double jdate)
 *
 * char *
 * julian2str(double jdate)
 *
 * char *
 * deg2dmsStr(double degrees)
 *	Convert degrees to dd°mm'ss.s
 *
 * char *
 * deg2dmStr(double degrees)
 *	Convert degrees to dd°mm.mm
 *
 * char *
 * hours2hmsStr(double hours)
 *	Convert hours to hh:mm:ss.s
 *
 * char *
 * hours2hmStr(double hours)
 *	Convert hours to hh:mm.mm
 *
 */


/* A word of explanation:  I defined a round-robin of four character
 * buffers.  Each function that returns a string value uses the next
 * buffer in turn.  Thus, you can call these functions four times
 * before any buffer becomes invalid.  This is very useful when you
 * use these functions inside a printf().
 */

#define	RVALN	4
static	char	rvals[RVALN][40] ;
static	int	rvalidx = 0 ;
#define	rval	(rvals[rvalidx])
#define	rvalNext()	(rvalidx = (rvalidx+1) % RVALN)

char *
julian2ymdStr(jdate)
	double	jdate ;
{
	int	y,m,d ;
static	char	*mnames[] = { "", "Jan", "Feb", "Mar", "Apr", "May",
		  "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"} ;

	rvalNext() ;
	julian2date(jdate, &y,&m,&d) ;
	sprintf(rval, "%d-%s-%d", d,mnames[m],y) ;
	return rval ;
}


char *
julian2hmsStr(jdate)
	double	jdate ;
{
	int	y,m,d, hh,mm ;
	double	ss ;
static	char	*mnames[] = { "", "Jan", "Feb", "Mar", "Apr", "May",
		  "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"} ;

	rvalNext() ;
	julian2time(jdate, &y,&m,&d, &hh,&mm,&ss) ;
	sprintf(rval, "%d:%2.2d:%2.2d.%d", hh,mm,(int)ss, (int)(ss*10+.5)%10) ;
	return rval ;
}


char *
julian2str(jdate)
	double	jdate ;
{
	int	y,m,d, hh,mm ;
	double	ss ;
static	char	*mnames[] = { "", "Jan", "Feb", "Mar", "Apr", "May",
		  "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"} ;

	rvalNext() ;
	julian2time(jdate, &y,&m,&d, &hh,&mm,&ss) ;
	sprintf(rval, "%d-%s-%d %d:%2.2d:%2.2d", d,mnames[m],y, hh,mm,(int)ss) ;
	return rval ;
}


char *
deg2dmsStr(degrees)
	double degrees ;
{
	int	d,m,s,ss ;

	rvalNext() ;
	d = degrees ;
	m = (int)(degrees * 60) % 60 ;
	s = (int)(degrees * 3600) % 60 ;
	ss = (int)(degrees * 36000 + .5) % 10 ;
	sprintf(rval, "%3.3d°%2.2d'%2.2d.%1.1d", d,m,s,ss) ;
	return rval ;
}


char *
deg2dmStr(degrees)
	double	degrees ;
{
	int	d,m,mm ;
	int	sign = 0 ;

	rvalNext() ;
	if( degrees < 0. ) {
	  sign = 1 ;
	  degrees = -degrees ;
	}
	d = degrees ;
	m = (int)(degrees * 60) % 60 ;
	mm = (int)(degrees * 6000 + .5) % 100 ;
	sprintf(rval, "%s%d°%2.2d.%2.2d", sign?"-":"", d,m,mm) ;
	return rval ;
}


char *
hours2hmsStr(hours)
	double hours ;
{
	int	h,m,s,ss ;

	rvalNext() ;
	h = hours ;
	m = (int)(hours * 60) % 60 ;
	s = (int)(hours * 3600) % 60 ;
	ss = (int)(hours * 36000 + .5) % 10 ;
	sprintf(rval, "%2d:%2.2d:%2.2d.%1.1d", h,m,s,ss) ;
	return rval ;
}


char *
hours2hmStr(hours)
	double	hours ;
{
	int	h,m,mm ;

	rvalNext() ;
	h = hours ;
	m = (int)(hours * 60) % 60 ;
	mm = (int)(hours * 6000 + .5) % 100 ;
	sprintf(rval, "%2d:%2.2d.%2.2d", h,m,mm) ;
	return rval ;
}
