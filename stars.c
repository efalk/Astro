#ifndef lint
static const char sccsid[] = "%Z%%M% %I% %E% SMI" ;
static const char rcsid[] = "$Id$" ;
#endif

	/* Read the Yale & PPM star catalogs */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>

#include "astro.h"

#define	PI	3.141592653589793

#define	JD2000	2451545.0

#define	MAX24	((1<<24)-1)	/* max 24-bit # */


/* 
 *   PPM (Positions & Proper Motion) catalog:
 * 
 *     record size = 19 bytes
 *       0 & 0xc0	type: 0=star, 1=star-like, 2=double
 *       0-3	PPM #
 * 
 *       3 & 0x80	1=SAO #, 0=HD #
 *       3-5	SAO or HD #
 * 
 *       6-8	RA, rads, 0 .. 2*PI mapped to 0 .. 1<<24
 *       9-11	Dec, rads, -PI/2 .. PI/2 mapped to 0 .. (1<<24)-1
 * 
 *       12	Magnitude, -2 .. 14 mapped as (Mag+2)*10
 * 
 *       13	spectrum type
 *       14	subclass
 * 
 *       15-16	PM RA, coded as (PMA*10000)+5000
 *       17-18	PM Dec, coded as (PMD*1000)+10000
 * 
 *       RA, Decl seem to be as of 2000 (jd 2451545.0)
 *       proper motion seems to be expressed as seconds/year
 */


/* conversions:
 *	There are 2^24 PPM units per 180 degrees of declination.
 *	This comes to (seconds) 16777216:648000 = 262144:10125
 *
 * You cannot convert directly with multiplication and division and
 * remain in the range of 32-bit integers.  We have to compromise.
 * Instead of asking the compiler to compute deg=ppm*180*60*60/2^24,
 * we ask for
 *	dec = ppm*180*60*60/2^24
 *	dec = ppm*648000/2^24
 *	dec = ppm*10125/2^18
 *	dec = ppm/64*10125/2^12
 *
 *	ra = ppm*360*60*60/2^24
 *	ra = ppm*1296000/2^24
 *	ra = ppm*10125/2^17
 *	ra = ppm/64*10125/2^11
 */

#define	ppm2dec(x)	((((x)/64)*10125)/4096-324000)
#define	ppm2ra(x)	((((x)/64)*10125)/2048)
#define	dec2ppm(d)	(((((d)+324000)/45)*0x40000)/225)

int
ReadPPMStars(maxmag, ra0,ra1, d0,d1, jd, filename,  rval)
	int	maxmag ;
	long	ra0,ra1, d0,d1 ;
	double	jd ;
	char	*filename ;
	PPMStar	**rval ;
{
	unsigned char	buf[19] ;
	FILE	*ifile ;
	int	type;
	long	ra,dec, mag, pma, pmd ;
	int	wrap ;
	float	minmag = 100. ;
	int	count = 0 ;
	int	nalloc = 256 ;
	int	n ;
	int	years = (float)(jd-JD2000)/365.24 ;
	PPMStar	*ptr ;

	if( filename == NULL )
	  filename = "ppm.xe" ;

	if( (ifile = fopen(filename, "r")) == NULL ) {
	  perror(filename) ;
	  return 0 ;
	}

	*rval = NULL ;
	ptr = (PPMStar *)malloc(nalloc*sizeof(PPMStar)) ;

	if( ptr == NULL )
	  return 0 ;

	n = 0 ;

	/* if minimum declination is much more than -90, it is worthwhile
	 * to execute a search for the first declination.  There are more
	 * effecient methods than binary search, but binary search is simpler
	 * and safer.
	 */

	if( d0 > -80 *60*60 )
	{
	  long	rec0,rec1,rec ;
	  u_long id, target = dec2ppm(d0) ;

	  rec0 = 0 ;
	  fseek(ifile, 0L, 2) ;
	  rec1 = ftell(ifile) / sizeof(buf) - 1 ;
	  while( rec0 < rec1 )
	  {
	    rec = (rec0+rec1)/2 ;
	    if( fseek(ifile, rec*sizeof(buf), 0) == -1 )
	      rec0 = rec1 = -1 ;
	    else if( fread(buf, sizeof(buf), 1, ifile) != 1 )
	      rec0 = rec1 = -1 ;
	    else
	    {
	      id = buf[9]<<16 | buf[10]<<8 | buf[11] ;
	      if( id >= target )
		rec1 = rec - 1 ;
	      else
		rec0 = rec + 1 ;
	    }
	  }
	  if( rec1 >= 0 )
	    fseek(ifile, rec1*sizeof(buf), 0) ;
	  else
	    rewind(ifile) ;
	}

	wrap = ra1 < ra0 ;
	if( wrap )
	  ra1 += 360*60*60 ;

	while( fread(buf, sizeof(buf), 1, ifile) == 1 )
	{
	  type = buf[0]>>6 ;

	  ra = ppm2ra(buf[6]<<16 | buf[7]<<8 | buf[8]) ;
	  dec = ppm2dec(buf[9]<<16 | buf[10]<<8 | buf[11]) ;
	  mag = buf[12]*10 - 200 ;
	  pma = (buf[15]<<8 | buf[16]) - 5000 ;		/* pma*10000 */
	  pmd = (buf[17]<<8 | buf[18]) - 10000 ;	/* pmd*1000 */

	  /* integrate proper motion */
	  ra += pma*years*15/10000 ;
	  dec += pmd*years/1000 ;

	  if( wrap && ra < ra0 )
	    ra += 360*60*60 ;

	  if( mag <= maxmag  && 
	      dec >= d0 && dec <= d1 && ra >= ra0 && ra <= ra1 )
	  {
	    if( mag < minmag ) minmag = mag ;
	    if( count >= nalloc ) {
	      nalloc *= 2 ;
	      ptr = (PPMStar *)realloc(ptr, nalloc*sizeof(PPMStar)) ;
	      if( ptr == NULL )
		return 0 ;
	    }
	    ptr[count].ra = ra ;
	    ptr[count].dec = dec ;
	    ptr[count].mag = mag ;
	    ptr[count].type[0] = 'S' ;
	    ptr[count].type[1] = type == 2 ? 'D' : 'S' ;
	    ptr[count].spec[0] = buf[13] ;
	    ptr[count].spec[1] = buf[14] ;
	    ++count ;
	  }
	}

	fprintf(stderr, "%d stars of magnitude %.1f to %.1f\n",
		count, .01*minmag, .01*maxmag) ;

	fclose(ifile) ;

	*rval = ptr ;
	return count ;
}
