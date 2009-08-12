
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

#include <stdio.h>

#define	PI	3.141592653589793

#define	JD2000	2451545.0

main()
{
	unsigned char buf[19] ;
	int	type, ppm, sao, hd ;
	double	ra,dec, mag, pma, pmd ;
	double	jd = JD2000 ;
	char	spect, class ;

	while( fread(buf, sizeof(buf), 1, stdin) == 1 )
	{

	  type = buf[0]>>6 ;
	  ppm = (buf[0]<<16 | buf[1]<<8 | buf[2]) & 0x3fffff ;

	  if( buf[3] & 0x80 ) {
	    sao = (buf[3]<<16 | buf[4]<<8 | buf[5]) & 0x7fffff ;
	    hd = 0 ;
	  } else {
	    sao = 0 ;
	    hd = buf[3]<<16 | buf[4]<<8 | buf[5] ;
	  }

	  ra = (buf[6]<<16 | buf[7]<<8 | buf[8]) * 24. / (1<<24) ;
	  dec = (buf[9]<<16 | buf[10]<<8 | buf[11]) * 180. / (1<<24) - 90. ;
	  mag = buf[12]/10. - 2. ;
	  spect = (char)buf[13] ;
	  class = (char)buf[14] ;
	  pma = ((buf[15]<<8 | buf[16]) - 5000.) / 10000. ;
	  pmd = ((buf[17]<<8 | buf[18]) - 10000.) / 1000. ;

	  /* integrate proper motion */
	  ra += pma/3600.*(jd-JD2000)/365.24 ;
	  dec += pmd/3600.*(jd-JD2000)/365.24 ;

	  printf("%c %6d %6d %6d %6.2f,%6.2f %4.1f %c%c %4.1f,%4.1f\n",
	    "STD"[type], ppm,sao,hd, ra,dec, mag, spect,class, pma,pmd) ;

	}

	exit(0) ;
}
