
	/* Read the Yale star catalog */

#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "astro.h"

static void rstrip(char *buf);


/**
 * PPM (Positions & Proper Motion) catalog:
 * See bsc5.readme for format.
 *
 * There are two files in the database, "bsc5.dat" which contains
 * the basic data, and "bsc5.notes" which contains metadata such as
 * the names of the stars. Both are sorted by the Yale catalog #
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

/**
 * Read the Yale Star Catalog (aka Bright Star Catalog) or a
 * portion thereof.
 *
 * @param maxmag - maximum magnitude
 * @param ra0,ra1 - right ascension bounds
 * @param d0,d1   - declination bounds
 * @param datfile - data filename, NULL defaults to "bsc5.dat"
 * @param notefile - notes filename, NULL defaults to "bsc5.notes"
 * @param rval     - returned array of YaleStar structs
 * @param recsize  - returned size of one struct
 *
 * @return number of structs returned in rval
 *
 * Only returns items inside the box defined by ra0,ra1, d0,d1 and
 * with magnitude maxmag or less.
 */
int
ReadYaleStars(float maxmag, double ra0, double ra1,
	double d0, double d1,
	const char *datfilename,
	const char *notefilename,
	YaleStar **rval, size_t *recsize)
{
	char	datbuf[300];
	char	notebuf[200];
	char	tmp[140];
	FILE	*dfile;		/* data file */
	FILE	*nfile;		/* notes file */
	long	idx, nidx, sao;
	double	ra,dec, mag;
	int	count = 0;
	int	nalloc = 256;
	long	d,h,m;
	double	s;
	int	wrap;
	YaleStar *ptr;

	if( datfilename == NULL )
	  datfilename = "bsc5.dat";
	if( notefilename == NULL )
	  notefilename = "bsc5.notes";

	if( (dfile = fopen(datfilename, "r")) == NULL ) {
	  perror(datfilename);
	  return 0;
	}

	if( (nfile = fopen(notefilename, "r")) == NULL ) {
	  perror(notefilename);
	  fclose(dfile);
	  return 0;
	}

	*rval = NULL;
	ptr = (YaleStar *)malloc(nalloc*sizeof(*ptr));

	if( ptr == NULL )
	  return 0;

	wrap = ra1 < ra0;
	if( wrap )
	    ra1 += 360;

	/* Sample line:
	 * 1       10        20        30        40        50        60        70        80        90       100       110       120       130       140       150       160       170
	 * |        |         |         |         |         |         |         |         |         |         |         |         |         |         |         |         |         |
	 *  372          BD+44  271   7647 37077    I                  011115.6+442231011705.1+445407127.70-17.73 6.34R                     K5                 +0.014-0.045      -052
	 */

	nidx = -1;

	while( fgets(datbuf, sizeof(datbuf), dfile) != NULL )
	{
	    idx = recLong(datbuf, 1,4);
	    sao = recLong(datbuf, 32,37);
	    h = recLong(datbuf, 76,77);
	    m = recLong(datbuf, 78,79);
	    s = recFloat(datbuf, 80,83);
	    ra = h*24*60 + m*60 + s;
	    ra = ra/3600.*(360/24);
	    d = recLong(datbuf, 84,86);
	    m = recLong(datbuf, 87,88);
	    s = recLong(datbuf, 89,90);
	    dec = d >= 0 ? d*3600. + m*60. + s : d*3600. - m*60. - s;;
	    dec /= 3600.;
	    mag = recFloat(datbuf, 103,107);

	    /* We have enough information to apply filter now */

	    if( wrap && ra < ra0 )
	      ra += 360*60*60;

	    if (dec < d0 || dec > d1 || ra < ra0 || ra > ra1 || mag > maxmag)
		continue;

	    /* Reserve space in the returned array */
	    if( count >= nalloc ) {
		nalloc *= 2;
		ptr = (YaleStar *)realloc(ptr, nalloc*sizeof(*ptr));
		if( ptr == NULL )
		  return 0;
	    }

	    ptr[count].s.ra = ra;
	    ptr[count].s.dec = dec;
	    ptr[count].s.mag = mag;
	    ptr[count].s.epoch = 2000;
	    ptr[count].s.pmr = recFloat(datbuf, 149,154);
	    ptr[count].s.pmd = recFloat(datbuf, 155,160);
	    strcpy(ptr[count].s.type, "SS");
	    recString(datbuf, 148,148, ptr[count].s.spec);
	    ptr[count].s.sao = recLong(datbuf, 32,37);
	    ptr[count].s.name = NULL;
	    ptr[count].cons[0] = '\0';
	    ptr[count].yale_cat = idx;

	    /* Fetch data from notes file, if any */
	    while (nidx < idx) {
		if (fgets(notebuf, sizeof(notebuf), nfile) == NULL) {
		    nidx = 99999999;
		    break;
		}
		nidx = recLong(notebuf, 2,5);
	    }
	    while (nidx == idx) {	/* Found some */
		recString(notebuf, 8,11, tmp);
		if (tmp[0] == 'N') {
		    recString(notebuf, 13,132, tmp);
		    rstrip(tmp);
		    ptr[count].s.name = strdup(tmp);
		}
		if (fgets(notebuf, sizeof(notebuf), nfile) == NULL) {
		    nidx = 99999999;
		    break;
		}
		nidx = recLong(notebuf, 2,5);
	    }

	    ++count;
	  }

	  fprintf(stderr, "%d stars up to magnitude %.1f\n", count, maxmag);

	  fclose(dfile);
	  fclose(nfile);

	  *rval = ptr;
	  *recsize = sizeof(*ptr);
	  return count;
}

static void
rstrip(char *buf)
{
    size_t len = strlen(buf);
    char *ptr = buf+len-1;
    while (ptr >= buf && isspace(*ptr))
	*ptr-- = '\0';
}

#ifdef STANDALONE

int
main()
{
    YaleStar *rval, *ptr;
    size_t recsize;
    int count;
    int i;

    count = ReadYaleStars(5.0, 0., 360., -90, 90.,
    	"../Yale/bsc5.dat", "../Yale/bsc5.notes",
	&rval, &recsize);

    ptr = rval;
    printf("%5s: %8s  %8s  %7s %7s %5s %8s %8s\n",
	"i", "Yale", "SAO", "RA", "decl", "mag", "pmr", "pmd");
    for (i=0; i < count; ++i) {
	printf("%5d: %8d: %8ld: %7.2f %7.2f %5.1f %8.5f %8.5f  %s\n",
	    i, ptr->yale_cat, ptr->s.sao,
	    ptr->s.ra, ptr->s.dec, ptr->s.mag,
	    ptr->s.pmr, ptr->s.pmd,
	    ptr->s.name != NULL ? ptr->s.name : "");
	++ptr;
    }

    return 0;
}
#endif	/* STANDALONE */
