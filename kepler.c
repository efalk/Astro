
#include <sys/types.h>
#include <time.h>
#include <math.h>

#include "astro.h"

/* Equation of Kepler, from Astronomical Formulae for Calculators,
 * by Jean Meeus, 4th edition, chapter 22.
 *
 * The equation of Kepler is
 *
 *	E = M + e*sin(E)
 *
 * where e is the eccentricity of the planet's orbit, M is the planet's
 * mean anomaly at a given instant, and E is the eccentric anomaly.
 *
 *
 * double
 * keplerE(double M, double e)
 *	Solve for E.  All values are in radians.
 */


#define	ACCURACY	(.000001*RAD)	/* desired accuracy */
#define	MAXLOOPS	40		/* give up after this many loops */

double
keplerE(M, e)
	double	M,e ;
{
	double	E ;	/* current guess */
	double	delta ;	/* new guess */
	int	i ;

	E = M ;
	for(i=0; i<MAXLOOPS; ++i)
	{
	  delta = (M + e*sin(E) - E) / (1 - e*cos(E)) ;
	  E += delta ;
	  if( delta >= -ACCURACY && delta <= ACCURACY )
	    break ;
	}
	return E ;
}
