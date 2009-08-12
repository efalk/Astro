/*
 *
 *General astronomical functions, from Astronomical Formulae for Calculators,
 *by by Jean Meeus, 4th edition.
 *
 *
 *General background:
 *
 *  There are several different coordinate systems used in astronomy, all
 *  different and all relatively unstable, because of precession, perturbation
 *  and so on.
 *
 *  Galactic Coordinates:
 *    This is a spherical coordinate system centered at the center of the
 *    galaxy.  The equator of this system is the mean galactic equator,
 *    and galactic "north" is in approximately the same direction as
 *    terrestrial north.  This coordinate system is not used here.
 *
 *  Ecliptic Coordinates:
 *    This is a spherical coordinate system either centered at the Earth
 *    ("geocentric") or centered at the Sun ("heliocentric").  The
 *    equator of this system is the plane of the Earth's orbit, with
 *    "north" 23.5 degrees off of the Earth's north.
 *
 *    Angles north and south are called latitude, and are measured in
 *    degrees, +north/-south.
 *
 *    Angles about the center are called longitude, and are measured in
 *    degrees.  Zero degrees longitude is at the first point in Aries,
 *    meaning the point where the Earth's equator crosses the ecliptic at
 *    the vernal equinox.  This point is no longer actually in Aries due
 *    to precession of Earth's rotation, which causes this point to
 *    continuously drift east about one degree/century.  Ecliptic
 *    longitudes increase to the east.
 *
 *    The name "ecliptic" comes from the fact that eclipses of the Moon
 *    or Sun can only happen when the moon is on the ecliptic.  By
 *    definition, the Sun is always on the ecliptic.
 *
 *  Celestial Coordinates:
 *    This is a spherical coordinate system centered at the Earth.  The
 *    equator of this system is the Earth's equator and the north pole
 *    of this system is the Earth's north pole.  Angle north and south
 *    are called declination, and are measured in degrees, +north/-south.
 *
 *    Angles about the equator are called right ascension, and are
 *    measured in hours, minutes and seconds.  0h right ascension is the
 *    first point in Aries, as described above.  Right ascension
 *    increases to the east.
 *
 *
 *  Rectangular Coordinates.
 *    Normally, coordinates are given in terms of lattitude, longitude
 *    (or declination and right ascension) and distance.  Occasionally
 *    rectangular coordinates are used.  In such cases, the axes are:
 *	X - towards vernal equinox = 0 degrees longitude
 *	Y - 90 degrees longitude
 *	Z - North
 *
 *
 *  Other Terms:
 *
 *    Meridian: line in the sky, directly over the observer's head,
 *	running north & south.
 *
 *    Zenith: point directly overhead.
 *
 *    Nadir: point directly down
 *
 *    Altitude: angle above the horizon
 *
 *    Azimuth: angle along the horizon, measured to the right of North
 *	(sometimes South, depending on author.)  0° = North, 90° = East,
 *	etc.
 *
 *    Hour Angle: an object's angle west of the meridian, measured in h:m:s.
 *	(= difference between object's RA and current sidereal time.)
 *	Hour Angle increases to the west.
 *	Hour Angle is measured in degrees for celestial navigation.
 *
 *    SHA:  Sidereal Hour Angle.  Object's angle west of Aries.
 *
 *    GHA:  Grenwich Hour Angle.  An object's angle west of Grenwich, England.
 *
 *    LHA:  Local Hour Angle.  An object's angle west of the local meridian.
 *
 *
 *  Time can be rather difficult to deal with, as the calendar has been
 *  changed several times over the centuries.  The current system is the
 *  Gregorian calendar.  The functions in dates.c convert civil
 *  calendars into Julian dates.  All work in the program is done in
 *  Julian time.
 *
 *  Julian time is simply days since noon, 1 Jan 4712 BC.  In the Julian
 *  system, the new day starts at noon (GMT) instead of midnight.
 *
 *  Sometimes, time is expressed in 36525-day centuries, often centuries
 *  since Jan 1, 1900.  The formula for this is
 *
 *		jd - 2415020
 *		------------
 *		    36525
 *
 *  Many measurements are made relative to the vernal equinox, which
 *  is the point where the ecliptic intersects the equator.  This
 *  point drifts constantly (TODO: which direction?) due to
 *  precession.  A full circle takes approximately 26,000 years.
 *  The North pole is currently moving towards Polaris, and should
 *  be within 28 minutes about the year 2102.  Vega should become the
 *  pole star around the year 14000.
 *
 *  Furthermore, there is a small elliptical oscillation in the Earth's
 *  axis caused by the Moon, called "nutation", which is 18.6 years in
 *  period and 9.2 seconds in amplitude.  This causes the vernal equinox
 *  to wobble back and forth.  Some measurements are relative to the
 *  exact vernal equinox, and some are relative to the average vernal
 *  equinox.
 *
 *  Often, to keep things simple, a "standard" equinox, typically from
 *  1950 or 2000, is used as the base point.
 *
 *
 *  Year:  These are the various years used in astronomy.
 *
 *    Tropical year: time for the earth to revolve from vernal equinox
 *	to vernal equinox.  Approx 365.24219 days, as of 1970.  This
 *	is constantly changing by very small amounts.
 *
 *    Julian year: 365.25 days.  Requires the addition of one leap day
 *	every four years.
 *
 *    Gregorian year: 365.2425 days.  Requires the addition of one leap
 *	day every four years except century years, except every fourth
 *	century.  Thus, 1900 was not a leap year, but 2000 will be.
 *	(Historical note:  John Hamilton Moore forgot this rule in
 *	his book "The Practical Navigator" and called 1800 a leap
 *	year.  This mistake caused the loss of at least one ship.
 *	Call it a Y1800 problem.)
 *
 *  Month:
 *    The Moon revolves around the Earth every 27 days (TODO: exact).
 *    This is a "sidereal month".
 *
 *    A solar month is the time from new moon to new moon.  Approximately
 *    29 days.
 *
 *  Day:
 *    The Earth rotates around its axis in 23:56:04.  (In this time, a
 *    star which was overhead will return to it's exact previous
 *    location.) However in this time, it also revolves around the Sun
 *    nearly one degree.  It takes another 3:56 minutes (average) to
 *    "catch up" with the Sun, for a total of (average) 24 hours.
 *
 *    The definition of a day can be somewhat elusive.  The Earth does
 *    not rotate at a constant rate -- it is slowing down, and with
 *    unpredictable irregularities.  Further, because the Earth's orbit
 *    is elliptical, the length of the day actually varies nearly a
 *    minute throughout the year.
 *
 *    Because the Earth's orbit around the Sun is elliptical, the Earth
 *    does not move around the Sun at a constant speed.  It speeds up as
 *    it gets near to the Sun, and slows down as it gets further away.
 *    As the Earth speeds up, the apparent length of the day increases,
 *    and as it slows down, the apparent length of the day decreases.
 *    The 24 hour day is only an average over the year.
 *
 *    Sidereal Day:  Time for the earth to rotate one full rotation
 *	relative to the vernal equinox.  23:56:4.099
 *
 *    Sidereal Time: hours since the vernal equinox passed overhead.
 *	(= right ascension of meridean = hour angle of vernal equinox)
 *
 *    Mean Sidereal Time: sidereal time based on mean equinox (corrected
 *	for nutation)
 *
 *    True Sidereal Time: sidereal time based on exact equinox.
 *
 *    Solar Time: hours since the Sun passed overhead, +12
 *	(= hour angle of Sun + 12)
 *
 *    Apparent Solar Time: solar time based on actual position of Sun.
 *
 *    Mean Solar Time: solar time based on imaginary average Sun.
 *
 *    The total difference between apparent solar time and mean solar time
 *    can be from zero to sixteen minutes.  The difference is computed
 *    with the "equation of time"
 *
 *    Ascending Node:  The point in a planet's orbit where it crosses the
 *    plane of Earth's orbit (the ecliptic) from south to north.
 *
 *    Descending Node:  The point where a planet passes the ecliptic
 *    from north to south.
 *
 *
 *
 *
 *  Orbits:
 *
 *    Planets orbit around the sun; moons orbit around their planets.
 *
 *    Orbits are not perfect circles (although Venus comes pretty close.)
 *    Instead, they are ellipses, with the primary at one focus.  The
 *    point in the orbit closest to the primary is perihelion and the
 *    furthest point is aphelion.  Not all orbits lie in the same plane;
 *    all planets have some amount of inclination relative to the
 *    ecliptic plane (Earth's orbital plane.)  There is no relationship
 *    between the angle of the orbit's major axis and the orbit's
 *    inclination.
 *
 *    If orbits were perfectly circular, orbits were not tilted relative
 *    to the ecliptic, perturbations caused by other planets did not
 *    exist, and the equinoxes did not precess, then the position of a
 *    planet relative to the Sun could be described as a point on a
 *    circle:
 *
 *	latitude = 0
 *	longitude = T*L1 + L0
 *
 *	(where T = time, L1 = angular rate, L0 = base angle)
 *
 *    In reality, the equation for the lat/lon of a planet is very
 *    complex, and varies from planet to planet.  The basic elements are:
 *
 *	Epoch	Base time.  Aka t0.
 *	a	semimajor axis: 1/2 length of major axis of planet's orbit.
 *		This is a constant for each planet.
 *	e	eccentricity: amount orbit varies from circular.  0 =
 *		perfect circle, 1=parabola.
 *	i	inclination of orbit relative to ecliptic
 *	omega	longitude of orbit's ascending node.  Aka W0.
 *	N0	mean motion.  Average speed of orbit, typically given
 *		as degrees/day, orbits/day, orbits/year, years/orbit, etc.
 *	M	Mean anomoly.  Angular distance from perihelion, as if
 *		the planet had a circular orbit.  0 = perihelion, 180 =
 *		aphelion.
 *	w	argument of perihelion.  The angle from the ascending
 *		node to the perihelion.
 *
 *
 *
 *	pi	longitude of perihelion = omega + w
 *	L	mean longitude of planet relative to Sun, as if the
 *		planet had a circular orbit.
 *
 *    Notes:
 *      "Ascending node" refers to the point where the planet crosses
 *      the ecliptic from south to north.
 *
 *      Inclination, longitude of ascending node define the plane of
 *      the orbit.
 *
 *      a,e,w define the shape of the orbit within that plane.  A defines
 *      the size, e defines the shape, w sets the angle of the major axis.
 *
 *      When describing satellite orbits, longitude of ascending node is
 *      replaced by Right Ascension of Ascending Node ("RAAN")
 *
 *      Excellent tutorial: http://www.amsat.org/amsat/keps/kepmodel.html
 */


#include <math.h>

#define	JD1900	2415020.0
#define	JD1950	2433282.423		/* standard equinox of 1950 */
#define	JD2000	2451545.0		/* Jan 0.5, 2000 */
#define	JDUnix	2440587.5		/* Unix epoch: Jan 0.0 1970 */

#define	EarthTilt	23.4392911	/* degrees, epoch JD2000 */

#define	RAD	(M_PI/180.)
#define	DEG	(180./M_PI)

#define	dsin(x)	(sin((x)*RAD))
#define	dcos(x)	(cos((x)*RAD))
#define	dtan(x)	(tan((x)*RAD))
#define	dasin(x) (DEG*asin(x))
#define	dacos(x) (DEG*acos(x))
#define	datan(x) (DEG*atan(x))
#define	datan2(x,y) (DEG*atan2(x,y))


	/* describe the current state of a planet.
	 * All angles in degrees.
	 * Exceptions:  For Sun, this structure is geocentric.
	 *	For earth satellites, this structure is geocentric,
	 *	distances in Earth diameters.
	 */

typedef	struct {
	  double date ;	/* julian date */
	  double L ;	/* mean longitude */
	  double dL ;	/* daily motion in longitude */
	  double e ;	/* eccentricity */
	  double i ;	/* inclination */
	  double om ;	/* longitude of ascending node */
	  double w ;	/* argument of perihelion */
	  double pi ;	/* longitude of perihelion */
	  double a ;	/* length of semi-major axis */
	  double ad ;	/* angular diameter from 1 AU, arc seconds */
	  double mag ;	/* magnitude at 1 AU */
	  double M ;	/* mean anomoly */
	  double v ;	/* true anomoly */
	  double lon ;	/* heliocentric longitude */
	  double lat ;	/* heliocentric lattitude */
	  double R ;	/* distance from Sun, AU */
	  double year ;	/* length of the year, in days */
	} PlanetState ;


	/* this structure describes items in the Yale Star Catalog */

typedef	struct {
	  long	ra ;		/* right ascension, seconds of arc, not time */
	  long	dec ;		/* declination, seconds of arc */
	  int	mag ;		/* magnitude * 100 */
	  char	type[2] ;	/* object type */
	  char	spec[2] ;	/* spectral type */
	  char	cons[3] ;	/* constellation */
	  char	*name ;		/* name, if any, or NULL */
	} YaleStar ;

	/* types */

/* CG-Glob.Cluster	CO-Open  Cluster	 GC-Galac.Cluster
 * GP-Sphere Galaxy	GS-Spiral Galaxy
 * ND-Difuse Nebula	NP-Planetary Nebula
 * P*-Planet
 * SS-star		SB-Binary Star		SD-Double Star
 * SV-Variable Star
 * VM-vector move	VS-vector draw (solid)	VD-vector draw (dotted)
 * VH-vector draw (hyphens = dashed)
 * I*-Invisible (for annotation)
 */

	/* this structure describes items in the PPM catalog */

typedef	struct {
	  long	ra ;		/* right ascension, seconds of arc, not time */
	  long	dec ;		/* declination, seconds of arc */
	  int	mag ;		/* magnitude * 100 */
	  char	type[2] ;	/* object type */
	  char	spec[2] ;	/* spectral type */
	} PPMStar ;



#ifdef	__STDC__

	/* time conversions */
extern	double	date2julian(int y, int m, int d) ;
extern	double	time2julian(int y, int m, int d, int hr, int mn, double s) ;
extern	void	julian2date(double jdate, int *y, int *m, int *d) ;
extern	void	julian2time(double jdate, int *y, int *m, int *d,
			int *hr, int *mn, double *s) ;
extern	int	date2yday(int y, int m, int d) ;
extern	void	yday2date(int y, int yday,  int *m, int *d) ;
extern	double	julian2sidereal(double jdate) ;
extern	double	julianTime2sidereal(double jdate, double time) ;
extern	double	time2sidereal(double jdate) ;
extern	double	siderealMean2Apparent(double jdate) ;
extern	double	unix2julian(time_t) ;
extern	double	jnow() ;

	/* coordinate conversions */
extern	void	equat2ecliptic(double decl, double RA,
			double *lat, double *lon, double jdate) ;
extern	void	ecliptic2equat(double lat, double lon,
			double *decl, double *RA, double jdate) ;
extern	void	equat2bearings(double decl, double RA,
			double lat, double lon, double *A, double *h,
			double jdate, double time) ;
extern	double	obliquity(double jdate) ;

	/* precession and nutation */
extern	void	precessionRate(double decl, double RA, double *dd, double *dr,
			double jdate) ;
extern	void	precession(double decl0, double RA0, double jdate0,
			   double *decl1, double *RA1, double jdate1) ;
extern	void	nutation(double *psi, double *eps, double jdate) ;

	/* Coordinates of the Sun */
extern	void	SunEcliptic(double jdate, double *lat,double *lon,double *rad);
extern	void	SunEquatorial(double jd, double *decl,double *RA,double *rad);

	/* Coordinates of the Planets (relative to the Sun) */
extern	void	Mercury(double date, PlanetState *p) ;
extern	void	Venus(double date, PlanetState *p) ;
extern	void	Earth(double date, PlanetState *p) ;
extern	void	Mars(double date, PlanetState *p) ;
extern	void	Jupiter(double date, PlanetState *p) ;
extern	void	Saturn(double date, PlanetState *p) ;
extern	void	Uranus(double date, PlanetState *p) ;
extern	void	Neptune(double date, PlanetState *p) ;
extern	void	Pluto(double date, PlanetState *p) ;

extern	void	MoonPrecise(double date, PlanetState *p) ;
extern	void	Moon(double date, PlanetState *p) ;

	/* star databases */
extern	int	ReadPPMStars(int maxmag, long ra0, long ra1, long d0, long d1,
			double jd, char *filename, PPMStar **rptr) ;

	/* navigation */

extern	double	ra2sha(double RA) ;
extern	double	sha2gha(double SHA, double jdate) ;
extern	double	gha2lha(double SHA, double jdate) ;


	/* utilities */
extern	double	hms2h(int h, int m, double s) ;
extern	void	h2hms(double hh, int *h, int *m, double *s) ;
extern	void	polar2rect(double lat,double lon,double R,
			double *X, double *Y, double *Z) ;
extern	void	rect2polar(double X, double Y, double Z, 
			double *lat, double *lon, double *R) ;
extern	void	deltaPolar(double lat1, double lon1, double r1,
			double lat2, double lon2, double r2,
			double *lat3, double *lon3, double *r3) ;
extern	double	keplerE(double M, double e) ;
extern	double	limitAngle(double angle) ;
extern	double	limitHour(double angle) ;

	/* I/O */

extern	char	*julian2ymdStr(double jdate) ;
extern	char	*julian2hmsStr(double jdate) ;
extern	char	*julian2str(double jdate) ;
extern	char	*deg2dmsStr(double jdate) ;
extern	char	*deg2dmStr(double jdate) ;
extern	char	*hours2hmsStr(double jdate) ;
extern	char	*hours2hmStr(double jdate) ;

#else
extern	double	date2julian() ;
extern	double	time2julian() ;
extern	void	julian2date() ;
extern	void	julian2time() ;
extern	int	date2yday() ;
extern	void	yday2date() ;
extern	double	julian2sidereal() ;
extern	double	julianTime2sidereal() ;
extern	double	siderealMean2Apparent() ;
extern	double	unix2julian() ;
extern	double	jnow() ;
extern	void	equat2ecliptic() ;
extern	void	ecliptic2equat() ;
extern	void	equat2bearings() ;
extern	void	precessionRate() ;
extern	void	precession() ;
extern	void	nutation() ;
extern	void	SunEcliptic() ;
extern	void	SunEquatorial() ;
extern	double	hms2h() ;
extern	void	h2hms() ;
extern	void	polar2rect() ;
extern	void	rect2polar() ;
extern	void	deltaPolar() ;
extern	double	keplerE() ;
extern	void	Mercury() ;
extern	void	Venus() ;
extern	void	Earth() ;
extern	void	Mars() ;
extern	void	Jupiter() ;
extern	void	Saturn() ;
extern	void	Uranus() ;
extern	void	Neptune() ;
extern	void	Pluto() ;
extern	void	MoonPrecise() ;
extern	void	Moon() ;
extern	int	ReadPPMStars() ;
#endif
