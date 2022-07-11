
Various astronomical utilities, by Edward Falk, based on the works of
Jean Meeus.

See `astro.h` for full notes

----

Functions defined in this library:

# Coords.c

Coordinate conversion routines. Based on Astronomical Formulae for
Calculators, by Jean Meeus, 4th edition, chapter 8.

## obliquity(date)

Find obliquity of eccliptic for a specified date.

## equat2ecliptic(double decl,RA, double \*lat, double \*lon, double jdate)
convert equatorial coordinates to ecliptic coordinates.

## ecliptic2equat(double lat,lon, double \*decl, double \*RA, double jdate)
convert ecliptic coordinates to equatorial

## equat2bearings(double decl,RA, double lat,lon, double A,h, double jdate)
convert equatorial coordinates to observer's local horizontal coords

# dates.c

Date conversion routines, from Astronomical Formulae for Calculators,
by Jean Meeus, 4th edition [1] and Astronimical Algorithms by Jean Meeus,
2nd edition [2].

## Overview:

In astronomy, dates are recorded as Julian dates, which are simply
days since noon, 1 Jan 4712 BC.  In the Julian system, the new day
starts at noon (GMT) instead of midnight.

Converting between Julian dates and standard (Gregorian) dates is
non-trivial, because of the various changes that have been made in
the civil calendar over the years.  For example, in September 1752,
11 days were removed from the calendar to make up for lack of leap year
adjustments (see cal(1) manpage).  More importantly, the modern
Gregorian calendar was created by Pope Gregory in 1582, which included
leap years.

Leap years themselves are fairly complex.  The length of the year
is not a whole number of days, but closer to 365.2424 days in
length (this varies, as the rate of the Earth's rotation is not
constant.)  To compensate, an extra day in inserted into the year
every fourth year, except for century years when it is not added,
except for fourth centuries when it is.  Thus, 1900 was not a leap
year, but 2000 was.

### Theory of Operations:
Meeus's book doesn't explain much, but this is how I think it
works.  The averge length of a month, ignoring February, is 30.60
days.  To get around the February problem, Meeus "rotates" the
year to put February at the end.  This is done by making January
and February the 13th and 14th months of the year and then subtracting
1 from the year.

Next, if the year is in the Gregorian calendar, the century #, a, is
computed from y/100.  The leap century # is a/4.  The number of skipped
leap years, b, is the number of centuries minus the number of leap
centuries, minus a fudge factor.

Finally, the Julian day is the year\*365.25 + month\*30.6001 + day +
julian(0,0,0) + b

BUG?  This code does not seem to allow for the September 1752 adjustment,
unless it's built into the Gregorian conversion.

### A few notes:

GMT is Greenwich Mean Time.

UT is Universal Time.
GMT and UT are both based on the rotation of the Earth, which is not
constant.  These are not uniform times.

Ephemeris time is based on the motion of the planets instead of the
rotation of the Earth.  It has been replaced by Dynamical Time.

Dynamical Time (TD) is based on atomic clocks.  There are actually
two Dynamical Times: Barycentric Dynamical Time, (TDB) based
on the center of mass of the solar system, and Terrestrial
Dynamical Time (TDT).  The difference is caused by relativistic effects
of the Earth's orbit.  Most of the time, you can ignore the difference.
The difference between TD and UT is determined by astronomical
observation.  As of 2004, TD - UT = 64.6 seconds.

Julian Days (JD) are days since the start of the year -4712.

Sidereal time (ST) measures the rotation of the Earth relative to
the vernal equinox rather than the Sun.  One sidereal day is 23 hours,
56 minutes, 4.091 seconds.  This is about .008 seconds shorter than
a day relative to the stars, due to the precession of the equinoxes.


## double date2julian(int y, int m, int d)
Convert (year, month, day) to a julian date

## double time2julian(int y, int m, int d, int h, int m, double s)

## void julian2date(double jdate, int \*y, int \*m, int \*d)
Convert a julian date to y/m/d

## void julian2time(double jdate, inty, int \*m, int \*d, int \*h, int \*m, double \*s)
Convert a julian date to y/m/d h:m:s GMT

## double unix2julian(time_t)
Convert unix time to julian

## double jnow()
Return the current julian date.

## double julian2sidereal(double jdate)
Return sidereal time in hours for a specific julian date
(Note: this is the sidereal time of midnight GMT.  I.e. draw a
line from the Sun through the center of the Earth and out
to the stars.  Wherever the line passes through the far side of
the Earth, that's where it's midnight.  Finding the sidereal
time at other longitudes will require further conversion.)
TODO: double-check this.  Where did I get it?

## double julianTime2sidereal(double jdate, double time)
Return sidereal time in hours for a specific julian date
and time other than 0h.  Time is specified in hours since midnight
TODO: double-check this.  Where did I get it?

## double time2sidereal(double jdate)
Return sidereal time at Grenwich

## double gmst2gast(double gmst, double JD)
Convert mean sidereal time to apparent sidereal time

## double siderealMean2Apparent(double jdate)
Return a correction factor to convert mean sidereal time in
hours to apparent sidereal time (correcting for nutation).
meanTime + siderealMean2Apparent() = apparentTime

## int date2yday(int y, int m, int d)
Obtain the year of the day from y/m/d

## void yday2date(int y, int yday,  intm, int \*d)
Obtain month and day from year and year of the day

# io.c

Assorted conversion utilities

## char \* julian2ymdStr(double jdate)
Convert Julian date to e.g. "11-Jul-2022"

## char \* julian2hmsStr(double jdate)
Convert Julian date to e.g. "16:04:16"

## char \* julian2str(double jdate)
Convert Julian date to e.g. "11-Jul-2022 16:04:16"

## char \* deg2dmsStr(double degrees)
Convert degrees to dd°mm'ss.s

## char \* deg2dmStr(double degrees)
Convert degrees to dd°mm.mm

## char \* hours2hmsStr(double hours)
Convert hours to hh:mm:ss.s

## char \* hours2hmStr(double hours)
Convert hours to hh:mm.mm

# kepler.c
Equation of Kepler, from Astronomical Formulae for Calculators,
by Jean Meeus, 4th edition, chapter 22

## double keplerE(double M, double e)
Solve for E.  All values are in radians.

# moon.c
Find the coordinates of the Moon, from Astronomical Formulae for
Calculators, by Jean Meeus, 4th edition, chapter 30.

## MoonPrecise(date, m)
Fill in the `PlanetState` structure **m** with the Moon's data for
the given Julian date.

## Moon(date, m)
Same as above, but lower precision.

# navigation.c
Assorted functions useful in celestial navigation

## double ra2sha(double ra)
convert Right Ascension (hours) to Sidereal Hour Angle (degrees)

## double sha2gha(double sha, double jdate)
convert Sidereal Hour Angle to Grenwich Hour Angle

## double gha2lha(double sha, double longitude)
convert Grenwich Hour Angle to Local Hour Angle

## double interpolate(double a, double b, double hours)
interpolate between a,b based on the fraction part of hours

## altaz(double lha, double decl, double lat, double \*alt, double \*az)
Compute altitude (Hc) and Azimuth (Z) from lha, declination, latitude

## double sext2obs(Hs, ie, h, T, P, HP)
Compute sextant altitude.

# planets.c
Find the coordinates of the Planets, from Astronomical Formulae for
Calculators, by Jean Meeus, 4th edition, chapter 23.

Meeus also gives formulae for computing the positions relative
to the standard ecliptics of 1950 and 2000, but I don't feel
like doing that right now.

Planetary orbits are described by these terms:

* L = mean longitude of planet
* a = semimajor axis of the orbit (a constant for each planet)
* e = eccentricity of orbit
* i = inclination on the plane of the ecliptic
* w = argument of perihelion
* omega = longitude of ascending node

* pi = longitude of perihelion = w + omega
* M = mean anomaly = L - pi
* q = perihelion distance = a\*(1-e)
* Q = aphelion distance = a\*(1+e)

Note: elements for longitude of perihelion taken from "ephem"
by Elwood Charles Downey.  These are probably just the sum of the
longitude of ascending node and argument of perihelion, but I don't
feel like double-checking right now.  Also, I'm not sure these will
work; I think he uses centuries since JD0 while Meeus uses
centuries since 1900

## void Mercury(double date, PlanetState \*p)
Fill in the `PlanetState` structure for the given Julian date.

## void Venus(double date, PlanetState \*p)

## void Earth(double date, PlanetState \*p)

## void Mars(double date, PlanetState \*p)

## void Jupiter(double date, PlanetState \*p)

## void Saturn(double date, PlanetState \*p)

## void Uranus(double date, PlanetState \*p)

## void Neptune(double date, PlanetState \*p)

# precession.c
Precession is the slow rotation of the Earth's axis of about 3 seconds
of right ascension per year.  One full rotation of the equinoxes takes
about 26,000 years.  Source: Meeus, chap. 14

References:

1. Astronomical Formulae for Calculators, * by Jean Meeus, 4th edition
1. Astronomical Algorithms by Jean Meeus, 2nd edition.

Nutation is the small elliptical wobble in the Earth's axis caused
by the Moon and other influences.  Nutation has a period of about 18.6
years, and an amplitude of about 9.2 seconds of arc.
Source: [1], chap 15, [2] chap 22

## precessionRate(double decl,RA, double \*dd, double \*dr, double jdate)
Accepts declination and right ascension in degrees and hours
respectively, and returns the annual differences, also in
degrees and hours

## precession(double decl0,RA0,jdate0, double \*decl1, double \*RA1,jdate1)
Accepts declination and right ascension for one date, and
returns declination and right ascension for another date.

## nutation(double \*psi, double \*eps, double jdate)
Given a julian date, return the nutation in longitude and
nutation in obliquity.  I have no idea what these really
mean, but they're used as inputs to other functions.
Return values are in seconds of arc.

# stars.c
Reads and returns the PPM (Positions and Proper Motions) star database.

## int ReadPPMStars(int maxmag, long ra0, long ra1, long d0, long d1, double jd, char \*filename, PPMStar \*\*rval)
Reads the PPM database. `maxmag, ra0, ra1, d0, d1` filter the data by limiting the maximum magnitude and specifying
a region of the sky. `jd` is the Julian date.

Allocates and returns an array of PPMStar structs. Return value is the number of stars in `rval`. `PPMStar` is
defined in `astro.h`

# sun.c
Find the coordinates of the Sun, from Astronomical Formulae for Calculators,
by Jean Meeus, 4th edition, chapter 18.

References:

1. Astronomical Formulae for Calculators, * by Jean Meeus, 4th edition
1. Astronomical Algorithms by Jean Meeus, 2nd edition.

Meeus also gives formulae for computing the Sun's position relative
to the standard ecliptic of 1950, but I don't feel like doing that
right now.

These functions return the mean position of the Sun.  I have not
written any code for the apparent position.

*Note:* Consider using VOSP87 instead.

## void SunEcliptic(double jdate, double \*lat, double \*lon, double \*rad)

Return lat,lon,radius from the date.  lat,lon are in degrees,
radius is in AU.  This is the Sun's coordinates relative to
the Earth in ecliptic coordinates.  Take the complement of
lat,lon to get the Earth's coordinates relative to the Sun.
Note: by definition, lat is always 0 for the mean ecliptic,
but I may change this function later to include perturbations.

## void      SunEquatorial(double jdate, double \*decl, double \*RA, double \*rad)

Return declination, right ascension and distance of the
Sun for a given date.  See [2], ch 25.  Accurate to about 0.01
degree.  If you need more accuracy, see [2], ch 26 or use the Vosp87
values. RA given in hours, declination given in degrees.

# utils.c

Various small utilities

## const char * convertHms(double hours)

convert hours to "hh:mm:ss.ssss"


## void h2hms(hours, h,m,s)

convert hours (or degrees) to hours:minutes:seconds

## void printHms(double hours)


## void polar2rect(lat, lon, R, X,Y,Z)
Convert spherical (e.g. ecliptic) coordinates to rectangular

## void rect2polar(X,Y,Z, lat,lon,R)
Convert rectangular coordinates to spherical (e.g. ecliptic)

## void deltaPolar(lat1,lon1,r1, lat2,lon2,r2, lat3,lon3,r3)

given polar coordinates of two objects, find the bearing
and distance of object2 relative to object1

## double limitAngle(double a)

Return a % 360

## double limitHour(double h)

Return h % 24


----

# ephem

Simple utility program to print out the positions of the Sun, Moon,
and planets (except Pluto). With no arguments, shows the current position. You can also
specify date and time. Run with `-h` for help.

Note that this program does not currently correct for speed-of-light delays.
That is, the output specifies where the various bodies are at the specified
time, but an observer on Earth is seeing where they were as much as six hours
ago.

# ppm

Demo program that shows how to read and parse the Positions & Proper Motion
star catalog.

