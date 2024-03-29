using MathNet.Numerics;

namespace AtmosTools;

public static class Geodesy
{
    public static double Haversine(double theta)
    {
        double x = Math.Sin(theta / 2.0);
        return x * x;
    }
    public static double CentralAngle(double lon0, double lat0, double lon1, double lat1)
    {
        // Returns the angle difference (central angle) between two locations on Earth, in radians
        double lon0Rad = Trig.DegreeToRadian(lon0);
        double lat0Rad = Trig.DegreeToRadian(lat0);
        double lon1Rad = Trig.DegreeToRadian(lon1);
        double lat1Rad = Trig.DegreeToRadian(lat1);
        double deltaLongitude = lon1Rad - lon0Rad;
        double deltaLatitude = lat1Rad - lat0Rad;
        double haversineDeltaLatitude = Haversine(deltaLatitude);
        double cosMeanLat = Math.Cos(0.5 * (lat0Rad + lat1Rad));
        double innerTerm = haversineDeltaLatitude + Haversine(deltaLongitude) * (cosMeanLat * cosMeanLat - haversineDeltaLatitude);
        double chordLength = 2.0 * Math.Sqrt(innerTerm);
        return 2.0 * Trig.Asin(chordLength / 2.0);
    }

    public static double GreatCircleDistance(double lon0, double lat0, double lon1, double lat1)
    {
        // Distance is in km
        return CentralAngle(lon0, lat0, lon1, lat1) * PhysConstants.EarthRadius;
    }

    public static (double[], double[]) GreatCircleWaypointsInaccurate(double lon0, double lat0, double lon1, double lat1,
        int nPoints)
    {
        /* Generate nPoints locations (vector of longitudes, vector of latitudes) corresponding to roughly
        equally-spaced waypoints along a great circle. nPoints must be >= 2 as the first and last waypoints are
        the start and end of the arc respectively. Likely to be very inaccurate. Uses the algorithm described at
        https://www.reddit.com/r/algorithms/comments/xljj7k/comment/ipjtx7o/?utm_source=share&utm_medium=web3x&utm_name=web3xcss&utm_term=1&utm_content=share_button
        */
        if (nPoints < 2) { throw new ArgumentException("Must request at least 2 waypoints"); }

        if (nPoints == 2)
        {
            return (new double[] { lon0, lon1 }, new double[] { lat0, lat1 });
        }

        double[] lons = new double[nPoints];
        double[] lats = new double[nPoints];
        lons[0] = lon0;
        lons[nPoints - 1] = lon1;
        lats[0] = lat0;
        lats[nPoints - 1] = lat1;
        
        // Convert to radians
        double lon0Rad = Trig.DegreeToRadian(lon0);
        double lat0Rad = Trig.DegreeToRadian(lat0);
        double lon1Rad = Trig.DegreeToRadian(lon1);
        double lat1Rad = Trig.DegreeToRadian(lat1);
        
        // Convert polar coordinates into cartesian on the unit sphere
        (double x0, double y0, double z0) = PolarToCartesian(lon0Rad, lat0Rad);
        (double x1, double y1, double z1) = PolarToCartesian(lon1Rad, lat1Rad);
        // Interpolate directly along the line between the two points
        double[] x = new double[nPoints];
        double[] y = new double[nPoints];
        double[] z = new double[nPoints];
        double[] xs = new double[nPoints];
        double[] ys = new double[nPoints];
        double[] zs = new double[nPoints];
        double dx = x1 - x0;
        double dy = y1 - y0;
        double dz = z1 - z0;
        for (int i = 0; i < nPoints; i++)
        {
            double frac = ((double)i / (double)(nPoints - 1));
            x[i] = x0 + frac * dx;
            y[i] = y0 + frac * dy;
            z[i] = z0 + frac * dz;
        }
        // Project back onto the sphere
        for (int i = 0; i < nPoints; i++)
        {
            double magnitude = Math.Sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
            xs[i] = x0 / magnitude;
            ys[i] = y0 / magnitude;
            zs[i] = z0 / magnitude;
        }
        
        for (int i = 0; i < nPoints; i++)
        {
            (double lonR, double latR) = CartesianToPolar(xs[i], ys[i], zs[i]);
            lons[i] = Trig.RadianToDegree(lonR);
            lats[i] = Trig.RadianToDegree(latR);
        }
        return (lons, lats);
    }
    
    public static (double[], double[], double[]) GreatCircleWaypoints(double lon0, double lat0, double lon1, double lat1,
        double distancePerWaypoint)
    {
        /* Generate nPoints locations (vector of longitudes, vector of latitudes) corresponding to exactly
        equally-spaced waypoints along a great circle. Largest source of inaccuracy is due to the treatment of Earth
        as a sphere rather than an ellipsoid.
        */
        // Get the great circle distance between the start and end. Since we will also need the
        // central angle, we break up the calculation here rather than using the GCD directly
        double totalAngle = CentralAngle(lon0, lat0, lon1, lat1);
        double startEndDistance = totalAngle * PhysConstants.EarthRadius;
        
        // Add one segment to deal with any leftover
        int nSegments = (int)Math.Floor(startEndDistance / distancePerWaypoint) + 1;
        int nPoints = nSegments + 1;

        double[] lons = new double[nPoints];
        double[] lats = new double[nPoints];
        double[] segmentLengths = new double[nSegments];
        lons[0] = lon0;
        lats[0] = lat0;
        lons[nPoints - 1] = lon1;
        lats[nPoints - 1] = lat1;
        for (int i = 0; i < (nSegments-1); i++)
        {
            segmentLengths[i] = distancePerWaypoint;
        }
        // Last segment is truncated
        segmentLengths[nSegments - 1] = startEndDistance - segmentLengths.Sum();
        
        // Convert to radians
        double lon0Rad = Trig.DegreeToRadian(lon0);
        double lat0Rad = Trig.DegreeToRadian(lat0);
        double lon1Rad = Trig.DegreeToRadian(lon1);
        double lat1Rad = Trig.DegreeToRadian(lat1);
        
        // Get the initial course
        double deltaLonRad = lon1Rad - lon0Rad;
        
        // These will be heavily reused
        double sinLat0 = Trig.Sin(lat0Rad);
        double cosLat0 = Trig.Cos(lat0Rad);
        
        double numerator = Trig.Cos(lat1Rad) * Trig.Sin(deltaLonRad);
        double denominator = cosLat0 * Trig.Sin(lat1Rad) -
                             sinLat0 * Trig.Cos(lat1Rad) * Trig.Cos(deltaLonRad);
        double initialCourse = Math.Atan2(numerator,denominator);
        // Alternative formula for the great circle distance
        //double gcd = Math.Atan2(Math.Sqrt(denominator * denominator + numerator * numerator),
        //                        Trig.Sin(lat0Rad) * Trig.Sin(lat1Rad) + Trig.Cos(lat0Rad) * Trig.Cos(lat1Rad) * Trig.Cos(deltaLonRad));
        
        // These will also be heavily used
        double sinCourse = Trig.Sin(initialCourse);
        double cosCourse = Trig.Cos(initialCourse);
        
        // Skip first and last (already known)
        for (int i = 1; i < (nPoints - 1); i++)
        {
            // Distance from first waypoint to this one
            double arcDistance = distancePerWaypoint * i;
            // Absolute angle to be covered
            double arcAngle = arcDistance / PhysConstants.EarthRadius;
            double sinArcAngle = Trig.Sin(arcAngle);
            double cosArcAngle = Trig.Cos(arcAngle);
            // This part seems fine
            double sinLat = sinLat0 * cosArcAngle + cosLat0 * sinArcAngle * cosCourse;
            double tdlNumerator = sinArcAngle * sinCourse;
            double tdlDenominator = cosLat0 * cosArcAngle - sinLat0 * sinArcAngle * cosCourse;
            double segmentDeltaLonRad = Math.Atan2(tdlNumerator, tdlDenominator);
            //Console.WriteLine($"Segment longitude change: {segmentDeltaLonRad * 180.0/Math.PI}");
            double newLon = Trig.RadianToDegree(lon0Rad + segmentDeltaLonRad);
            double newLat = Trig.RadianToDegree(Trig.Asin(sinLat));
            // The % operator is remainder, not modulo, in C# so cannot be used here
            while (newLon < -180.0)
            {
                newLon += 360.0;
            }

            while (newLon >= 180.0)
            {
                newLon -= 360.0;
            }
            lons[i] = newLon;
            lats[i] = newLat;
        }
        
        return (lons, lats, segmentLengths);
    }
    
    public static (double, double, double) PolarToCartesian(double lon, double lat)
    {
        // Assumes lon and lat are in radians; converts for the unit sphere
        return (Math.Cos(lon) * Math.Cos(lat), Math.Sin(lat), Math.Sin(lon) * Math.Cos(lat));
    }

    public static (double, double) CartesianToPolar(double x, double y, double z)
    {
        // Take x, y, and z for the unit sphere and return lon and lat in radians
        return (Math.Atan(z / x), Math.Asin(y));
    }
}