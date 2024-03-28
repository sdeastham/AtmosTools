using MathNet.Numerics;

namespace AtmosTools;

public static class Distances
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
}