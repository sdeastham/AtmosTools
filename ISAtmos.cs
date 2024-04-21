using MathNet.Numerics;
using MathNet.Numerics.Interpolation;

namespace AtmosTools;

public static class ISAtmos
{
    // Reference altitudes in km
    private static readonly double[] ReferenceAltitudes = [-0.1, 0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.8520, 1.0e6];
    // Reference lapse rates in K/km
    private static readonly double[] ReferenceLapseRates = [0.0, -6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0];
    private const double ReferenceSurfaceTemperature = 288.15;

    private static double[] _referenceTemperatures;
    // Standard gases:                N2,      O2,      Ar,     CO2,      Neon,   He,     Kr,    Xe,     CH4,      H2
    private static double[] _gasMw = [28.0134, 31.9988, 39.948, 44.00995, 20.183, 4.0026, 83.80, 131.30, 16.04303, 2.01594];
    private static double[] _gasFrac;
    private static double _gasPremultiplier;
    private static IInterpolation _temperatureInterp;
    private static IInterpolation _pressureInterp;

    static ISAtmos()
    {
        int nReference = ReferenceAltitudes.Length;
        int nDeltas = nReference - 1;
        
        // Derive the temperature at each reference altitude (linear interpolation between)
        double[] temperatures = new double[nReference];
        temperatures[0] = ReferenceSurfaceTemperature;
        for (int i = 0; i < nDeltas; i++)
        {
            double heightDelta = ReferenceAltitudes[i + 1] - ReferenceAltitudes[i];
            double temperatureDelta = ReferenceLapseRates[i] * heightDelta;
            temperatures[i + 1] = temperatures[i] + temperatureDelta;
        }
        
        // Create the reference gas composition
        double[] gasFrac = [0.78084, 0.209476, 0.00934, 0.000314, 0.00001818, 0.00000524, 0.00000114, 0.000000087, 0.000002, 0.0000005];
        // Normalize
        double gasCorrection = 1.0 / gasFrac.Sum();
        for (int i = 0; i < gasFrac.Length; i++)
        {
            gasFrac[i] = gasFrac[i] * gasCorrection;
        }
        
        // Copy in for storage
        _gasFrac = gasFrac;
        _referenceTemperatures = temperatures;
        
        // Create interpolation for temperature
        _temperatureInterp = Interpolate.Linear(ReferenceAltitudes, _referenceTemperatures);

        double[] gasProduct = new double[gasFrac.Length];
        for (int i = 0; i < gasFrac.Length; i++)
        {
            gasProduct[i] = gasFrac[i] * _gasMw[i];
        }
        _gasPremultiplier = gasProduct.Sum() * 1.0e-3 * PhysConstants.G0 / PhysConstants.RGasUniversal;
        
        // Now an interpolation for pressure
        // Create an array of altitudes in 1 m increments from -1 km to +82 km, then represent as km
        double[] testAltitudes = Enumerable.Range(-1000, 82000).ToArray().Select(x => (double)x * 0.001).ToArray();
        double[] referencePressures = AltitudeToPressure(testAltitudes);
        _pressureInterp = Interpolate.Linear(testAltitudes, referencePressures);
    }
    
    public static double[] AltitudeToTemperature(double[] altitudes)
    {
        // Input: vector of altitudes in km
        // Output: vector of temperatures in K
        int nInput = altitudes.Length;
        double[] outputTemperatures = new double[nInput];
        for (int i = 0; i < nInput; i++)
        {
            outputTemperatures[i] = _temperatureInterp.Interpolate(altitudes[i]);
        }
        return outputTemperatures;
    }
    public static double[] AltitudeToPressure(double[] altitudes)
    {
        // Input: vector of altitudes in km
        // Output: vector of pressures in Pa
        int nInput = altitudes.Length;
        double[] outputPressures = new double[nInput];
        
        // Sort the data using algorithm from https://stackoverflow.com/a/1760202
        var sorted = altitudes
            .Select((x, i) => new KeyValuePair<double, int>(x, i))
            .OrderBy(x => x.Key)
            .ToList();
        double[] sortedAltitudes = sorted.Select(x => x.Key).ToArray();
        int[] sortingIndices = sorted.Select(x => x.Value).ToArray();
        // Now, altitudes[sortingIndices[i]] = sortedAltitudes[i]
        
        // First - determine the temperature at each requested altitude
        double[] sortedTemperatures = AltitudeToTemperature(sortedAltitudes);
        double[] outputTemperatures = new double[nInput];
        for (int i = 0; i < nInput; i++)
        {
            outputTemperatures[sortingIndices[i]] = sortedTemperatures[i];
        }

        int iBase = 0;
        int iCeiling = 1;
        double zBase_m = ReferenceAltitudes[0] * 1000.0;
        double zCeiling_m = ReferenceAltitudes[1] * 1000.0;
        double temperatureBase = _referenceTemperatures[0];
        double pressureBase = 101325 * Math.Exp(-1.0 * _gasPremultiplier * zBase_m / temperatureBase);
        double localLapseRate = 0.0; // Lapse rate in K/m
        // Iterate through the sorted values
        for (int i = 0; i < nInput; i++)
        {
            double zCurrent_m = sortedAltitudes[i] * 1000.0;
            while (zCurrent_m > zCeiling_m)
            {
                double pressureNew;
                if (Math.Abs(localLapseRate) > 0)
                {
                    pressureNew = pressureBase * Math.Pow(_referenceTemperatures[iCeiling] / _referenceTemperatures[iBase],
                        _gasPremultiplier / (-1.0 * localLapseRate));
                }
                else
                {
                    pressureNew = pressureBase * Math.Exp(_gasPremultiplier * (zBase_m - zCeiling_m) / temperatureBase);
                }
                pressureBase = pressureNew;
                iBase = iCeiling;
                iCeiling++;
                zBase_m = zCeiling_m;
                zCeiling_m = ReferenceAltitudes[iCeiling] * 1000.0;
                temperatureBase = _referenceTemperatures[iBase];
                localLapseRate = ReferenceLapseRates[iBase] / 1000.0;
            }

            double pressure;
            if (Math.Abs(localLapseRate) > 0)
            {
                pressure = pressureBase * Math.Pow(sortedTemperatures[i] / _referenceTemperatures[iBase],
                    _gasPremultiplier / (-1.0 * localLapseRate));
            }
            else
            {
                pressure = pressureBase * Math.Exp(_gasPremultiplier * (zBase_m - zCeiling_m) / temperatureBase);
            }
            // Put the value back in the output
            outputPressures[sortingIndices[i]] = pressure;
        }
        return outputPressures;
    }

    public static double AltitudeToPressure(double altitude)
    {
        // Input: altitude in km
        // Output: pressure in Pa
        return AltitudeToPressure([altitude])[0];
    }
    
    public static double AltitudeToTemperature(double altitude)
    {
        // Input: altitude in km
        // Output: pressure in Pa
        return AltitudeToTemperature([altitude])[0];
    }
}