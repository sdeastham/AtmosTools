namespace AtmosTools;

public static class IdealGases
{
    public static readonly double GammaAir = 1.4;
    public static readonly double CpAir = 1005.0; // J/kg/K
    public static readonly double MolarMassAir = 28.97e-3; // kg/mol

    // Precomputed factor for convenience
    private static readonly double MachFactor = GammaAir * PhysConstants.RGasUniversal / MolarMassAir; 
    public static double SpeedOfSound(double temperature)
    {
        return Math.Sqrt(temperature * MachFactor);
    }
}