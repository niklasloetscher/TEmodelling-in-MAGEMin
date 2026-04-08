module LSMUtils

# using ..LatticeStrainModels
# using ..PartitionCoefficients

export Eu_ratio_Burnham, Kd_Eu, optical_basicity

# Eu-ratio
function Eu_ratio_Burnham(logfO2::Float64, T::Float64, Λ::Float64)
    """
        Burnham et al. 2015

        Inputs:
            - logfO2... Oxygen fugacity in base 10 log units
            - T... Temperature in °C
            - Λ... Optical basicity (Duffy, 1993)
    """
    burnham_ratio = 1/(1+10^(-0.25*logfO2 - 6410/(T+273.15) -14.2Λ + 10.1)) #Eu3+/∑Eu
    ratio = (1-burnham_ratio)/burnham_ratio # Convert to Eu2+/Eu3+
    return ratio
end

# Eu Kd
function Kd_Eu(Kd_Eu2::Float64, Kd_Eu3::Float64, ratio::Float64)
    """
        Inputs:
            - Kd_Eu2... Mineral-melt partition coefficient for Eu 2+
            - Kd_Eu3... Mineral-melt partition coefficient for Eu 3+
            - ratio... Eu2+/Eu3+
    """
    return (Kd_Eu2 * ratio + Kd_Eu3) / (ratio + 1)
end

# Optical basicity 
function optical_basicity(bulk)
    """
        Input:
            - bulk... bulk melt molar composition
    """
    O_units = [2, 3, 1, 1, 1, 1, 1, 2, 1, 3, 1]
    # "SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"
    lc98 = [0.48, 0.6, 1.00, 0.78, 1.00, 1.40, 1.15, 0.75, 0.0, 0.70, 0.40]

    Λ = sum(bulk.*O_units.*lc98)/sum(bulk.*O_units)

    return Λ
end




end # Module