module PartitionCoefficientsV2

using ..LatticeStrainModels
using ..LSMUtils

export get_D, set_lsm_authors!

using CSV, DataFrames

# Gas constant
const R = 8.3144626  # J/mol/K

# Ionic radii from Shannon(1976) in °A
const _SHANNON = let
    path = joinpath(@__DIR__, "..", "Data", "Ionic radii", "shannon.csv")
    df = CSV.read(path, DataFrame)
    # Build lookup: (element, charge, coordination) => ionicradius
    Dict{Tuple{String,Int,Int}, Float64}(
        (row.element, row.charge, row.coordination) => row.ionicradius
        for row in eachrow(df)
        if !ismissing(row.ionicradius)
    )
end

# Phase aliases for some phases in MAGEMin based on composition
const PHASE_ALIASES = Dict{String,String}(
    "Na-cpx" => "cpx",
    "pig"    => "cpx",
    "pl"     => "fsp",
    "afs"    => "fsp",
    "hem"    => "ilm",
    "gl"     => "amp",
    "act"    => "amp",
    "cumm"   => "amp",
    "tr"     => "amp",
    "cm"     => "spl",
    "usp"    => "spl",
    "mgt"    => "spl",
	"pat"	 => "mu",
	"sill"   => "ky",
	"and"    => "ky"
)

# Phases for which LSMs are implemented in LatticeStrainModels.jl
const LSM_PHASES = Set{String}(model.mineral for model in values(LatticeStrainModels._REGISTRY))

# LSM_PHASES are only implemented for canonical phase names, so I need a function tha maps phase name to canonical phase name based on PHASE_ALIASES
# This function needs to map e.g. "Na-cpx" to "cpx", but also handle the case where the phase is simply named "cpx" in whic case it also maps to "cpx"
canonical(phase::String) = get(PHASE_ALIASES, phase, phase)

function _get_ri(element::Symbol, charge::Int, coordination::Int)
    key = (String(element), charge, coordination)
    haskey(_SHANNON, key) || error("No ionic radius for element=$(element), charge=$(charge), coordination=$(coordination)")
    return _SHANNON[key]
end

# Fixed partition coefficients for non-LSM minerals
# Dict: mineral => (element => D)
const FIXED_D = Dict{String, Dict{Symbol, Float64}}(
    "ilm" => Dict(
        :La  => 0.002,
        :Ce  => 0.003,
        :Pr  => 0.0035,
        :Nd  => 0.004,
        :Sm  => 0.011,
		:Eu  => 0.010, # Table 1 Bédard 2006
        :Gd  => 0.007,
        :Tb  => 0.013,
        :Dy  => 0.019,
        :Y   => 0.021,
        :Ho  => 0.029,
        :Er  => 0.036,
        :Tm  => 0.054,
        :Yb  => 0.071,
        :Lu  => 0.112,
        :Sr  => 0.004,
    ), # Supplement Table S4 Soderman & Weller 2026
    "bi" => Dict(
        :La  => 0.61,
        :Ce  => 0.32,
        :Pr  => 0.29,
        :Nd  => 0.18,
        :Sm  => 0.05,
		:Eu  => 0.031, # Table 1 Bédard 2006
        :Gd  => 0.04,
        :Tb  => 0.06,
        :Dy  => 0.29,
        :Y   => 1.00,
        :Ho  => 0.07,
        :Er  => 0.10,
        :Tm  => 0.29,
        :Yb  => 0.30,
        :Lu  => 0.33,
        :Sr  => 0.60,
    ), # Supplement Table S4 Soderman & Weller 2026
	"mu" => Dict(
        :La  => 0.61,
        :Ce  => 0.32,
        :Pr  => 0.29,
        :Nd  => 0.18,
        :Sm  => 0.05,
		:Eu  => 0.031,
        :Gd  => 0.04,
        :Tb  => 0.06,
        :Dy  => 0.29,
        :Y   => 1.00,
        :Ho  => 0.07,
        :Er  => 0.10,
        :Tm  => 0.29,
        :Yb  => 0.30,
        :Lu  => 0.33,
        :Sr  => 0.60,
    ), # THIS IS JUST A COPY OF "bi"!
	"chl" => Dict(
        :La  => 0.61,
        :Ce  => 0.32,
        :Pr  => 0.29,
        :Nd  => 0.18,
        :Sm  => 0.05,
		:Eu  => 0.031,
        :Gd  => 0.04,
        :Tb  => 0.06,
        :Dy  => 0.29,
        :Y   => 1.00,
        :Ho  => 0.07,
        :Er  => 0.10,
        :Tm  => 0.29,
        :Yb  => 0.30,
        :Lu  => 0.33,
        :Sr  => 0.60,
    ), # THIS IS JUST A COPY OF "bi"!
    "fl" => Dict(
        :La  => 1.84,
        :Ce  => 2.11,
        :Pr  => 2.14,
        :Nd  => 2.16,
        :Sm  => 1.68,
		:Eu  => 1.57, # 0.5*(Sm+Gd) 
        :Gd  => 1.46,
        :Tb  => 1.42,
        :Dy  => 1.26,
        :Y   => 1.11,
        :Ho  => 0.95,
        :Er  => 0.88,
        :Tm  => 0.81,
        :Yb  => 0.74,
        :Lu  => 0.54,
        :Sr  => 0.24,
    ), # Supplement Table S4 Soderman & Weller 2026
    "ru" => Dict(
        :La  => 0.0001,
        :Ce  => 0.0001,
        :Pr  => 0.0001,
        :Nd  => 1e-5,
        :Sm  => 0.0018,
		:Eu  => 0.00037, # Table 1 Bédard 2006
        :Gd  => 0.0007,
        :Tb  => 1e-5,
        :Dy  => 0.010,
        :Y   => 0.007,
        :Ho  => 1e-5,
        :Er  => 0.012,
        :Tm  => 1e-5,
        :Yb  => 0.016,
        :Lu  => 0.018,
        :Sr  => 0.003,
    ), # Supplement Table S4 Soderman & Weller 2026
	"sph" => Dict(
        :La  => 4.73,
        :Ce  => 7.57,
        :Pr  => 9.0,
        :Nd  => 12.4,
        :Sm  => 14.0,
		:Eu  => 13.8,
        :Gd  => 11.9,
        :Tb  => 10.0,
        :Dy  => 8.27,
        :Y   => 5.42,
        :Ho  => 5.5,
        :Er  => 5.54,
        :Tm  => 4.00,
        :Yb  => 3.02,
        :Lu  => 2.00,
        :Sr  => 2.68,
    ), # Table 1 Bédard 2006
    "q" => Dict(
        :La  => 0.016,
        :Ce  => 0.014,
        :Pr  => 1e-5,
        :Nd  => 0.016,
        :Sm  => 0.014,
		:Eu => 1e-5, # Arbitrary
        :Gd  => 1e-5,
        :Tb  => 0.017,
        :Dy  => 0.015,
        :Y   => 1e-5,
        :Ho  => 1e-5,
        :Er  => 1e-5,
        :Tm  => 1e-5,
        :Yb  => 0.017,
        :Lu  => 0.012,
        :Sr  => 1e-5,
    ), # Supplement Table S4 Soderman & Weller 2026
	"cd" => Dict(
        :La  => 0.016,
        :Ce  => 0.014,
        :Pr  => 1e-5,
        :Nd  => 0.016,
        :Sm  => 0.014,
		:Eu => 1e-5, # Arbitrary
        :Gd  => 1e-5,
        :Tb  => 0.017,
        :Dy  => 0.015,
        :Y   => 1e-5,
        :Ho  => 1e-5,
        :Er  => 1e-5,
        :Tm  => 1e-5,
        :Yb  => 0.017,
        :Lu  => 0.012,
        :Sr  => 1e-5,
    ), # THIS IS JUST A COPY OF "q"!
	"ne" => Dict(
        :La  => 0.016,
        :Ce  => 0.014,
        :Pr  => 1e-5,
        :Nd  => 0.016,
        :Sm  => 0.014,
		:Eu => 1e-5, # Arbitrary
        :Gd  => 1e-5,
        :Tb  => 0.017,
        :Dy  => 0.015,
        :Y   => 1e-5,
        :Ho  => 1e-5,
        :Er  => 1e-5,
        :Tm  => 1e-5,
        :Yb  => 0.017,
        :Lu  => 0.012,
        :Sr  => 1e-5,
    ), # THIS IS JUST A COPY OF "q"!
	"ky" => Dict(
        :La  => 0.016,
        :Ce  => 0.014,
        :Pr  => 1e-5,
        :Nd  => 0.016,
        :Sm  => 0.014,
		:Eu => 1e-5, # Arbitrary
        :Gd  => 1e-5,
        :Tb  => 0.017,
        :Dy  => 0.015,
        :Y   => 1e-5,
        :Ho  => 1e-5,
        :Er  => 1e-5,
        :Tm  => 1e-5,
        :Yb  => 0.017,
        :Lu  => 0.012,
        :Sr  => 1e-5,
    ), # THIS IS JUST A COPY OF "q"!
		"ky" => Dict(
        :La  => 0.016,
        :Ce  => 0.014,
        :Pr  => 1e-5,
        :Nd  => 0.016,
        :Sm  => 0.014,
		:Eu => 1e-5, # Arbitrary
        :Gd  => 1e-5,
        :Tb  => 0.017,
        :Dy  => 0.015,
        :Y   => 1e-5,
        :Ho  => 1e-5,
        :Er  => 1e-5,
        :Tm  => 1e-5,
        :Yb  => 0.017,
        :Lu  => 0.012,
        :Sr  => 1e-5,
    ), # THIS IS JUST A COPY OF "q"!
	"ep" => Dict(
        :La  => 2.05,
        :Ce  => 2.44,
        :Pr  => 2.86,
        :Nd  => 3.34,
        :Sm  => 4.22,
		:Eu => 3.78,
        :Gd  => 4.67,
        :Tb  => 4.67,
        :Dy  => 4.50,
        :Y   => 4.30,
        :Ho  => 4.18,
        :Er  => 3.78,
        :Tm  => 3.36,
        :Yb  => 2.96,
        :Lu  => 2.59,
        :Sr  => 2.0,
    ) # Table 1 Bédard 2006
)


# Default author assignments per (mineral, charge)
# Modify via set_lsm_authors!()
const LSM_AUTHORS = Dict{Tuple{String,Int}, String}(
    ("fsp", 3) => "sun17_3",
    ("fsp", 2) => "sun17_2",
    ("cpx", 3) => "sl12",
	("cpx", 2) => "wb14cpx", # Based on "sl12"!
    ("opx", 3) => "bed25",
	("opx", 2) => "wb14opx", # Based on "bed25"!
    ("amp", 3) => "shi17melt",
    ("g", 3) => "mk20",
	("g", 2) => "wb14g", # Based on "mk20"!
    ("spl", 3) => "sie20",
    ("ol", 3) => "sl13ol",
	("ol", 2) => "wb14ol", # Based on "bed05"!
	("ap", 3) => "jir25",
	("zr", 3) => "str23"	
)

# Function to set alternative author for a (mineral, charge) pair
function set_lsm_authors!(mineral::String, charge::Int, authors::String)
    """
    Input:
        mineral: Mineral name (e.g. "cpx")
        charge: Charge state (e.g. 3 for trivalent)
        authors: Author string corresponding to a model in LatticeStrainModels._REGISTRY (e.g. "bed14M1")
    """
    LSM_AUTHORS[(mineral, charge)] = authors
end


# Function to get partition coeffficients for given mineral, list of elements and their cation charge.
function get_D(mineral::String, elements::Vector{Symbol}, charges::Vector{Int}, params::Dict{Symbol, Float64}) :: Dict{Symbol, Float64}
    """
    Input:
        mineral: Mineral name (e.g. 'cpx')
        elements: Vector of element symbols (e.g. [:La, :Ce, :Nd])
        charges: Vector of cation charges corresponding to the elements (e.g. [3, 3, 3])
        params: Dict of required parameters for the LSM model (e.g. T, logfO2, Λ, compositional, etc.)

    Output: Dict{Symbol, Float64} mapping element symbol to partition coefficient D for the given mineral.

    """

    # Internal constant to adjust for charge mismatch after Wood & Blundy (2014)
    A = 28.0e3

    # Canonical mineral name
    mineral = canonical(mineral)
    # Ensure elements are Symbols (accept string or symbol input)
    # Check if canonical(mineral) is in LSM_PHASES
    hasLSM = mineral in LSM_PHASES

    # If there is no LSM, just return the fixed D values from FIXED_D for canonical(mineral). Throw an error if either the mineral is not implemented in FIXED_D or if any of the elements are not implemented for that mineral.
    if !hasLSM
        haskey(FIXED_D, mineral) || error("No fixed D values for mineral=$(canonical(mineral)). Available: $(keys(FIXED_D))")
        D_dict = FIXED_D[mineral]
        # Collect results but provide diagnostics if a key is missing (helps catch Symbol vs String mismatches)
        if all(haskey(D_dict, el) for el in elements)
            return Dict(el => D_dict[el] for el in elements)
        else
            # Diagnostic output to help debug key/type mismatches
            println("[PartitionCoefficientsV2 DEBUG] mineral=", mineral)
            println("[PartitionCoefficientsV2 DEBUG] FIXED_D keys: ", collect(keys(D_dict)))
            println("[PartitionCoefficientsV2 DEBUG] elements types: ", typeof.(elements))
            println("[PartitionCoefficientsV2 DEBUG] elements values: ", elements)
            # Raise the original informative error for the first missing element
            for el in elements
                if !haskey(D_dict, el)
                    error("No fixed D value for element=$(el) in mineral=$(canonical(mineral)). Available: $(keys(D_dict))")
                end
            end
        end
        #return Dict(Symbol(el) => get(D_dict, Symbol(el), error("No fixed D value for element=$(Symbol(el)) in mineral=$(canonical(mineral)). Available: $(keys(D_dict))")) for el in elements)
    # If there is an LSM, start by reserving a D_dict of length elements
    else
        D_dict = Dict{Symbol, Float64}()
        # Iterate over elements
        for (i, el) in enumerate(elements)
            # Get the charge of the element
            charge = charges[i]
            dcharge = 0
            closest_charge = charge
            # Get the LSM_AUTHORS for this mineral and charge, if there is none for this particular charge, then get the one which has the closest charge state.

            # Case 1) There is an exact match for (mineral, charge) in LSM_AUTHORS
            if haskey(LSM_AUTHORS, (mineral, charge))
                authors = LSM_AUTHORS[(mineral, charge)]
            else # Case 2) There is no exact match, so we need to find the closest charge state for this mineral in LSM_AUTHORS
                available_charges = [c for (m, c) in keys(LSM_AUTHORS) if m == mineral]
                if isempty(available_charges)
                    error("No LSM authors found for mineral=$(mineral). Available: $(keys(LSM_AUTHORS))")
                end
                closest_charge = available_charges[argmin(abs.(available_charges .- charge))]
                dcharge = abs(closest_charge - charge)
                authors = LSM_AUTHORS[(mineral, closest_charge)]
                #@warn "No LSM authors found for mineral=$(mineral) and charge=$(charge). Using authors=$(authors) for closest charge=$(closest_charge) instead."
            end

            # Eu needs to be handled specially because it is in a mixed valence state
            isEu = el == :Eu

            if isEu
                # Compute Eu ratio 
                eu_ratio = LSMUtils.Eu_ratio_Burnham(params[:logfO2], params[:T], params[:Λ])
                # Compute both D for Eu2+ and Eu3+; :Eu is given with charge 3
                D_Eu3 = _LSM_D(mineral, el, 3, authors, params)
                # Need to check whether there is a model for 2+ again here
                if haskey(LSM_AUTHORS, (mineral, 2))
                    D_Eu2 = _LSM_D(mineral, el, 2, LSM_AUTHORS[(mineral, 2)], params)
                else
                    # Find model with closest charge state to 2
                    available_charges = [c for (m, c) in keys(LSM_AUTHORS) if m == mineral]
                    if isempty(available_charges)
                        error("No LSM authors found for mineral=$(mineral). Available: $(keys(LSM_AUTHORS))")
                    end
                    closest_chargeEu = available_charges[argmin(abs.(available_charges .- 2))]
                    dchargeEu = abs(closest_chargeEu - 2)
                    authorsEu = LSM_AUTHORS[(mineral, closest_chargeEu)]
                    D_Eu2 = _LSM_D(mineral, el, 2, authorsEu, params) * exp(-A * dchargeEu^2 / (R * params[:T])) # Adjust for charge mismatch using Wood & Blundy (2014)
                end
                # Compute Eu D
                D_dict[el] = LSMUtils.D_Eu(D_Eu2, D_Eu3, eu_ratio)
            else
                # Compute _LSM_D(mineral, charge, authors, el) and add it to D_dict
                if dcharge == 0
                    D_dict[el] = _LSM_D(mineral, el, charge, authors, params)
                else
                    D_dict[el] = _LSM_D(mineral, el, charge, authors, params) * exp(-A * dcharge^2 / (R * params[:T])) # Adjust for charge mismatch using Wood & Blundy (2014)
                end
            end

        end
    end

    return D_dict

end

function _LSM_D(mineral::String, element::Symbol, elcharge::Int, authors::String, params::Dict{Symbol, Float64}) :: Float64
    """
        Input:
            mineral: Mineral name (e.g. "cpx")
            element: Element symbol (e.g. :La)
            elcharge: Element charge state for getting ri, can be different from the charge for which LSM of authors is defined (e.g. 3 for trivalent)
            authors: Author string corresponding to a model in LatticeStrainModels._REGISTRY (e.g. "bed14M1")
            params: Dict of required parameters for the LSM model (e.g. T, logfO2, Λ, compositional, etc.)
    
        Output: Partition coefficient D for the given mineral, element and charge state computed using the Lattice Strain Model with the parameters from LatticeStrainModels._REGISTRY corresponding to (mineral, charge, authors).

    """

    NA = 6.02214076e23 # Avogadro's number

    # Multisite model: weighted average of site-specific D values
    if LatticeStrainModels.has_multisite_model(mineral, elcharge, authors)
        variants = LatticeStrainModels.get_multisite_variants(mineral, elcharge, authors)
        variants === nothing && error("No multisite variants found for mineral=$(mineral), charge=$(elcharge), authors=$(authors)")

        total_weight = sum(weight for (_, weight) in variants)
        total_weight > 0.0 || error("Invalid multisite weights for mineral=$(mineral), charge=$(elcharge), authors=$(authors): sum(weights) must be > 0")

        return sum(weight * _LSM_D(mineral, element, elcharge, variant_authors, params) for (variant_authors, weight) in variants) / total_weight
    end

    # Single-site model

    # Blundy & Wood model
    # D = D0 * exp(-σ/(R*T))
    # \sigma = 4*π*E*NA*(r0/2*(r0-ri)^2 - 1/3*(r0-ri)^3))

    # For ri, need to get the ionic radius for this element with given charge and coordination number

    lsm = get_lsm(authors=authors) # The model charge is defined by the registered model for authors; elcharge is only used for the ionic radius.

    r0 = lsm.r0 isa Function ? lsm.r0(params) : lsm.r0
    E  = lsm.E  isa Function ? lsm.E(params)  : lsm.E
    D0 = lsm.D0 isa Function ? lsm.D0(params) : lsm.D0

    ri = _get_ri(element, elcharge, lsm.coordination)

    # Convert to SI units
    r0_si = r0 * 1e-10 # Angstrom to meters
    ri_si = ri * 1e-10 # Angstrom to meters
    E_si  = E * 1e9 # GPa to Pa
    T_K   = params[:T] + 273.15 # Celsius to Kelvin

    dr = r0_si - ri_si
    strain = -4π * E_si * NA * (r0_si / 2 * dr^2 - dr^3 / 3)
    return D0 * exp(strain / (R * T_K))

    
end

end # module