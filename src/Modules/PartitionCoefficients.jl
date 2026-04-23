module PartitionCoefficients

using ..LatticeStrainModels
using ..LSMUtils

export get_D

using CSV, DataFrames

const R = 8.314  # J/mol/K

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
	"pat"	 => "mu"
)

const LSM_PHASES = Set{String}(["ol", "cpx", "fsp", "amp", "spl", "opx", "g", "ap", "zr"])

canonical(phase::String) = get(PHASE_ALIASES, phase, phase)

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
    ("opx", 3) => "bed25",
    ("amp", 3) => "shi17melt",
    ("g", 3) => "mk20",
    ("spl", 3) => "sie20",
    ("ol", 3) => "bed05",
	("ap", 3) => "jir25",
	("zr", 3) => "str23"	
)

function _get_ri(element::Symbol, charge::Int, coordination::Int)
    key = (string(element), charge, coordination)
    haskey(_SHANNON, key) || error("No ionic radius for $element, charge=$charge, coordination=$coordination")
    return _SHANNON[key]
end

function set_lsm_authors!(mineral::String, charge::Int, authors::String)
    LSM_AUTHORS[(mineral, charge)] = authors
end


function get_D(
    mineral::String,
    elements::Vector{Symbol},
    charges::Vector{Int},
    params::Dict;
    authors_map::Union{Dict{Tuple{String,Int},String}, Nothing} = nothing,
    A::Float64 = 28.0
) :: Dict{Symbol, Float64}
    cmineral = canonical(mineral)

    # Split Eu from the rest
    eu_idx    = findall(==(:Eu), elements)
    other_idx = findall(!=(:Eu), elements)

    other_elements = elements[other_idx]
    other_charges  = charges[other_idx]

    result = Dict{Symbol, Float64}()

    # Standard path for all non-Eu elements
    if !isempty(other_elements)
        if cmineral in LSM_PHASES
            merge!(result, _D_lsm(cmineral, other_elements, other_charges, params, authors_map, A))
        else
            merge!(result, Dict(el => _D_fixed(cmineral, el) for el in other_elements))
        end
    end

    # Eu path
    if !isempty(eu_idx)
        result[:Eu] = _D_Eu(cmineral, params, authors_map, A)
    end

    return result
end


function _D_Eu(mineral, params, authors_map, A)
    # Eu3+ — standard LSM or fixed D path
    Kd_Eu3 = if mineral in LSM_PHASES
        _D_lsm(mineral, [:Eu], [3], params, authors_map, A)[:Eu]
    else
        _D_fixed(mineral, :Eu)
    end

    # Eu2+ — LSM or fixed D with charge=2
    Kd_Eu2 = if mineral in LSM_PHASES
        _D_lsm(mineral, [:Eu], [2], params, authors_map, A)[:Eu]
    else
        _D_fixed(mineral, :Eu)   # no charge distinction for fixed D
    end

    # Eu2+/Eu3+; requires :logfO2, :T, :Λ in params
    ratio = LSMUtils.Eu_ratio_Burnham(
        Float64(params[:logfO2]),
        Float64(params[:T]),
        Float64(params[:Λ])
    )

    return LSMUtils.Kd_Eu(Kd_Eu2, Kd_Eu3, ratio)
end

function _D_lsm(mineral, elements, charges, params, authors_map, A)
    T = Float64(params[:T])
    length(elements) == length(charges) || error("elements and charges must have the same length")

    # Resolve and run one LSM per requested charge, then recombine.
    result = Dict{Symbol, Float64}()
    for charge in unique(charges)
        idx = findall(==(charge), charges)
        group_elements = elements[idx]
        group_charges  = charges[idx]

        authors, lsm_charge = _resolve_authors_and_charge(mineral, charge, authors_map)
        lsm = LatticeStrainModels.get_lsm(mineral=mineral, charge=lsm_charge, authors=authors)

        r0 = lsm.r0 isa Function ? lsm.r0(params) : lsm.r0
        E  = lsm.E  isa Function ? lsm.E(params)  : lsm.E
        D0 = lsm.D0 isa Function ? lsm.D0(params) : lsm.D0

        merge!(result, _lattice_strain_D(r0, E, D0, group_elements, group_charges, lsm_charge, lsm.coordination, T; A=A))
    end

    return result
end


function _resolve_authors_and_charge(mineral, target_charge, authors_map)
    # Prefer exact (mineral, charge); otherwise use closest available charge.
    target_charge = Int(target_charge)

    map = authors_map !== nothing ? authors_map : LSM_AUTHORS
    key = (mineral, target_charge)

    if haskey(map, key)
        return map[key], target_charge
    else
        # Find closest available charge for this mineral
        available = [(k[2], v) for (k, v) in map if k[1] == mineral]
        isempty(available) && error("No LSM author entry for mineral: $mineral")
        closest_charge, closest_authors = argmin(x -> abs(x[1] - target_charge), available)
        #@warn "No LSM for ($mineral, charge=$target_charge) — using charge=$closest_charge"
        return closest_authors, closest_charge
    end
end

function _lattice_strain_D(
    r0::Float64, E::Float64, D0::Float64,
    elements::Vector{Symbol}, charges::Vector{Int},
    lsm_charge::Int, coordination::Int,
    T::Float64;
    A::Float64 = 28.0
)  # default A — adjust or pass in
							
	"""
		r0... Ideal cation radius for this site [°A]
		E... Effective Young's modulus [GPa]
		D0... Zero-strain partition coefficient
		elements... Elements for which Ds are to be returned
		charges... Charges of elements
		lsm_charge... Charge for this phase
		T... Temeperature [°C]
		
		Opt:
		A... Charge factor [kJ]
	"""

	length(elements) == length(charges) || error("elements and charges must have the same length")

    # Convert to SI
    r0_si = r0 * 1e-10
    E_si  = E  * 1e9
    T_K   = T  + 273.15
    A_J   = A  * 1e3
    NA    = 6.02214e23

    result = Dict{Symbol, Float64}()
    for (el, z) in zip(elements, charges)
        ri      = _get_ri(el, z, coordination)
        ri_si   = ri * 1e-10
        dcharge = z - lsm_charge

        dr     = ri_si - r0_si
        strain = -4π * E_si * NA * (r0_si/2 * dr^2 + dr^3/3)
        D      = D0 * exp(strain / (R * T_K))

        if dcharge != 0
            D *= exp(-A_J * abs(dcharge) / (R * T_K))
        end

        result[el] = D
    end
    return result
end

function _D_fixed(mineral, element)
    haskey(FIXED_D, mineral) || error("No fixed-D data for mineral: $mineral")
    d = FIXED_D[mineral]
    haskey(d, element) || error("No fixed D for $element in $mineral")
    return d[element]
end

lsm_phases() = unique(first.(keys(LatticeStrainModels._REGISTRY)))

end # module
