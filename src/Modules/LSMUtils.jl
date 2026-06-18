module LSMUtils

# using ..LatticeStrainModels
# using ..PartitionCoefficients

export Eu_ratio_Burnham, Kd_Eu, optical_basicity, get_params_for_phase, calc_comp_fo2

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
function D_Eu(Kd_Eu2::Float64, Kd_Eu3::Float64, ratio::Float64)
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


function _oxide_idx(oxides::Vector{String}, name::String)
    return findfirst(==(name), oxides)
end

function _phase_index(phase::String, phases::Vector{String}, solidx)
    sym = Symbol(phase)
    if haskey(solidx, sym)
        return Int(solidx[sym])
    end
    idx = findfirst(==(phase), phases)
    return idx === nothing ? nothing : Int(idx)
end

function _safe_comp_wt(solvecs, idx::Int, oxides::Vector{String}, oxide::String)
    oid = _oxide_idx(oxides, oxide)
    oid === nothing && return 0.0
    return Float64(solvecs[idx].Comp_wt[oid])
end

function get_params_for_phase(
    phase::String,
    out,
    out_sat,
    T::Float64,
    P::Float64;
    phases::Union{Nothing, Vector{String}}=nothing,
    melt_comp_wt_override=nothing,
    bulk_override=nothing
)::Dict{Symbol, Float64}
    phases_local = phases === nothing ? out.ph : phases
    solvecs = out.SS_vec
    solidx = out.SS_syms
    oxides = out.oxides
    bulk = bulk_override === nothing ? Vector{Float64}(out.bulk_M) : Vector{Float64}(bulk_override)

    melt_comp_wt = if melt_comp_wt_override === nothing
        mc = Vector{Float64}(out.bulk_M_wt)
        s = sum(mc)
        s > 0.0 ? mc ./ s : mc
    else
        Vector{Float64}(melt_comp_wt_override)
    end

    logfO2 = Float64(out.fO2)

    params = Dict{Symbol, Float64}()

    params[:T] = T
    params[:P] = P
    params[:logfO2] = logfO2
    params[:Λ] = optical_basicity(bulk)

    # Defaults for all phase-specific fields expected by LSM models.
    params[:XAn] = 0.0
    params[:fspXNa] = 0.0
    params[:fspXCa] = 0.0
    params[:cpxXMgM2] = 0.0
    params[:cpxXAlM1] = 0.0
    params[:cpxXAlT] = 0.0
    params[:cpxCaO] = 0.0
    params[:opxMgO] = 0.0
    params[:ampXFmM4] = 0.0
    params[:gtFeO] = 0.0
    params[:gtMgO] = 0.0
    params[:olMgO] = 0.0
	params[:olFo] = 0.0

    # Melt/system fields (shared across phases at a given step).
    iSiO2 = _oxide_idx(oxides, "SiO2")
    iTiO2 = _oxide_idx(oxides, "TiO2")
    iCaO  = _oxide_idx(oxides, "CaO")
    iMgO  = _oxide_idx(oxides, "MgO")
    iFeO  = _oxide_idx(oxides, "FeO")
    iH2O  = _oxide_idx(oxides, "H2O")
    iO    = _oxide_idx(oxides, "O")

    meltXSi  = iSiO2 === nothing ? 0.0 : Float64(bulk[iSiO2])
    meltXTi  = iTiO2 === nothing ? 0.0 : Float64(bulk[iTiO2])
    meltXCa  = iCaO  === nothing ? 0.0 : Float64(bulk[iCaO])
    meltXMg  = iMgO  === nothing ? 0.0 : Float64(bulk[iMgO])
    meltXFe  = iFeO  === nothing ? 0.0 : Float64(bulk[iFeO])
    meltXH2O = iH2O  === nothing ? 0.0 : Float64(bulk[iH2O])
    meltXO   = iO    === nothing ? 0.0 : Float64(bulk[iO])

    meltP2O5 = 0.0
    if hasproperty(out_sat, :Cliq)
        cliq = getproperty(out_sat, :Cliq)
        if length(cliq) >= 2
            meltP2O5 = Float64(cliq[2]) * 1e-4
        end
    end

    meltMgN = (meltXMg + meltXFe) > 0.0 ? meltXMg / (meltXMg + meltXFe) : 0.0
    meltχH2O = sum(bulk) > 0.0 ? meltXH2O / sum(bulk) : 0.0
    meltH2O = isempty(melt_comp_wt) ? 0.0 : Float64(last(melt_comp_wt))

    meltMgO = iMgO === nothing ? 0.0 : Float64(melt_comp_wt[iMgO])
    meltCaO = iCaO === nothing ? 0.0 : Float64(melt_comp_wt[iCaO])
    meltFeO = iFeO === nothing ? 0.0 : Float64(melt_comp_wt[iFeO])
    meltSiO2 = iSiO2 === nothing ? 0.0 : Float64(melt_comp_wt[iSiO2])

    anh_norm = iH2O === nothing ? sum(bulk) : sum(bulk) - Float64(bulk[iH2O])
    anh_norm_wt = iH2O === nothing ? sum(melt_comp_wt) : sum(melt_comp_wt) - Float64(melt_comp_wt[iH2O])
    anh_norm = anh_norm == 0.0 ? 1.0 : anh_norm
    anh_norm_wt = anh_norm_wt == 0.0 ? 1.0 : anh_norm_wt

    params[:meltXSi] = meltXSi / anh_norm
    params[:meltXTi] = meltXTi / anh_norm
    params[:meltXCa] = meltXCa / anh_norm
    params[:meltXMg] = meltXMg / anh_norm
    params[:meltXO] = meltXO / anh_norm
    params[:meltMgNhydr] = meltMgN
    params[:meltXH2O] = meltXH2O
    params[:meltP2O5] = meltP2O5
    params[:meltχH2O] = meltχH2O
    params[:meltH2O] = meltH2O
    params[:meltMgO] = meltMgO / anh_norm_wt
    params[:meltCaO] = meltCaO / anh_norm_wt
    params[:meltFeOhydr] = meltFeO
    params[:meltMgOhydr] = meltMgO
    params[:meltSiO2hydr] = meltSiO2

    idx = _phase_index(phase, phases_local, solidx)
    idx === nothing && return params
	
	
	params[:meltXCatAl] = 0.0
	if "liq" in phases_local
		meltidx = _phase_index("liq", phases_local, solidx)
		params[:meltXCatAl] = solvecs[meltidx].Comp_apfu[2]
	end

    # Feldspar aliases: compute from the actual phase instance currently requested.
    if phase in ("pl", "afs", "fsp")
        XCaA = Float64(solvecs[idx].siteFractions[2])
        XNaA = Float64(solvecs[idx].siteFractions[1])
        denom = XCaA + XNaA
        params[:XAn] = denom > 0.0 ? XCaA / denom : 0.0
        params[:fspXNa] = Float64(solvecs[idx].Comp_apfu[7])
        params[:fspXCa] = Float64(solvecs[idx].Comp_apfu[3])
    end

    # Clinopyroxene aliases.
    if phase in ("cpx", "pig", "Na-cpx")
        params[:cpxXMgM2] = Float64(solvecs[idx].siteFractions[7])
        params[:cpxXAlM1] = Float64(solvecs[idx].siteFractions[3])
        params[:cpxXAlT] = Float64(solvecs[idx].siteFractions[13])
        cpxH2O = _safe_comp_wt(solvecs, idx, oxides, "H2O")
        cpxCaO = _safe_comp_wt(solvecs, idx, oxides, "CaO")
        params[:cpxCaO] = cpxCaO * (sum(solvecs[idx].Comp_wt) - cpxH2O)
    end

    if phase == "opx"
        opxH2O = _safe_comp_wt(solvecs, idx, oxides, "H2O")
        opxMgO = _safe_comp_wt(solvecs, idx, oxides, "MgO")
        params[:opxMgO] = opxMgO * (sum(solvecs[idx].Comp_wt) - opxH2O)
    end

    if phase in ("amp", "gl", "act", "cumm", "tr")
        XFeM4 = Float64(solvecs[idx].siteFractions[13])
        XMgM4 = Float64(solvecs[idx].siteFractions[12])
        params[:ampXFmM4] = XFeM4 + XMgM4
    end

    if phase == "g"
        gtH2O = _safe_comp_wt(solvecs, idx, oxides, "H2O")
        gtFeO = _safe_comp_wt(solvecs, idx, oxides, "FeO")
        gtMgO = _safe_comp_wt(solvecs, idx, oxides, "MgO")
        params[:gtFeO] = gtFeO * (sum(solvecs[idx].Comp_wt) - gtH2O)
        params[:gtMgO] = gtMgO * (sum(solvecs[idx].Comp_wt) - gtH2O)
    end

    if phase == "ol"
        olH2O = _safe_comp_wt(solvecs, idx, oxides, "H2O")
        olMgO = _safe_comp_wt(solvecs, idx, oxides, "MgO")
		olXMgO = Float64(solvecs[idx].Comp[4])
		olXFeO = Float64(solvecs[idx].Comp[5])
        params[:olMgO] = olMgO * (sum(solvecs[idx].Comp_wt) - olH2O)
		params[:olFo] = 100.0*olXFeO/(olXMgO+olXFeO)
    end

    return params
end


# Getting Fe2O3, FeO speciation after Kress & Carmichael (1991)
function calc_comp_fo2(bulkin, t, p, delta)
    # Calculates new composition with the fo2 (and thus Fe valence) taken into account for a new temp or pressure
    # T in C, P in MPa, bulkin as array of wt% corresponding to:
    # [SiO2, Al2O3, CaO, MgO, FeO, K2O, Na2O, TiO2, O, Cr2O3, H2O]
    # Returns new_comp as array of wt% corresponding to:
    # [SiO2, Al2O3, CaO, MgO, FeO, Fe2O3, K2O, Na2O, TiO2, Cr2O3, H2O]

    t = t + 273.15

    logfo2FMQ = 8.58 - 25050 / t
    logfo2 = delta + logfo2FMQ
    fo2 = 10^logfo2

    mm = Dict("SiO2" => 60.08, "TiO2" => 79.866, "Al2O3" => 101.96, "Fe2O3" => 159.69,
              "FeO" => 71.844, "MnO" => 70.9374, "MgO" => 40.3044, "CaO" => 56.0774,
              "Na2O" => 61.9789, "K2O" => 94.2, "P2O5" => 283.89, "H2O" => 18.01528,
              "CO2" => 44.01, "S" => 32.065, "Cl" => 35.453, "F" => 18.9984)

    p = p * 1e6
    param = [0.196, 11492.0, -6.675, -2.243, -1.828, 3.201, 5.854, 6.215, -3.36, 1673.0, -7.01e-7, -1.54e-10, 3.85e-17]

    # Build comp_dict from fixed input oxide order, ignoring "O"
    input_oxides = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"]
    comp_dict = Dict(input_oxides[i] => bulkin[i] for i in eachindex(input_oxides) if input_oxides[i] != "O")

    # Ensure Fe2O3 is present (may not be in input)
    get!(comp_dict, "Fe2O3", 0.0)

    # Divide composition by molar mass to get moles (only for keys present in both comp_dict and mm)
    mol = Dict(k => comp_dict[k] / mm[k] for k in keys(comp_dict) if haskey(mm, k))

    mol_sum = sum(values(mol))
    mol_frac = Dict(k => v / mol_sum for (k, v) in mol)

    # Use Kress & Carmichael Eq 7 to calculate Fe2O3/FeO ratio
    Fe2O3FeO = exp((param[1] * log(fo2)) + (param[2] / t) + param[3] +
                   ((param[4] * get(mol_frac, "Al2O3", 0.0)) + (param[5] * get(mol_frac, "FeO", 0.0)) +
                    (param[6] * get(mol_frac, "CaO", 0.0)) + (param[7] * get(mol_frac, "Na2O", 0.0)) +
                    (param[8] * get(mol_frac, "K2O", 0.0))) +
                   param[9] * (1 - param[10] / t - log(t / param[10]) +
                   param[11] * (p / t) + param[12] * ((t - param[10]) * (p / t) +
                   param[13] * (p^2) / t))
               )

    FeOtotal = ((comp_dict["FeO"] + comp_dict["Fe2O3"] / 1.1111) / mm["FeO"]) / mol_sum
    XFe2O3 = Fe2O3FeO * FeOtotal / (2 * Fe2O3FeO + 1)
    XFeO = FeOtotal / (1 + 2 * Fe2O3FeO)
    molFe3 = XFe2O3 * 2
    molFe2 = XFeO
    molFe3_Fetot_ratio = molFe3 / (molFe2 + molFe3)
    Fe3_Fe2 = molFe3 / molFe2

    mol_frac["FeO"] = XFeO
    mol_frac["Fe2O3"] = XFe2O3

    xmw = Dict(k => mol_frac[k] * mm[k] for k in keys(mol_frac) if haskey(mm, k))
    xmw_sum = sum(values(xmw))
    wt = Dict(k => v * (100 / xmw_sum) for (k, v) in xmw)

    # Build output array in fixed order, replacing O with Fe2O3
    output_oxides = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "Fe2O3", "K2O", "Na2O", "TiO2", "Cr2O3", "H2O"]
    new_comp = [get(wt, ox, get(comp_dict, ox, 0.0)) for ox in output_oxides]

    return new_comp, Fe3_Fe2
end


end # Module