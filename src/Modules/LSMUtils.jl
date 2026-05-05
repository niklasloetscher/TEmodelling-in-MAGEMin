module LSMUtils

# using ..LatticeStrainModels
# using ..PartitionCoefficients

export Eu_ratio_Burnham, Kd_Eu, optical_basicity, get_params_for_phase

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


end # Module