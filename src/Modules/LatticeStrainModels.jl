module LatticeStrainModels

export LatticeStrainModel, get_lsm

const R = 8.314  # J/mol/K

struct LatticeStrainModel
    mineral::String
    charge::Int
	coordination::Int
    authors::String
    r0::Union{Float64, Function}
    E::Union{Float64, Function}
    D0::Union{Float64, Function}
end


""" Feldspar; 3+; Bédard (2023) """
const _fsp3_bed23 = LatticeStrainModel(
    "fsp", 3, 8, "bed23",
    params -> -0.067*params[:params[:XAn]] + 1.331,
    params -> begin
		if params[:params[:meltH2O]] < 1.0
			E = 1.2808*params[:XAn] + 3.4892
		elseif params[:meltH2O] >= 1.0 & params[:meltH2O] < 5.0
			E = 1.247*params[:XAn] + 3.4238
		else
			E = 1.2675*params[:XAn] + 3.3288
		end
		return exp(E)
	end,
    params -> exp(-2.8328*params[:XAn] - 0.0538)
)



""" Feldspar; 3+; Sun et al. (2017) """
const _fsp3_sun17 = LatticeStrainModel(
    "fsp", 3, 8, "sun17_3",
    1.179,
    196.0,
    params -> exp(16.05 - (19.45 + 1.17*(params[:P]/10)^2) / (R*(params[:T]+273.15)) * 1e4 - 5.17*params[:fspXCa]^2)
)


""" Feldspar; 2+; Sun et al. (2017) """
const _fsp2_sun17 = LatticeStrainModel(
    "fsp", 2, 8, "sun17_2",
    params -> 0.075*params[:fspXNa]+1.189,
    params -> -487*(0.075*params[:fspXNa]+1.189) + 719,
    params -> exp(2.39*params[:fspXNa]^2+(6910-2542*(params[:P]/10)^2)/(R*(params[:T]+273.15)))
)


""" Clinopyroxene; 3+; Sun & Liang (2012) """
const _cpx3_sl12 = LatticeStrainModel(
    "cpx", 3, 8, "sl12",
    params -> -0.212*params[:cpxXMgM2] - 0.104*params[:cpxXAlM1] + 1.066,
    params -> 1000.0*(2.27*(-0.212*params[:cpxXMgM2] - 0.104*params[:cpxXAlM1] + 1.066)-2.00),
    params -> exp(-0.91*params[:meltXH2O] + 1.98*params[:cpxXMgM2] + 4.37*params[:cpxXAlT] + 71900/(R*(params[:T]+273.15)) - 7.14)
)

""" Clinopyroxene M1; 3+; Bédard (2014) """
const _cpx3_bed14M1 = LatticeStrainModel(
    "cpx", 3, 8, "bed14M1",
    params -> begin
		if (18e4/(params[:T]+273.15)-13.0) < -0.17112
			r0 = -0.0040115*(18e4/(params[:T]+273.15)-13.0)^2-0.019187*(18e4/(params[:T]+273.15)-13.0)+0.74173
		else
			r0 = 0.744896
		end
		return r0
	end,
    params -> begin
		if (18e4/(params[:T]+273.15)-13.0) < -0.15489
			E = 47.214*(18e4/(params[:T]+273.15)-13.0)^2-135.08*(18e4/(params[:T]+273.15)-13.0)+694.21
		else
			E = -0.86319*(18e4/(params[:T]+273.15)-13.0)^2+31.794*(18e4/(params[:T]+273.15)-13.0)+721.21
		end
		return E
	end,
    params -> exp(1.2849*(18e4/(params[:T]+273.15)-13.0)+1.3636)
)

""" Clinopyroxene M2; 3+; Bédard (2014) """
const _cpx3_bed14M2 = LatticeStrainModel(
    "cpx", 3, 8, "bed14M2",
	params -> 0.00020599*(18e4/(params[:T]+273.15)-13.0)^3-0.0017107*(18e4/(params[:T]+273.15)-13.0)^2+0.013443*(18e4/(params[:T]+273.15)-13.0)+1.0401,
    params -> begin
		if (18e4/(params[:T]+273.15)-13.0) < -3.0176
			E = -1.4171*(18e4/(params[:T]+273.15)-13.0)^2 + 1.2721*(18e4/(params[:T]+273.15)-13.0) + 694.21
		elseif (18e4/(params[:T]+273.15)-13.0) >= -3.0176 & (18e4/(params[:T]+273.15)-13.0) < -1.95982
			E = 8.5143*(18e4/(params[:T]+273.15)-13.0)^2 + 71.2*(18e4/(params[:T]+273.15)-13.0) + 429.54
		else
			E = -0.60535*(18e4/(params[:T]+273.15)-13.0)^2 + 13.104*(18e4/(params[:T]+273.15)-13.0) + 350.71
		end
		return E
	end,
    params -> exp(-0.0020082*(18e4/(params[:T]+273.15)-13.0)^3+0.02549*(18e4/(params[:T]+273.15)-13.0)^2+0.86352*(18e4/(params[:T]+273.15)-13.0)+0.17177)
)

""" Orthopyroxene; 3; Bédard (2025) """
const _opx3_bed25 = LatticeStrainModel(
    "opx", 3, 6, "bed25",
    0.81,
    params -> begin
		if (params[:meltMgO]*100.0) > 3.58
			E = 78.495 + 26.867 * log(params[:meltMgO]*100.0)
		else
			E = 103.564 + 6.194 * log(params[:meltMgO]*100.0)
		end
		return E
	end,
    params -> begin
		if (params[:meltMgO]*100.0) > 3.4
			D0 = exp(0.5345 - 0.713 * log(params[:meltMgO]*100.0))
		else
			D0 = exp(0.6559 - 0.8031 * log(params[:meltMgO]*100.0))
		end
		return D0
	end
)

""" Amphibole; 3; Shimizu et al. (2017) mineral """
const _amp3_shi17mineral = LatticeStrainModel(
    "amp", 3, 8, "shi17mineral",
    params -> -0.039*params[:ampXFmM4] + 1.043,
    337.0,
    params -> exp(-2.95*params[:ampXK] - 1.83*params[:ampXNa] - 0.35*params[:ampXMg] - 1.52*params[:ampXTi] + 72700/(R*(params[:T]+273.15)) - 4.21)
)

""" Amphibole; 3; Shimizu et al. (2017) melt """
const _amp3_shi17melt = LatticeStrainModel(
    "amp", 3, 8, "shi17melt",
    params -> -0.048*params[:ampXFmM4] + 1.045,
    341.0,
    params -> exp(-3.08)*params[:meltXSi]^0.74*params[:meltXTi]^(-0.33)*params[:meltXCa]^(-0.84)
)

""" Garnet; 3; Meltzer & Kessel (2020) """
const _g3_mk20 = LatticeStrainModel(
    "g", 3, 8, "mk20",
    params -> 1.083 - 9.027e-5*(params[:T]+273.15) - 7.865e-4*params[:gtMgO]/params[:meltMgOhydr],
    params -> 350.0*params[:meltχH2O] - 542.0*params[:meltMgNhydr] + 1854.0*params[:meltFeOhydr]/params[:meltSiO2hydr] + 485.0,
    params -> 7.2*params[:gtFeO]/params[:meltFeOhydr]
)

""" Garnet; 3; Sun & Liang (2013) """
const _g3_sl13 = LatticeStrainModel(
    "g", 3, 8, "sl13",
    params -> 0.155*params[:XCa] + 0.78,
    params -> 1000.0*(2.29*(0.155*params[:XCa] + 0.78)-1.62),
    params -> exp(-1.02*params[:gtXCa] + (91700.0 - 91.34*(params[:P]/10.0)*(38.0-(params[:P]/10.0)))/(R*(params[:T]+273.15)) - 2.05)
)

""" Spinel/Magentite; Sievwright et al. (2020) """
const _spl3_sie20 = LatticeStrainModel(
    "spl", 3, 6, "sie20",
    0.63,
    340.7,
    0.94
)

""" Olivine; 3; Bédard (2005) """
const _ol3_bed05 = LatticeStrainModel(
    "ol", 3, 6, "bed05",
    0.807,
    params -> begin
		if (params[:meltMgO]*100.0) > 1.5
			E = 5.346736 * params[:meltMgO]*100.0 + 170.196
		else
			E = 38.85714 * params[:meltMgO]*100.0 + 120.6667
		end
		return E 
	end,
    params -> exp(-0.6126633*log(params[:meltMgO]*100.0) - 0.02713243)
)

""" Apatite; 3; Jirku et al. (2025) """
const _ap3_jir25 = LatticeStrainModel(
    "ap", 3, 9, "jir25",
    params -> 7.191e-5*(params[:T]+273.15) + 3.723e-3*params[:meltXP2O5] - 15.6e-4*params[:meltXCa] - 1.047,
    params -> 0.366*(params[:T]+273.15) + 6.112*params[:meltXP2O5] - 7.388*params[:meltXCa] - 135.8,
    params -> exp(5.16*params[:meltXCa]^(-0.416))
)

""" Zircon; 3; Streicher et al. (2023) """
const _zr3_str23 = LatticeStrainModel(
    "zr", 3, 8, "str23",
    0.94, # Average of Q1,Q2 in the paper
    572.5, # Average of Q1,Q2 in the paper
    params -> exp(13594.0/(params[:T]+273.15) - 7.1266)
)


const _REGISTRY = Dict{Tuple{String,Int,String}, LatticeStrainModel}(
	("fsp", 3, "bed23") => _fsp3_bed23,
    ("fsp", 3, "sun17_3") => _fsp3_sun17,
	("fsp", 2, "sun17_2") => _fsp2_sun17,
	("cpx", 3, "sl12") => _cpx3_sl12,
	("cpx", 3, "bed14M1") => _cpx3_bed14M1,
	("cpx", 3, "bed14M2") => _cpx3_bed14M2,
	("opx", 3, "bed25") => _opx3_bed25,
	("amp", 3, "shi17mineral") => _amp3_shi17mineral,
	("amp", 3, "shi17melt") => _amp3_shi17melt,
	("g", 3, "mk20") => _g3_mk20,
	("g", 3, "sl13") => _g3_sl13,
	("spl", 3, "sie20") => _spl3_sie20,
	("ol", 3, "bed05") => _ol3_bed05,
	("ap", 3, "jir25") => _ap3_jir25,
	("zr", 3, "str23") => _zr3_str23
)

function get_lsm(; mineral::String, charge::Int, authors::String)
    key = (mineral, charge, authors)
    haskey(_REGISTRY, key) || error(
        "No model for mineral=$(mineral), charge=$(charge), authors=$(authors). " *
        "Available: $(collect(keys(_REGISTRY)))"
    )
    return _REGISTRY[key]
end

end  # module