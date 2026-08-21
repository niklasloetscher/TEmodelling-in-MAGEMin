module LSMPlotting

using CairoMakie

const PHASE_COLORS = Dict{String, String}([("amp", "#777E49"),("cpx", "#B7E892"), ("opx", "#4682B4"), ("fsp", "#A9A9A9"), ("ilm", "#708090"), ("mu", "#F0E8BB"), ("pat", "#F0E8BB"), ("nph", "#EB7347"), ("spl", "#CE2029"), ("q", "#28282D"),
        ("bi", "#2F1B12"), ("chl", "#7CFC00"), ("sph", "#942222"), ("ru", "#FF1649"), ("ep", "#08F26E"), ("g", "#7B1B38"), ("ol", "#FF7415"), ("fl", "#DDEEFA"), ("ne", "#FFFFFF"), ("fper", "#FFFFFF"), ("cd", "#007FFF"), ("ap", "#FF00FF"), ("zr", "#CC7722"), ("liq", "#FFFFFF")])
		
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

const CHONDRITIC_TE = Dict{String, Float64}([])

const colorbarpalette = :buda

canonical(phase::String) = get(PHASE_ALIASES, phase, phase)
		
		
function phasefracplot(phasefracs::Vector{Dict{String, Float64}}, xvals::Union{Vector{Float64}, Vector{Int}}, xvalname::String; noliq=true)
"""
A stacked area plot of phasefractions vs. xvals. Phasefracs contains a series of dicts which in turn contain phases as keys and their fractions as values. noliq=True means that 'liq' is to be skipped as a phase.
"""

    fig = Figure(size = (900, 600), fontsize = 14)
    ax  = Axis(fig[1, 1], xlabel = xvalname, ylabel = "Phase fraction")

    # 1. Collect canonical phase names, collating aliases
    all_phases = String[]
    for d in phasefracs
        for k in keys(d)
            (noliq && k == "liq") && continue
            cph = canonical(k)
            cph ∉ all_phases && push!(all_phases, cph)
        end
    end

    # 2. Build matrix, accumulating aliases into canonical column
    n = length(xvals)
    m = length(all_phases)
    phase_idx = Dict(ph => j for (j, ph) in enumerate(all_phases))
    fracs = zeros(Float64, n, m)
    for (i, d) in enumerate(phasefracs)
        for (rawph, val) in d
            (noliq && rawph == "liq") && continue
            j = get(phase_idx, canonical(rawph), 0)
            j == 0 && continue
            fracs[i, j] += val
        end
    end

    # 3. Compute cumulative sums for stacked bands
    lower = zeros(Float64, n, m)
    upper = zeros(Float64, n, m)
    for j in 1:m
        j > 1 && (lower[:, j] .= upper[:, j-1])
        upper[:, j] .= lower[:, j] .+ fracs[:, j]
    end

    # 4. Draw one band per canonical phase
    for j in 1:m
        band!(ax, xvals, lower[:, j], upper[:, j];
            color = PHASE_COLORS[all_phases[j]],
            label = all_phases[j]
        )
    end

    axislegend(ax; position = :rt)
    ylims!(ax, 0.0, 1.0)
    return fig

end

function phaselineplot(phasefracs::Vector{Dict{String, Float64}}, phases::Vector{String}, xvals::Union{Vector{Float64}, Vector{Int}}, xvalname::String)
"""
A plot containing line traces of phase fractions of selected phases vs. xvals
"""

	fig = Figure(size = (900, 600), fontsize = 14)
    ax  = Axis(fig[1, 1], xlabel = xvalname, ylabel = "Phase fraction")

    # 1. Build a matrix: rows = xvals, cols = phases  (missing → NaN)
    n = length(xvals)
    m = length(phases)

    # Keep NaN if a requested phase is completely absent at a step
    fracs = fill(NaN, n, m)

    # Robust lookup: accept canonical names, aliases, and non-aliased raw names
    phase_idx = Dict{String, Int}()
    for (j, ph) in enumerate(phases)
        phase_idx[ph] = j
        phase_idx[canonical(ph)] = j
    end

    for (i, d) in enumerate(phasefracs)
        for (rawph, val) in d
            cph = canonical(rawph)  # alias -> canonical; non-aliased stays unchanged

            j = get(phase_idx, cph, 0)
            j == 0 && (j = get(phase_idx, rawph, 0))
            j == 0 && continue

            if isnan(fracs[i, j])
                fracs[i, j] = 0.0
            end
            fracs[i, j] += val  # sum canonical + all aliases into one requested phase
        end
    end
	
	
	for (j, ph) in enumerate(phases)
		lines!(ax, xvals, fracs[:,j], color=PHASE_COLORS[canonical(ph)], label = canonical(ph))
	end
	
	axislegend(ax; position = :rt)
	
	return fig

end

function cumulphaseplot(phasefracs::Vector{Dict{String, Float64}}, xvals::Union{Vector{Float64}, Vector{Int}}, xvalname::String)
"""
A stacked area plot of cumulative phase fractions vs. xvals. Phasefracs contains a series of dicts which in turn contain phases as keys and their fractions as values.
"""

    fig = Figure(size = (900, 600), fontsize = 14)
    ax  = Axis(fig[1, 1], xlabel = xvalname, ylabel = "Cumulative phase fraction")

    # 1. Collect canonical phase names, collating aliases
    n = length(xvals)
    meltfracs = zeros(Float64, n)
    all_phases = String[]
    for (i, d) in enumerate(phasefracs)
        for (k, v) in d
            if k == "liq"
                meltfracs[i] = v
            else
                cph = canonical(k)
                cph ∉ all_phases && push!(all_phases, cph)
            end
        end
    end

    # 2. Build matrix, accumulating aliases into canonical column
    m = length(all_phases)
    phase_idx = Dict(ph => j for (j, ph) in enumerate(all_phases))
    fracs = zeros(Float64, n, m)
    cumulative_sums = zeros(Float64, m)
    for (i, d) in enumerate(phasefracs)
        for (rawph, val) in d
            rawph == "liq" && continue
            j = get(phase_idx, canonical(rawph), 0)
            j == 0 && continue
            cumulative_sums[j] += val * (1.0 - meltfracs[i])
            fracs[i, j] = cumulative_sums[j]
        end
        # Assign cumulative sums for phases not present in current dict
        for j in 1:m
            if fracs[i, j] == 0.0
                fracs[i, j] = cumulative_sums[j]
            end
        end
    end

    # 3. Compute cumulative sums for the stacked bands
    lower = zeros(Float64, n, m)
    upper = zeros(Float64, n, m)
    for j in 1:m
        j > 1 && (lower[:, j] .= upper[:, j-1])
        upper[:, j] .= lower[:, j] .+ fracs[:, j]
    end

    # 4. Draw one band per phase
    for j in 1:m
        band!(ax, xvals, lower[:, j], upper[:, j];
            color = PHASE_COLORS[all_phases[j]],
            label = all_phases[j]
        )
    end

    axislegend(ax; position = :rt)
    #ylims!(ax, 0.0, 1.0)
    return fig
	
	
end

function TEplot(TEout::Vector{Dict}, phase::String, elements::Vector{String},
                colorvals::Union{Vector{Float64}, Vector{Int}}, colorvalname::String)
    """
    A plot containing spidergrams of selected TEs of this phase.
    """
    fig = Figure(size = (600, 900), fontsize = 14)
    ax  = Axis(fig[1, 1],
               ylabel = "$(phase) / Chondritic",
               xticks = (1:length(elements), elements),
               xticklabelrotation = π/3)

    # Normalisation reference for selected elements
    chondrite_vals = [CHONDRITIC_TE[el] for el in elements]

    mincval = minimum(colorvals)
    maxcval = maximum(colorvals)
    cmap    = cgrad(colorbarpalette)      # materialise the gradient

    xvals = 1:length(elements)

    # 1. Extract data for each step, normalise, and plot
    for (i, stepdict) in enumerate(TEout)
        # stepdict is Dict(phase => Dict(element => concentration))
        phasedict = get(stepdict, phase, nothing)
        phasedict === nothing && continue   # skip if phase absent in this step

        yvals = [get(phasedict, el, NaN) for el in elements]
        yvals ./= chondrite_vals            # normalise to CI chondrite

        # Map colorval to [0, 1] for palette lookup
        t     = mincval == maxcval ? 0.5 :
                (colorvals[i] - mincval) / (maxcval - mincval)
        color = cmap[t]

        lines!(ax, xvals, yvals; color = color, linewidth = 1.5)
        scatter!(ax, xvals, yvals; color = color, markersize = 6)
    end

    # 2. Add a colorbar
    Colorbar(fig[1, 2];
             colormap = colorbarpalette,
             limits   = (mincval, maxcval),
             label    = colorvalname)

    return fig
	
end

#=function liqlinedesc(TEout::IMPLEMENT, MEout::IMPLEMENT, elements::Vector{String}, xvals::Union{Vector{Float64}, Vector{Int}}, xvalname::String)
"""
Plots of elements vs. xvals, where elements are either major oxides in MEout or trace elements in TEout.
"""



end=#

#=function cumulinedesc(MEout::IMPLEMENT, elements::Vector{String}, xvals::Union{Vector{Float64}, Vector{Int}}, xvalname::String)
"""
Plot of elements vs. xvals where elements are either major oxides in MEout
"""
end=#


		
end # Module