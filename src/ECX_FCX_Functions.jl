""" Functions I am using in forward and inverse modelling of equilibrium and fractional crystallization. """



# Equilibrium crystallization

function ECX_run_buffered(bulk_0::Vector{Float64}, TE_0::Vector{Float64}, P::Vector{Float64}, T::Vector{Float64}, boff::Vector{Float64})

    """
     Input:

        bulk_0: Vector of initial bulk composition 
                Input format, mol% or wt%, as well as oxides need to be defined outside of this function
        
        TE_0: Vector of initial trace element concentrations in ppm
              Elements and their charges need to be defined outside of this function

        P: Vector of pressures for each step in kbar

        T: Vector of temperatures for each step in °C

        boff: Vector of buffer offsets for each step in log units
              Buffer needs to be defined outside of this function

    Output:

        Out_XY: Vector of MAGEMin outputs from single point minimization for each step

        Out_Sat: Vector of outputs from MAGEMin saturation model for each step

        Out_TE: Vector of dictionaries of trace element concentrations for each step
                Out_TE[Int(step)][String("phase")][Symbol(:element)] gives the concentration of that element in that phase for that step

        Out_Ds: Vector of dictionaries ofpartition coefficients for each step
                Out_Ds[Int(step)][String("phase")][Symbol(:element)] gives the concentration of that element in that phase for that step
                
        Out_EuR: Vector of Eu/R ratio for each step
    
    """

    # Globals defined in other cells; copy into local bindings for safe, fast access
    sys_in_local    = sys_in
    oxides_local    = oxides
    sat_C0s_local   = sat_C0s
    sat_KdsDB_local = sat_KdsDB
    elements_local  = elements
    charges_local   = charges
    buffer_local    = buffer


    # Initalize MAGEMin
    data = Initialize_MAGEMin("ig", buffer=buffer_local, verbose=false)


    # Initialize melt fraction
    meltF = 1.0
    # Initialize counter
    np = 0


    # Initialize output storage
    Out_XY = Vector{out_struct}(undef, nsteps)
    Out_Sat = Vector{out_TE_struct}(undef, nsteps)


    # Out_TE should be a vector of dictionaries, where each dictionary corresponds to a step. Each step-dictionary is to be structured like Dict("phase" => Dict(:element => value))
    Out_TE = Vector{Dict{String,Dict{Symbol,Float64}}}(undef, nsteps)


    # Similar to Out_TE 
    Out_Ds = Vector{Dict{String,Dict{Symbol,Float64}}}(undef, nsteps)

    Out_EuR = zeros(Float64, nsteps)

    param_dicts = [Dict{Symbol,Float64}() for _ in 1:nsteps]

    while meltF > 0.0 && np < nsteps

        np += 1

        # Reset initial bulk composition -> No melt sequestration
        bulk = deepcopy(bulk_0)
        sys_in_local = sys_in

        Out_TE[np] = Dict{String,Dict{Symbol,Float64}}()
        Out_Ds[np] = Dict{String,Dict{Symbol,Float64}}()

        # If input bulk composition is in wt.%, convert to mol.%
        if sys_in_local == "wt"
            bulk, oxides_local = convertBulk4MAGEMin(bulk, oxides_local, "wt", "ig")
            sys_in_local = "mol"
        end

        # Routine from https://computationalthermodynamics.github.io/MAGEMin_C.jl/dev/MAGEMin_C/saturation_models#E.4-Saturation-models-and-bulk-correction
        tol = 1e-6
        res = 1.0
        n0 = 0.0
        ite = 0

        while res > tol && ite < 32
            #println(bulk)
            out = single_point_minimization(P[np], T[np], data, X=bulk, Xoxides=oxides_local, B=boff[np], sys_in=sys_in_local, name_solvus=true)
            Out_XY[np] = deepcopy(out) # Store output

            Out_Sat[np] = TE_prediction(out, sat_C0s_local, sat_KdsDB_local, "ig";
                ZrSat_model="CB",
                P2O5Sat_model="Klein26")


            # Reset P_kbar 
            #out = @set out.P_kbar = pkbar
            bulk .-= Out_Sat[np].bulk_cor_mol
            res = abs(n0 - vec_norm(Out_Sat[np].bulk_cor_mol))
            n0 = vec_norm(Out_Sat[np].bulk_cor_mol)

            ite += 1
            if ite == 32
                @warn "Saturation model did not converge in 32 iterations, residual is $res"
            end

        end

        out = Out_XY[np] # Get the output of the single point minimization for this step

        out_sat = Out_Sat[np] # Get the output of the saturation model for this step

        # Reset oxides
        if np == 1
            oxides_local = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]
        end


        # If any of the entries in bulk are NaN, set them to the original value from bulk_0
        for i in eachindex(bulk)
            if isnan(bulk[i])
                bulk[i] = bulk_0[i]
            else
                # Nothing to do
            end
        end


        meltF = out.frac_M
        println("Step $np: Melt Fraction = $meltF")
        # Break loop if melt fraction is zero or negative (which can happen due to numerical issues)
        if meltF <= 0.0
            println("Melt fraction is zero or negative, stopping iteration.")
            break
        end
        meltF_wt = out.frac_M_wt

        outte = Dict{String,Dict{Symbol,Float64}}()
        out_d = Dict{String,Dict{Symbol,Float64}}()

        phases = copy(out.ph)

        # Work on a local copy of phase fractions so the output `out` is unchanged
        phase_fracs_wt = copy(out.ph_frac_wt)

        # Add 'zr' and 'ap' to phases if zrc_wt > 0 and fapt_wt > 0, respectively
        # Append fractions to the local `phase_fracs_wt` to keep alignment
        if out_sat.zrc_wt > 0.0 && !("zr" in phases)
            push!(phases, "zr")
            push!(phase_fracs_wt, out_sat.zrc_wt)
        end
        if out_sat.fapt_wt > 0.0 && !("ap" in phases)
            push!(phases, "ap")
            push!(phase_fracs_wt, out_sat.fapt_wt)
        end

        for phase in phases
            params = LSMUtils.get_params_for_phase(phase, out, out_sat, T[np], P[np])
            # Skip buffer and melt
            if phase in ["liq", buffer]
                out_d[phase] = Dict(el => 0.0 for el in elements)
            else
                Dphase = PartitionCoefficientsV2.get_D(phase, elements, charges, params)
                out_d[phase] = Dphase
                Out_TE[np][phase] = Dict(el => TE_0[i] * Dphase[el] for (i, el) in enumerate(elements))
            end
        end

        params = LSMUtils.get_params_for_phase("liq", out, out_sat, T[np], P[np])

        Out_EuR[np] = LSMUtils.Eu_ratio_Burnham(params[:logfO2], params[:T], params[:Λ])

        # Need to normalize phase_fracs_wt to 1 while omitting "liq" and "buffer" fractions, it can be done for all fractions as "liq" and "buffer" have zero Kds and thus do not contribute to the sum
        if "liq" in phases && buffer in phases
            phase_fracs_wt_used = phase_fracs_wt ./ (1.0 - phase_fracs_wt[findfirst(==(buffer), phases)] - phase_fracs_wt[findfirst(==("liq"), phases)])
        else
            phase_fracs_wt_used = phase_fracs_wt ./ sum(phase_fracs_wt)
        end

        has_other_phases = any(phase -> !(phase in ["liq", buffer]), phases)

        if has_other_phases
            Kd = Dict(el => sum((phase_fracs_wt_used[i] * get(out_d[phases[i]], el, 0.0) for i in eachindex(phases)); init=0.0) for el in elements)

            # Add Kd to out_d with key "bulk"
            out_d["bulk"] = Kd

            Out_TE[np]["liq"] = Dict(el => TE_0[i] / (Kd[el] + meltF_wt * (1.0 - Kd[el])) for (i, el) in enumerate(elements))
        else
            out_d["bulk"] = Dict(el => NaN for el in elements)

            Out_TE[np]["liq"] = Dict(el => TE_0[i] for (i, el) in enumerate(elements))
        end

        Out_Ds[np] = out_d

    end

    return Out_XY[1:np], Out_Sat[1:np], Out_TE[1:np], Out_Ds[1:np], Out_EuR[1:np], param_dicts[1:np]

end


# Fractional crystallization