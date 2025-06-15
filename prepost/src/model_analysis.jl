module BoModels

export sim_params, overshoot_model, dhr_model, hyst_model, model_summary, intrinsic_bond_num

using ..CurvePlots: _read_vf, _probe2df, FluidisedBed, _calc_vel!

using DataFrames
using CSV
using Statistics

function sim_params(flbed::FluidisedBed; 
    p_diameter::Float64,
    rho_p::Real,
    cg_factor::Real,
    poisson_ratio::Real,
    youngs_modulus::Real,
    ced::Union{Real, Nothing} = nothing
)

    if !isnothing(cg_factor)
        flbed.𝜌_p = rho_p / cg_factor
        flbed.p_diameter = p_diameter * cg_factor
    else
        flbed.𝜌_p = rho_p
        flbed.p_diameter = p_diameter
    end

    flbed.poisson_ratio = poisson_ratio
    flbed.youngs_modulus = youngs_modulus
    flbed.ced = ced

    println("FluidisedBed parameters set")
end

function _store_pressure(flbed)

    if (flbed._df_store != ['z', "pressure"]) || (flbed._timeser_df == nothing)
        flbed._timeser_df = _probe2df(flbed, true, 'z', nothing)
        _calc_vel!(flbed, flbed._timeser_df)
    end

    pressure_gdf = groupby(flbed._timeser_df, [:direction, :v_z])
    cols_to_average = valuecols(pressure_gdf)
    vel_plot_df = combine(pressure_gdf, cols_to_average .=> mean)

    vel_up = vel_up = filter(:direction => in(["up", "max"]), vel_plot_df)
    sort!(vel_up, :v_z)

    vel_down = filter(:direction => in(["down", "max"]), vel_plot_df)
    sort!(vel_down, :v_z)

    println("vel_up: ", vel_up)

    flbed._model_store["p_up"] = vel_up[!, :probe_0_mean]
    flbed._model_store["p_down"] = vel_down[!, :probe_0_mean]

    max_idx = argmax(flbed._model_store["p_up"])
    flbed._model_store["v_z_up"] = vel_up.v_z
    flbed._model_store["v_z_down"] = vel_down.v_z

    flbed._model_store["u_mf"] = [vel_up[max_idx, :v_z]]

end

function _store_vf(flbed::FluidisedBed)

    if !isnothing(flbed.void_frac_path)
        post_dir = flbed.void_frac_path
    else
        post_dir = "CFD/postProcessing/cuttingPlane/"
    end

    if (flbed._df_store != ['y', "pressure"]) || (flbed._timeser_df == nothing)
        flbed._timeser_df = _read_vf(flbed, post_dir, 'y')
        _calc_vel!(flbed, flbed._timeser_df)
    end

    vf_gdf = groupby(flbed._timeser_df, [:direction, :v_z])
    cols_to_average = valuecols(vf_gdf)
    vf_plot_df = combine(vf_gdf, cols_to_average .=> mean)

    vel_up = filter(:direction => in(["up", "max"]), vf_plot_df)
    sort!(vel_up, :v_z)

    vel_down = filter(:direction => in(["down", "max"]), vf_plot_df)
    sort!(vel_down, :v_z)

    flbed._model_store["vf_up"] = vel_up[!, :void_fraction_mean]
    flbed._model_store["vf_down"] = vel_down[!, :void_fraction_mean]
end

function _store_contactn(flbed::FluidisedBed)
    if !isfile("DEM/post/collisions.csv")
        throw(LoadError("collisions.csv not found. Check source code for details."))
    end

    contact_df = DataFrame(CSV.File("DEM/post/collisions.csv"))
    contact_df[!, :time] = contact_df[!, :time] .- minimum(contact_df[!, :time])

    _calc_vel!(flbed, contact_df)
    contact_gdf = groupby(contact_df, [:direction, :v_z])
    cols_to_average = valuecols(contact_gdf)
    contact_plot_df = combine(contact_gdf, cols_to_average .=> mean)

    ncontact_up = filter(:direction => in(["up", "max"]), contact_plot_df)
    sort!(ncontact_up, :v_z)
    ncontact_down = filter(:direction => in(["down", "max"]), contact_plot_df)
    sort!(ncontact_down, :v_z)



    flbed._model_store["ncontact_up"] = ncontact_up[!, :n_contact_mean] ./ ncontact_up[!, :n_atoms_mean]
    flbed._model_store["ncontact_down"] = ncontact_down[!, :n_contact_mean] ./ ncontact_down[!, :n_atoms_mean]
end

function _store_data(flbed::FluidisedBed)

    if !haskey(flbed._model_store, "p_up") || !haskey(flbed._model_store, "p_down")
        println("Storing pressure data...")
        _store_pressure(flbed)
    end

    if !haskey(flbed._model_store, "vf_up") || !haskey(flbed._model_store, "vf_down")
        println("Storing void fraction data...")
        _store_vf(flbed)
    end

    if !haskey(flbed._model_store, "ncontact_up") || !haskey(flbed._model_store, "ncontact_down")
        println("Storing contact number data...")
        _store_contactn(flbed)
    end

end

function _validate_params_exist(flbed::FluidisedBed)

    if isnothing(flbed.𝜌_p) || isnothing(flbed.p_diameter)
        throw(ArgumentError("FluidisedBed parameters not set. Use sim_params() to set them."))
    end
end

function overshoot_model(flbed::FluidisedBed)

    _validate_params_exist(flbed)

    p_1 = maximum(flbed._model_store["p_up"])
    p_ss = flbed._model_store["p_down"][end]
    p_over = p_1 - p_ss

    println("p_1: ", p_1, " p_ss: ", p_ss, " p_over: ", p_over)
    
    avg_contactn = mean(
        vcat(flbed._model_store["ncontact_up"], flbed._model_store["ncontact_down"])
    )

    println("avg_contactn: ", avg_contactn)

    avg_vf = mean(
        vcat(flbed._model_store["vf_up"], flbed._model_store["vf_down"])
    )

    return (6 * p_over)/(avg_contactn^2 * (1-avg_vf) * flbed.p_diameter * flbed.𝜌_p * 9.81)
end

function dhr_model(flbed::FluidisedBed)

    _validate_params_exist(flbed)

    println("U_mf: ", flbed._model_store["u_mf"])
    println("v_z_up: ", flbed._model_store["v_z_up"])

    umf_idx_up = findall(x -> x == flbed._model_store["u_mf"][1], flbed._model_store["v_z_up"])
    umf_idx_down = findall(x -> x == flbed._model_store["u_mf"][1], flbed._model_store["v_z_down"])

    println("umf_idx_up: ", umf_idx_up, " umf_idx_down: ", umf_idx_down)
    

    ε1 = flbed._model_store["vf_up"][umf_idx_up][1]
    ε2 = flbed._model_store["vf_down"][umf_idx_down][1]

    println("ε1: ", ε1, " ε2: ", ε2)

    return (ε2/ε1)^3 * (1-ε2)/(1-ε2) - 1
end

function hyst_model(flbed::FluidisedBed)

    _validate_params_exist(flbed)

    umf_idx_up = findall(x -> x == flbed._model_store["u_mf"][1], flbed._model_store["v_z_up"])
    umf_idx_down = findall(x -> x == flbed._model_store["u_mf"][1], flbed._model_store["v_z_down"])

    println("p_up:" , flbed._model_store["p_up"])

    p_1 = maximum(flbed._model_store["p_up"])
    p_ss = flbed._model_store["p_up"][end]
    p_2 = flbed._model_store["p_down"][umf_idx_down][1]

    Δk = (flbed._model_store["ncontact_up"][umf_idx_up]
        - flbed._model_store["ncontact_down"][umf_idx_down])

    return (p_1 - p_2) / (p_ss * Δk)
end

function intrinsic_bond_num(flbed::FluidisedBed)

    _validate_params_exist(flbed)
    
    E_eq = 0.5 * flbed.youngs_modulus / (1 - flbed.poisson_ratio)^2
    R_eq = flbed.p_diameter / 4
    p_mass = flbed.𝜌_p * flbed.p_diameter^3 / 6

    numerator = 27 * (π^3) * (flbed.ced^3) * (R_eq^2)
    denominator = 64 * p_mass * 9.81 * (E_eq^2)

    return numerator / denominator
end
    

function model_summary(flbed::FluidisedBed)
    
    _store_data(flbed)

    overshoot = overshoot_model(flbed)
    dhr = dhr_model(flbed)
    hyst = hyst_model(flbed)
    intr = intrinsic_bond_num(flbed)

    return Dict(
        "overshoot_model" => overshoot,
        "dhr_model" => dhr,
        "hyst_model" => hyst,
        "intrinsic_bo" => intr
    )
end


end