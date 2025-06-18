module BoModels

export sim_params, overshoot_model, dhr_model, hyst_model, model_summary, intrinsic_bond_num

using ..CurvePlots: _read_vf, _probe2df, FluidisedBed, _calc_vel!

using DataFrames
using CSV
using Statistics

"""
    sim_params(flbed::FluidisedBed; 
        p_diameter::Float64,
        rho_p::Real,
        cg_factor::Real,
        poisson_ratio::Real,
        youngs_modulus::Real,
        ced::Union{Real, Nothing} = nothing
    )

Sets simulation parameters for a fluidised bed model.

# Arguments
- `flbed::FluidisedBed`: The fluidised bed object to be modified
- `p_diameter::Float64`: Particle diameter
- `rho_p::Real`: Particle density
- `cg_factor::Real`: Coarse-graining factor for particle scaling
- `poisson_ratio::Real`: Poisson's ratio for particle material
- `youngs_modulus::Real`: Young's modulus for particle material
- `ced::Union{Real, Nothing}=nothing`: Coefficient of elastic damping (optional)

# Effects
Updates the properties of the fluidised bed object, including:
- Particle density (scaled by coarse-graining factor if provided)
- Particle diameter (scaled by coarse-graining factor if provided)
- Material properties (Poisson's ratio, Young's modulus)
- Coefficient of elastic damping (if provided)

Prints a confirmation message after parameters are set.
"""
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

"""
    _store_pressure(flbed)

Calculate and store pressure-related data for a fluidized bed.

This function processes pressure measurements from the fluidized bed simulation
to determine key parameters such as minimum fluidization velocity (u_mf).

# Arguments
- `flbed`: A fluidized bed object containing simulation data

# Operation
1. Ensures the appropriate data frame with pressure and position data is available
2. Groups the data by direction and vertical velocity
3. Calculates average pressure values for each group
4. Separates data into upward and downward directions
5. Stores the processed pressure data in the model_store dictionary:
   - `p_up`: Average pressure values for upward flow
   - `p_down`: Average pressure values for downward flow
   - `v_z_up`: Vertical velocities for upward flow
   - `v_z_down`: Vertical velocities for downward flow
   - `u_mf`: Minimum fluidization velocity (determined from maximum pressure in upward flow)

# Notes
- Relies on helper functions `_probe2df` and `_calc_vel!`
- Modifies the `flbed._model_store` and potentially `flbed._timeser_df` objects
"""
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

"""
    _store_vf(flbed::FluidisedBed)

Store void fraction data for a fluidised bed simulation.

This function reads and processes void fraction data from simulation results.
It performs the following steps:
1. Determines the post-processing directory from `flbed.void_frac_path` or uses a default path
2. Reads void fraction data along the y-direction if not already loaded
3. Calculates velocity information for the data
4. Groups the data by direction and vertical velocity
5. Computes mean values for all value columns
6. Separates and sorts data for upward and downward flow
7. Stores the processed void fraction data in the fluidised bed model's storage dictionary

The function stores two key datasets in `flbed._model_store`:
- "vf_up": Mean void fraction for upward flow
- "vf_down": Mean void fraction for downward flow
"""
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

"""
    _store_contactn(flbed::FluidisedBed)

Store the normalized contact number data for a fluidised bed simulation.

This function processes collision data from 'DEM/post/collisions.csv' to calculate 
and store the average number of contacts per particle during upward and downward 
bed movement phases.

# Process:
1. Reads collision data from CSV file
2. Normalizes time values relative to the minimum time
3. Calculates velocities using `_calc_vel!`
4. Groups data by direction and vertical velocity
5. Computes mean values for all numeric columns
6. Separates data for upward and downward movement
7. Stores normalized contact numbers (contacts per particle) in the model store

# Arguments
- `flbed::FluidisedBed`: The fluidised bed object containing simulation data

# Throws
- `LoadError`: If the collisions.csv file is not found

# Notes
- Results are stored in `flbed._model_store` under keys "ncontact_up" and "ncontact_down"
- The function requires `_calc_vel!` to be defined elsewhere to calculate velocities
"""
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

"""
    _store_data(flbed::FluidisedBed)

Internal function that ensures required data is stored in the model for analysis.

This function checks if specific data keys exist in the `_model_store` dictionary
of the `FluidisedBed` object and calls the appropriate storage functions if they don't.

It manages storage of:
- Pressure data (via `_store_pressure`)
- Void fraction data (via `_store_vf`)
- Contact number data (via `_store_contactn`)

Each data type is only stored if it doesn't already exist in the model store.
"""
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

"""
    _validate_params_exist(flbed::FluidisedBed)

Validates that essential parameters of a `FluidisedBed` object are set.

# Arguments
- `flbed::FluidisedBed`: The fluidised bed object to validate.

# Throws
- `ArgumentError`: If the particle density (`𝜌_p`) or particle diameter (`p_diameter`) are not set.

# Details
This internal function is used to ensure that essential parameters are properly 
defined before performing calculations. It checks if `𝜌_p` or `p_diameter` are 
`nothing` and throws an error if they are, suggesting to use the `sim_params()` 
function to set the required parameters.
"""
function _validate_params_exist(flbed::FluidisedBed)

    if isnothing(flbed.𝜌_p) || isnothing(flbed.p_diameter)
        throw(ArgumentError("FluidisedBed parameters not set. Use sim_params() to set them."))
    end
end

"""
    overshoot_model(flbed::FluidisedBed) -> Float64

Calculate the model parameter from pressure overshoot in a fluidised bed.

This function computes a dimensionless parameter based on the pressure overshoot 
observed in a fluidisation-defluidisation cycle. It uses the maximum pressure during 
fluidisation (`p_1`), the steady-state pressure (`p_ss`), the average number of contacts, 
and the average void fraction.

# Arguments
- `flbed::FluidisedBed`: A fluidised bed object containing simulation data and parameters

# Returns
- A dimensionless parameter calculated as `(6 * p_over)/(avg_contactn^2 * (1-avg_vf) * flbed.p_diameter * flbed.𝜌_p * 9.81)`
  where `p_over = p_1 - p_ss`

# Notes
- Validates that required parameters exist in the fluidised bed object before calculation
- Uses both the fluidisation (up) and defluidisation (down) data for averaging contact numbers and void fractions
"""
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

"""
    dhr_model(flbed::FluidisedBed)

Calculate the Richardson-Zaki bed expansion parameter using the Davidson-Harrison-Richardson model.

This function computes the expansion parameter from void fraction data at minimum fluidization velocity
during both the upward (increasing velocity) and downward (decreasing velocity) phases of fluidization.

The calculation follows the relation:
(ε₂/ε₁)³ * (1-ε₁)/(1-ε₂) - 1

where:
- ε₁ is the void fraction at minimum fluidization during upward flow
- ε₂ is the void fraction at minimum fluidization during downward flow

# Arguments
- `flbed::FluidisedBed`: A fluidised bed object containing simulation data with:
  - `_model_store["u_mf"]`: Minimum fluidization velocity
  - `_model_store["v_z_up"]`: Upward flow velocities
  - `_model_store["v_z_down"]`: Downward flow velocities
  - `_model_store["vf_up"]`: Void fractions for upward flow
  - `_model_store["vf_down"]`: Void fractions for downward flow

# Returns
The calculated expansion parameter value.

# Notes
- The function validates that required parameters exist in the `flbed` object before calculation.
- The function prints intermediate values for debugging purposes.
"""
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

"""
    hyst_model(flbed::FluidisedBed) -> Float64

Calculate the hysteresis model parameter based on pressure differences in a fluidised bed.

This function analyzes the hysteresis behavior by comparing pressure and contact data
between the up-fluidisation and down-fluidisation processes.

# Arguments
- `flbed::FluidisedBed`: A fluidised bed object containing model data in `_model_store`

# Returns
The hysteresis model parameter calculated as (p₁ - p₂)/(pₛₛ * Δk), where:
- p₁: maximum pressure during up-fluidisation
- p₂: pressure at minimum fluidisation velocity during down-fluidisation
- pₛₛ: steady-state pressure at the end of up-fluidisation
- Δk: difference in number of contacts between up and down fluidisation at minimum fluidisation velocity

# Implementation Details
The function first validates that required parameters exist in the fluidised bed object,
then identifies indices corresponding to minimum fluidisation velocity in both up and down
fluidisation data, and finally calculates the hysteresis parameter.
"""
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

"""
    intrinsic_bond_num(flbed::FluidisedBed)

Calculate the intrinsic bond number for a fluidised bed system.

The intrinsic bond number is a dimensionless parameter that quantifies the ratio 
of cohesive forces to gravitational forces acting on particles in a fluidised bed.

# Arguments
- `flbed::FluidisedBed`: A fluidised bed object containing necessary parameters

# Required Parameters
The function validates that the following parameters exist in `flbed`:
- `youngs_modulus`: Young's modulus of the particle material
- `poisson_ratio`: Poisson's ratio of the particle material
- `p_diameter`: Particle diameter
- `𝜌_p`: Particle density
- `ced`: Cohesive energy density

# Returns
- The intrinsic bond number (dimensionless)

# Formula
The intrinsic bond number is calculated as:
    
    Bon = [27π³(ced)³(R_eq)²] / [64 × p_mass × g × (E_eq)²]

where:
- E_eq = 0.5 × youngs_modulus / (1 - poisson_ratio)²
- R_eq = p_diameter / 4
- p_mass = 𝜌_p × p_diameter³ / 6
- g = 9.81 m/s²
"""
function intrinsic_bond_num(flbed::FluidisedBed)

    _validate_params_exist(flbed)
    
    E_eq = 0.5 * flbed.youngs_modulus / (1 - flbed.poisson_ratio)^2
    R_eq = flbed.p_diameter / 4
    p_mass = flbed.𝜌_p * flbed.p_diameter^3 / 6

    numerator = 27 * (π^3) * (flbed.ced^3) * (R_eq^2)
    denominator = 64 * p_mass * 9.81 * (E_eq^2)

    return numerator / denominator
end
    

"""
    model_summary(flbed::FluidisedBed)

Calculate and return a summary of various fluidisation models for a given fluidised bed.

This function evaluates several models to characterize the behavior of the fluidised bed:
- Overshoot model
- Dimensionless height ratio (DHR) model
- Hysteresis model
- Intrinsic bond number

Before performing the calculations, it stores the current state of the fluidised bed
using the `_store_data` function.

# Arguments
- `flbed::FluidisedBed`: The fluidised bed object to analyze.

# Returns
- `Dict`: A dictionary containing the results of each model:
  - `"overshoot_model"`: Results from the overshoot model
  - `"dhr_model"`: Results from the dimensionless height ratio model
  - `"hyst_model"`: Results from the hysteresis model
  - `"intrinsic_bo"`: Calculated intrinsic bond number

# See also
- [`overshoot_model`](@ref)
- [`dhr_model`](@ref)
- [`hyst_model`](@ref)
- [`intrinsic_bond_num`](@ref)
"""
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