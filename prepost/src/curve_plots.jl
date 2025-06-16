module CurvePlots
"""
Module for plotting pressure and void fraction data from fluidised bed simulations.
This module provides functionality to read pressure and void fraction data from simulation results,
and plot them against time or velocity. It supports both probe-based and slice-based data extraction.
"""

using VTKDataIO
using StatsBase
using Statistics
using DataFrames
using CSV
using Plots

export FluidisedBed,
       plot_pressure,
       plot_voidfrac,
       _read_vf,
       _probe2df

       
"""
    FluidisedBed

A mutable struct representing a fluidised bed simulation configuration and results.

# Fields
- `presure_path::String`: Path to the pressure data files.
- `n_probes::Int`: Number of pressure probes used in the simulation.
- `dump2csv::Bool`: Flag indicating whether to dump data to CSV files.
- `velcfg_path::String`: Path to the velocity configuration files.
- `plots_dir::String`: Directory where plots will be saved.
- `t::Vector{Float32}`: Time values for simulation data.

# Internals
- `v_z::Vector{Float32}`: Vertical velocity values.
- `void_frac_path::Union{String, Nothing}`: Path to void fraction data, if available.
- `_timeser_df::Union{DataFrame, Nothing}`: DataFrame for time series data.
- `_df_store::Union{Array, Nothing}`: Storage for dataframes.
- `p_diameter::Union{Real, Nothing}`: Particle diameter.
- `𝜌_p::Union{Real, Nothing}`: Particle density.
- `poisson_ratio::Union{Real, Nothing}`: Poisson ratio for particles.
- `youngs_modulus::Union{Real, Nothing}`: Young's modulus for particles.
- `ced::Union{Real, Nothing}`: Coefficient of elastic damping.
- `_model_store`: Dictionary storing model data as arrays of Float32.
"""
@kwdef mutable struct FluidisedBed
    presure_path::String
    n_probes::Int
    dump2csv::Bool
    velcfg_path::String
    plots_dir::String
    
    t::Vector{Float32} = Vector{Float32}()
    v_z::Vector{Float32} = Vector{Float32}()
    void_frac_path::Union{String, Nothing} = nothing
    
    _timeser_df::Union{DataFrame, Nothing} = nothing
    _df_store::Union{Array, Nothing} = nothing
    p_diameter::Union{Real, Nothing} = nothing
    𝜌_p::Union{Real, Nothing} = nothing
    poisson_ratio::Union{Real, Nothing} = nothing
    youngs_modulus::Union{Real, Nothing} = nothing
    ced::Union{Real, Nothing} = nothing
    _model_store = Dict{String, Array{Float32}}()
end

"""
    _find_cdfmedian(x::Vector{Float32})

Calculate the median value of a vector using its cumulative distribution function (CDF).

This function computes the median by:
1. Finding unique values and sorting them
2. Counting occurrences of each unique value
3. Computing the cumulative distribution function
4. Finding the value corresponding to CDF = 0.5

If no exact value exists where CDF = 0.5, the function performs linear interpolation
between the closest values below and above the median.

# Arguments
- `x::Vector{Float32}`: Input vector of floating-point values

# Returns
- `Float32`: The median value of the input vector

# Note
This is an internal function as indicated by the leading underscore.
"""
function _find_cdfmedian(x::Vector{Float32})
    x_unique = sort(unique(x))
    x_counts = values(countmap(x))
    csum = accumulate(+, x_counts)
    cdf = csum / csum[end]

    median_idx = findfirst(==(0.5), cdf)

    if isnothing(median_idx)
        upper_idx = findfirst(>=(0.5), cdf)
        lower_idx = findlast(<=(0.5), cdf)
        
        if isnothing(upper_idx)
            return x_unique[lower_idx]
        elseif isnothing(lower_idx)
            return x_unique[upper_idx]
        else
            intrpl8(x1, p1, x2, p2, p) = x1 + (x2 - x1) * (p - p1) / (p2 - p1)

            median_val = intrpl8(
                x_unique[lower_idx], cdf[lower_idx],
                x_unique[upper_idx], cdf[upper_idx], 0.5
            )
            return median_val
        end
    else
        return x_unique[median_idx]
    end
end

"""
    _read_probetxt(flbed::FluidisedBed)

Reads velocity probe data from a text file specified in `flbed.velcfg_path` and populates 
the time series (`flbed.t`) and vertical velocity (`flbed.v_z`) arrays in the FluidisedBed object.

The function processes each line of the file by:
1. Removing parentheses from the line
2. Splitting the line into components
3. Extracting the time value (first component) and vertical velocity (last component)
4. Adding these values to the corresponding arrays in the FluidisedBed object

# Arguments
- `flbed::FluidisedBed`: A FluidisedBed object with a valid `velcfg_path` and initialized `t` and `v_z` arrays
"""
function _read_probetxt(flbed::FluidisedBed)

    open(flbed.velcfg_path, "r") do file
        for line in eachline(file)
            line = replace(line, "(" => "")
            line = replace(line, ")" => "")
            line_splt = split(line, " ")

            t_val = parse(Float32, line_splt[1])
            push!(flbed.t, t_val)
            v_z_val = parse(Float32, line_splt[end])
            push!(flbed.v_z, v_z_val)
        end
    end
end

"""
    _calc_vel!(flbed::FluidisedBed, time_df::DataFrame)

Calculate velocities for a fluidised bed and add them to the provided DataFrame.

This function performs the following operations:
- Reads probe data if not already loaded
- Identifies constant velocity segments in the data
- Determines the time range for maximum velocity
- Maps velocities to corresponding time points in the provided DataFrame
- Adds a 'direction' column based on time position relative to maximum velocity period
- Removes rows with missing velocity values

# Arguments
- `flbed::FluidisedBed`: A fluidised bed object containing time and velocity data
- `time_df::DataFrame`: DataFrame with a 'time' column to which velocity data will be added

# Side Effects
- Modifies `time_df` by adding 'v_z' and 'direction' columns
- Removes rows with missing values using `dropmissing!`

# Notes
- The 'direction' column categorizes time points as "up", "max", or "down" based on their 
  relation to the maximum velocity time range
- If velocity data isn't already loaded in `flbed`, it's read from probe text files
"""
function _calc_vel!(flbed::FluidisedBed, time_df::DataFrame)
    function _map_direction(x)
        if x < max_vel_t1
            return "up"
        elseif (x >= max_vel_t1) & (x <= max_vel_t2)
            return "max"
        else
            return "down"
        end
    end

    bounds = []
    vel = []
    max_vel_t1 = 0.0
    max_vel_t2 = 0.0

    if !(length(flbed.t) > 0) && !(length(flbed.v_z) > 0)
        _read_probetxt(flbed)
    end
    
    for i in 1:(length(flbed.t)-1)
        if flbed.v_z[i] == flbed.v_z[i+1]
            push!(bounds, [flbed.t[i], flbed.t[i+1]])
            push!(vel, flbed.v_z[i])

            if flbed.v_z[i] == maximum(flbed.v_z)
                max_vel_t1 = flbed.t[i]
                max_vel_t2 = flbed.t[i+1]
            end
        else
        end
    end

    lb = [b[1] for b in bounds]
    ub = [b[2] for b in bounds]

    vz_arr = Vector{Union{Float32, Missing}}(undef, length(time_df.time))
    fill!(vz_arr, 0.0f0)

    for i in 1:length(bounds)
        mask = (time_df.time .>= lb[i]) .& (time_df.time .<= ub[i])
        vz_arr[mask] .= vel[i]
        if i < length(bounds)
            gap_mask = (time_df.time .>= ub[i]) .& (time_df.time .<= lb[i+1])
            vz_arr[gap_mask] .= missing
        end
    end
    time_df.v_z = vz_arr
    time_df.direction = map(_map_direction, time_df.time)
    dropmissing!(time_df)
end

"""
    _probe2df(flbed::FluidisedBed, use_slices::Bool, slice_dirn::Char, y_agg)

Converts pressure probe data from a fluidized bed simulation to a DataFrame.

# Arguments
- `flbed::FluidisedBed`: A FluidisedBed object containing metadata about the simulation.
- `use_slices::Bool`: If true, read data from VTK slice files; if false, read from CSV file.
- `slice_dirn::Char`: Direction of slices ('z' or 'y') when `use_slices` is true.
- `y_agg`: Aggregation method for y-normal slices. Options: "cdf_median", "mean", "median".

# Returns
- `DataFrame`: Contains time series pressure data for each probe.

# Behavior
- For z-normal slices: Reads VTK files for each probe at each time step and computes mean pressure.
- For y-normal slices: Reads a single VTK file and applies the specified aggregation.
- For non-slice data: Reads directly from CSV file.
- Always sorts the resulting DataFrame by time.
- Optionally writes the DataFrame to CSV if `flbed.dump2csv` is true.

# Throws
- `LoadError`: If VTK files are missing or pressure data cannot be read.
- `ErrorException`: If an unknown slice direction or aggregation method is specified.
"""
function _probe2df(flbed::FluidisedBed, use_slices::Bool, slice_dirn::Char, y_agg)

    headers = ["time"]
    for i = 0:flbed.n_probes-1
        push!(headers, "probe_$(i)")
    end

    if use_slices
        times = readdir(flbed.presure_path)
        n_times = length(times)
        time_idx = 1
        pressure_arr = zeros(Float32, n_times, flbed.n_probes)

        while time_idx <= n_times
            if slice_dirn == 'z'

                for i = 0:flbed.n_probes-1
                    dir = times[time_idx]
                    pdata_path = joinpath(
                        flbed.presure_path, dir, "p_zNormal$i.vtk"
                    )
                    if !isfile(pdata_path)
                        throw(LoadError("VTK file not found: $pdata_path"))
                    end

                    println(pdata_path)
                    p_data = read_vtk(pdata_path)
                    try
                        p_arr = p_data.point_data["p"]
                    catch e
                        throw(LoadError("Error reading pressure data from $pdata_path: $e"))
                    end
                    pressure_arr[time_idx, i+1] = mean(p_arr)
                end

            elseif slice_dirn == 'y'
                p_data = read_vtk(
                    joinpath(flbed.presure_path, "p_yNormal.vtk")
                )
                p_arr = p_data.cell_data["p"]

                if y_agg == "cdf_median"
                    p_arr = median(p_arr)
                    p_arr = Float32.(p_arr)
                elseif y_agg == "mean"
                    p_arr = mean(p_arr)
                elseif y_agg == "median"
                    p_arr = median(p_arr)
                else
                    error("Unknown aggregation method: $y_agg")
                end

                pressure_df = DataFrame(
                    time=parse.(Float32, times),
                    pressure=p_arr
                )

            else
                error("Unknown slice direction: $slice_dirn",
                    " (expected 'z' or 'y')")
            end
            time_idx += 1
        end
        df_dict = Dict{String, Vector{Float32}}("time" => parse.(Float32, times))
        
        for i in 0:flbed.n_probes-1
            df_dict["probe_$(i)"] = pressure_arr[:, i+1]
        end

        pressure_df = DataFrame(df_dict)
    else

        pressure_df = CSV.read(
            flbed.presure_path, DataFrame, skipto=8,
            delim=' ', ignorerepeated=true, header=headers
        )
    end

    sort!(pressure_df, :time)

    if flbed.dump2csv
        CSV.write(joinpath(flbed.plots_dir, "pressure.csv"), pressure_df)
    end

    return pressure_df
end

"""
    plot_pressure(
        flbed::FluidisedBed;
        x_var::String="velocity",
        png_name=nothing,
        use_slices::Bool=true,
        slice_dirn::Char='z',
        y_agg=nothing
    )

Generate plots of pressure data from a fluidised bed simulation.

# Arguments
- `flbed::FluidisedBed`: The fluidised bed object containing simulation data.
- `x_var::String="velocity"`: The x-axis variable, either "time" or "velocity".
- `png_name=nothing`: Custom name for the output PNG file. If `nothing`, a default name is generated.
- `use_slices::Bool=true`: If `true`, use slice data; otherwise, use probe data.
- `slice_dirn::Char='z'`: Direction for slicing ('z' supported; 'y' not fully supported).
- `y_agg=nothing`: Aggregation function for y-direction data if needed.

# Description
This function creates pressure plots for fluidised bed simulations with two main options:
- Time series plots: When `x_var="time"`, plots pressure vs. time for each probe.
- Velocity plots: When `x_var="velocity"`, plots pressure vs. velocity showing hysteresis 
  between increasing and decreasing velocity conditions.

The function handles data preparation, grouping, and visualization with appropriate labels 
and styling. Results are saved as PNG files in the `plots_dir` directory of the fluidised bed object.

# Notes
- Y-direction slicing is currently not fully supported.
- For velocity plots, separate curves are drawn for increasing and decreasing velocity.
- The function caches processed data in `flbed._timeser_df` to avoid redundant calculations.
"""
function plot_pressure(
    flbed::FluidisedBed;
    x_var::String="velocity",
    png_name=nothing,
    use_slices::Bool=true,
    slice_dirn::Char='z',
    y_agg=nothing)

    if slice_dirn == "y"
        @warn("Plotting pressure in y direction is not supported yet.")
    end

    plot_suffix = use_slices ? "_slice" : "_probes"

    if (isnothing(flbed._timeser_df)) || (flbed._df_store != [slice_dirn, "pressure"])
        flbed._timeser_df = _probe2df(flbed, use_slices, slice_dirn, y_agg)
        flbed._df_store = [slice_dirn, "pressure"]
    end

    if x_var == "time"
        x_data = flbed._timeser_df.time
        if "direction" in names(flbed._timeser_df)
            y_df = flbed._timeser_df[!, Not([:time, :direction, :v_z])]
        else
            y_df = flbed._timeser_df[!, Not(:time)]
        end

        println(y_df)
        y_data = Matrix(y_df)

        label_names = ["Probe $(i)" for i in 0:flbed.n_probes-1]
        label_matrix = reshape(label_names, 1, length(label_names))

        plot(
            x_data, y_data, label=label_matrix,
            xlabel="Time (s)", ylabel="Pressure (Pa)",
            title="Pressure vs Time"
        )

    elseif x_var == "velocity"
        _calc_vel!(flbed, flbed._timeser_df)
        pressure_gdf = groupby(flbed._timeser_df, [:direction, :v_z])
        cols_to_average = valuecols(pressure_gdf)
        vel_plot_df = combine(pressure_gdf, cols_to_average .=> mean)

        vel_up = filter(:direction => in(["up", "max"]), vel_plot_df)
        sort!(vel_up, :v_z)

        vel_down = filter(:direction => in(["down", "max"]), vel_plot_df)
        sort!(vel_down, :v_z)

        println(vel_up)
        println(vel_down)

        if slice_dirn == 'z'
            vel_data_up = vel_up.v_z
            p_data_up = Matrix(
                vel_up[!, Not([:v_z, :direction])]
            )
            vel_down_data = vel_down.v_z
            p_data_down = Matrix(
                vel_down[!, Not([:v_z, :direction])]
            )
            label_names = Array{String}(undef, flbed.n_probes)
            for i in 0:flbed.n_probes-1
                label_names[i+1] = "Probe $(i)"
            end

            for i in 0:flbed.n_probes-1
                plot!(
                    vel_data_up, p_data_up[:, i+1],
                    label="$(label_names[i+1]) (Increasing Velocity)", xlabel="V    elocity (m/s)", dashed=false, color=i+1, marker=:circle
                )
                plot!(
                    vel_down_data, p_data_down[:, i+1],
                    label="$(label_names[i+1]) (Decreasing Velocity)", xlabel="Velocity (m/s)", linestyle=:dash, color=i+1, marker=:circle
                )
            end

        else
            plot(
                vel_up.v_z, vel_down.pressure,
                label="Increasing Velocity",
                color=:blue, xlabel="Velocity (m/s)",
                ylabel="Pressure (Pa)", dashed=false, marker=:circle
            )
            plot!(
                vel_down.v_z, vel_down.pressure,
                label="Decreasing Velocity",
                color=:red, xlabel="Velocity (m/s)",
                ylabel="Pressure (Pa)", marker=:circle, dashed=true
            )
        end
    end
    if isnothing(png_name)
        png_name = "pressure_$(x_var)_$(slice_dirn)_$(plot_suffix)"
    elseif !endswith(png_name, ".png")
        png_name *= ".png"
    end

    savefig(joinpath(flbed.plots_dir, "pressure_$(x_var)_$(plot_suffix).png"))
end


"""
    _read_vf(flbed::FluidisedBed, post_dir::String, slice_dirn::Char)

Reads void fraction data from VTK files in the specified post-processing directory.

# Arguments
- `flbed::FluidisedBed`: A fluidised bed object containing simulation parameters, including `n_probes`.
- `post_dir::String`: Path to the directory containing post-processing data organized in time-stamped folders.
- `slice_dirn::Char`: Direction of the slice, either 'z' or 'y'.
  - 'z': Reads void fraction data from multiple probes along the z-axis.
  - 'y': Reads void fraction data from a single y-normal slice and calculates the CDF median.

# Returns
- `DataFrame`: A sorted DataFrame containing time and void fraction data. 
  - For 'z' direction: columns include 'time' and 'probe_0', 'probe_1', etc.
  - For 'y' direction: columns include 'time' and 'void_fraction'.

# Throws
- `LoadError`: If VTK files are not found or void fraction data cannot be read.
- `ErrorException`: If an unsupported slice direction is specified.
"""
function _read_vf(flbed::FluidisedBed, post_dir::String, slice_dirn::Char)
    
    times = readdir(post_dir)
    n_times = length(times)

    if slice_dirn == 'y'
        voidfrac_arr = zeros(Float32, n_times)
    elseif slice_dirn == 'z'
        voidfrac_arr = zeros(Float32, n_times, flbed.n_probes)
    else
        error("Unknown slice direction: $slice_dirn (expected 'z' or 'y')")
    end

    time_idx = 1
    while time_idx <= n_times
        dir = times[time_idx]
        if slice_dirn == 'z'
            for i = 0:flbed.n_probes-1
                vfdata_path = joinpath(
                    post_dir, dir, "voidfraction_zNormal$i.vtk"
                )
                if !isfile(vfdata_path)
                    throw(LoadError("VTK file not found: $vfdata_path"))
                end

                vf_data = read_vtk(vfdata_path)
                try
                    vf_arr = Float32.(vf_data.point_data["voidfraction"])
                catch e
                    throw(LoadError("Error reading void fraction data from $vfdata_path: $e"))
                end
                voidfrac_arr[time_idx, i+1] = mean(vf_arr)
            end

        elseif slice_dirn == 'y'
            vf_data = read_vtk(
                joinpath(post_dir, dir, "voidfraction_yNormal.vtk")
            )

            vf_arr = Float32.(vf_data.point_data["voidfraction"])
            voidfrac_arr[time_idx] = _find_cdfmedian(vf_arr)
            println("Void fraction at time $(times[time_idx]): $(voidfrac_arr[time_idx])")

        else
            error("Unknown slice direction: $slice_dirn (expected 'z')")
        end
        time_idx += 1
    end

    if slice_dirn == 'z'
        vdf_dict = Dict{String, Vector{Float32}}("time" => parse.(Float32, times))
        for i in 0:flbed.n_probes-1
            vdf_dict["probe_$(i)"] = voidfrac_arr[:, i+1]
        end
        voidfrac_df = DataFrame(vdf_dict)
    else
        voidfrac_df = DataFrame(
            time=parse.(Float32, times),
            void_fraction=voidfrac_arr
        )
    end
    sort!(voidfrac_df, :time)
    return voidfrac_df
end

"""
    plot_voidfrac(
        flbed::FluidisedBed;
        slice_dirn::Char='y',
        x_var::String="velocity",
        png_name=nothing,
    )

Generate plots of void fraction data for a fluidised bed simulation.

# Arguments
- `flbed::FluidisedBed`: The fluidised bed object containing simulation data.
- `slice_dirn::Char='y'`: The direction of the cutting plane ('x', 'y', or 'z').
- `x_var::String="velocity"`: The x-axis variable for plotting. Options are "time" or "velocity".
- `png_name=nothing`: Custom name for the output PNG file. If `nothing`, a default name will be generated.

# Details
- When `x_var="time"`, generates a plot of void fraction vs time.
- When `x_var="velocity"`, generates plots showing hysteresis effects with separate curves for increasing and decreasing velocity.
- For `slice_dirn='z'`, multiple probe points are plotted separately.
- Saves the plot to `flbed.plots_dir` and optionally exports the data to a CSV file if `flbed.dump2csv` is `true`.

# Note
The function checks if the void fraction data has already been loaded and stored in `flbed._timeser_df`
to avoid redundant reading of data.
"""
function plot_voidfrac(
    flbed::FluidisedBed;
    slice_dirn::Char='y',
    x_var::String="velocity",
    png_name=nothing,
)   

    if flbed.void_frac_path == nothing
        post_dir = "CFD/postProcessing/cuttingPlane/"
    else
        post_dir = flbed.void_frac_path
    end

    if flbed._df_store != [slice_dirn, "void_fraction"]
        flbed._df_store = [slice_dirn, "void_fraction"]
        flbed._timeser_df = _read_vf(flbed, post_dir, slice_dirn)
    end
    println(flbed._timeser_df)
    
    if x_var == "time"
        x_data = flbed._timeser_df.time
        y_data = Matrix(flbed._timeser_df[!, Not(:time)])
        if slice_dirn == 'z'
            label_names = ["Probe $(i)" for i in 0:flbed.n_probes-1]
        else
            label_names = ["Void Fraction"]
        end
        label_matrix = reshape(label_names, 1, length(label_names))

        plot(
            x_data, y_data, label=label_matrix,
            xlabel="Time (s)", ylabel="Void Fraction",
            title="Void Fraction vs Time"
        )
    elseif x_var == "velocity"
        _calc_vel!(flbed, flbed._timeser_df)
        voidfrac_gdf = groupby(flbed._timeser_df, [:direction, :v_z])
        cols_to_average = valuecols(voidfrac_gdf)
        vel_plot_df = combine(voidfrac_gdf, cols_to_average .=> mean)

        vel_up = filter(:direction => in(["up", "max"]), vel_plot_df)
        sort!(vel_up, :v_z)
        vel_down = filter(:direction => in(["down", "max"]), vel_plot_df)
        sort!(vel_down, :v_z)

        if slice_dirn == 'z'
            vel_data_up = vel_up.v_z
            vf_data_up = Matrix(
                vel_up[!, Not([:v_z, :direction])]
            )
            vel_down_data = vel_down.v_z
            vf_data_down = Matrix(
                vel_down[!, Not([:v_z, :direction])]
            )
            label_names = Array{String}(undef, flbed.n_probes)
            for i in 0:flbed.n_probes-1
                label_names[i+1] = "Probe $(i)"
            end

            for i in 0:flbed.n_probes-1
                plot!(
                    vel_data_up, vf_data_up[:, i+1],
                    label="$(label_names[i+1]) (Increasing Velocity)", xlabel="Velocity (m/s)", dashed=false, color=i+1, marker=:circle
                )
                plot!(
                    vel_down_data, vf_data_down[:, i+1],
                    label="$(label_names[i+1]) (Decreasing Velocity)", xlabel="Velocity (m/s)", linestyle=:dash, color=i+1, marker=:circle
                )
            end
        else
            plot(
                vel_up.v_z, vel_up.void_fraction_mean,
                label="Increasing Velocity",
                color=:blue, xlabel="Velocity (m/s)",
                ylabel="Void Fraction", dashed=false, marker=:circle
            )
            plot!(
                vel_down.v_z, vel_down.void_fraction_mean,
                label="Decreasing Velocity",
                color=:red, xlabel="Velocity (m/s)",
                ylabel="Void Fraction", marker=:circle, linestyle=:dash
            )
        end
    end
    if isnothing(png_name)
        png_name = "void_fraction_$(x_var)_$(slice_dirn)"
    elseif !endswith(png_name, ".png")
        png_name *= ".png"
    end
    savefig(joinpath(flbed.plots_dir, png_name))
    if flbed.dump2csv
        CSV.write(joinpath(flbed.plots_dir, "void_fraction.csv"), flbed._timeser_df)
    end
end
end
