module CurvePlots

using VTKDataIO
using StatsBase
using Statistics
using DataFrames
using CSV
using Plots

export FluidisedBed, plot_pressure

@kwdef mutable struct FluidisedBed
    presure_path::String
    n_probes::Int
    dump2csv::Bool
    velcfg_path::String
    plots_dir::String
    t::Vector{Float32} = Vector{Float32}()
    v_z::Vector{Float32} = Vector{Float32}()
    
    _timeser_df::Union{DataFrame, Nothing} = nothing
    _df_store::Union{Array, Nothing} = nothing
end

function _find_cdfmedian(x::Vector{Float32})
    x_unique = unique(x)
    x_counts = [count(==(element), x) for element in unique(x)]
    csum = cumsum(x_counts)
    cdf = csum / csum[end]

    pct_50 = percentile(cdf, 0.5)
    median_idx = findfirst(==(pct_50), cdf)

    return x_unique[median_idx]
end

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

    vz_arr = zeros(Float32, length(time_df.time))

    for i in 1:length(bounds)
        mask = (time_df.time .>= lb[i]) .& (time_df.time .<= ub[i])
        vz_arr[mask] .= vel[i]
        if i < length(bounds)
            gap_mask = (time_df.time .>= ub[i]) .& (time_df.time .<= lb[i+1])
            vz_arr[gap_mask] .= NaN
        end
    end
    time_df.v_z = vz_arr
    time_df.direction = map(_map_direction, time_df.time)
end

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
                        throw("VTK file not found: $pdata_path")
                    end

                    println(pdata_path)
                    p_data = read_vtk(pdata_path)
                    try
                        p_arr = p_data.point_data["p"]
                    catch e
                        throw("Error reading pressure data from $pdata_path: $e")
                    end
                    pressure_arr[time_idx, i+1] = mean(p_arr)
                end

            elseif slice_dirn == 'y'
                p_data = read_vtk(
                    joinpath(flbed.presure_path, "p_yNormal.vtk")
                )
                p_arr = p_data.cell_data["p"]

                if y_agg == "cdf_median"
                    p_arr = _find_cdfmedian(p_arr)
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
end
