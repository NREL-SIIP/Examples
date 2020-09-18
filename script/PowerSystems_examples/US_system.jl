# # Creating a `System` representing the entire U.S.

# **Originally Contributed by**: Clayton Barrows

# ## Introduction

# This example demonstrates how to assemble a `System` representing the entire U.S. using
# [PowerSystems.jl](https://github.com/NREL-SIIP/powersystems.jl) and the data assembled by
# [Xu, et. al.](https://arxiv.org/abs/2002.06155). We'll use the same tabular data parsing
# capability [demonstrated on the RTS-GMLC dataset](../../notebook/PowerSystems_examples/parse_tabulardata.ipynb).

# ### Environment
# This notebook depends on the SIIPExamples.jl environment which is loaded by default

#nb using Pkg
#nb Pkg.status()

# ### Dependencies
using SIIPExamples
using PowerSystems
using TimeSeries
using Dates
using TimeZones
using DataFrames
using CSV

# ### Fetch Data
# PowerSystems.jl links to some test data that is suitable for this example.
# Let's download the test data
@info "downloading data..."
datadir = joinpath(dirname(dirname(pathof(SIIPExamples))), "US-System")
siip_data = joinpath(datadir, "SIIP")
if !isdir(datadir)
    mkdir(datadir)
    mkdir(siip_data)
    tempfilename =
        download("https://zenodo.org/record/3753177/files/USATestSystem.zip?download=1")
    SIIPExamples.unzip(SIIPExamples.os, tempfilename, datadir)
end

config_dir = joinpath(
    dirname(dirname(pathof(SIIPExamples))),
    "script",
    "PowerSystems_examples",
    "US_config",
)

# ### Data Formatting
# This is a big dataset. Typically one would only want to include one of the interconnects
# available. Lets use Texas to start. You can set `interconnect = nothing` if you want everything.
interconnect = "Texas"
timezone = FixedTimeZone("UTC-6")
initial_time = ZonedDateTime(DateTime("2016-01-01T00:00:00"), timezone)

# There are a few minor incompatibilities between the data and the supported tabular data
# format. We can resolve those here.
#
# First, PowerSystems.jl only supports parsing piecewise linear generator costs from tabular
# data. So, we can sample the quadratic polynomial cost curves and provide PWL points.
@info "formatting data ..."
!isnothing(interconnect) && @info "filtering data to include $interconnect ..."
gen = DataFrame(CSV.File(joinpath(datadir, "plant.csv")))
!isnothing(interconnect) && filter!(row -> row[:interconnect] == interconnect, gen)
gencost = DataFrame(CSV.File(joinpath(datadir, "gencost.csv")))
gen = innerjoin(gen, gencost, on = :plant_id, makeunique = true, validate = (false, false))

function make_pwl(gen::DataFrame, traunches = 2)
    output_pct_cols = ["output_point_" * string(i) for i = 0:traunches]
    hr_cols = ["heat_rate_incr_" * string(i) for i = 1:traunches]
    pushfirst!(hr_cols, "heat_rate_avg_0")
    pwl = DataFrame(repeat([Float64], 6), Symbol.(vcat(output_pct_cols, hr_cols)))
    for row in eachrow(gen)
        traunch_len = (1.0 - row.Pmin / row.Pmax) / traunches
        pct = [row.Pmin / row.Pmax + i * traunch_len for i = 0:traunches]
        #c(pct) = pct * row.Pmax * (row.GenIOB + row.GenIOC^2 + row.GenIOD^3)
        c(pct) = pct * row.Pmax * (row.c1 + row.c2^2) + row.c0 #this formats the "c" columns to hack the heat rate parser in PSY
        hr = [c(pct[1])]
        [push!(hr, c(pct[i+1]) - hr[i]) for i = 1:traunches]
        push!(pwl, vcat(pct, hr))
    end
    return hcat(gen, pwl)
end

gen = make_pwl(gen);

gen[!, "fuel_price"] .= 1000.0;  #this formats the "c" columns to hack the heat rate parser in PSY

# There are some incomplete aspects of this dataset. Here, I've assigned some approximate
# minimum up/down times, and some minor adjustments to categories. There are better
# ways to do this, but this works for this script...
gen[:, :unit_type] .= "OT"
gen[:, :min_up_time] .= 0.0
gen[:, :min_down_time] .= 0.0
gen[:, :ramp_30] .= gen[:, :ramp_30] ./ 30.0 # we need ramp rates in MW/min
[
    gen[gen.type.=="wind", col] .= ["Wind", 0.0, 0.0][ix]
    for (ix, col) in enumerate([:unit_type, :min_up_time, :min_down_time])
]
[
    gen[gen.type.=="solar", col] .= ["PV", 0.0, 0.0][ix]
    for (ix, col) in enumerate([:unit_type, :min_up_time, :min_down_time])
]
[
    gen[gen.type.=="hydro", col] .= ["HY", 0.0, 0.0][ix]
    for (ix, col) in enumerate([:unit_type, :min_up_time, :min_down_time])
]
[
    gen[gen.type.=="ng", col] .= [4.5, 8][ix]
    for (ix, col) in enumerate([:min_up_time, :min_down_time])
]
[
    gen[gen.type.=="coal", col] .= [24, 48][ix]
    for (ix, col) in enumerate([:min_up_time, :min_down_time])
]
[
    gen[gen.type.=="nuclear", col] .= [72, 72][ix]
    for (ix, col) in enumerate([:min_up_time, :min_down_time])
]

# At the moment, PowerSimulations can't do unit commitment with generators that have Pmin = 0.0
idx_zero_pmin = [
    g.type in ["ng", "coal", "hydro", "nuclear"] && g.Pmin <= 0
    for g in eachrow(gen[:, [:type, :Pmin]])
]
gen[idx_zero_pmin, :Pmin] = gen[idx_zero_pmin, :Pmax] .* 0.05

gen[:, :name] = "gen" .* string.(gen.plant_id)
CSV.write(joinpath(siip_data, "gen.csv"), gen)

# Let's also merge the zone.csv with the bus.csv and identify bus types
bus = DataFrame(CSV.File(joinpath(datadir, "bus.csv")))
!isnothing(interconnect) && filter!(row -> row[:interconnect] == interconnect, bus)
zone = DataFrame(CSV.File(joinpath(datadir, "zone.csv")))
bus = leftjoin(bus, zone, on = :zone_id)
int2bustype(b) =
    replace(split(string(PowerSystems.BusTypes.BusType(b)), ".")[end], "]" => "")
bus.bustype = int2bustype.(bus.type)
bus.name = "bus" .* string.(bus.bus_id)
CSV.write(joinpath(siip_data, "bus.csv"), bus)

# We need branch names as strings
branch = DataFrame(CSV.File(joinpath(datadir, "branch.csv")))
branch = join(
    branch,
    DataFrames.rename!(bus[:, [:bus_id, :baseKV]], [:from_bus_id, :from_baseKV]),
    on = :from_bus_id,
    kind = :left,
)
branch = join(
    branch,
    DataFrames.rename!(bus[:, [:bus_id, :baseKV]], [:to_bus_id, :to_baseKV]),
    on = :to_bus_id,
    kind = :left,
)
!isnothing(interconnect) && filter!(row -> row[:interconnect] == interconnect, branch)
branch.name = "branch" .* string.(branch.branch_id)
branch.tr_ratio = branch.from_baseKV ./ branch.to_baseKV
CSV.write(joinpath(siip_data, "branch.csv"), branch)

# The PowerSystems parser expects the files to be named a certain way.
# And, we need a "control_mode" column in dc-line data
dcbranch = DataFrame(CSV.File(joinpath(datadir, "dcline.csv")))
!isnothing(interconnect) && filter!(row -> row[:from_bus_id] in bus.bus_id, dcbranch)
!isnothing(interconnect) && filter!(row -> row[:to_bus_id] in bus.bus_id, dcbranch)
dcbranch.name = "dcbranch" .* string.(dcbranch.dcline_id)
CSV.write(joinpath(siip_data, "dc_branch.csv"), dcbranch)

# ### We need to create a reference for where to get timeseries data for each component.
timeseries = []
ts_csv = ["wind", "solar", "hydro", "demand"]
plant_ids = Symbol.(string.(gen.plant_id))
for f in ts_csv
    @info "formatting $f.csv ..."
    csvpath = joinpath(siip_data, f * ".csv")
    csv = DataFrame(CSV.File(joinpath(datadir, f * ".csv")))
    (category, name_prefix, label) = f == "demand" ? ("Area", "", "get_max_active_power") :
        ("Generator", "gen", "get_max_active_power")
    if !(:DateTime in names(csv))
        DataFrames.rename!(
            csv,
            (names(csv)[occursin.("UTC", String.(names(csv)))][1] => :DateTime),
        )
        #The timeseries data is in UTC, this converts it to a fixed UTC offset
        csv.DateTime =
            ZonedDateTime.(
                DateTime.(csv.DateTime, "yyyy-mm-dd HH:MM:SS"),
                timezone,
                from_utc = true,
            )
        delete!(csv, csv.DateTime .< initial_time)
        csv.DateTime = Dates.format.(csv.DateTime, "yyyy-mm-ddTHH:MM:SS")
    end
    device_names = f == "demand" ? unique(bus.zone_name) : gen.name
    for id in names(csv)
        colname = id
        if f == "demand"
            if Symbol(id) in Symbol.(zone.zone_id)
                colname = Symbol(zone[Symbol.(zone.zone_id).==Symbol(id), :zone_name][1])
                DataFrames.rename!(csv, (id => colname))
            end
            sf = sum(bus[string.(bus.zone_id).==id, :Pd])
        else
            if Symbol(id) in plant_ids
                colname = Symbol(gen[Symbol.(gen.plant_id).==Symbol(id), :name][1])
                DataFrames.rename!(csv, (id => colname))
            end
            sf = maximum(csv[:, colname]) == 0.0 ? 1.0 : "Max"
        end
        if String(colname) in device_names
            push!(
                timeseries,
                Dict(
                    "simulation" => "DA",
                    "category" => category,
                    "type" => "Deterministic",
                    "label" => replace(label, "get_"=>""),
                    "component_name" => String(colname),
                    "scaling_factor_multiplier" => label,
                    "scaling_factor_multiplier_module" => "PowerSystems",
                    "normalization_factor" => sf,
                    "data_file" => csvpath,
                ),
            )
        end
    end
    CSV.write(csvpath, csv)
end

timeseries_pointers = joinpath(siip_data, "timeseries_pointers.json")
open(timeseries_pointers, "w") do io
    JSON2.write(io, timeseries)
end

# ### The tabular data format relies on a folder containing `*.csv` files and `.yaml` files
# describing the column names of each file in PowerSystems terms, and the PowerSystems
# data type that should be created for each generator type. The respective "us_decriptors.yaml"
# and "US_generator_mapping.yaml" files have already been tailored to this dataset.
@info "parsing csv files..."
rawsys = PowerSystems.PowerSystemTableData(
    siip_data,
    100.0,
    joinpath(config_dir, "us_descriptors.yaml"),
    generator_mapping_file = joinpath(config_dir, "US_generator_mapping.yaml"),
)

# ### Create a `System`
# Next, we'll create a `System` from the `rawsys` data. Since a `System` is predicated on a
# forecast resolution and the `rawsys` data includes both 5-minute and 1-hour resolution
# time series, we also need to specify which time series we want to include in the `System`.
# The `forecast_resolution` kwarg filters to only include time series with a matching resolution.

@info "creating System"
sys = System(rawsys; config_path = joinpath(config_dir, "us_system_validation.json"));
sys

# This all took reasonably long, so we can save our `System` using the serialization
# capability included with PowerSystems.jl:
to_json(sys, joinpath(siip_data, "sys.json"), force = true)
