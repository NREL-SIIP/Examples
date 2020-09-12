# This script is intended to be a dependency for the hydropower example.
using PowerSystems
using Statistics
using Dates
using TimeSeries

# We can use some of the simple data that's been assembled for testing PowerSimulations.
PowerSystems.download(PowerSystems.TestData; branch = "master") # *note* add `force=true` to get a fresh copy
data_dir = joinpath(dirname(dirname(pathof(PowerSystems))), "data");

rawsys = PowerSystems.PowerSystemTableData(
    joinpath(data_dir, "5-bus-hydro"),
    100.0,
    joinpath(data_dir, "5-bus-hydro", "user_descriptors.yaml");
    generator_mapping_file = joinpath(data_dir, "5-bus-hydro", "generator_mapping.yaml"),
)

demo_sys = System(
    rawsys,
    timeseries_metadata_file = joinpath(
        data_dir,
        "5-bus-hydro",
        "timeseries_pointers_da.json",
    ),
    time_series_in_memory = true,
)

demo_rt_sys = System(
    rawsys,
    timeseries_metadata_file = joinpath(
        data_dir,
        "forecasts",
        "5bus_ts",
        "7day",
        "timeseries_pointers_rt_7day.json",
    ),
    time_series_in_memory = true,
)

demo_da_sys =  System(
    rawsys,
    timeseries_metadata_file = joinpath(
        data_dir,
        "forecasts",
        "5bus_ts",
        "7day",
        "timeseries_pointers_da_7day.json"),
    time_series_in_memory = true,
    )

demo_wk_sys =  System(
    rawsys,
    timeseries_metadata_file = joinpath(
        data_dir,
        "forecasts",
        "5bus_ts",
        "7day",
        "timeseries_pointers_wk_7day.json"),
    time_series_in_memory = true,
    )
### helper to convert ts to daily resolution
#=
n = 7# number of days
# And an hourly system with longer time scales
MultiDay = collect(
    get_forecast_initial_times(demo_da_sys)[1]:Hour(24):(get_forecast_initial_times(
        demo_da_sys,
    )[1]+Day(n-1)),
);

fname = joinpath(data_dir,"forecasts","5bus_ts","7day","wk","wk_wind.csv")
df = CSV.read(fname, DataFrame)
ta = []
for col in names(df)
    if col != "DateTime"
        fc = TimeArray(df.DateTime, df[:,col])
        fc_values = [mean(values(fc[Date.(timestamp(fc)).==t])) for t in MultiDay]

        wk_fc = TimeArray(MultiDay, fc_values)
        push!(ta,wk_fc)
    end
end
table = reduce(hcat, ta)
TimeSeries.rename!(table, Symbol.(names(df)[names(df) .!= "DateTime"]))
CSV.write(fname, table)
=#