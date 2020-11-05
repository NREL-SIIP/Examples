using Pkg
Pkg.status()

using SIIPExamples;
using PowerSystems;
using D3TypeTrees;

print_struct(Bus)

TypeTree(Forecast)

print_struct(Deterministic)

print_struct(System)

BASE_DIR = abspath(joinpath(dirname(Base.find_package("PowerSystems")), ".."))
include(joinpath(BASE_DIR, "test", "data_5bus_pu.jl")) #.jl file containing 5-bus system data
nodes_5 = nodes5() # function to create 5-bus buses

sys = System(
    100.0,
    nodes_5,
    vcat(thermal_generators5(nodes_5), renewable_generators5(nodes_5)),
    loads5(nodes_5),
    branches5(nodes_5),
)

@show get_component(Bus, sys, "nodeA")
@show get_component(Line, sys, "1")

get_components(Bus, sys) |> collect

get_components(Branch, sys) |> collect

bus1 = get_component(Bus, sys, "nodeA")
@show get_name(bus1);
@show get_magnitude(bus1);

loads = collect(get_components(PowerLoad, sys))
for (l, ts) in zip(loads, load_timeseries_DA[2])
    add_time_series!(sys, l, Deterministic("activepower", Dict(TimeSeries.timestamp(load_timeseries_DA[2][1])[1] => ts)))
end

get_forecast_initial_times(sys)

@show labels = IS.get_time_series_names(Deterministic, loads[1])

@show initial_times = IS.get_initial_timestamp(get_time_series(Deterministic, loads[1], labels[1]));

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

