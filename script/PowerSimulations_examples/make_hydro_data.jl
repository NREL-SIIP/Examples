# This script is intended to be a dependency for the hydropower example. #src
using InfrastructureSystems #src
const IS = InfrastructureSystems #src
using PowerSystems #src
const PSY = PowerSystems #src
using PowerSimulations #src
const PSI = PowerSimulations #src
 #src
# We can use some of the simple data that's been assembled for testing PowerSimulations. #src
include(joinpath(pathof(PSI), "../../test/test_utils/get_test_data.jl")) #src
 #src
# Additionally, let's add two hydro generators. One of each type supported by PowerSystems. #src
hydro_generators5(nodes5) = [ #src
    HydroFix("HydroFix", #src
        true, #src
        nodes5[2], #src
        0.0, #src
        0.0, #src
        TechHydro(0.600, #src
            PowerSystems.HY, #src
            (min = 0.0, max = 60.0), #src
            (min = 0.0, max = 60.0), #src
            nothing, #src
            nothing, #src
        ), #src
    ), #src
    HydroDispatch("HydroDispatch", #src
        true, #src
        nodes5[3], #src
        0.0, #src
        0.0, #src
        TechHydro(0.600, #src
            PowerSystems.HY, #src
            (min = 0.0, max = 60.0), #src
            (min = 0.0, max = 60.0), #src
            (up = 10.0, down = 10.0), #src
            nothing, #src
        ), #src
        TwoPartCost(15.0, 0.0), #src
        10.0, #src
        2.0, #src
        5.0, #src
    ), #src
]; #src
 #src
# We can add some random time series information too. #src
 #src
hydro_timeseries_DA = [ #src
    [TimeSeries.TimeArray(DayAhead, wind_ts_DA)], #src
    [TimeSeries.TimeArray(DayAhead + Day(1), wind_ts_DA)], #src
]; #src
 #src
 #src
hydro_timeseries_RT = [ #src
    [TimeArray(RealTime, repeat(wind_ts_DA, inner = 12))], #src
    [TimeArray(RealTime + Day(1), repeat(wind_ts_DA, inner = 12))], #src
]; #src
 #src
# Now we can create a system with hourly resolution and add forecast data to it. #src
 #src
c_sys5_hy = System(nodes, #src
    vcat(thermal_generators5_uc_testing(nodes), #src
        hydro_generators5(nodes), #src
        renewable_generators5(nodes), #src
    ), #src
    loads5(nodes), #src
    branches5(nodes), #src
    nothing, #src
    100.0, #src
    nothing, #src
    nothing, #src
) #src
 #src
for t in 1:2 #src
    for (ix, l) in enumerate(get_components(PowerLoad, c_sys5_hy)) #src
        add_forecast!(c_sys5_hy, l, #src
            Deterministic("get_maxactivepower", load_timeseries_DA[t][ix]), #src
        ) #src
    end #src
    for (ix, h) in enumerate(get_components(HydroDispatch, c_sys5_hy)) #src
        add_forecast!(c_sys5_hy, #src
            h, #src
            Deterministic("get_rating", hydro_timeseries_DA[t][ix]), #src
        ) #src
    end #src
    for (ix, h) in enumerate(get_components(HydroDispatch, c_sys5_hy)) #src
        add_forecast!(c_sys5_hy, #src
            h, #src
            Deterministic("get_storage_capacity", hydro_timeseries_DA[t][ix]), #src
        ) #src
    end #src
    for (ix, h) in enumerate(get_components(HydroDispatch, c_sys5_hy)) #src
        add_forecast!(c_sys5_hy, #src
            h, #src
            Deterministic("get_inflow", hydro_timeseries_DA[t][ix]), #src
        ) #src
    end #src
    for (ix, h) in enumerate(get_components(HydroFix, c_sys5_hy)) #src
        add_forecast!(c_sys5_hy, h, Deterministic("get_rating", hydro_timeseries_DA[t][ix])) #src
    end #src
    for (ix, r) in enumerate(get_components(RenewableGen, c_sys5_hy)) #src
        add_forecast!(c_sys5_hy, r, Deterministic("get_rating", ren_timeseries_DA[t][ix])) #src
    end #src
    for (ix, i) in enumerate(get_components(InterruptibleLoad, c_sys5_hy)) #src
        add_forecast!(c_sys5_hy, #src
            i, #src
            Deterministic("get_maxactivepower", Iload_timeseries_DA[t][ix]), #src
        ) #src
    end #src
end #src
 #src
# And we can make system with 5-minute resolution #src
 #src
c_sys5_hy_ed = System(nodes, #src
    vcat(thermal_generators5_uc_testing(nodes), #src
        hydro_generators5(nodes), #src
        renewable_generators5(nodes), #src
    ), #src
    vcat(loads5(nodes), interruptible(nodes)), #src
    branches5(nodes), #src
    nothing, #src
    100.0, #src
    nothing, #src
    nothing, #src
) #src
for t in 1:2 #src
    for (ix, l) in enumerate(get_components(PowerLoad, c_sys5_hy_ed)) #src
        ta = load_timeseries_DA[t][ix] #src
        for i in 1:length(ta) # loop over hours #src
            ini_time = timestamp(ta[i]) # get the hour #src
            data = when(load_timeseries_RT[t][ix], hour, hour(ini_time[1])) # get the subset ts for that hour #src
            add_forecast!(c_sys5_hy_ed, l, Deterministic("get_maxactivepower", data)) #src
        end #src
    end #src
    for (ix, l) in enumerate(get_components(HydroDispatch, c_sys5_hy_ed)) #src
        ta = hydro_timeseries_DA[t][ix] #src
        for i in 1:length(ta) # loop over hours #src
            ini_time = timestamp(ta[i]) # get the hour #src
            data = when(hydro_timeseries_RT[t][ix], hour, hour(ini_time[1])) # get the subset ts for that hour #src
            add_forecast!(c_sys5_hy_ed, l, Deterministic("get_rating", data)) #src
        end #src
    end #src
    for (ix, l) in enumerate(get_components(RenewableGen, c_sys5_hy_ed)) #src
        ta = load_timeseries_DA[t][ix] #src
        for i in 1:length(ta) # loop over hours #src
            ini_time = timestamp(ta[i]) # get the hour #src
            data = when(load_timeseries_RT[t][ix], hour, hour(ini_time[1])) # get the subset ts for that hour #src
            add_forecast!(c_sys5_hy_ed, l, Deterministic("get_rating", data)) #src
        end #src
    end #src
    for (ix, l) in enumerate(get_components(HydroDispatch, c_sys5_hy_ed)) #src
        ta = hydro_timeseries_DA[t][ix] #src
        for i in 1:length(ta) # loop over hours #src
            ini_time = timestamp(ta[i]) # get the hour #src
            data = when(hydro_timeseries_RT[t][ix], hour, hour(ini_time[1])) # get the subset ts for that hour #src
            add_forecast!(c_sys5_hy_ed, l, Deterministic("get_storage_capacity", data)) #src
        end #src
    end #src
    for (ix, l) in enumerate(get_components(InterruptibleLoad, c_sys5_hy_ed)) #src
        ta = load_timeseries_DA[t][ix] #src
        for i in 1:length(ta) # loop over hours #src
            ini_time = timestamp(ta[i]) # get the hour #src
            data = when(load_timeseries_RT[t][ix], hour, hour(ini_time[1])) # get the subset ts for that hour #src
            add_forecast!(c_sys5_hy_ed, l, Deterministic("get_maxactivepower", data)) #src
        end #src
    end #src
    for (ix, l) in enumerate(get_components(HydroFix, c_sys5_hy_ed)) #src
        ta = hydro_timeseries_DA[t][ix] #src
        for i in 1:length(ta) # loop over hours #src
            ini_time = timestamp(ta[i]) # get the hour #src
            data = when(hydro_timeseries_RT[t][ix], hour, hour(ini_time[1])) # get the subset ts for that hour #src
            add_forecast!(c_sys5_hy_ed, l, Deterministic("get_rating", data)) #src
        end #src
    end #src
end #src
 #src