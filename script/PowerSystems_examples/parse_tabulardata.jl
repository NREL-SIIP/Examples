# # Parsing Tabular Data

# **Originally Contributed by**: Clayton Barrows

# ## Introduction

# An example of how to parse tabular files (CSV) files similar to the format established in
# the [RTS-GMLC](github.com/gridmod/rts-gmlc/RTS_Data/SourceData) and create a `System` using
# [PowerSystems.jl](github.com/NREL-SIIP/PowerSystems.jl)

# ### Environment
# This notebook depends on the SIIPExamples.jl environment which is loaded by default

#nb using Pkg
#nb Pkg.status()

# ### Dependencies
using SIIPExamples
using PowerSystems
using TimeSeries
using Dates

# ### Fetch Data
# PowerSystems.jl links to some test data that is suitable for this example.
# Let's download the test data
PowerSystems.download(PowerSystems.TestData; branch = "master") # *note* add `force=true` to get a fresh copy
base_dir = dirname(dirname(pathof(PowerSystems)));

# ### The tabular data format relies on a folder containing `*.csv` files and a `user_descriptors.yaml` file
# First, we'll read the tabular data
RTS_GMLC_DIR = joinpath(base_dir, "data", "RTS_GMLC");
rawsys = PowerSystems.PowerSystemTableData(
    RTS_GMLC_DIR,
    100.0,
    joinpath(RTS_GMLC_DIR, "user_descriptors.yaml"),
    timeseries_metadata_file = joinpath(RTS_GMLC_DIR, "timeseries_pointers.json"),
    generator_mapping_file = joinpath(RTS_GMLC_DIR, "generator_mapping.yaml"),
)

# ### Create a `System`
# Next, we'll create a `System` from the `rawsys` data. Since a `System` is predicated on a
# forecast resolution and the `rawsys` data includes both 5-minute and 1-hour resolution
# time series, we also need to specify which time series we want to include in the `System`.
# The `time_series_resolution` kwarg filters to only include time series with a matching resolution.

sys = System(rawsys; time_series_resolution = Dates.Hour(1));
sys
