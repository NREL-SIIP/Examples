# # Hydropower Simulations with [PowerSimulations.jl](https://github.com/NREL-SIIP/PowerSimulations.jl)

# **Originally Contributed by**: Clayton Barrows and Sourabh Dalvi

# ## Introduction

# PowerSimulations.jl supports simulations that consist of sequential optimization problems
# where results from previous problems inform subsequent problems in a variety of ways.
# This example demonstrates a few of the options for modeling hydropower generation.

# ## Dependencies
using SIIPExamples
pkgpath = dirname(dirname(pathof(SIIPExamples)))

# ### Modeling Packages
using PowerSystems
using PowerSimulations
const PSI = PowerSimulations
using D3TypeTrees

# ### Data management packages
using Dates
using DataFrames

# ### Optimization packages
using Cbc # solver
solver = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 1, "ratioGap" => 0.05)

# ### Data
# PowerSystems.jl links to some meaningless test data that is suitable for this example.
# The [make_hydro_data.jl](../../script/PowerSimulations_examples/make_hydro_data.jl) script
# loads three systems suitable for the examples here.

include(joinpath(pkgpath, "script", "PowerSimulations_examples", "make_hydro_data.jl"))

# This line just overloads JuMP printing to remove double underscores added by PowerSimulations.jl
PSI.JuMP._wrap_in_math_mode(str) = "\$\$ $(replace(str, "__"=>"")) \$\$"

# ## Two PowerSimulations features determine hydropower representation.
# There are two principal ways that we can customize hydropower representation in
# PowerSimulations. First, we can play with the formulation applied to hydropower generators
# using the `DeviceModel`. We can also adjust how simulations are configured to represent
# different decision making processes and the information flow between those processes.

# ## Hydropower `DeviceModel`s

# First, the assignment of device formulations to particular device types gives us control
# over the representation of devices. This is accomplished by defining `DeviceModel`
# instances. For hydro power representations, we have two available generator types in
# PowerSystems:

#nb TypeTree(PowerSystems.HydroGen)

# And in PowerSimulations, we have several available formulations that can be applied to
# the hydropower generation devices:

#nb TypeTree(PSI.AbstractHydroFormulation, scopesep = "\n", init_expand = 5)

# ### `FixedOutput`

# Let's see what some of the different combinations create. First, let's apply the
# `FixedOutput` formulation to the `HydroEnergyReservoir` generators.
#  - The `FixedOutput` formulation just acts
# like a load subtractor, forcing the system to accept it's generation.

devices = Dict{Symbol,DeviceModel}(:Hyd1 => DeviceModel(HydroDispatch, FixedOutput));

template = PSI.OperationsProblemTemplate(CopperPlatePowerModel, devices, Dict(), Dict());

op_problem = PSI.OperationsProblem(GenericOpProblem, template, demo_sys, horizon = 2)

# Now we can see the resulting JuMP model:

op_problem.psi_container.JuMPmodel

# The two constrants are the power balance constarints in each of the two periods. Since
# we haven't added a load or any coontrollable generation to this model, it is infeasible.
# But you can see that the values correspond to the `get_max_active_power` forecast for the
# `HydroDispatch2` generator.

get_forecast_values(
    Deterministic,
    get_component(HydroDispatch, demo_sys, "HydroDispatch2"),
    DateTime(today()),
    "get_max_active_power",
    2,
)

# ### `HydroDispatchRunOfRiver`

# Next, let's apply the `HydroDispatchRunOfRiver` formulation to the `HydroEnergyReservoir`
# generators.
#  - The `HydroDispatchRunOfRiver` formulation represents the the energy flowing out of
# a reservoir. The model can choose to produce power with that energy or just let it spill by.

devices = Dict{Symbol,DeviceModel}(
    :Hyd1 => DeviceModel(HydroEnergyReservoir, HydroDispatchRunOfRiver),
);

template = PSI.OperationsProblemTemplate(CopperPlatePowerModel, devices, Dict(), Dict());

op_problem = PSI.OperationsProblem(GenericOpProblem, template, demo_sys, horizon = 2)

op_problem.psi_container.JuMPmodel

# The first two constraints are the power balance constraints that require the generation
# from the controllable `HydroEnergyReservoir` generators to be equal to the (non-existant)
# load. The 3rd through 6th constraints limit the output of the `HydroEnergyReservoir`
# generators to the limit defined by the `get_max_active_pwoer` forecast.
# And the remainider of the constraints are the lower and upper bounds of
# the `HydroEnergyReservoir` operating ranges.

# ### `HydroDispatchReservoirStorage`

# Next, let's apply the `HydroDispatchReservoirStorage` formulation to the `HydroEnergyReservoir` generators.
devices = Dict{Symbol,DeviceModel}(
    :Hyd1 => DeviceModel(HydroEnergyReservoir, HydroDispatchReservoirStorage),
);

template = PSI.OperationsProblemTemplate(CopperPlatePowerModel, devices, Dict(), Dict());

op_problem = PSI.OperationsProblem(GenericOpProblem, template, demo_sys, horizon = 2)

# And, the resulting JuMP model:

op_problem.psi_container.JuMPmodel

# This time, the first and third constraints represent the energy balance in each `HydroEnergyReservoir`
# reservoirs in the first time period, accounting for  the energy released through spillage
# `Sp` and power production `P`, and the inflow energy and the initial stored energy are rerpresented
# on the right hand side from the corresponding `get_inflow` and `get_storage_capacity` forecasts, reespectively.

get_forecast_values(
    Deterministic,
    get_component(HydroEnergyReservoir, demo_sys, "HydroDispatch3"),
    DateTime(today()),
    "get_storage_capacity",
    2,
)

get_forecast_values(
    Deterministic,
    get_component(HydroEnergyReservoir, demo_sys, "HydroDispatch3"),
    DateTime(today()),
    "get_inflow",
    2,
)

# The second and fourth constraints represent the energy balance in the reservoirs in subsequent
# periods. Constraints 5 and 6 are the power balance constraints, followed by lower bounds of
# zero for power (`P`), energy (`E`), and spillage (`Sp`) variables, and upper boounds for
# `P` and `E` defined by the `active_power_limits.max` and `storage_capacity` parameters.

# ### `HydroDispatchReservoirBudget`

# Finally, let's apply the `HydroDispatchReservoirBudget` formulation to the `HydroEnergyReservoir` generators.

devices = Dict{Symbol,DeviceModel}(
    :Hyd1 => DeviceModel(HydroEnergyReservoir, HydroDispatchReservoirBudget),
);

template = PSI.OperationsProblemTemplate(CopperPlatePowerModel, devices, Dict(), Dict());

op_problem = PSI.OperationsProblem(GenericOpProblem, template, demo_sys, horizon = 2)

#-

op_problem.psi_container.JuMPmodel

# Again, the first and second constraints represent the energy balance. The third and fourth
# constraints represent the energy budget for each `HydroEnergyReservoir` generator over the
# optimization horizon. The right hand side of these constraints correspond to the sums of
# the values of the "get_hydro_budget" forecast multiplied by the `storage_capacity`:
g = get_component(HydroEnergyReservoir, demo_sys, "HydroDispatch3")
get_forecast(
    Deterministic,
    g,
    DateTime(today()),
    "get_hydro_budget",
    2,
).data .* get_storage_capacity(g) |> sum

# The remainder of the constraints correspond to the lower and upper bounds of each generator in each time period.

# ### Multi-Stage `SimulationSequence`
# The purpose of a multi-stage simulation is to represent scheduling decisions consistently
# with the time scales that govern different elements of power systems.

# #### Multi-Day to Daily Simulation:
# In the multi-day model, we'll use a really simple representation of all system devices
# so that we can maintain computational tractability while getting an estimate of system
# requirements/capabilities.

devices = Dict(
    :Generators => DeviceModel(ThermalStandard, ThermalBasicUnitCommitment),
    :Loads => DeviceModel(PowerLoad, StaticPowerLoad),
    :HydroEnergyReservoir =>
        DeviceModel(HydroEnergyReservoir, HydroDispatchReservoirBudget),
)
template_md = OperationsProblemTemplate(CopperPlatePowerModel, devices, Dict(), Dict());

# For the daily model, we can increase the modeling detail since we'll be solving shorter
# problems.

devices = Dict(
    :Generators => DeviceModel(ThermalStandard, ThermalBasicUnitCommitment),
    :Loads => DeviceModel(PowerLoad, StaticPowerLoad),
    :HydroEnergyReservoir =>
        DeviceModel(HydroEnergyReservoir, HydroDispatchReservoirBudget),
)
template_da = OperationsProblemTemplate(CopperPlatePowerModel, devices, Dict(), Dict());

op_problem = PSI.OperationsProblem(GenericOpProblem, template_md, demo_wk_sys, horizon = 2)

#-

stages_definition = Dict(
    "MD" => Stage(
        GenericOpProblem,
        template_md,
        demo_wk_sys,
        solver,
        system_to_file = false,
    ),
    "DA" => Stage(
        GenericOpProblem,
        template_da,
        demo_da_sys,
        solver,
        system_to_file = false,
    ),
)

# This builds the sequence and passes the the energy dispatch schedule for the `HydroEnergyReservoir`
# generator from the "MD" stage to the "DA" stage in the form of an energy limit over the
# synchronized periods.

sequence = SimulationSequence(
    step_resolution = Hour(48),
    order = Dict(1 => "MD", 2 => "DA"),
    feedforward_chronologies = Dict(("MD" => "DA") => Synchronize(periods = 2)),
    horizons = Dict("MD" => 2, "DA" => 24),
    intervals = Dict("MD" => (Hour(48), Consecutive()), "DA" => (Hour(24), Consecutive())),
    feedforward = Dict(
        ("DA", :devices, :HydroEnergyReservoir) => IntegralLimitFF(
            variable_source_stage = PSI.ACTIVE_POWER,
            affected_variables = [PSI.ACTIVE_POWER],
        ),
    ),
    #cache = Dict(("MD", "DA") => StoredEnergy(HydroEnergyReservoir, PSI.ENERGY)),
    ini_cond_chronology = InterStageChronology(),
);

#-

file_path = tempdir()

sim = Simulation(
    name = "hydro",
    steps = 1,
    stages = stages_definition,
    stages_sequence = sequence,
    simulation_folder = file_path,
)

build!(sim)

# We can look at the "MD" Model

sim.stages["MD"].internal.psi_container.JuMPmodel

# And we can look at the "DA" model

sim.stages["DA"].internal.psi_container.JuMPmodel

# And we can execute the simulation by running the following command
# ```julia
# sim_results = execute!(sim)
# ```
#-

# #### 3-Stage Simulation:

stages_definition = Dict(
    "MD" => Stage(
        GenericOpProblem,
        template_md,
        demo_wk_sys,
        solver,
        system_to_file = false,
    ),
    "DA" => Stage(
        GenericOpProblem,
        template_da,
        demo_da_sys,
        solver,
        system_to_file = false,
    ),
    "ED" => Stage(
        GenericOpProblem,
        template_da,
        demo_rt_sys,
        solver,
        system_to_file = false,
    ),
)

sequence = SimulationSequence(
    step_resolution = Hour(72),
    order = Dict(1 => "MD", 2 => "DA", 3 => "ED"),
    feedforward_chronologies = Dict(
        ("MD" => "DA") => Synchronize(periods = 1),
        ("DA" => "ED") => Synchronize(periods = 24),
    ),
    intervals = Dict(
        "MD" => (Hour(72), Consecutive()),
        "DA" => (Hour(24), Consecutive()),
        "ED" => (Hour(1), Consecutive()),
    ),
    horizons = Dict("MD" => 3, "DA" => 24, "ED" => 12),
    feedforward = Dict(
        ("DA", :devices, :HydroEnergyReservoir) => IntegralLimitFF(
            variable_source_stage = PSI.ACTIVE_POWER,
            affected_variables = [PSI.ACTIVE_POWER],
        ),
        ("ED", :devices, :HydroEnergyReservoir) => IntegralLimitFF(
            variable_source_stage = PSI.ACTIVE_POWER,
            affected_variables = [PSI.ACTIVE_POWER],
        ),
    ),
    #cache = Dict(("MD", "DA") => StoredEnergy(HydroEnergyReservoir, PSI.ENERGY)),
    ini_cond_chronology = InterStageChronology(),
);

#-

sim = Simulation(
    name = "hydro",
    steps = 1,
    stages = stages_definition,
    stages_sequence = sequence,
    simulation_folder = file_path,
)

#-

build!(sim)

# We can look at the "MD" Model

sim.stages["MD"].internal.psi_container.JuMPmodel

# And we can look at the "DA" model

sim.stages["DA"].internal.psi_container.JuMPmodel

# And we can look at the "ED" model

sim.stages["ED"].internal.psi_container.JuMPmodel

# And we can execute the simulation by running the following command
# ```julia
# sim_results = execute!(sim)
# ```
#-
