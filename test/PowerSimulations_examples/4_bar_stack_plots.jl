using SIIPExamples #for path locations
pkgpath = dirname(dirname(pathof(SIIPExamples)))
using PowerSystems #to load results
using PowerSimulations #to load results
using PowerGraphics

simulation_folder = joinpath(pkgpath, "RTS-GMLC-master", "rts-test")
simulation_folder =
    joinpath(simulation_folder, "$(maximum(parse.(Int64,readdir(simulation_folder))))")
res = load_simulation_results(simulation_folder, "UC")

gr() # loads the GR backend
bar_plot(res)

plotlyjs()

bar_plot(res)

stack_plot(res, [Symbol("P__ThermalStandard"), Symbol("P__RenewableDispatch")], load = true)

uc_sys = System(joinpath(
    simulation_folder,
    "models_json",
    "stage_UC_model",
    "Stage1_sys_data.json",
))

fuel_plot(res, uc_sys, load = true)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

