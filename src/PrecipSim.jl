module PrecipSim

# using Modules
using Parameters
using eFEM
using FileIO: save

using Test

# include Files
include("parameter_definition.jl")
include("problem_structs.jl")
include("mesh_functions.jl")
include("output_functions.jl")
include("driver_functions.jl")
include("matrix_functions.jl")
include("solution_methods.jl")

function runTest()
    p = precipParam()
    runSimulation(p)
end

# export functions
export precipParam, runSimulation

###
export runTest, exportParam_test          # JUST FOR TESTING.
export QSSRstep_test
###

end # module
