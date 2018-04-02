module FDFD2d

( export InputStruct, ChannelStruct, processInputs, updateInputs!, eigCF, eigKL,
eigNKL, scatter, smatrix, analyze_output, analyze_input)

using NLsolve
using SpecialFunctions
using Formatting
using Interpolations

################################################################################
include("FDFD2d_core.jl")
include("FDFD2d_eigensolvers.jl")
include("FDFD2d_scattering.jl")
include("FDFD2d_scattering_analysis.jl")
include("FDFD2d_parallel.jl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end #end of module FDFD2d
