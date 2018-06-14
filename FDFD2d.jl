module FDFD2d

using Formatting
using Interpolations
using SpecialFunctions

################################################################################
include("FDFD2d_defaults.jl")
include("FDFD2d_structures.jl")
include("FDFD2d_core_numerics.jl")
include("FDFD2d_scattering.jl")
include("FDFD2d_synthesis.jl")
include("FDFD2d_smatrix.jl")
include("FDFD2d_analysis.jl")
include("FDFD2d_eigensolvers.jl")
include("FDFD2d_parallel.jl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end #end of module FDFD2d
