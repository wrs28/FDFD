module FDFD1d

using Formatting
using Interpolations

################################################################################
include("FDFD1d_defaults.jl")
include("FDFD1d_structures.jl")
include("FDFD1d_core_numerics.jl")
include("FDFD1d_scattering.jl")
include("FDFD1d_synthesis.jl")
include("FDFD1d_smatrix.jl")
include("FDFD1d_analysis.jl")
include("FDFD1d_eigensolvers.jl")
include("FDFD1d_parallel.jl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end #end of module FDFD21
