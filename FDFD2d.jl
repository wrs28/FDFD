module FDFD2d

( export InputStruct, ChannelStruct, process_input, update_input!, eig_cf, eig_k, eig_kl,
eig_nkl, eig_nklp)#, scattering, smatrix, analyzeOutput, analyzeInput)

using NLsolve
using SpecialFunctions
using Formatting
using Interpolations

################################################################################
include("FDFD2d_defaults.jl")
include("FDFD2d_core_auxilliaries.jl")
include("FDFD2d_core_numerics.jl")
# include("FDFD2d_eigensolvers.jl")
include("FDFD2d_scattering.jl")
include("FDFD2d_incident_fields.jl")
# include("FDFD2d_scattering_auxilliaries.jl")
# include("FDFD2d_parallel.jl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end #end of module FDFD2d
