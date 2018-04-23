################################################################################
#### CONSTANT FLUX EIGENVALUE SOLVERS
################################################################################

"""
η,u = eig_cf(inputs, k; nCF=1, F=[1.], η_init=0., u_init=[], truncate=false)

    Compute CF state according to bc set in inputs.
"""
function eig_cf(inputs::InputStruct, k::Union{Complex128,Float64,Int}; nCF::Int=1,
    F::Array{Float64,1}=[1.], η_init::Complex128=complex(0.),
    u_init::Array{Complex128,1}=Complex128[],
    truncate::Bool=false)::Tuple{Array{Complex128,1},Array{Complex128,2}}

    k = complex(float(k))

    k²= k^2

    ∇² = laplacian(k, inputs)

    N = inputs.N_ext; ε_sm = inputs.ε_sm; F_sm = inputs.F_sm
    ɛk² = sparse(1:prod(N), 1:prod(N), ɛ_sm[:]*k²       , prod(N), prod(N), +)
    sF  = sparse(1:prod(N), 1:prod(N), sign.(F.*F_sm[:]), prod(N), prod(N), +)
    FF  = sparse(1:prod(N), 1:prod(N), abs.(F.*F_sm[:]) , prod(N), prod(N), +)

    if isempty(u_init)
        (η,u,nconv,niter,nmult,resid) = eigs(-sF*(∇²+ɛk²)./k²,FF, which = :LM,
            nev = nCF, sigma = η_init)
    elseif !isempty(u_init)
        (η,u,nconv,niter,nmult,resid) = eigs(-sF*(∇²+ɛk²)./k²,FF, which = :LM,
            nev = nCF, sigma = η_init, v0 = u_init)
    end

    for ii = 1:nCF
        u[:,ii] = u[:,ii]/sqrt( trapz( u[:,ii].*F.*F_sm[:].*u[:,ii], inputs.dx̄ ) )
    end

    if truncate
        return η,u[inputs.x̄_inds,:]
    else
        return η,u
    end
end # end of function eig_cf for complex
function eig_cf(inputs::InputStruct, k::Union{Complex128,Float64,Int}, k_type::String; nCF::Int=1,
    F::Array{Float64,1}=[1.], η_init::Complex128=complex(0.),
    u_init::Array{Complex128,1}=Complex128[],
    truncate::Bool=false)::Tuple{Array{Complex128,1},Array{Complex128,2}}

    inputs1 = set_bc(inputs,k_type)

    η, u = eig_cf(inputs1, complex(k); nCF=nCF, F=F, η_init=η_init, u_init=u_init,
                    truncate=truncate)
    return η, u
end # end of function eig_cf with modified boundary
