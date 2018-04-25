export eig_k, eig_cf

################################################################################
#### CONSTANT FLUX EIGENVALUE SOLVERS
################################################################################
"""
η,u = eig_cf(input, k; ncf=1, F=[1.], η_init=0., u_init=[], truncate=false)
η,u = eig_cf(input, k, k_type; ncf=1, F=[1.], η_init=0., u_init=[], truncate=false)

    Compute CF state according to bc set in input or string k_type.
"""
function eig_cf(input::InputStruct, k::Union{Complex128,Float64,Int}; ncf::Int=1,
    F::Array{Float64,1}=[1.], η_init::Union{Complex128,Float64,Int}=complex(0.),
    u_init::Array{Complex128,1}=Complex128[],
    truncate::Bool=false)::Tuple{Array{Complex128,1},Array{Complex128,2}}

    k = complex(float(k))

    k²= k^2

    ∇² = laplacian(input, k)

    N = input.dis.N_PML; ε = input.sys.ε_PML; F1 = input.sys.F_PML
    ɛk² = sparse(1:prod(N), 1:prod(N), ɛ[:]*k²       , prod(N), prod(N), +)
    sF  = sparse(1:prod(N), 1:prod(N), sign.(F.*F1[:]), prod(N), prod(N), +)
    FF  = sparse(1:prod(N), 1:prod(N), abs.(F.*F1[:]) , prod(N), prod(N), +)

    if isempty(u_init)
        (η,u,nconv,niter,nmult,resid) = eigs(-sF*(∇²+ɛk²)./k²,FF, which = :LM,
            nev = ncf, sigma = η_init)
    elseif !isempty(u_init)
        (η,u,nconv,niter,nmult,resid) = eigs(-sF*(∇²+ɛk²)./k²,FF, which = :LM,
            nev = ncf, sigma = η_init, v0 = u_init)
    end

    for ii = 1:ncf
        u[:,ii] = u[:,ii]/sqrt( trapz( u[:,ii].*F.*F_sm[:].*u[:,ii], input.dx ) )
    end

    if truncate
        return η,u[input.dis.xy_inds,:]
    else
        return η,u
    end
end # end of function eig_cf for complex
function eig_cf(input::InputStruct, k::Union{Complex128,Float64,Int}, k_type::String; ncf::Int=1,
    F::Array{Float64,1}=[1.], η_init::Union{Complex128,Float64,Int}=complex(0.),
    u_init::Array{Complex128,1}=Complex128[],
    truncate::Bool=false, direction::Array{Int,1}=[1,0])::
    Tuple{Array{Complex128,1},Array{Complex128,2}}

    input, bc_original = set_bc(input, k_type, direction)

    η, u = eig_cf(input1, k; ncf=ncf, F=F, η_init=η_init, u_init=u_init,
                    truncate=truncate)

    reset_bc!(input, bc_original)
    return η, u
end # end of function eig_cf with modified boundary

function eig_k(input::InputStruct, k::Union{Complex128,Float64,Int};
    isLinear::Bool=true, nk::Int=1, F::Array{Float64,1}=[1.], truncate::Bool=false,
    ψ_init::Array{Complex128,1}=Complex128[])::Tuple{Array{Complex128,1},Array{Complex128,2}}

    K = complex(float(k))

    if isLinear
        k,ψ = eig_kl(input, K; nk=nk, F=F, truncate=truncate, ψ_init=ψ_init)
    else
        k,ψ = eig_knl(input, K; nk=nk, F=F, truncate=truncate, ψ_init=ψ_init)
    end

    return k,ψ
end # end of function eig_k (wrapper for linear or nonlinear cf root finder)
function eig_k(input::InputStruct, k::Union{Complex128,Float64,Int}, k_type::String;
    isLinear::Bool=true, nk::Int=1, F::Array{Float64,1}=[1.], truncate::Bool=false,
    ψ_init::Array{Complex128,1}=Complex128[], direction::Array{Int,1}=[1,0])::
    Tuple{Array{Complex128,1},Array{Complex128,2}}

    input, bc_original = set_bc(input,k_type, direction)

    k, ψ = eig_k(input, k; isLinear=isLinear, nk=nk, F=F, truncate=truncate, ψ_init=ψ_init)

    reset_bc!(input, bc_original)
    return k, ψ
end # end of function eig_k with modified boundary

################################################################################
####### LINEAR EIGENVALUE SOLVERS
################################################################################

"""
k,ψ =  eig_kl(input, k; nk=1, F=[1.], truncate=false, ψ_init=[])
k,ψ =  eig_kl(input, k, k_type; nk=1, F=[1.], truncate=false, ψ_init=[])

    Compute eigenfrequencies k w/o line pulling, according to boundary conditions
    set in input.bc or determined by k_type

    Set D or F to zero to get passive cavity.

    k = k[# of poles]

    ψ = ψ[cavity size, # of poles]
"""
function eig_kl(input::InputStruct, K::Union{Complex128,Float64,Int}; nk::Int=1,
    F::Array{Float64,1}=[1.], truncate::Bool=false, ψ_init::Array{Complex128,1}=Complex128[])::
    Tuple{Array{Complex128,1},Array{Complex128,2}}

    k = complex(float(K))

    ∇² = laplacian(input, k)

    N = prod(input.dis.N_PML); ε_sm = input.sys.ε_PML; D₀ = input.tls.D₀; F_sm = input.sys.F_PML
    ɛ⁻¹ = sparse(1:N, 1:N, 1./(ɛ_sm[:]-1im*D₀.*F.*F_sm[:]), N, N, +)
    if isempty(ψ_init)
        (k²,ψ,nconv,niter,nmult,resid) = eigs(-ɛ⁻¹*∇²,which = :LM, nev = nk,
            sigma = k^2)
    else
        (k²,ψ,nconv,niter,nmult,resid) = eigs(-ɛ⁻¹*∇²,which = :LM, nev = nk,
            sigma = k^2, v0 = ψ_init)
    end

    inds = input.dis.xy_inds
    if length(F) == 1
        F_temp = F
    else
        F_temp = F[inds]
    end

    for ii = 1:nk
        ψ[:,ii] = ψ[:,ii]/sqrt( trapz( (ɛ_sm[inds]-1im*D₀.*F_temp[:].*F_sm[inds]).*abs2.(ψ[inds,ii]), input.dis.dx) )
    end

    if truncate
        return sqrt.(k²),ψ[inds,:]
    else
        return sqrt.(k²),ψ
    end
end # end of function eig_kl (main)

################################################################################
#### NONLINEAR EIGENVALUE SOLVERS
################################################################################

"""
k, u, η, conv = computeK_NL1(input, k_init; F=[1.], dispOpt=false,η_init = 0., u_init = [],
    k_avoid = 0., tol = .5, max_count = 15, max_iter = 50)

    compute eigenfrequency k with dispersion via root-finding of η(k) - D₀γ(k).

    max_count: the max number of CF states to compute at a time.
    max_iter: maximum number of iterations to include in nonlinear solve

k = computeK_NL2(input, kc, Radii; nk=3, Nq=100, F=[1.], R_min=.01, rank_tol=1e-8)

    Compute eigenfrequency with dispersion, using contour integration. BC's set
    by input.bc

    Contour is centered on kc, Radii = (x-diameter, y-diameter).

    nk is an upper bound on the number of eigenfrequencies contained in the contour.

    Nq is the number of contour quadrature points.
"""
function eig_knl(input::InputStruct, k_init::Union{Complex128,Float64,Int};
    F::Array{Float64,1}=[1.], dispOpt::Bool=false, η_init::Complex128=complex(0.),
    u_init::Array{Complex128,1}=Complex128[], k_avoid::Array{Complex128,1}=Complex128[0],
    tol::Float64=.5, max_count::Int = 15, max_iter::Int=50)::
    Tuple{Array{Complex128,1},Array{Complex128,2},Array{Complex128,1},Bool}

    ηt, ut = eig_cf(input, k_init; ncf = 1, η_init=η_init, u_init = u_init)

    γ⟂ = input.γ⟂; k₀ = input.k₀; D₀ = input.D₀; dx̄ = input.dx̄
    function f!(z, fvec)

        k = z[1]+1im*z[2]

        flag = true
        count = 1
        M = 1
        ind = Int

        while flag

            η_temp, u_temp = eig_cf(input, k; ncf = M, F=F, η_init=ηt[1], u_init = ut[:,1])

            overlap = zeros(Float64,M)
            for i in 1:M
                overlap[i] = abs(trapz(ut[:,1].*input.F_sm[:].*F.*u_temp[:,i],dx̄))
            end

            max_overlap, ind = findmax(overlap)

            if (max_overlap > (1-tol))
                flag = false
                ηt[1] = η_temp[ind]
                ut[:,1] = u_temp[:,ind]
            elseif  (count < max_count) & (max_overlap ≤ (1-tol))
                M += 1
            else
                flag = false
                ηt[1] = η_temp[ind]
                ut[:,1] = u_temp[:,ind]
                println("Warning: overlap less than tolerance, skipped to neighboring TCF.")
            end

            count += 1

        end

        fvec[1] = real((ηt[1]-D₀*γ(k,input))/prod(k-k_avoid))
        fvec[2] = imag((ηt[1]-D₀*γ(k,input))/prod(k-k_avoid))

    end

    z = nlsolve(f!,[real(k_init),imag(k_init)]; iterations = max_iter, show_trace = dispOpt)
    k = z.zero[1]+1im*z.zero[2]
    conv = converged(z)

    ηt[1] = D₀*γ(k,input)
    η, ψ = eig_cf(input, k; ncf = 1, η_init=ηt[1], F=F, u_init = ut[:,1])

    return [k], ψ, η, conv
end # end of function eig_knl, nonlinear cf rootfinder
function eig_knl(input::InputStruct, k_init::Union{Complex128,Float64,Int}, k_type::String;
    F::Array{Float64,1}=[1.], dispOpt::Bool=false, η_init::Complex128=complex(0.),
    u_init::Array{Complex128,1}=Complex128[], k_avoid::Array{Complex128,1}=Complex128[0],
    tol::Float64=.5, max_count::Int=15, max_iter::Int=50, direction::Array{Int,1}=[1,0])::
    Tuple{Array{Complex128,1},Array{Complex128,2},Array{Complex128,1},Bool}

    input1 = set_bc(input,k_type)

    k, ψ, η, conv = eigKNL(input1, k_init; F=F, dispOpt=dispOpt, η_init=η_init,
        u_init=u_init, k_avoid=k_avoid, tol=tol, max_count=max_count, max_iter=max_iter)
end # end of function eig_knl, nonlinear cf rootfinder with modified boundaries


function eig_knl(input::InputStruct, kc::Union{Complex128,Float64,Int}, Radii::Tuple{Float64,Float64};
    nk::Int=3, Nq::Int=100, F::Array{Float64,1}=[1.], R_min::Float64=.01,
    rank_tol::Float64=1e-8)::Array{Complex{Float64},1}

    k = complex(float(kc))

    ∇² = laplacian(input, k)

    N_ext = prod(input.N_ext)
    A  = zeros(Complex128,N_ext,nk)
    A₀ = zeros(Complex128,N_ext,nk)
    A₁ = zeros(Complex128,N_ext,nk)
    M = rand(N_ext,nk)

    ϕ = 2π*(0:1/Nq:(1-1/Nq))
    Ω = k + Radii[1]*cos.(ϕ) + 1im*Radii[2]*sin.(ϕ)

    θ =  angle(input.k₀ - 1im*input.γ⟂-k)
    flag = abs(input.k₀ - 1im*input.γ⟂-k) < rad(Radii[1],Radii[2],θ)

    if flag
        AA  = zeros(Complex128,N_ext,nk)
        AA₀ = zeros(Complex128,N_ext,nk)
        AA₁ = zeros(Complex128,N_ext,nk)
        RR = 2*R_min
        ΩΩ = k₀-1im*γ⟂ + (RR/2)*cos.(ϕ) + 1im*(RR/2)*sin.(ϕ)
    end

    ε_sm = input.ε_sm; F_sm = input.F_sm; D₀ = input.D₀; dx̄ = input.dx̄
    for i in 1:Nq

        k′ = Ω[i]
        k′² = k′^2

        if (i > 1) & (i < Nq)
            dk′ = (Ω[i+1] - Ω[i-1]  )/2
        elseif i == Nq
            dk′ = (Ω[1]   - Ω[end-1])/2
        elseif i == 1
            dk′ = (Ω[2]   - Ω[end]  )/2
        end

        ɛk′² = sparse(1:N_ext, 1:N_ext, ɛ_sm[:]*k′², N_ext, N_ext, +)
        χk′² = sparse(1:N_ext, 1:N_ext, D₀*γ(k′,input)*F.*F_sm[:]*k′², N_ext, N_ext, +)

        A = (∇²+ɛk′²+χk′²)\M
        A₀ += A*dk′/(2π*1im)
        A₁ += A*k′*dk′/(2π*1im)

        if flag
            kk′ = ΩΩ[i]
            kk′² = kk′^2
            if (i > 1) & (i < Nq)
                dkk′ = (ΩΩ[i+1] - ΩΩ[i-1]  )/2
            elseif i == Nq
                dkk′ = (ΩΩ[1]   - ΩΩ[end-1])/2
            elseif i == 1
                dkk′ = (ΩΩ[2]   - ΩΩ[end]  )/2
            end
            ɛkk′² = sparse(1:N_ext, 1:N_ext, ɛ_sm[:]*kk′², N_ext, N_ext, +)
            χkk′² = sparse(1:N_ext, 1:N_ext, D₀*γ(kk′,input)*F.*F_sm[:]*kk′², N_ext, N_ext, +)

            AA = (∇²+ɛkk′²+χkk′²)\M
            AA₀ += AA*dkk′/(2π*1im)
            AA₁ += AA*kk′*dkk′/(2π*1im)

       end

    end

    if flag
        A₀ = A₀-AA₀
        A₁ = A₁-AA₁
    end

    P = svdfact(A₀,thin = true)
    temp = find(P[:S] .< rank_tol)
    if isempty(temp)
        println("Error. Need more nevals")
        return Complex128[NaN]
    else
        k = temp[1]-1
    end

    B = (P[:U][:,1:k])'*A₁*(P[:Vt][1:k,:])'*diagm(1./P[:S][1:k])

    D,V = eig(B)

    return D
end # end of function eig_knl, contour integral
function eig_knl(input::InputStruct, kc::Union{Complex128,Float64,Int}, Radii::Tuple{Float64,Float64},
    k_type::String;
    nk::Int=3, Nq::Int=100, F::Array{Float64,1}=[1.], R_min::Float64=.01,
    rank_tol::Float64=1e-8)::Array{Complex{Float64},1}

    input1 = set_bc(input,k_type)

    k = eig_knl(input, kc, Radii; nk=nk, Nq=Nq, F=F, R_min=R_min, rank_tol=rank_tol)
end # end of function eig_knl, contour integral with modified boundaries

################################################################################
#### AUXILLIARIES
################################################################################


"""
rad(a,b,θ) is the complex radius as a function of angle for an ellipse with
    semi-major axis a and semi-minor axis b.

    For use in contour eigensolver.
"""
function rad(a::Float64, b::Float64, θ::Float64)::Float64
    return b./sqrt.(sin(θ).^2+(b/a)^2.*cos(θ).^2)
end

function set_bc(input::InputStruct, k_type::String, direction::Array{Int,1})::
    Tuple{InputStruct, Array{String, 1}}

    bc_original = input.bnd.bc

    if k_type ∈ ["Pole","pole","P","p"]
        for i in 1:4
            if input.bnd.bc[i] == "I"
                input.bnd.bc[i] = "O"
            end
        end
    elseif k_type ∈ ["Zero","zero","Z","z"]
        for i in 1:4
            if input.bnd.bc[i] == "O"
                input.bnd.bc[i] = "I"
            end
        end
    elseif k_type ∈ ["UZR","uzr","U","u"]
        if (direction[1]==+1) && (input.bnd.bc[1:2] !== ["I", "O"])
            input.bnd.bc[1] = "I"
            input.bnd.bc[2] = "O"
        elseif (direction[1]==-1) && (input1.bc[1:2] !== ["O", "I"])
            input.bnd.bc[1] = "O"
            input.bnd.bc[2] = "I"
        end
        if (direction[2]==+1) && (input1.bc[3:4] !== ["I", "O"])
            input.bnd.bc[3] = "I"
            input.bnd.bc[4] = "O"
        elseif (direction[2]==-1) && (input1.bc[3:4] !== ["O", "I"])
            input.bnd.bc[3] = "O"
            input.bnd.bc[4] = "I"
        end
    end

    input.bnd.bc_sig = prod(input.bnd.bc)

    return input, bc_original
end

function reset_bc!(input::InputStruct,bc_original::Array{String,1})::InputStruct

    input.bnd.bc = bc_original
    input.bnd.bc_sig = prod(bc_original)

    return input
end
