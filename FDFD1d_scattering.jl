export scattering

################################################################################
####### SCATTERING
################################################################################
"""
ψ, H = scattering(input, k, a; isNonLinear=false, dispOpt=false, ψ_init=[],
    F=[1.], truncate=false, H=[], fileName="", ftol=2e-8, iter=750)

    Solves inhomogeneous problem with source defined in incident wave file.

    k is the injection frequency.

    A is a factorized wave operator.
"""
function scattering(input::InputStruct, K::Union{Complex128,Float64,Int},
    A::Union{Array{Complex128,1},Array{Float64,1},Array{Int,1}};
    isNonLinear::Bool=false, dispOpt::Bool = false, ψ_init::Array{Complex128,1}=Complex128[],
    F::Array{Float64,1}=[1.], truncate::Bool=false,
    H::Base.SparseArrays.UMFPACK.UmfpackLU=lufact(speye(1,1)), fileName::String = "",
    ftol::Float64=2e-8, iter::Int=150)::
    Tuple{Array{Complex128,1},Array{Complex128,1},Base.SparseArrays.UMFPACK.UmfpackLU}

    k = complex.(float.(K))
    a = complex.(float.(A))

    bc_original = set_bc!(input)

    if !isNonLinear
        ψ, φ, H = scattering_l(input, k, a; H=H, F=F)
    elseif isNonLinear
        ψ, φ, H = scattering_nl(input, k, a; H=H, F=F, dispOpt=dispOpt, ψ_init=ψ_init, ftol=ftol, iter=iter)
    end

    if !isempty(fileName)
        if truncate
            fid = open(fileName,"w")
            serialize(fid, (ψ[input.dis.xy_inds], input, k) )
            close(fid)
        else
            fid = open(fileName,"w")
            serialize(fid, (ψ, input, k) )
            close(fid)
        end
    end

    reset_bc!(input, bc_original)

    if truncate
        return ψ[input.dis.xy_inds], φ[input.dis.xy_inds], H
    else
        return ψ, φ, H
    end
end # end of function scattering

"""
ψ, A = scattering_l(input, k, a; H=[], F=[1.])

    Solves linear inhomogeneous problem.

    k is the frequency.

    a is a vector of amplitudes for each channel given in input.

    H is a factorized wave operator (useful for repeated calculation)
"""
function scattering_l(input::InputStruct, k::Complex128, a::Array{Complex128,1};
    H::Base.SparseArrays.UMFPACK.UmfpackLU=lufact(speye(1,1)), F::Array{Float64,1}=[1.])::
    Tuple{Array{Complex128,1},Array{Complex128,1},Base.SparseArrays.UMFPACK.UmfpackLU}

    k²= k^2

    bc_original = set_bc!(input)

    ∇² = laplacian(input, k)

    j, φ = synthesize_source(input, k, a)

    if (H.m == H.n == 1)
        N = prod(input.dis.N_PML)
        εt = input.sys.ε_PML
        Ft = input.sys.F_PML

        ɛk² = sparse(1:N, 1:N, ɛt[:]*k², N, N, +)
        χk² = sparse(1:N, 1:N, input.tls.D₀*γ(input,k)*F.*Ft[:]*k², N, N, +)

        H = factorize(∇²+ɛk²+χk²)
    end

    ψ = H\j

    reset_bc!(input, bc_original)

    return ψ, φ, H
end # end of function scattering_l

"""
j, ∇² = synthesize_source(input, k, a)

    Synthesizes source j for scattering calculation.

    input defines channels

    k is frequency

    a is amplitude vector.

    j is source

    ∇² is laplacian, returned because it was computed along the way.
"""
function synthesize_source(input::InputStruct, k::Complex128, a::Array{Complex128,1})::
    Tuple{Array{Complex128,1},Array{Complex128,1}}

    N = prod(input.dis.N_PML)
    φ₊ = zeros(Complex128,N)
    φ₋ = zeros(Complex128,N)
    φt₊ = zeros(Complex128,N)
    φt₋ = zeros(Complex128,N)

    for m in 1:length(input.sct.channels)
        φt₊, φt₋ = incident_mode(input, k, m)
        φ₊ += a[m]*φt₊
        φ₋ += a[m]*φt₋
    end

    φ = (input.sys.ε_PML[:] .!== input.sct.ε₀_PML[:]).*(φ₊+φ₋) + (input.sys.ε_PML[:] .== input.sct.ε₀_PML[:]).*φ₋

    k² = k^2
    k²δ = k²*(input.sys.ε_PML[:]-input.sct.ε₀_PML[:])
    j = -k²δ.*φ

    return j, φ
end # end of function synthesize_source


################################################################################
### NONLINEAR SOLVERS
################################################################################
"""
ψ, ϕ, A, input = scattering_nl(input, k; dispOpt=false, ψ_init=[],
    F=[1.], A=[], fileName="", ftol=2e-8, iter=750)

    Solves inhomogeneous problem with source defined in incident wave file.

    k is the injection frequency.

    A is a factorized wave operator.
"""
function scattering_nl(input1::InputStruct, k::Complex128;
    dispOpt::Bool = false, ψ_init::Array{Complex128,1}=Complex128[], F::Array{Float64,1}=[1.],
    A::Base.SparseArrays.UMFPACK.UmfpackLU=lufact(speye(1,1)),
    ftol::Float64=2e-8, iter::Int=150)::Tuple{Array{Complex128,1},Array{Complex128,1},Base.SparseArrays.UMFPACK.UmfpackLU,InputStruct}

        input = open_to_pml(input1)

        j, ∇², φ₊, φ₋ = createJ(input, k, m)

        N = prod(input.N_ext); ε_sm = input.ε_sm; k²= k^2;
        D₀ = input.D₀; F_sm = input.F_sm
        ɛk² = sparse(1:N, 1:N, ɛ_sm[:]*k², N, N, +)
        χk² = sparse(1:N, 1:N, D₀*γ(k,input)*F.*F_sm[:]*k², N, N, +)

        f!(Ψ, fvec) = scattering_residual(Ψ, fvec, j, ∇², εk², k, input)
        jac!(Ψ, jacarray) = scattering_jacobian(Ψ, jacarray, j, ∇², εk², k, input)
        df = DifferentiableSparseMultivariateFunction(f!, jac!)

        Ψ_init = Array{Float64}(2*length(j))
        Ψ_init[1:length(j)]     = real(ψ_init)
        Ψ_init[length(j)+1:2*length(j)] = imag(ψ_init)

        z = nlsolve(df, Ψ_init, show_trace=dispOpt, ftol=ftol, iterations=iter)

        if converged(z)
            ψ = z.zero[1:length(j)] + 1im*z.zero[length(j)+1:2*length(j)]
        else
            ψ = NaN*ψ
            println("Warning, solve_scattered did not converge. Returning NaN.")
        end

        if !isempty(fileName)
            if truncate
                fid = open(fileName,"w")
                serialize(fid, ((ψ-φ₊)[input.x̄_inds], φ₊[input.x̄_inds], input, k) )
                close(fid)
            else
                fid = open(fileName,"w")
                serialize(fid, ((ψ-φ₊), φ₊, input, k) )
                close(fid)
            end
        end

        if truncate
            return (ψ-φ₊)[input.x̄_inds], φ₊[input.x̄_inds], A, input
        else
            return ψ-φ₊, φ₊, A, input
        end
end # end of function computePsi_nonlinear

################################################################################
###  NONLINEAR AUXILLIARIES
################################################################################
"""
scattering_residual(Ψ, fvec, j, ∇², εk², k, input)
"""
function scattering_residual(Ψ::Array{Float64,1}, fvec::Array{Float64,1}, j,
                                    ∇², εk², k::Complex128, input::InputStruct)

    ψ = similar(j,Complex128)
    ind_r = 1:length(j)
    ind_i = length(j)+1:2*length(j)
    ψ = Ψ[ind_r] + 1im*Ψ[ind_i]
    temp = (∇²+ɛk²+χ(ψ,k,input)*k^2)*ψ - j
    fvec[ind_r] = real(temp)
    fvec[ind_i] = imag(temp)
end

"""
scattering_jacobian(Ψ, jacarray, j, ∇², εk², k, input)
"""
function scattering_jacobian(Ψ::Array{Float64,1}, jacarray, j, ∇², εk², k::Complex128, input::InputStruct)
    ψ = similar(j,Complex128)
    ind_r = 1:length(j)
    ind_i = length(j)+1:2*length(j)
    ψ = Ψ[ind_r] + 1im*Ψ[ind_i]
    temp = ∇²+ɛk²+χ(ψ,k,input)*k^2
    tempr = similar(temp,Float64)
    tempi = similar(temp,Float64)
    tr = nonzeros(tempr)
    ti = nonzeros(tempi)
    tr[:] = real((nonzeros(temp)))
    ti[:] = imag((nonzeros(temp)))
    tempj = [tempr+χʳʳ′(ψ,k,input) -tempi+χʳⁱ′(ψ,k,input); tempi+χⁱʳ′(ψ,k,input) tempr+χⁱⁱ′(ψ,k,input)]
    jacarray[:,:] = tempj[:,:]
end

"""
χ(ψ, k, input)
"""
function χ(ψ::Array{Complex128,1}, k::Complex128, input::InputStruct)::SparseMatrixCSC{Complex{Float64},Int64}
    N = prod(input.N_ext)
    h = hb(ψ,k,input)
    V = input.F_sm[:].*γ(k,input)*input.D₀./(1+h)
    return sparse(1:N, 1:N, V, N, N, +)
end

"""
"""
function hb(ψ::Array{Complex128,1},k::Complex128,input::InputStruct)::Array{Float64,1}
    N = prod(input.N_ext)
    h = abs2.(γ.(k, input)*ψ)
    return h
end

"""
χʳʳ′(Ψ, k, input)
"""
function χʳʳ′(ψ::Array{Complex128,1}, k::Complex128, input::InputStruct)::SparseMatrixCSC{Float64,Int64}
    N = prod(input.N_ext); D₀ = input.D₀; F_sm = input.F_sm; γt = γ(k, input)
    h = hb(ψ,k,input)
    V = -2D₀.*F_sm[:].*abs2(γt).*real.(γt.*ψ).*real.(ψ)./(1+h).^2
    return sparse(1:N,1:N,V,N,N,+)
end

"""
χⁱʳ′(Ψ, k, input)
"""
function χⁱʳ′(ψ::Array{Complex128,1}, k::Complex128, input::InputStruct)::SparseMatrixCSC{Float64,Int64}
    N = prod(input.N_ext); D₀ = input.D₀; F_sm = input.F_sm; γt = γ(k, input)
    h = hb(ψ,k,input)
    V = -2D₀.*F_sm[:].*abs2(γt).*imag.(γt.*ψ).*real.(ψ)./(1+h).^2
    return sparse(1:N,1:N,V,N,N,+)
end

"""
χʳⁱ′(Ψ, k, input)
"""
function χʳⁱ′(ψ::Array{Complex128,1}, k::Complex128, input::InputStruct)::SparseMatrixCSC{Float64,Int64}
    N = prod(input.N_ext); D₀ = input.D₀; F_sm = input.F_sm; γt = γ(k, input)
    h = hb(ψ,k,input)
    V = -2D₀.*F_sm[:].*abs2(γt).*real.(γt.*ψ).*imag.(ψ)./(1+h).^2
    return sparse(1:N,1:N,V,N,N,+)
end

"""
χⁱⁱ′(Ψ, k, input)
"""
function χⁱⁱ′(ψ::Array{Complex128,1}, k::Complex128, input::InputStruct)::SparseMatrixCSC{Float64,Int64}
    N = prod(input.N_ext); D₀ = input.D₀; F_sm = input.F_sm; γt = γ(k, input)
    h = hb(ψ,k,input)
    V = -2D₀.*F_sm[:].*abs2(γt).*imag.(γt.*ψ).*imag.(ψ)./(1+h).^2
    return sparse(1:N,1:N,V,N,N,+)
end
