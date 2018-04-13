"""
ψ, A = scattering_l(input, k, a; A=[], F=[1.])

    Solves linear inhomogeneous problem.

    k is the frequency.

    a is a vector of amplitudes for each channel given in input.

    A is a factorized wave operator (useful for repeated calculation)
"""
function scattering_l(input::InputStruct, K::Union{Complex128,Float64,Int},
    A::Union{Array{Complex128,1},Array{Float64,1},Array{Int,1}};
    H::Base.SparseArrays.UMFPACK.UmfpackLU=lufact(speye(1,1)), F::Array{Float64,1}=[1.])::
    Tuple{Array{Complex128,1},Base.SparseArrays.UMFPACK.UmfpackLU,InputStruct}

    k = complex(float(K))
    k²= k^2
    a = complex(float(A))

    j, ∇² = synthesize_source(input, k, a)

    if (H.m == H.n == 1)
        N = prod(input.N_sct)
        ε_sct = input.ε_sct
        F_sct = input.F_sct

        ɛk² = sparse(1:N, 1:N, ɛ_sct[:]*k², N, N, +)
        χk² = sparse(1:N, 1:N, input.D₀*γ(input,k)*F.*F_sct[:]*k², N, N, +)

        H = factorize(∇²+ɛk²+χk²)
    end

    ψ = H\j

    return ψ, H, input
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
    Tuple{Array{Complex128,1},SparseMatrixCSC{Complex128,Int64}}

    N = prod(inputs.N_sct)
    φ₊ = zeros(Complex128,N)
    φ₋ = zeros(Complex128,N)
    φt₊ = zeros(Complex128,N)
    φt₋ = zeros(Complex128,N)

    M₊, M₋ = source_mask(input)

    for m in 1:length(input.channels)
        φt₊, φt₋ = incident_mode(input, k, m)
        φ₊ += a[m]*φt₊
        φ₋ += a[m]*φt₋
    end

    ∇² = laplacian(input, k; PML=true)

    k² = k^2
    ɛk² = sparse(1:N, 1:N, input.ɛ_wgs[:]*k², N, N, +)

    j = (∇²+ɛk²)*(M₊.*φ₊ + M₋.*φ₋)

    return j, ∇²
end # end of function synthesize_source


"""
M₊, M₋ = source_mask(input)

    M₊ is the mask for the incident field, corresponds to a circle or radius ∂S
        if length(∂S) = 1, and a rectangle with boundaries at the elements of ∂S
        otherwise
    M₋ is the mask for the outgoing field, corresponds to a rectangle at PML
        boundary
"""
function source_mask(input::InputStruct)::Tuple{Array{Bool,1},Array{Bool,1}}

    ∂S = input.∂S
    ∂R = input.∂R

    if length(∂S)>1
        M₊ = (∂S[1] .≤ input.x̄_sct[1] .≤ ∂S[2]) .& (∂S[3] .≤ input.x̄_sct[2] .≤ ∂S[4])
    elseif length(∂S)==1
        r = sqrt.(inputs.x̄_sct[1].^2 + inputs.x̄_sct[2].^2)
        M₊ = r .≤ ∂S[1]
    end

    M₋ = (∂R[1] .≤ input.x̄_sct[1] .≤ ∂R[2]) .& (∂R[3] .≤ input.x̄_sct[2] .≤ ∂R[4])

    return M₊, M₋
end # end of function source_mask



















################################################################################
### INHOMOGENEOUS SOLVERS
################################################################################

"""
ψ, A, input = scattering(input, k; isNonLinear=false, dispOpt=false, ψ_init=[],
    F=[1.], truncate=false, A=[], fileName="", ftol=2e-8, iter=750)

    Solves inhomogeneous problem with source defined in incident wave file.

    k is the injection frequency.

    A is a factorized wave operator.
"""
function scattering(input::InputStruct, k::Union{Complex128,Float64,Int};
    isNonLinear::Bool=false, dispOpt::Bool = false, ψ_init::Array{Complex128,1}=Complex128[],
    F::Array{Float64,1}=[1.], truncate::Bool=false,
    A::Base.SparseArrays.UMFPACK.UmfpackLU=lufact(speye(1,1)), fileName::String = "",
    ftol::Float64=2e-8, iter::Int=150)::
    Tuple{Array{Complex128,1},Base.SparseArrays.UMFPACK.UmfpackLU,InputStruct}

    if !isNonLinear
        ψ, A, input = scattering_l(input, k; A=A, F=F)
    elseif isNonLinear
        ψ, A, input = scattering_nl(input, k; A=A, F=F, dispOpt=dispOpt, ψ_init=ψ_init, ftol=ftol, iter=iter)
    end

    if !isempty(fileName)
        if truncate
            fid = open(fileName,"w")
            serialize(fid, (ψ[input.x̄_inds], input, k) )
            close(fid)
        else
            fid = open(fileName,"w")
            serialize(fid, (ψ, input, k) )
            close(fid)
        end
    end

    if truncate
        return ψ[input.x̄_inds], A, input
    else
        return ψ, A, input
    end
end # end of function scattering





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
### S-MATRIX ROUTINES
################################################################################


"""
S =  smatrix(input, k; isNonLinear=false, F=[1.], dispOpt = true,
                                    N=1, N_Type="D", ψ_init = [], fileName = "")

    N is the number of steps to go from D0 = 0 to given D0
"""
function smatrix(input::InputStruct, k::Union{Complex128,Float64,Int};
    channels::Array{Int,1}=Array(1:length(input.channels)),
    isNonLinear::Bool=false, F::Array{Float64,1}=[1.], dispOpt::Bool=true,
    fileName::String = "", N::Int=1, N_Type::String="D",
    ψ_init::Array{Complex128,1}=Complex128[], parallel::Bool=false)::
    Array{Complex128,4}

    K = [complex(1.0*k)]
    S = smatrix(input, K; channels=channels, isNonLinear=isNonLinear, F=F,
                dispOpt=dispOpt, fileName=fileName, N=N, N_type=N_type,
                ψ_init=ψ_init, parallel=parallel)

    return S
end # end of fuction computeS
function smatrix(input::InputStruct, k::Array{Float64,1};
    channels::Array{Int,1}=Array(1:length(input.channels)), isNonLinear::Bool=false,
    F::Array{Float64,1}=[1.], dispOpt::Bool=true, fileName::String = "",
    N::Int=1, N_Type::String="D", ψ_init::Array{Complex128,1}=Complex128[],
    parallel::Bool=false)::Array{Complex128,4}

    S = smatrix(input, complex(k); channels=channels, isNonLinear=isNonLinear, F=F,
                    dispOpt=dispOpt, fileName=fileName, N=N, N_Type=N_tyoe,
                    ψ_init=ψ_init, parallel=parallel)

    return S
end # end of fuction computeS
function smatrix(input::InputStruct, k::Array{Complex128,1};
    channels::Array{Int,1}=Array(1:length(input.channels)), isNonLinear::Bool=false,
    F::Array{Float64,1}=[1.], dispOpt::Bool=true, fileName::String = "",
    N::Int=1, N_Type::String="D", ψ_init::Array{Complex128,1}=Complex128[],
    parallel::Bool=false)::Array{Complex128,4}

    if !isNonLinear & parallel
        S = smatrixL(input, k, parallel; channels=channels, F=F, dispOpt=dispOpt,
                        fileName=fileName)
    elseif !isNonLinear & !parallel
        S = smatrixL(input, k; channels=channels, F=F, dispOpt=dispOpt,
                        fileName=fileName)
    elseif isNonLinear & parallel
        S = smatrixNL(input, k, parallel; channels=channels, N=N, N_Type=N_Type,
                        isNonLinear=isNonLinear, F=F, dispOpt=dispOpt,
                        ψ_init=ψ_init, fileName=fileName)
    elseif isNonLinear & !parallel
        S = smatrixNL(input, k; channels=channels, N=N, N_Type=N_Type,
                        isNonLinear=isNonLinear, F=F, dispOpt=dispOpt,
                        ψ_init=ψ_init, fileName=fileName)
    end

    return S
end # end of fuction computeS



"""
S =  smatrixL(input; F=[1.], dispOpt=true, fileName = "")
"""
function smatrixL(input::InputStruct, k::Array{Complex128,1};
    channels::Array{Int,1}=Array(1:length(input.channels)), F::Array{Float64,1}=[1.],
    dispOpt::Bool=true, fileName::String = "")::Array{Complex128,4}

    nc = length(channels)
    nk = length(k)

    M = length(input.channels)
    S = NaN*ones(Complex128,nk,length(channels),M,1)
    a_original = input.a
    a = zeros(Complex128,M)
    ψ = Array{Complex128}(1)

    for ii in 1:nk
        if (ii/1 == round(ii/1)) & dispOpt
            if typeof(k)<:Real
                printfmtln("Solving for frequency {1:d} of {2:d}, ω = {3:2.3f}.",ii,nk,k[ii])
            else
                printfmtln("Solving for frequency {1:d} of {2:d}, ω = {3:2.3f}{4:+2.3f}i.",ii,nk,real(k[ii]),imag(k[ii]))
            end
        end

        ζ = lufact(speye(1,1))

        for m in 1:nc
            a = 0*a
            a[channels[m]] = 1.
            updateInputs!(input, :a, a)
            ψ, ζ, inputs_s = scattering_l(input, k[ii]; A=ζ)
            for m′ in 1:M
                 cm = analyze_output(inputs_s, k[ii], ψ, m′)
                 S[ii,m,m′,1] = cm
            end
        end

        updateInputs!(input, :a, a_original)

        if !isempty(fileName)
            fid = open(fileName,"w")
            serialize(fid,(S,input,ii))
            close(fid)
        end
    end

    return S
end # end of fuction computeS_linear


"""
S =  computeS(input::InputStruct; N=10, N_Type="D", isNonLinear=false, F=1.,
    dispOpt = true, ψ_init = [], fileName = "")

    N is the number of steps to go from D0 = 0 to given D0
"""
function smatrixNL(input1::InputStruct, k::Array{Complex128,1};
    N::Int=1, N_Type::String="D", isNonLinear::Bool=false, F::Array{Float64,1}=[1.],
    dispOpt::Bool=true, ψ_init::Array{Complex128,1}=Complex128[],
    fileName::String = "")::Array{Complex128,4}

    if !isempty(input1.bc ∩ ["o", "open", "pml_in"])
        input = deepcopy(input1)
        for i in 1:4
            if input.bc[i] in ["o", "open", "pml_in"]
                input.bc[i] = "pml_out"
                updateInputs!(input, :bc, input.bc)
            end
        end
    else
        input = input1
    end

    M = length(input.channels)
    D₀ = input.D₀
    A = input.a

    ψ₋ = Array{Complex128}(prod(input.N_ext))
    ψ  = Array{Complex128}(prod(input.N_ext))

    # if isempty(ψ_init) && (N>1)
    S = NaN*ones(Complex128,length(input.k),N,M,M)
    # else
        # S = NaN*ones(Complex128,length(input.k),1,2,2)
    # end

    if N > 1
        D = linspace(0,D₀,N)
        A_vec = linspace(0.001,1,N)
    else
        D = [input.D₀]
        A_vec = [1]
    end

    for ii in 1:length(input.k)

        k = input.k[ii]

        if (ii/1 == round(ii/1)) & dispOpt
            if typeof(k)<:Real
                printfmtln("Solving for frequency {1:d} of {2:d}, ω = {3:2.3f}.",ii,length(input.k),k)
            else
                printfmtln("Solving for frequency {1:d} of {2:d}, ω = {3:2.3f}{4:+2.3f}i.",ii,length(input.k),real(k),imag(k))
            end
        end

        if isempty(ψ_init)

            for j in 1:N

                if N_Type == "D"
                    updateInputs!(input, :D₀, D[j])
                elseif N_Type == "A"
                    updateInputs!(input, :a, A_vec[j]*A)
                end

                if isNonLinear
                    if isnan(ψ[1]) | j==1
                        ψ, W = computePsi(input,k,isNonLinear = true, F = F)
                    else
                        ψ, W = computePsi(input,k,isNonLinear = true, F = F, ψ_init = ψ)
                    end
                else
                    ψ = 0.
                end

                # compute linearized S-matrix
                ζ = lufact(speye(1,1)) # initialize lufact variable (speeds up repeated inhomogeneous solves)
                for m1 in 1:M
                    # set flux-normalized input amplitudes to unity
                    at = zeros(Complex128,M)
                    at[m1] = 1.
                    updateInputs!(input,:a,at)
                    # solve for scattered field
                    ψ₋, dummy1, ζ, inputs_s = computePsi(input, k; isNonLinear=false, F=F./(1+abs2.(γ(k,input)*ψ)), A=ζ)
                    # analyze into channels
                    println("here 1")
                    for m2 in m1:M
                        dt=0
                        if input.channelBoundaries[m2] in [1,2]
                            t = inputs_s.x₂_ext
                            u = inputs_s.x₁_inds
                            dt = input.dx̄[2]
                        else
                            t = inputs_s.x₁_ext
                            u = inputs_s.x₂_inds
                            dt = input.dx̄[1]
                        end
                        println("here 2")
                        ϕ = zeros(Complex128,length(t))
                        for q in 1:length(t)
                             ϕ[q] = conj(input.incidentWave(k, m2, input.∂R[input.channelBoundaries[m2]], t[ii], input.∂R, input.geometry, input.geoParams)[1])
                        end
                        println("here 3")
                        P = 0
                        if input.channelBoundaries[m2] in [1,3]
                            println("here 4")
                            println(size(ψ₋))
                            println(size(t))
                            println(size(u))
                            P = reshape(ψ₋,:,length(t))[u[1],:]
                            println("here 5")
                        else
                            println("here 6")
                            P = reshape(ψ₋,:,length(t))[u[end],:]
                        end
                        println("here 7")
                        println(size(P))
                        println(size(S))
                        println(j)
                        println(ii)
                        println(size(ϕ))
                        S[ii,j,m1,m2] = sum(P.*ϕ)*dt
                        S[ii,j,m2,m1] = S[ii,j,m1,m2]
                    end
                end
    println("here 8")
                updateInputs!(input,:D₀, D₀)
                updateInputs!(input,:a , A )

                if !isempty(fileName)
                    fid = open(fileName,"w")
                    serialize(fid,(input,D,S,ii,j))
                    close(fid)
                end

            end

        else
            if isNonLinear
                ψ = scattering(input,k,isNonLinear = true, F = F, ψ_init = ψ_init)
            else
                ψ = 0.
            end

            input.a = [1.0,0.0]
            ψ₊, W = scattering(input,k,isNonLinear = false, F = F./(1+abs2.(γ(input,k)*ψ)))

            input.a = [0.0,1.0]
            ψ₋, W = scattering(input,k,isNonLinear = false, F = F./(1+abs2.(γ(input,k)*ψ)), A = W)

            S[1,1,ii,1] = ψ₊[x_inds[1]]*exp(+1im*input.dx*k)
            S[2,1,ii,1] = ψ₊[x_inds[end]]
            S[1,2,ii,1] = ψ₋[x_inds[1]]
            S[2,2,ii,1] = ψ₋[x_inds[end]]*exp(-1im*input.dx*k)

            input.D₀ = D₀
            input.a = A

            if !isempty(fileName)
                foo2(fid) = serialize(fid,(input,D,S,ii,1))
                open(foo2,fileName,"w")
            end

        end

    end
    println("here 9")
    return S
end # end of fuction computeS


################################################################################
###  NONLINEAR INHOMOGENEOUS AUXILLIARIES
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
