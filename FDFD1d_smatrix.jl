export smatrix, smatrix_p

################################################################################
### S-MATRIX MAIN ROUTINES
################################################################################
"""
S =  smatrix(input, k; channels, isNonLinear, F, dispOpt, fileName, N, N_type, ψ_init = [])
    N is the number of steps to go from D0 = 0 to given D0
"""
function smatrix(input::InputStruct, k::Union{Complex128,Float64,Int};
    channels::Array{Int,1}=Array(1:length(input.sct.channels)),
    isNonLinear::Bool=false, F::Array{Float64,1}=[1.], dispOpt::Bool=true,
    fileName::String = "", N::Int=1, N_type::String="D",
    ψ_init::Array{Complex128,1}=Complex128[])::Array{Complex128}

    K = [k]
    S = smatrix(input, K; channels=channels, isNonLinear=isNonLinear, F=F,
                dispOpt=dispOpt, fileName=fileName, N=N, N_type=N_type,
                ψ_init=ψ_init)

    return S
end
function smatrix(input::InputStruct, k::Union{
    StepRangeLen{Complex{Float64},Base.TwicePrecision{Complex{Float64}},Base.TwicePrecision{Float64}},
    StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
    StepRange{Int64,Int64} };
    channels::Array{Int,1}=Array(1:length(input.sct.channels)),
    isNonLinear::Bool=false, F::Array{Float64,1}=[1.], dispOpt::Bool=true,
    fileName::String = "", N::Int=1, N_type::String="D",
    ψ_init::Array{Complex128,1}=Complex128[])::Array{Complex128}

    K = Array(complex.(float.(k)))
    S = smatrix(input, K; channels=channels, isNonLinear=isNonLinear, F=F,
            dispOpt=dispOpt, fileName=fileName, N=N, N_type=N_type,
            ψ_init=ψ_init)

            return S
end
function smatrix(input::InputStruct, k::Union{Array{Complex128,1},Array{Float64,1},Array{Int,1}};
    channels::Array{Int,1}=Array(1:length(input.sct.channels)), isNonLinear::Bool=false,
    F::Array{Float64,1}=[1.], dispOpt::Bool=true, fileName::String = "",
    N::Int=1, N_type::String="D", ψ_init::Array{Complex128,1}=Complex128[])::Array{Complex128}

    if !isNonLinear
        S = smatrix_l(input, complex(k); channels=channels, F=F, dispOpt=dispOpt,
                        fileName=fileName)
    else
        S = smatrix_nl(input, complex(k); channels=channels, N=N, N_type=N_type,
                        isNonLinear=isNonLinear, F=F, dispOpt=dispOpt,
                        ψ_init=ψ_init, fileName=fileName)
    end

    return S
end # end of fuction smatrix

"""
S = smatrix_p(input, k; channels, isNonLinear, F, dispOpt, fileName, N, N_type, ψ_init, num_blocks)
    PARALLEL
"""
function smatrix_p(input::InputStruct, k::Union{Complex128,Float64,Int};
    isNonLinear::Bool=false, F::Array{Float64,1}=[1.], dispOpt::Bool=true,
    fileName::String = "", N::Int=1, N_type::String="D",
    ψ_init::Array{Complex128,1}=Complex128[], num_blocks::Int=3)::
    Tuple{SharedArray{Complex128},Channel}

    K = [k]
    S, r = smatrix_p(input, K; isNonLinear=isNonLinear, F=F,
        dispOpt=dispOpt, fileName=fileName, N=N, N_type=N_type,
        ψ_init=ψ_init, num_blocks=num_blocks)

    return S, r
end
function smatrix_p(input::InputStruct, k::Union{
    StepRangeLen{Complex{Float64},Base.TwicePrecision{Complex{Float64}},Base.TwicePrecision{Float64}},
    StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
    StepRange{Int64,Int64} };
    isNonLinear::Bool=false, F::Array{Float64,1}=[1.], dispOpt::Bool=true,
    fileName::String = "", N::Int=1, N_type::String="D",
    ψ_init::Array{Complex128,1}=Complex128[], num_blocks::Int=3)::
    Tuple{SharedArray{Complex128},Channel}

    K = Array(complex.(float.(k)))
    S, r = smatrix_p(input, K; isNonLinear=isNonLinear, F=F,
        dispOpt=dispOpt, fileName=fileName, N=N, N_type=N_type,
        ψ_init=ψ_init, num_blocks=num_blocks)

    return S, r
end
function smatrix_p(input::InputStruct, k::Union{Array{Complex128,1},Array{Float64,1},Array{Int,1}};
    isNonLinear::Bool=false,
    F::Array{Float64,1}=[1.], dispOpt::Bool=true, fileName::String = "", N::Int=1,
    N_type::String="D", ψ_init::Array{Complex128,1}=Complex128[], num_blocks::Int=3)::
    Tuple{SharedArray{Complex128},Channel}

    if !isNonLinear
        S, r = smatrix_lp(input, complex.(float.(k)); F=F, dispOpt=dispOpt,
                        fileName=fileName, num_blocks=num_blocks)
    else
        S, r = smatrix_nlp(input, complex.(float.(k)); N=N, N_type=N_type,
                        isNonLinear=isNonLinear, F=F, dispOpt=dispOpt,
                        ψ_init=ψ_init, fileName=fileName, num_blocks=num_blocks)
    end

    return S, r
end # end of fuction smatrix_p

################################################################################
### S-MATRIX LINEAR ROUTINES
################################################################################

"""
S =  smatrix_l(input; F=[1.], dispOpt=true, fileName = "")
"""
function smatrix_l(input::InputStruct, k::Array{Complex128,1};
    channels::Array{Int,1}=Array(1:length(input.sct.channels)), F::Array{Float64,1}=[1.],
    dispOpt::Bool=true, fileName::String = "")::Array{Complex128,3}

    nc = length(channels)
    nk = length(k)

    M = length(input.sct.channels)
    S = NaN*ones(Complex128,nk,length(channels),M)
    a = zeros(Complex128,M)
    ψ = Array{Complex128}(1)

    ∂ = input.sys.∂C(input.sys.geoParams)
    ∂R = input.bnd.∂R

    for ii in 1:nk
        if (ii/1 == round(ii/1)) & dispOpt
            if imag(k[ii]) < 5e-3
                printfmt("Solving for frequency {1:d} of {2:d}, ω = {3:2.3f}\n",ii,nk,real(k[ii]))
            else
                printfmt("Solving for frequency {1:d} of {2:d}, ω = {3:2.3f}{4:+2.3f}i\n",ii,nk,real(k[ii]),imag(k[ii]))
            end
        end

        ζ = lufact(speye(1,1))

        for m in 1:nc
            a = 0*a
            a[channels[m]] = 1.
            ψ, φ, ζ = scattering_l(input, k[ii], a; H=ζ)
            for m′ in 1:M

                if input.sct.channels[m′].side in ["L", "l", "left", "Left"]
                    phs = exp(-1im*k[ii]*abs(∂[1]-∂R[1]))
                else
                    phs = exp(-1im*k[ii]*abs(∂[end]-∂R[end]))
                end

                bc_sig = input.bnd.bc_sig
                if bc_sig in ["dO","Od"]
                    cm = analyze_output(input, k[ii], ψ + φ, m′)
                elseif bc_sig == "OO"
                    if input.sct.channels[m′].side == input.sct.channels[m].side
                        cm = analyze_output(input, k[ii], ψ, m′)
                    else
                        cm = analyze_output(input, k[ii], ψ + φ, m′)
                    end
                end

                S[ii,m,m′] = phs*cm
            end
        end

        if !isempty(fileName)
            fid = open(fileName,"w")
            serialize(fid,(S,input,ii))
            close(fid)
        end
    end

    return S
end # end of fuction smatrix_l

"""
smatrix_lp(input, k; channels, F, disOpt, fileName, num_blocks)
defaults to running S only on workers, not on head node
"""
function smatrix_lp(input::InputStruct, k::Array{Complex128,1};
    F::Array{Float64,1}=[1.], dispOpt::Bool = true, fileName::String = "",
    num_blocks::Int=3)::Tuple{SharedArray{Complex128,3}, Channel}

    M = length(input.sct.channels)

    if isempty(fileName)
        S = SharedArray{Complex128}((length(k),M,M), pids=workers())
    else
        S = SharedArray{Complex128}(abspath(fileName),(length(k),M,M), pids=workers(), mode="w+")
    end

    for i in 1:length(S)
        S[i]=1im*NaN
    end

    P = procs(S)

    nc = length(input.sct.channels)
    nk = length(k)
    M = num_blocks*length(P)

    mchunks = min(max(1,floor(Int,sqrt(M)*nk/nc)), nk, M)
    nchunks = min(max(1,floor(Int,M/mchunks)), nc, M)

    if dispOpt
        printfmt("{1:d} blocks: {2:d} frequency, {3:d} channel\n", M, mchunks, nchunks)
    end

    r = Channel(M)
    for idx in 1:M

        p = P[ mod(idx-1, length(P)) + 1]

        if idx == 0 || idx > nchunks*mchunks # This worker is not assigned a piece
            a_inds = Array(1:0)
            k_inds = Array(1:0)
        else
            ncsplits = [round(Int, s) for s in linspace(0, nc, nchunks+1)]
            nksplits = [round(Int, s) for s in linspace(0, nk, mchunks+1)]
            a_idx, k_idx = ind2sub(zeros(Int,nchunks,num_blocks*mchunks),idx)

            a_inds = Array(ncsplits[a_idx]+1:ncsplits[a_idx+1])
            k_inds = Array(nksplits[k_idx]+1:nksplits[k_idx+1])
        end

        @async put!(r, remotecall_fetch(smatrix_lp!, p, S, input, k, k_inds, a_inds;
                F=F, dispOpt=dispOpt) )
    end

    return S, r
end # end of function smatrix_lp

"""
smatrix_lp!
"""
function smatrix_lp!(S::SharedArray{Complex128,3}, input::InputStruct, k::Array{Complex128,1},
    k_inds::Array{Int,1}, a_inds::Array{Int,1};
    F::Array{Float64,1}=[1.], dispOpt::Bool=true)::SharedArray{Complex128,3}

    S[k_inds,a_inds,:] = smatrix_l(input, k[k_inds]; channels=a_inds, F=F, dispOpt=dispOpt)

    return S
end # end of function smatrix_lp!


################################################################################
### S-MATRIX NONLINEAR ROUTINES
################################################################################
"""
S =  smatrix_nl(input::InputStruct; N=10, N_type="D", isNonLinear=false, F=1.,
    dispOpt = true, ψ_init = [], fileName = "")

    N is the number of steps to go from D0 = 0 to given D0
"""
function smatrix_nl(input1::InputStruct, k::Array{Complex128,1};
    N::Int=1, N_type::String="D", isNonLinear::Bool=false, F::Array{Float64,1}=[1.],
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

                if N_type == "D"
                    updateInputs!(input, :D₀, D[j])
                elseif N_type == "A"
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
end # end of fuction smatrix_nl
