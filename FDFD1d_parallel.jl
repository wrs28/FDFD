################################################################################
### Linear Eigenvalue, Embarrassingly Parallel over Parameter Space
################################################################################
export eig_kl

function eig_kl(input::InputStruct, k::Union{Complex128,Float64,Int},
    fields::Array{Symbol,1}, field_inds::Array{Int,1}, field_vals::Array{Array{Float64,1},1};
    nk::Int=1, F::Array{Float64,1}=[1.], truncate::Bool=false,
    ψ_init::Array{Complex128,1}=Complex128[], fileName::String="")::Tuple{SharedArray,Channel}

    dims = tuple(nk, length.(field_vals)...)
    if isempty(fileName)
        K = SharedArray{Complex128}(dims)
    else
        K = SharedArray{Complex128}(abspath(fileName), dims)
    end

    r = Channel(length(procs(K)))
    for p in procs(K)
        @async put!(r, remotecall_fetch(eig_klp!, p, K, input, complex(float(k)), fields, field_inds,
                                field_vals, nk, F, truncate, ψ_init))
    end

    return K,r
end  # end of function eig_kl, parallel computation over parameters
function eig_kl(input::InputStruct, k::Union{Complex128,Float64,Int}, k_type::String,
    fields::Array{Symbol,1}, field_inds::Array{Int,1}, field_vals::Array{Array{Float64,1},1};
    nk::Int=1, F::Array{Float64,1}=[1.], truncate::Bool=false,
    ψ_init::Array{Complex128,1}=Complex128[], direction::Array{Int,1}=[1,0], fileName::String="")::Tuple{SharedArray,Channel}

    input1 = set_bc(input, k_type, direction)

    K,r = eig_kl(input1, k, fields, field_inds, field_vals; nk=nk, F=F,
                truncate=truncate, ψ_init=ψ_init, fileName=fileName)
end # end of function eig_kl, parallel computation over parameters with modified boundaries
function eig_klp!(K::SharedArray, input::InputStruct, k::Complex128,
    fields::Array{Symbol,1}, field_inds::Array{Int,1}, field_vals::Array{Array{Float64,1},1},
    nk::Int, F::Array{Float64,1}, truncate::Bool, ψ_init::Array{Complex128,1})

    inds = p_range(K)
    subs = ind2sub(size(K)[2:end],inds)

    f = fieldnames(input)
    top_field = copy(fields)
    for i in 1:length(f)
        g = fieldnames(getfield(input,f[i]))
        for m in 1:length(fields)
            if fields[m] in g
                top_field[m] = f[i]
            end
        end
    end

    for i in 1:length(inds)
        for f in 1:length(fields)
            if !isempty(size( getfield( getfield(input,top_field[f]) , fields[f]) ))
                vals_temp = getfield( getfield(input, top_field[f]) , fields[f])
                vals_temp[field_inds[f]] = field_vals[f][subs[f][i]]
                input = update_input(input,fields[f],vals_temp)
            else
                input = update_input(input,fields[f],field_vals[f][subs[f][i]])
            end
        end

        K[:,[subs[j][i] for j in 1:length(subs)]...], ψ = eig_kl(input, k; nk=nk, F=F, truncate=truncate, ψ_init=ψ_init)
    end

    return
end # end of function eig_klp! for parallel computation over parameters
"""
inds = p_range(q)
inds = p_range(q, dim)

    q is a shared array.
    returns the indices to be computed on this process.
"""
function p_range(q::SharedArray)::Array{Int,1}
    idx = indexpids(q)
    if idx == 0 # This worker is not assigned a piece
        return 1:0
    end
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in linspace(0, prod(size(q)[2:end]) , nchunks+1)]
    splits[idx]+1:splits[idx+1]
end


################################################################################
### Linear Eigenvalue, Bootstrapped Parallel over Parameter Space
################################################################################

function eig_kl(input::InputStruct, k::Array{Complex128,1},
    fields::Array{Symbol,1}, field_inds::Array{Int,1}, field_vals::Array{Array{Float64,1},1};
    F::Array{Float64,1}=[1.], truncate::Bool=true,
    ψ_init::Array{Complex128,1}=Complex128[], dispOpt::Bool=true, fileName::String="")::
    Tuple{SharedArray,Channel}

    nk = length(k)
    dims = tuple(nk, length.(field_vals)...)
    if isempty(fileName)
        K = SharedArray{Complex128}(dims)
    else
        K = SharedArray{Complex128}(abspath(fileName),dims)
    end
    ψ = Complex128[]

    f = fieldnames(input)
    top_field = copy(fields)
    for i in 1:length(f)
        g = fieldnames(getfield(input,f[i]))
        for m in 1:length(fields)
            if fields[m] in g
                top_field[m] = f[i]
            end
        end
    end

    input1 = deepcopy(input)
    for i in 1:nk
            if !isempty(size( getfield( getfield(input1,top_field[1]) , fields[1]) ))
                vals_temp = getfield( getfield(input1, top_field[1]) , fields[1])
                vals_temp[field_inds[1]] = field_vals[1][1]
                input1 = update_input(input1,fields[1],vals_temp)
            else
                input1 = update_input(input1,fields[1],field_vals[1][1])
            end
                k_temp, ψ_temp = eig_k(input1, k[i]; nk=1, F=F, truncate=truncate, ψ_init=ψ_init)
                K[i,1,ones(Int64,ndims(K)-2)...] = k_temp[1]
                ψ = ψ_temp[:,1]
    end

    r = Channel(length(procs(K)))
    for d in 2:ndims(K)
        if dispOpt
            println("Computing dimension $(d-1)")
        end
        if d < ndims(K)
            @sync for p in procs(K)
                @async remotecall_fetch(eig_klp!, p, K, input, fields, field_inds,
                                        field_vals, d, F, truncate, ψ_init)
            end
        else

            for p in procs(K)
                @async put!(r,remotecall_fetch(eig_klp!, p, K, input, fields, field_inds,
                                        field_vals, d, F, truncate, ψ_init))
            end
        end
    end

    return K, r
end # end of function eig_kl, bootstrapped parallel computation over parameters
function eig_kl(input::InputStruct, k::Array{Complex128,1}, k_type::String,
    fields::Array{Symbol,1}, field_inds::Array{Int,1}, field_vals::Array{Array{Float64,1},1};
    F::Array{Float64,1}=[1.], truncate::Bool=true,
    ψ_init::Array{Complex128,1}=Complex128[], dispOpt::Bool=true,  direction::Array{Int,1}=[1,0],
    fileName::String="")::Tuple{SharedArray,Channel}

    input, bc_original = set_bc(input,k_type, direction)

    K, r = eig_kl(input, k, fields, field_inds, field_vals; F=F,
        truncate=truncate, ψ_init=ψ_init, dispOpt=dispOpt, fileName=fileName)

    reset_bc!(input, bc_original)

    return K, r
end # end of function eig_kl, bootstrapped parallel computation over parameters with modified boundaries
function eig_klp!(K::SharedArray, input::InputStruct,
    fields::Array{Symbol,1}, field_inds::Array{Int,1}, field_vals::Array{Array{Float64,1},1},
    dim::Int64, F::Array{Float64,1}, truncate::Bool, ψ_init::Array{Complex128,1})

    inds = p_range(K,dim)
    subs = ind2sub(size(K)[1:dim-1],inds)

    f = fieldnames(input)
    top_field = copy(fields)
    for i in 1:length(f)
        g = fieldnames(getfield(input,f[i]))
        for m in 1:length(fields)
            if fields[m] in g
                top_field[m] = f[i]
            end
        end
    end

    for d in 2:size(K,dim)
        for i in 1:length(inds)
            for f in 1:length(fields)
                if f < dim-1
                    val_ind = subs[1+f][i]
                elseif f == dim-1
                    val_ind = d
                else
                    val_ind = 1
                end
                if !isempty(size( getfield( getfield(input, top_field[f]) , fields[f]) ))
                    vals_temp = getfield(getfield(input,top_field[f]),fields[f])
                    vals_temp[field_inds[f]] = field_vals[f][val_ind]
                    input = update_input(input,fields[f],vals_temp)
                else
                    input = update_input(input,fields[f],field_vals[f][val_ind])
                end
            end
            k_temp, ψ = eig_kl(input, K[[subs[j][i] for j in 1:length(subs)]..., d-1, ones(Int64,ndims(K)-dim)...]; nk=1, F=F, truncate=truncate, ψ_init=ψ_init)
            K[[subs[j][i] for j in 1:length(subs)]..., d, ones(Int64,ndims(K)-dim)...] = k_temp[1]
        end
    end

    return
end # end of function eig_klp! for bootstrapped parallel computation over parameters
"""
inds = p_range(q)
inds = p_range(q, dim)

    q is a shared array.
    returns the indices to be computed on this process.
"""
function p_range(q::SharedArray, dim::Int64)::Array{Int,1}
    idx = indexpids(q)
    if idx == 0 # This worker is not assigned a piece
        return 1:0
    end
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in linspace(0, prod(size(q)[1:dim-1]) , nchunks+1)]
    splits[idx]+1:splits[idx+1]
end


################################################################################
### NL2 Eigenvalue routines, parallel quadrature
################################################################################


"""
k = eigKNL(input, kc, Radii; nk=3, Nq=100, F=[1.], R_min=.01, rank_tol=1e-8)

    Compute eigenfrequency with dispersion, using contour integration. BC's set
    by input.bc

    Contour is centered on kc, Radii = (x-diameter, y-diameter).

    nk is an upper bound on the number of eigenfrequencies contained in the contour.

    Nq is the number of contour quadrature points.

    Parallelizes quadrature.
"""
function eig_knlp(input::InputStruct, kc::Union{Complex128,Float64,Int},
    Radii::Tuple{Float64,Float64};
    nk::Int=3, Nq::Int=100, F::Array{Float64,1}=[1.],
    R_min::Float64=.01, rank_tol::Float64=1e-8)::Array{Complex128,1}

    N_ext = prod(input.N_ext); ε_sm = input.ε_sm; F_sm = input.F_sm
    D₀ = input.D₀; γ⟂ = input.γ⟂; k₀ = input.k₀

    k = Complex128(kc)

    ∇² = laplacian(k,input)

    M = rand(N_ext,nk)
    ϕ = 2π*(0:1/Nq:(1-1/Nq))
    Ω = k + Radii[1]*cos.(ϕ) + 1im*Radii[2]*sin.(ϕ)

    θ =  angle(k₀-1im*γ⟂-k)
    flag = abs(k₀-1im*γ⟂-k) < rad(Radii[1],Radii[2],θ)
    if flag
        RR = 2*R_min
        ΩΩ = k₀-1im*γ⟂ + (RR/2)*cos.(ϕ) + 1im*(RR/2)*sin.(ϕ)
    end

    AA = @parallel (+) for i in 1:Nq

        k′ = Ω[i]
        k′² = k′^2

        if (i > 1) & (i < Nq)
            dk′ = (Ω[i+1]-Ω[i-1]  )/2
        elseif i == Nq
            dk′ = (Ω[1]  -Ω[end-1])/2
        elseif i == 1
            dk′ = (Ω[2]  -Ω[end]  )/2
        end

        ɛk′² = sparse(1:N_ext, 1:N_ext, ɛ_sm[:]*k′²            , N_ext, N_ext, +)
        χk′² = sparse(1:N_ext, 1:N_ext, D₀*(γ⟂/(k′-k₀+1im*γ⟂))*F.*F_sm[:]*k′², N_ext, N_ext, +)

        A  = (∇²+ɛk′²+χk′²)\M
        A₀ = A*dk′/(2π*1im)
        A₁ = A*k′*dk′/(2π*1im)

        if flag
            kk′ = ΩΩ[i]
            kk′² = kk′^2
            if (i > 1) & (i < Nq)
                dkk′ = (ΩΩ[i+1]-ΩΩ[i-1]  )/2
            elseif i == Nq
                dkk′ = (ΩΩ[1]  -ΩΩ[end-1])/2
            elseif i == 1
                dkk′ = (ΩΩ[2]  -ΩΩ[end]  )/2
            end
            ɛkk′² = sparse(1:N_ext, 1:N_ext, ɛ_sm[:]*kk′², N_ext, N_ext, +)
            χkk′² = sparse(1:N_ext, 1:N_ext, D₀*(γ⟂/(kk′-k₀+1im*γ⟂))*F.*F_sm[:]*kk′², N_ext, N_ext, +)

            AA  = (∇²+ɛkk′²+χkk′²)\M
            AA₀ = AA*dkk′/(2π*1im)
            AA₁ = AA*kk′*dkk′/(2π*1im)

            A₀ = A₀-AA₀
            A₁ = A₁-AA₁
        end

        [A₀ A₁]

    end

    A₀ = AA[:,1:nk]
    A₁ = AA[:,nk + (1:nk)]

    P = svdfact(A₀,thin = true)
    temp = find(P[:S] .< rank_tol)
    if isempty(temp)
        println("Error. Need more nevals")
        return [NaN]
    else
        k = temp[1]-1
    end

    B = (P[:U][:,1:k])'*A₁*(P[:Vt][1:k,:])'*diagm(1./P[:S][1:k])

    D,V = eig(B)

    return D
end # end of eig_knlp, contour integral parallel
function eig_knlp(input::InputStruct, kc::Union{Complex128,Float64,Int},
    Radii::Tuple{Float64,Float64}, k_type::String;
    nk::Int=3, Nq::Int=100, F::Array{Float64,1}=[1.],
    R_min::Float64=.01, rank_tol::Float64=1e-8, direction::Array{Int,1}=[1,0])::Array{Complex128,1}

    input1 = set_bc(input,k_type, direction)

    k = eig_knlp(input1, kc, Radii; nk=nk, Nq=Nq, F=F, R_min=R_min, rank_tol=rank_tol)
end # end of eig_knlp, contour integral parallel with modified boundary
