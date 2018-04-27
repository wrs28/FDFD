export analyze_input, analyze_output
"""
s = analyze_output(input, k, ψ, m)
    s is the output coefficient in the mth channel

    S is constructed from s for unit input on each channel
"""
function analyze_output(input::InputStruct, K::Union{Complex128,Float64,Int},
    ψ::Array{Complex{Float64},1}, m::Int)::Complex128

    bc_sig = input.bnd.bc_sig
    k = complex(float(K))

    if bc_sig in ["Oddd", "Odnn", "Oddn", "Odnd"]
        x = input.dis.xy[1] - input.bnd.∂R[2]
        y = input.dis.xy_PML[2]
        φy = quasi_1d_transverse_y.(input,m,y)
        kᵤ = quasi_1d_transverse_y(input,m)
        kₓ = sqrt(k^2 - kᵤ^2)
        φ = +sqrt(1/kₓ)*exp(+1im*kₓ*x)*φy
        P = reshape(ψ[input.x̄_inds],input.N[1],:)[1,:]
        cm = sum(φ.*P)*input.dx̄[2]
        bm = -input.a[m]
    elseif bc_sig in ["dOdd", "dOnn", "dOdn", "dOnd"]
        x = input.x₁[end] - input.∂R[1]
        y = input.x₂_ext
        φy = quasi_1d_transverse_y.(input,m,y)
        kᵤ = quasi_1d_transverse_y(input,m)
        kₓ = sqrt(k^2 - kᵤ^2)
        φ = +sqrt(1/kₓ)*exp(-1im*kₓ*x)*φy
        P = reshape(ψ[input.x̄_inds],input.N[1],:)[end,:]
        cm = sum(φ.*P)*input.dx̄[2]
        bm = -input.a[m]
    elseif (bc_sig in ["OOOO", "IIII"]) && (!isempty(input.wgs.dir))
        if (input.wgs.dir[input.sct.channels[m].wg] in ["x", "X"])
            kₓ, φy = wg_transverse_y(input, k, m)
            φy, φy2 = incident_mode(input, k, m)
            if input.sct.channels[m].side in ["l", "L", "left", "Left"]
                x = input.dis.xy[1][1]
                phs = exp(+1im*kₓ*x)
                P = real(kₓ)*reshape(ψ[input.dis.xy_inds],input.dis.N[1],:)[1,:]
            elseif input.sct.channels[m].side in ["r", "R", "right", "Right"]
                x = input.dis.xy[1][end]
                phs = exp(-1im*kₓ*x)
                P = real(kₓ)*reshape(ψ[input.dis.xy_inds],input.dis.N[1],:)[end,:]
            end
            φ = phs*reshape(φy[input.dis.xy_inds],input.dis.N[1],:)[1,:]*sqrt(1/real(kₓ))
        elseif input.wgs.dir[input.sct.channels[m].wg] in ["y", "Y"]
            error("Haven't written vertical waveguide code yet.")
        end

        cm = sum(conj(φ).*P)*input.dis.dx[2]
    elseif (bc_sig in ["OOOO", "IIII"])
        cm = analyze_into_angular_momentum(input, k, ψ, m, "out")
    end

    return cm
end

"""
cm = analyze_input(input, k, ψ, m)

    cm is the input power for a CPA input, that is, one where there is no output
"""
function analyze_input(input::InputStruct, K::Union{Complex128,Float64,Int},
    ψ::Array{Complex{Float64},1}, m::Int)::Complex128

    bc_sig = input.bnd.bc_sig
    k = complex(float(K))

    if bc_sig in ["Oddd", "Odnn", "Oddn", "Odnd"]
        error("Haven't written input analyzer for one-sided input in a waveguide")
    elseif bc_sig in ["dOdd", "dOnn", "dOdn", "dOnd"]
        error("Haven't written input analyzer for one-sided input in a waveguide")
    elseif (bc_sig in ["OOOO", "IIII"]) && (!isempty(input.wgs.dir))
        if (input.wgs.dir[input.sct.channels[m].wg] in ["x", "X"])
            kₓ, φy = wg_transverse_y(input, k, m)
            φy, φy2 = incident_mode(input, k, m)
            if input.sct.channels[m].side in ["l", "L", "left", "Left"]
                x = input.dis.xy[1][1]
                phs = exp(-1im*kₓ*x)
                P = real(kₓ)*reshape(ψ[input.dis.xy_inds],input.dis.N[1],:)[1,:]
            elseif input.sct.channels[m].side in ["r", "R", "right", "Right"]
                x = input.dis.xy[1][end]
                phs = exp(+1im*kₓ*x)
                P = real(kₓ)*reshape(ψ[input.dis.xy_inds],input.dis.N[1],:)[end,:]
            end
            φ = reshape(φy[input.dis.xy_inds],input.dis.N[1],:)[1,:]*sqrt(1/real(kₓ))
        elseif input.wgs.dir[input.sct.channels[m].wg] in ["y", "Y"]
            error("Haven't written vertical waveguide code yet.")
        end

        cm = sum(conj(φ).*P)*input.dis.dx[2]
    elseif (bc_sig in ["OOOO", "IIII"])
        cm = analyze_into_angular_momentum(input, k, ψ, m, "in")
    end

    return cm
end

################################################################################
### Analyzer Subroutines
################################################################################
"""
analyze_into_angular_momentum(input, k, ψ, m, direction)
"""
function analyze_into_angular_momentum(input::InputStruct, k::Complex128,
    ψ::Array{Complex{Float64},1}, m::Int, direction::String)::Complex128

    nθ = Int(5e3)+1
    θ = linspace(0,2π,nθ)
    dθ = θ[2]-θ[1]

    # R is radius at which to interpolate
    R = (findmin(abs.(input.∂R))[1] + findmin(abs.(input.∂S))[1])/2 # FIX THIS

    # interpolate wavefunction at r=R, result is P(θ)
    p = interpolate(reshape(ψ,input.N_ext[1],:), BSpline(Linear()), OnGrid() )
    X = R*cos.(θ[1:end-1])
    Y = R*sin.(θ[1:end-1])
    X_int = input.N_ext[1]*(X-input.∂R_ext[1])/(input.∂R_ext[2]-input.∂R_ext[1])
    Y_int = input.N_ext[2]*(Y-input.∂R_ext[3])/(input.∂R_ext[4]-input.∂R_ext[3])
    P = [p[X_int[ii],Y_int[ii]] for ii in 1:(nθ-1)]

    q = input.channels[m].tqn

    if direction == "in"
        cm = sum(exp.(-1im*q*θ[1:end-1]).*P)*dθ./(π*hankelh2(q,k*R))
    elseif direction == "out"
        cm = sum(exp.(-1im*q*θ[1:end-1]).*P)*dθ./(π*hankelh1(q,k*R))
    else
        error("Invalid direction.")
    end

    return cm
end


"""
surface_K(ψ, pos, coord, inputs)
"""
function surface_K(ψ::Array{Complex128,Int}, pos::Float64, coord::String, inputs::InputStruct)::Array{Float64,Int}

    if coord ∈ ["x","X"]
        x = inputs.x₁
        N = inputs.N[1]
        dx = inputs.dx̄[1]
    elseif coord ∈ ["y","Y"]
        x = inputs.x₂
        N = inputs.N[2]
        dx = inputs.dx̄[2]
    end

    ind = findmin(abs.(x-pos))[2]

    if abs(x[ind+1]-pos)>abs(x[ind-1]-pos)
        ind1 = ind-1
        ind2 = ind
    else
        ind1 = ind
        ind2 = ind+1
    end

    w1 = (pos-x[ind1])/(x[ind2]-x[ind1])
    w = [w1, 1-w1]

    K = zeros(Float64,N,size(ψ,2))
    for i in 1:size(ψ,2)
        dΨ = zeros(Complex128,N,2)
        if coord=="x"
            Ψ = reshape(ψ[inputs.x̄_inds,i],N,:)[min(ind1-1,ind2-1):max(ind1+1,ind2+1),:]
            dΨ[:,1] = Ψ[3,:]-Ψ[1,:]
            dΨ[:,2] = Ψ[4,:]-Ψ[2,:]
            Ψ = transpose(Ψ)

        elseif coord=="y"
            Ψ = reshape(ψ[inputs.x̄_inds,i],:,N)[:,min(ind1-1,ind2-1):max(ind1+1,ind2+1)]
            return min(ind1-1,ind2-1):max(ind1+1,ind2+1)
            dΨ[:,1] = Ψ[:,3]-Ψ[:,1]
            dΨ[:,2] = Ψ[:,4]-Ψ[:,2]
        end
        K[:,i] = imag((conj(Ψ[:,2:3])*w).*(dΨ*w/dx))
    end

    return K
end # end of surface_K



"""
surface_flux(ψ, pos, coord, inputs)
"""
function surface_flux(ψ::Array{Complex128,Int}, pos::Float64, coord::String, inputs::InputStruct)::Array{Float64,Int}

    if coord ∈ ["x","X"]
        dx = inputs.dx̄[1]
    elseif coord ∈ ["y","Y"]
        dx = inputs.dx̄[2]
    end

    return sum(surface_K(ψ, pos::Float64, coord::String, inputs))*dx
end



"""
compute_loss(ψ,k,inputs)
"""
function compute_loss(ψ,k,inputs)

    loss = zeros(Float64,length(k))

    for i in 1:length(k)
        Ψ = reshape(ψ[inputs.x̄_inds,i],inputs.N[1],:)
        ε = reshape(inputs.ε_sm[inputs.x̄_inds],inputs.N[1],:)
        F = reshape(inputs.F_sm[inputs.x̄_inds],inputs.N[1],:)
        Ψ = (Ψ[1:end-1,:]+Ψ[2:end,:])/2
        Ψ = (Ψ[:,1:end-1]+Ψ[:,2:end])/2
        F = (F[1:end-1,:]+F[2:end,:])/2
        F = (F[:,1:end-1]+F[:,2:end])/2
        ε = (ε[1:end-1,:]+ε[2:end,:])/2
        ε = (ε[:,1:end-1]+ε[:,2:end])/2
        loss[i] = sum(abs2.(Ψ).*imag((ε-1im*inputs.D₀.*F)*k[i]^2))*prod(inputs.dx̄)
    end

    return loss
end
