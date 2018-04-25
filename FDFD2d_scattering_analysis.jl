"""
s = analyze_output(input, k, ψ, m)
    s is the output coefficient in the mth channel

    S is constructed from s for unit input on each channel
"""
function analyze_output(input::InputStruct, k::Complex128,
    ψ::Array{Complex{Float64},1}, m::Int)::Complex128

    bc_sig = input.bnd.bc_sig

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
            if input.sct.channels[m].side in ["l", "L", "left", "Left"]
                x = input.dis.xy[1][1] - input.bnd.∂R[1]
                phs = exp.(+1im*kₓ*x)
                P = reshape(ψ,input.dis.N_PML[1],:)[input.dis.dN_PML[1]+input.dis.N[1]-1,:]
                ε = input.wgs.ε_PML[input.sct.channels[m].wg]
            elseif input.sct.channels[m].side in ["r", "R", "right", "Right"]
                x = input.dis.xy[1][end] - input.bnd.∂R[2]
                phs = exp.(-1im*kₓ*x)
                P = reshape(ψ,input.dis.N_PML[1],:)[input.dis.dN_PML[1]+2,:]
                ε = input.wgs.ε_PML[input.sct.channels[m].wg]
            end
            φ = reshape(φy,input.dis.N_PML[1],:)[input.dis.dN_PML[1]+2,:]
        elseif input.wgs.dir[input.sct.channels[m].wg] in ["y", "Y"]
            error("Haven't written vertical waveguide code yet.")
        end

        cm = sqrt(kₓ)*phs*sum(φ.*ε.*P)*input.dis.dx[2]
    elseif (bc_sig in ["OOOO", "IIII"])
        cm = analyze_into_angular_momentum(input, k, ψ, m, "out")
    end

    return cm
end

"""
cm = analyze_input(input, k, ψ, m)

    cm is the input power for a CPA input, that is, one where there is no output
"""
function analyze_input(input::InputStruct, k::Complex128,
    ψ::Array{Complex{Float64},1}, m::Int)::Complex128

    bc_sig = input.bnd.bc_sig

    if bc_sig in ["Oddd", "Odnn", "Oddn", "Odnd"]
        error("Haven't written input analyzer for one-sided input in a waveguide")
    elseif bc_sig in ["dOdd", "dOnn", "dOdn", "dOnd"]
        error("Haven't written input analyzer for one-sided input in a waveguide")
    elseif (bc_sig in ["OOOO", "IIII"]) && (!isempty(input.wgs.dir))
        if (input.wgs.dir[input.sct.channels[m].wg] in ["x", "X"])
            kₓ, φy = wg_transverse_y(input, k, m)
            if input.sct.channels[m].side in ["l", "L", "left", "Left"]
                x = input.dis.xy[1][1] - input.bnd.∂R[1]
                phs = exp.(-1im*kₓ*x)
                P = reshape(ψ,input.dis.N_PML[1],:)[input.dis.dN_PML[1]+input.dis.N[1]-1,:]
                ε = input.wgs.ε_PML[input.sct.channels[m].wg]
            elseif input.sct.channels[m].side in ["r", "R", "right", "Right"]
                x = input.dis.xy[1][end] - input.bnd.∂R[2]
                phs = exp.(+1im*kₓ*x)
                P = reshape(ψ,input.dis.N_PML[1],:)[input.dis.dN_PML[1]+2,:]
                ε = input.wgs.ε_PML[input.sct.channels[m].wg]
            end
            φ = reshape(φy,input.dis.N_PML[1],:)[input.dis.dN_PML[1]+2,:]
        elseif input.wgs.dir[input.channels[m].wg] in ["y", "Y"]
            error("Haven't written vertical waveguide code yet.")
        end

        cm = sqrt(kₓ)*phs*sum(φ.*ε.*P)*input.dis.dx[2]
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
