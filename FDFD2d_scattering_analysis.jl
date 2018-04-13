











"""
s = analyze_output(inputs, k, ψ, m)
    s is the output coefficient in the mth channel

    S is constructed from s for unit inputs on each channel
"""
function analyze_output(inputs::InputStruct, k::Complex128,
    ψ::Array{Complex{Float64},1}, m::Int)::Complex128

    bc_sig = inputs.bc_sig

    if bc_sig in ["Oddd", "Odnn", "Oddn", "Odnd"]
        x = inputs.x₁[1] - inputs.∂R[2]
        y = inputs.x₂_ext
        φy = quasi_1d_transverse_y.(inputs,m,y)
        kᵤ = quasi_1d_transverse_y(inputs,m)
        kₓ = sqrt(k^2 - kᵤ^2)
        φ = +sqrt(1/kₓ)*exp(+1im*kₓ*x)*φy
        P = reshape(ψ[inputs.x̄_inds],inputs.N[1],:)[1,:]
        cm = sum(φ.*P)*inputs.dx̄[2]
        bm = -inputs.a[m]
    elseif bc_sig in ["dOdd", "dOnn", "dOdn", "dOnd"]
        x = inputs.x₁[end] - inputs.∂R[1]
        y = inputs.x₂_ext
        φy = quasi_1d_transverse_y.(inputs,m,y)
        kᵤ = quasi_1d_transverse_y(inputs,m)
        kₓ = sqrt(k^2 - kᵤ^2)
        φ = +sqrt(1/kₓ)*exp(-1im*kₓ*x)*φy
        P = reshape(ψ[inputs.x̄_inds],inputs.N[1],:)[end,:]
        cm = sum(φ.*P)*inputs.dx̄[2]
        bm = -inputs.a[m]
    elseif (bc_sig in ["OOOO", "IIII"]) && (!isempty(inputs.wgd))
        if (inputs.wgd[inputs.channels[m].wg] in ["x", "X"])
            kₓ, φy = wg_transverse_y(inputs, k, m)
            if inputs.channels[m].side in ["l", "L", "left", "Left"]
                x = inputs.x₁[1] - inputs.∂R[1]
                phs = exp.(+1im*kₓ*x)
                xb = inputs.x₁[1] - inputs.∂R[2] # ballistic
                phsb = exp(-1im*kₓ*xb)
                P = reshape(ψ[inputs.x̄_inds],inputs.N[1],:)[1,:]
                ε = inputs.ε_sm[1,inputs.x₂_inds]
            elseif inputs.channels[m].side in ["r", "R", "right", "Right"]
                x = inputs.x₁[end] - inputs.∂R[2]
                phs = exp.(-1im*kₓ*x)
                xb = inputs.x₁[end] - inputs.∂R[1] # ballistic
                phsb = exp(+1im*kₓ*xb)
                P = reshape(ψ[inputs.x̄_inds],inputs.N[1],:)[end,:]
                ε = inputs.ε_sm[end,inputs.x₂_inds]
            end
            φ = reshape(φy[inputs.x̄_inds],inputs.N[1],:)[1,:]
        elseif inputs.channels[m].wgd in ["y", "Y"]
            error("Haven't written vertical waveguide code yet.")
        end

        wg_bool = [inputs.channels[q].wg for q in 1:length(inputs.channels)] .== inputs.channels[m].wg
        tqn_bool = [inputs.channels[q].tqn for q in 1:length(inputs.channels)] .== inputs.channels[m].tqn
        side_bool = [inputs.channels[q].side for q in 1:length(inputs.channels)] .== inputs.channels[m].side
        wg_ind = find(wg_bool .& tqn_bool .& side_bool)
        wg_bal_ind = find(wg_bool .& tqn_bool .& .!side_bool)
        if (length(wg_ind)>1) | (length(wg_bal_ind)>1)
            error("Channels not uniquely defined.")
        end

        cm = sqrt(kₓ)*phs*sum(φ.*ε.*P)*inputs.dx̄[2]
        bm = inputs.a[wg_bal_ind[1]]*phsb
    elseif (bc_sig in ["OOOO", "IIII"])
        cm = analyze_into_angular_momentum(inputs, k, ψ, m, "out")
    end

    return cm
end


"""
cm = analyze_input(inputs, k, ψ, m)

    cm is the input power for a CPA input, that is, one where there is no output
"""
function analyze_input(inputs::InputStruct, k::Complex128,
    ψ::Array{Complex{Float64},1}, m::Int)::Complex128

    bc_sig = inputs.bc_sig

    if bc_sig in ["Oddd", "Odnn", "Oddn", "Odnd"]
        error("Haven't written input analyzer for one-sided input in a waveguide")
    elseif bc_sig in ["dOdd", "dOnn", "dOdn", "dOnd"]
        error("Haven't written input analyzer for one-sided input in a waveguide")
    elseif (bc_sig in ["OOOO", "IIII"]) && (!isempty(inputs.wgd))
        if (inputs.wgd[inputs.channels[m].wg] in ["x", "X"])
            kₓ, φy = wg_transverse_y(inputs, k, m)
            if inputs.channels[m].side in ["l", "L", "left", "Left"]
                x = inputs.x₁[1] - inputs.∂R[1]
                phs = exp.(-1im*kₓ*x)
                P = reshape(ψ[inputs.x̄_inds],inputs.N[1],:)[1,:]
                ε = inputs.ε_sm[1,inputs.x₂_inds]
            elseif inputs.channels[m].side in ["r", "R", "right", "Right"]
                x = inputs.x₁[end] - inputs.∂R[2]
                phs = exp.(+1im*kₓ*x)
                P = reshape(ψ[inputs.x̄_inds],inputs.N[1],:)[end,:]
                ε = inputs.ε_sm[end,inputs.x₂_inds]
            end
            φ = reshape(φy[inputs.x̄_inds],inputs.N[1],:)[1,:]
        elseif inputs.channels[m].wgd in ["y", "Y"]
            error("Haven't written vertical waveguide code yet.")
        end

        wg_bool = [inputs.channels[q].wg for q in 1:length(inputs.channels)] .== inputs.channels[m].wg
        tqn_bool = [inputs.channels[q].tqn for q in 1:length(inputs.channels)] .== inputs.channels[m].tqn
        side_bool = [inputs.channels[q].side for q in 1:length(inputs.channels)] .!== inputs.channels[m].side
        wg_ind = find(wg_bool .& tqn_bool .& side_bool)
        if length(wg_ind) > 1
            error("Channels not uniquely defined.")
        end

        cm = sqrt(kₓ)*phs*sum(φ.*ε.*P)*inputs.dx̄[2]

    elseif (bc_sig in ["OOOO", "IIII"])
        cm = analyze_into_angular_momentum(inputs, k, ψ, m, "in")
    end

    return cm
end


################################################################################
### Analyzer Subroutines
################################################################################

"""
analyze_into_angular_momentum(inputs, k, ψ, m, direction)
"""
function analyze_into_angular_momentum(inputs::InputStruct, k::Complex128,
    ψ::Array{Complex{Float64},1}, m::Int, direction::String)::Complex128

    nθ = Int(5e3)+1
    θ = linspace(0,2π,nθ)
    dθ = θ[2]-θ[1]

    # R is radius at which to interpolate
    R = (findmin(abs.(inputs.∂R))[1] + findmin(abs.(inputs.∂S))[1])/2 # FIX THIS

    # interpolate wavefunction at r=R, result is P(θ)
    p = interpolate(reshape(ψ,inputs.N_ext[1],:), BSpline(Linear()), OnGrid() )
    X = R*cos.(θ[1:end-1])
    Y = R*sin.(θ[1:end-1])
    X_int = inputs.N_ext[1]*(X-inputs.∂R_ext[1])/(inputs.∂R_ext[2]-inputs.∂R_ext[1])
    Y_int = inputs.N_ext[2]*(Y-inputs.∂R_ext[3])/(inputs.∂R_ext[4]-inputs.∂R_ext[3])
    P = [p[X_int[ii],Y_int[ii]] for ii in 1:(nθ-1)]

    q = inputs.channels[m].tqn

    if direction == "in"
        cm = sum(exp.(-1im*q*θ[1:end-1]).*P)*dθ./(π*hankelh2(q,k*R))
    elseif direction == "out"
        cm = sum(exp.(-1im*q*θ[1:end-1]).*P)*dθ./(π*hankelh1(q,k*R))
    else
        error("Invalid direction.")
    end

    return cm
end
