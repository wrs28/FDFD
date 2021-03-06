"""
incident_mode(input, k, m)

    k frequency

    m is channel number, as defined in input

    φ₊ is the incident field, except in the 2-dim case where it is the superposition
    of incident and outgoing fields in the absence of scattering.

    φ₋ is the outgoing field in the absence of scattering
"""
function incident_mode(input::InputStruct, k::Complex128, m::Int)::
    Tuple{Array{Complex128,1},Array{Complex128,1}}

    bc_original = set_bc!(input)

    N = prod(input.dis.N_PML)
    φ = zeros(Complex128, N)

    x = input.dis.XY_PML[1][:]
    y = input.dis.XY_PML[2][:]
    M = (input.bnd.∂R[1] .≤ x .≤ input.bnd.∂R[2]) .& (input.bnd.∂R[3] .≤ y .≤ input.bnd.∂R[4])

    bc_sig = input.bnd.bc_sig
    if bc_sig in ["Oddd", "Odnn", "Oddn", "Odnd"]
        # metallic waveguide in x-direction, open on left
        y = input.dis.x_PML[2]
        x = input.dis.x_PML[1] - input.bnd.∂R[2]
        φy = quasi_1d_transverse_y.(input,m,y)
        kᵤ = quasi_1d_transverse_y(input,m)
        kₓ = sqrt(k^2 - kᵤ^2)
        L = input.ℓ[1]
        φ₊ = +exp(2L*imag(kₓ))*sqrt(1/real(kₓ))*exp(+1im*kₓ*x).*φy
        φ₋ = -exp(2L*imag(kₓ))*sqrt(1/real(kₓ))*exp(-1im*kₓ*x).*φy
    elseif bc_sig in ["dOdd", "dOnn", "dOdn", "dOnd"]
        # metallic waveguide in x-direction, open on right
        y = input.x̄_sct[2]
        x = input.x̄_sct[1] - input.∂R[1]
        L = input.ℓ[1]
        φy = quasi_1d_transverse_y.(input,m,y)
        kᵤ = quasi_1d_transverse_y(input,m)
        kₓ = sqrt(k^2 - kᵤ^2)
        φ₊ = +exp(2L*imag(kₓ))*sqrt(1/real(kₓ))*exp(-1im*kₓ*x).*φy
        φ₋ = -exp(2L*imag(kₓ))*sqrt(1/real(kₓ))*exp(+1im*kₓ*x).*φy
    elseif (bc_sig in ["OOOO", "IIII"]) && (!isempty(input.wgs.dir))
        # dielectric waveguide
        kₓ, φy = wg_transverse_y(input, k, m)
        if input.sct.channels[m].side in ["l", "L", "left", "Left"]
            φ₊ = sqrt(1/real(kₓ))*exp.(+1im*kₓ*x).*φy
        elseif input.sct.channels[m].side in ["r", "R", "right", "Right"]
            φ₊ = sqrt(1/real(kₓ))*exp.(-1im*kₓ*x).*φy
        end
        φ₋ = copy(φ₊)
    elseif (bc_sig in ["OOOO", "IIII"])
        r = sqrt.(x.^2 + y.^2)
        θ = atan2.(y,x)
        q = input.sct.channels[m].tqn
        φ₊ = exp.(1im*q*θ).*besselj.(q,k*r)
        φ₋ = exp.(1im*q*θ).*hankelh1.(q,k*r + k*(r.==0).*eps.(r) )/2
    end



    reset_bc!(input,bc_original)

    return φ₊, φ₋.*M
end # end of function incident_mode

"""
kᵤ, φᵤ = wg_transverse_y(input, k, m)

    k is the frequency
    m is the channel number as defined in input.

    kᵤ is the transverse eigenvalue
    φᵤ is the transverse mode defined over whole domain
"""
function wg_transverse_y(input::InputStruct, k::Complex128, m::Int)::
    Tuple{Complex128, Array{Complex128,1}}

    bc_original = set_bc!(input)

    n = input.sct.channels[m].wg
    q = input.sct.channels[m].tqn

    wg_pos_ind = 3
    ind = findmin(input.wgs.pos[n]-input.dis.xy[2])[2] + wg_pos_ind

    N = input.dis.N_PML[2]
    k² = k^2
    ε = input.wgs.ε_PML[n]
    εk² = sparse(1:N, 1:N, ε[:]*k², N, N, +)
    ∇₁², ∇₂² = laplacians(input,k)

    nev, k_scale = wg_transverse_y_params(q)

    nev = 4 + 2*q
    kₓ²,φ = eigs(∇₂²+εk², nev=nev, sigma=k_scale*k², which = :LM)
    perm1 = sortperm(kₓ²; by = x -> imag(sqrt.(x)), rev=false)
    perm2 = sortperm(kₓ²[perm1[1:q]]; by = x -> real(sqrt.(x)), rev=true)
    perm = perm1[(1:q)[perm2[1]]]
    φ_temp = φ[:,perm]
    φ_temp = φ_temp*( conj(φ_temp[ind])/abs(φ_temp[ind]) ) #makes field positive at wg_pos_ind
    φ_temp = φ_temp./sqrt(sum(abs2,φ_temp)*input.dis.dx[2])

    φy = repeat(transpose(φ_temp); outer=(input.dis.N_PML[1],1))[:]

    reset_bc!(input, bc_original)

    return (sqrt.(kₓ²[perm]), φy)
end # end of function wg_transverse_y




















"""
kₓ = quasi_1d_transverse_x(inputs, m)
    OR
φ = quasi_1d_transverse_x(inputs, m, x)
"""
function quasi_1d_transverse_x(inputs::InputStruct, m::Int)::Float64

    ℓ = inputs.ℓ[1]
    bc = inputs.bc
    q = inputs.channels[m].tqn

    if bc[1:2] == ["n", "n"]
        kₓ = (q-1)*π/ℓ
    elseif bc[1:2] == ["n", "d"]
        kₓ = (q-1/2)*π/ℓ
    elseif bc[1:2] == ["d", "n"]
        kₓ = (q-1/2)*π/ℓ
    elseif bc[1:2] == ["d", "d"]
        kₓ = q*π/ℓ
    end

    return kₓ
end
function quasi_1d_transverse_x(inputs::InputStruct, m::Int, x::Float64)::Float64

    ℓ = inputs.ℓ[1]
    bc = inputs.bc
    kₓ = quasi_1d_transverse_y(inputs,m)

    if bc[1:2] == ["n", "n"]
        φ = sqrt(2/ℓ)*cos.(kₓ*(x-inputs.∂R[1]))
    elseif bc[1:2] == ["n", "d"]
        φ = sqrt(2/ℓ)*cos.(kₓ*(x-inputs.∂R[1]))
    elseif bc[1:2] == ["d", "n"]
        φ = sqrt(2/ℓ)*sin.(kₓ*(x-inputs.∂R[1]))
    elseif bc[1:2] == ["d", "d"]
        φ = sqrt(2/ℓ)*sin.(kₓ*(x-inputs.∂R[1]))
    end

    return φ
end


"""
kᵤ = quasi_1d_transverse_y(inputs, m)
    OR
φ = quasi_1d_transverse_y(inputs, m, y)
"""
function quasi_1d_transverse_y(inputs::InputStruct, m::Int)::Float64

    ℓ = inputs.ℓ[2]
    bc = inputs.bc
    q = inputs.channels[m].tqn

    if bc[3:4] == ["n", "n"]
        kᵤ = (q-1)*π/ℓ
    elseif bc[3:4] == ["n", "d"]
        kᵤ = (q-1/2)*π/ℓ
    elseif bc[3:4] == ["d", "n"]
        kᵤ = (q-1/2)*π/ℓ
    elseif bc[3:4] == ["d", "d"]
        kᵤ = q*π/ℓ
    end

    return kᵤ
end
function quasi_1d_transverse_y(inputs::InputStruct, m::Int, y::Float64)::Float64

    ℓ = inputs.ℓ[2]
    bc = inputs.bc
    kᵤ = quasi_1d_transverse_y(inputs,m)

    if bc[3:4] == ["n", "n"]
        φ = sqrt(2/ℓ)*cos.(kᵤ*(y-inputs.∂R[3]))
    elseif bc[3:4] == ["n", "d"]
        φ = sqrt(2/ℓ)*cos.(kᵤ*(y-inputs.∂R[3]))
    elseif bc[3:4] == ["d", "n"]
        φ = sqrt(2/ℓ)*sin.(kᵤ*(y-inputs.∂R[3]))
    elseif bc[3:4] == ["d", "d"]
        φ = sqrt(2/ℓ)*sin.(kᵤ*(y-inputs.∂R[3]))
    end

    return φ
end
