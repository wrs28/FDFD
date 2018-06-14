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

    N = input.dis.N_PML
    φ = zeros(Complex128, N)

    x = input.dis.x_PML
    M = input.bnd.∂R[1] .≤ x .≤ input.bnd.∂R[2]

    bc_sig = input.bnd.bc_sig
    if bc_sig == "Od"
        # metallic waveguide in x-direction, open on left
        x = input.dis.x_PML[1] - input.bnd.∂R[2]
        L = input.ℓ
        φ₊ = +exp(2L*imag(k))*sqrt(1/real(k))*exp(+1im*k*x)
        φ₋ = -exp(2L*imag(k))*sqrt(1/real(k))*exp(-1im*k*x)
    elseif bc_sig == "dO"
        # metallic waveguide in x-direction, open on right
        x = input.x̄_sct[1] - input.∂R[1]
        L = input.ℓ[1]
        φy = quasi_1d_transverse_y.(input,m,y)
        kᵤ = quasi_1d_transverse_y(input,m)
        kₓ = sqrt(k^2 - kᵤ^2)
        φ₊ = +exp(2L*imag(kₓ))*sqrt(1/real(kₓ))*exp(-1im*kₓ*x).*φy
        φ₋ = -exp(2L*imag(kₓ))*sqrt(1/real(kₓ))*exp(+1im*kₓ*x).*φy
    elseif bc_sig in ["OO", "II"]
        ∂ = input.sys.∂C(input.sys.geoParams)
        if input.sct.channels[m].side in ["l", "L", "left", "Left"]
            φ₊ = sqrt(1/real(k))*exp.(+1im*k*(x-∂[1]))
        elseif input.sct.channels[m].side in ["r", "R", "right", "Right"]
            φ₊ = sqrt(1/real(k))*exp.(-1im*k*(x-∂[end]))
        end
        φ₋ = copy(φ₊)
    end

    reset_bc!(input,bc_original)

    return φ₊, φ₋.*M
end # end of function incident_mode
