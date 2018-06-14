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

    if bc_sig == "Od"
        x = input.dis.x[1] - input.bnd.∂R[2]
        φ = +sqrt(1/k)*exp(+1im*k*x)
        P = reshape(ψ[input.x̄_inds],input.N[1],:)[1,:]
        cm = sum(φ.*P)*input.dx̄[2]
        bm = -input.a[m]
    elseif bc_sig == "dO"
        x = input.x₁[end] - input.∂R[1]
        φ = +sqrt(1/k)*exp(-1im*k*x)
        P = reshape(ψ[input.x̄_inds],input.N[1],:)[end,:]
        cm = sum(φ.*P)*input.dx̄[2]
        bm = -input.a[m]
    elseif bc_sig in ["OO", "II"]
        cm = analyze_into_waveguides(input, k, ψ, m, "out")
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

    if bc_sig == "Od"
        error("Haven't written input analyzer for one-sided input in a waveguide")
    elseif bc_sig == "dO"
        error("Haven't written input analyzer for one-sided input in a waveguide")
    elseif bc_sig in ["OO", "II"]
        cm = analyze_into_waveguides(input, k, ψ, m, "in")
    end

    return cm
end

################################################################################
### Analyzer Subroutines
################################################################################

"""
analyze_into_waveguides(input, k, ψ, m, direction)
"""
function analyze_into_waveguides(input::InputStruct, k::Complex128,
    ψ::Array{Complex{Float64},1}, m::Int, direction::String)::Complex128

    if input.sct.channels[m].side in ["l", "L", "left", "Left"]
        p = ψ[input.dis.x_inds[1]]
        # x = input.bnd.∂R[1]
        # if direction == "in"
        #     phs = exp(-1im*k*(x-input.bnd.∂S[2]))
        # elseif direction == "out"
        #     phs = exp(+1im*k*(x-input.bnd.∂S[2]))
        # end
    elseif input.sct.channels[m].side in ["r", "R", "right", "Right"]
        p = ψ[input.dis.x_inds[end]]
        # x = input.bnd.∂R[2]
        # if direction == "in"
        #     phs = exp(+1im*k*(x-input.bnd.∂S[1]))
        # elseif direction == "out"
        #     phs = exp(-1im*k*(x-input.bnd.∂S[1]))
        # end
    end

    cm = sqrt(real(k))*p#analysis_interpolate(input, x, ψ)

    return cm
end

"""
P = analysis_interpolate(input, X, ψ)
"""
function analysis_interpolate(input::InputStruct, X::Float64,
    ψ::Array{Complex128,1})::Complex128

    p = interpolate(ψ, BSpline(Quadratic(Reflect())), OnGrid())
    Ax = (input.dis.N_PML-1)/(input.dis.x_PML[end]-input.dis.x_PML[1])
    Bx = (input.dis.N_PML*input.dis.x_PML[1]-input.dis.x_PML[end])/(input.dis.N_PML-1)
    X_int = Ax*(X-Bx)
    P = p[X_int]

    return P
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
