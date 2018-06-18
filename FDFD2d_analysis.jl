export analyze_input, analyze_output, surface_flux, compute_loss
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
        cm = analyze_into_waveguides(input, k, ψ, m, "out")
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
        cm = analyze_into_waveguides(input, k, ψ, m, "in")
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

    nθ = analysis_quadrature_defaults()
    θ = linspace(0,2π,nθ)
    dθ = θ[2]-θ[1]

    # R is radius at which to interpolate
    R = (findmin(abs.(input.bnd.∂R))[1] + findmin(abs.(input.sct.∂S))[1])/2 # FIX THIS
    X = R*cos.(θ[1:end-1])
    Y = R*sin.(θ[1:end-1])

    # interpolate wavefunction at r=R, result is P(θ)
    P = analysis_interpolate(input, X, Y, ψ)

    q = input.sct.channels[m].tqn

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
analyze_into_waveguides(input, k, ψ, m, direction)
"""
function analyze_into_waveguides(input::InputStruct, k::Complex128,
    ψ::Array{Complex{Float64},1}, m::Int, direction::String)::Complex128

    if (input.wgs.dir[input.sct.channels[m].wg] in ["x", "X"])
        kₓ, φy = wg_transverse_y(input, k, m)
        ny = analysis_quadrature_defaults()
        Y = Array(linspace(input.bnd.∂R[3], input.bnd.∂R[4], ny))
        dy = Y[2]-Y[1]
        if input.sct.channels[m].side in ["l", "L", "left", "Left"]
            x = input.sct.∂S[1]
            if direction == "in"
                phs = exp(-1im*kₓ*x)
            elseif direction == "out"
                phs = exp(+1im*kₓ*x)
            end
        elseif input.sct.channels[m].side in ["r", "R", "right", "Right"]
            x = input.sct.∂S[2]
            if direction == "in"
                phs = exp(+1im*kₓ*x)
            elseif direction == "out"
                phs = exp(-1im*kₓ*x)
            end
        end
        X = x*ones(Y)

        P = sqrt(real(kₓ))*analysis_interpolate(input, X, Y, conj(φy).*ψ)
        # φ = reshape(φy[input.dis.xy_inds],input.dis.N[1],:)[1,:]*sqrt(1/real(kₓ))
# println(size(φ))
# println(size(P))
        cm = sum(P)*dy
    elseif input.wgs.dir[input.sct.channels[m].wg] in ["y", "Y"]
        error("Haven't written vertical waveguide code yet.")
    end

    return cm
end

"""
P = analysis_interpolate(input, X, Y, ψ)
"""
function analysis_interpolate(input::InputStruct, X::Array{Float64,1}, Y::Array{Float64,1},
    ψ::Array{Complex128,1})::Array{Complex128,1}

    p = interpolate(reshape(ψ,input.dis.N_PML[1],:), BSpline(Linear()), OnGrid())
    # p = interpolate(reshape(ψ,input.dis.N_PML[1],:), BSpline(Quadratic(Reflect())), OnGrid())
    Ax = (input.dis.N_PML[1]-1)/(input.dis.xy_PML[1][end]-input.dis.xy_PML[1][1])
    Ay = (input.dis.N_PML[2]-1)/(input.dis.xy_PML[2][end]-input.dis.xy_PML[2][1])
    Bx = (input.dis.N_PML[1]*input.dis.xy_PML[1][1]-input.dis.xy_PML[1][end])/(input.dis.N_PML[1]-1)
    By = (input.dis.N_PML[2]*input.dis.xy_PML[2][1]-input.dis.xy_PML[2][end])/(input.dis.N_PML[2]-1)
    X_int = Ax*(X-Bx)
    Y_int = Ay*(Y-By)
    P = [p[X_int[ii],Y_int[ii]] for ii in 1:length(X)]

    return P
end

"""
surface_flux(input,ψ)
"""
function surface_flux(input::InputStruct, ψ::Array{Complex128,1}; nx::Int=20, ny::Int=20)::
    Tuple{Float64,Tuple{Float64,Float64,Float64,Float64}}

    flux, (left, right, bottom, top) = surface_flux(input, hcat(ψ,); nx=nx, ny=ny)

    return flux[1], (left[1],right[1],bottom[1],top[1])
end
function surface_flux(input::InputStruct,Ψ::Array{Complex128,2}; nx::Int=20, ny::Int=20)::
    Tuple{Array{Float64,1},Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1}}}

    flux = zeros(Complex128,size(Ψ,2))
    top = copy(flux)
    bottom = copy(flux)
    left = copy(flux)
    right = copy(flux)

    for i in 1:size(Ψ,2)
        ψ = reshape(Ψ[input.dis.xy_inds,i],input.dis.N[1],:)
        dψdx = (ψ[[nx+1,end-nx+1],:] - ψ[[nx-1, end-nx-1],:]) /2input.dis.dx[1]
        dψdy = (ψ[:,[ny+1,end-ny+1]] - ψ[:,[ny-1, end-ny-1]]) /2input.dis.dx[2]
        ψx = ψ[[nx,end-nx],:]
        ψy = ψ[:,[ny,end-ny]]
        kx = -real((conj(ψx).*dψdx - ψx.*conj(dψdx))/2im)
        ky = -real((conj(ψy).*dψdy - ψy.*conj(dψdy))/2im)
        top[i]    = (4*sum(ky[1:2:end-1,end])*input.dis.dx[1] + 2*sum(ky[2:2:end,end])*input.dis.dx[1])/3
        bottom[i] = (4*sum(-ky[1:2:end-1,1 ])*input.dis.dx[1] + 2*sum(-ky[2:2:end,1 ])*input.dis.dx[1])/3
        right[i]  = (4*sum(kx[end,1:2:end-1])*input.dis.dx[2] + 2*sum(kx[end,2:2:end])*input.dis.dx[2])/3
        left[i] =   (4*sum(-kx[1 ,1:2:end-1])*input.dis.dx[2] + 2*sum(-kx[1 ,2:2:end])*input.dis.dx[2])/3
        flux[i] = top[i] + bottom[i] + right[i] + left[i]
    end

    return flux,(left,right,bottom,top)
end # end of surface_flux


"""
compute_loss(ψ,k,inputs)
"""
function compute_loss(input::InputStruct, k::Union{Complex128,Float64}, ψ::Array{Complex128,1})::Float64

    loss = compute_loss(input,[complex(k)],hcat(ψ,))

    return loss[1]
end
function compute_loss(input::InputStruct,k::Array{Complex128,1},ψ::Array{Complex128,2})::Array{Float64,1}

    loss = zeros(Float64,length(k))

    for i in 1:length(k)

        Ψ = reshape(ψ[input.dis.xy_inds,i],input.dis.N[1],:)
        A1 = imag((input.sys.ε_sm-1im*input.tls.D₀.*input.sys.F_sm)*k[i]^2).*abs2.(Ψ)
        A1 = sum(4A1[1:2:end-1,:] + 2A1[2:2:end,:],1)/3
        loss[i] = (sum(4A1[1:2:end-1] + 2A1[2:2:end])/3)[1]*prod(input.dis.dx)
    end

    return loss
end
