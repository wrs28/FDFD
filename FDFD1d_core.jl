################################################################################
##### SYSTEM DEFINITION
################################################################################
(export SystemStruct, BoundaryStruct, DiscretizationStruct, ScatteringStruct,
            TwoLevelSystemStruct, ChannelStruct, WaveguideStruct, InputStruct)

"""
sys = SystemStruct(geometry, geoParams, n₁_vals, n₁_inds, n₂_vals, n₂_inds,
    F_vals, F_inds, n, ε, F, ε_sm, F_sm, regions)
sys = SystemStruct(geometry, geoParams, n₁_vals, n₁_inds, n₂_vals, n₂_inds, F_vals, F_inds)
"""
mutable struct SystemStruct
    geometry::Function
    geoParams::Array{Float64,1}
    n₁_vals::Array{Float64,1}
    n₁_inds::Array{Int,1}
    n₂_vals::Array{Float64,1}
    n₂_inds::Array{Int,1}
    F_vals::Array{Float64,1}
    F_inds::Array{Int,1}
    n::Array{Complex128,1}
    ε::Array{Complex128,1}
    F::Array{Float64,1}
    ε_sm::Array{Complex128,2}
    F_sm::Array{Float64,2}
    ε_PML::Array{Complex128,2}
    F_PML::Array{Float64,2}
    regions::Array{Int,2}

    function SystemStruct(geometry::Function, geoParams::Array{Float64,1},
        n₁_vals::Array{Float64,1}, n₁_inds::Array{Int,1}, n₂_vals::Array{Float64,1},
        n₂_inds::Array{Int,1}, F_vals::Array{Float64,1}, F_inds::Array{Int,1})::SystemStruct
        n = n₁_vals[n₁_inds]+1.0im*n₂_vals[n₂_inds]
        ε = n.^2
        F = F_vals[F_inds]
        ε_sm = Complex128[NaN NaN]
        F_sm = Float64[NaN NaN]
        ε_PML = Complex128[NaN NaN]
        F_PML = Float64[NaN NaN]
        regions = ones(Int,0,0)
        return new(geometry, geoParams, n₁_vals, n₁_inds, n₂_vals, n₂_inds,
                    F_vals, F_inds, n, ε, F, ε_sm, F_sm, ε_PML, F_PML, regions)
    end
    function SystemStruct(sys::SystemStruct)::SystemStruct
        return SystemStruct(sys.geometry, sys.geoParams, sys.n₁_vals, sys.n₁_inds,
                            sys.n₂_vals, sys.n₂_inds, sys.F_vals, sys.F_inds)
    end
end # end of struct SystemStruct

"""
bnd = BoundaryStruct(∂R, bc)
bnd = BoundaryStruct(∂R, bc, bk, input_modes)
"""
mutable struct BoundaryStruct
    ∂R::Array{Float64,1}
    ∂R_PML::Array{Float64,1}
    ℓ::Float64
    ℓ_PML::Float64
    bc::Array{String,1}
    bc_sig::String
    bk::Complex{Float64}
    incident_modes::Array{Int,1},1

    function BoundaryStruct(∂R::Array{Float64,1}, bc::Array{String,1})::BoundaryStruct
        ℓ = ∂R[2]-∂R[1]
        ∂R_PML = ∂R
        ℓ_PML = ∂R_PML[2]-∂R_PML[1]
        fix_bc!(bc)
        bc_sig = prod(bc)
        bk = complex(float(0))
        incident_modes = Int[]
        return new(∂R, ∂R_PML, ℓ, ℓ_PML, bc, bc_sig, bk, incident_modes)
    end
    function BoundaryStruct(∂R::Array{Float64,1}, bc::Array{String,1}, bk::Complex128,
        incident_modes::Array{Int,1})::BoundaryStruct
        x = BoundaryStruct(∂R, bc)
        x.bk = bk
        x.incident_modes = incident_modes
        return x
    end
    function BoundaryStruct(bnd::BoundaryStruct)::BoundaryStruct
        x = BoundaryStruct(bnd.∂R, bnd.bc)
        return x
    end
end # end of struct BoundaryStruct

"""
dis = DiscretizationStruct(N, sub_pixel_num)
"""
mutable struct DiscretizationStruct
    N::Int
    N_PML::Int
    dN_PML::Int

    dx::Float64

    x::Array{Float64,1}
    x_PML::Array{Float64,1}

    x_inds::Array{Int,1}

    sub_pixel_num::Int

    function DiscretizationStruct(N::Int, sub_pixel_num::Int)
        return new(N, N, 0, NaN, [NaN], [NaN], [0], sub_pixel_num)
    end
    function DiscretizationStruct(dis::DiscretizationStruct)
        return DiscretizationStruct(dis.N, dis.sub_pixel_num)
    end
end

"""
channel = ChannelStruct()
channel = ChannelStruct(tqn, wg, side)
"""
mutable struct ChannelStruct
    wg::Int         # waveguide number
    tqn::Int        # transverse quantum number
    side::String    # side

    function ChannelStruct()::ChannelStruct
        return new(0, 0, "")
    end
    function ChannelStruct(wg::Int, tqn::Int, side::String)::ChannelStruct
        return  new(wg, tqn, side)
    end
end # end of struct ChannelStruct

"""
sct = ScatteringStruct(∂S, channels, a)
"""
mutable struct ScatteringStruct
    ∂S::Array{Float64,1}
    channels::Array{ChannelStruct,1}

    function ScatteringStruct(∂S::Array{Float64,1}, channels::Array{ChannelStruct,1})::ScattteringStruct
        return new(∂S,channels)
    end
    function ScatteringStruct(sct::ScatteringStruct)::ScatteringStruct
        return ScatteringStruct(sct.∂S, sct.channels)
    end
end # end of struct ScatteringStruct

"""
tls = TwoLevelSystemStruct()
tls = TwoLevelSystemStruct(D₀, k₀, γ⟂)
"""
mutable struct TwoLevelSystemStruct
    D₀::Float64
    k₀::Float64
    γ⟂::Float64
    ω₀::Float64

    function TwoLevelSystemStruct()::TwoLevelSystemStruct
        D₀ = 0.
        k₀ = 10.
        γ⟂ = 1e8
        return new(D₀, k₀, γ⟂, k₀)
    end
    function TwoLevelSystemStruct(D₀::Float64,k₀::Float64,γ⟂::Float64)::TwoLevelSystemStruct
        return new(D₀, k₀, γ⟂, k₀)
    end
    function TwoLevelSystemStruct(tls::TwoLevelSystemStruct)::TwoLevelSystemStruct
        return TwoLevelSystemStruct(tls.D₀,tls.k₀,tls.γ⟂)
    end
end # end of struct TwoLevelSystemStruct

"""
input = InputStruct(system, boundary, discretization, scattering, tls, channels, waveguides)
"""
mutable struct InputStruct
    sys::SystemStruct
    bnd::BoundaryStruct
    dis::DiscretizationStruct
    sct::ScatteringStruct
    tls::TwoLevelSystemStruct

    function InputStruct(sys::SystemStruct, bnd::BoundaryStruct,
        dis::DiscretizationStruct, sct::ScatteringStruct,
        tls::TwoLevelSystemStruct)::InputStruct

        # define discretization parameters and PML parameters
        N = dis.N
        ∂R = bnd.∂R
        bc = bnd.bc
        ℓ = bnd.ℓ

        dx = ℓ[1]/N[1]
        x = ∂R[1] + dx[1]/2 + dx[1]*(0:N[1]-1)
        X = repeat(x[1];outer=(1,N[2]))

        extinction, change_per_site, power_law, α_imag = PML_params()
        dN_PML = Int[0, 0]
        for i in 1:2
            if bc[i] in ["O", "I"]
                dN_PML[i] = ceil(Int,(power_law+1)*log(extinction)/change_per_site)
            elseif bc[i] == "d"
                dN_PML[i] = 0
            elseif bc[i] == "n"
                dN_PML[i] = 0
            elseif bc[i] == "o"
                dN_PML[i] = 0
            elseif bc[i] == "p"
                dN_PML[i] = 0
            end
        end
        N_PML = N[1] + dN_PML[1] + dN_PML[2]
        ∂R_PML = [∂R[1] - dx[1]*(dN_PML[1]) - dx[1]/2,
                ∂R[2] + dx[1]*(dN_PML[2]) + dx[1]/2]
        ℓ_PML = ∂R_PML[2]-∂R_PML[1]
        x_PML = ∂R_PML[1] + dx[1]/2 + dx[1]*(0:N_PML[1]-1)
        X_PML = repeat(x_PML[1]; outer=(1,N_PML[2]))
        #xy_inds = reshape(kron(dN_PML[1]+(1:N[1]), ones(dN_PML[3]+(1:N[2]))), N[1],N[2])[:]
        # @time xy_inds = ( (A,B) -> sub2ind((N_PML...),A,B)).( repeat(dN_PML[1]+(1:N[1]); inner=(1,N[2])) , repeat(transpose(dN_PML[3]+(1:N[2])); outer=(N[1],1)))[:]
        x_inds = zeros(Int,N[1])
        for i in 1:N[1]
            xy_inds[i] = sub2ind((N_PML...), dN_PML[1]+(1:N[1])[i])
        end

        bnd.∂R_PML = ∂R_PML
        bnd.ℓ = ℓ
        bnd.ℓ_PML = ℓ_PML

        dis.N_PML = N_PML
        dis.dN_PML = dN_PML
        dis.x = x
        dis.X = X
        dis.dx = dx
        dis.x_PML = x_PML
        dis.X_PML = X_PML
        dis.x_inds = x_inds[:]

        # define system parameters
        n₁_vals = sys.n₁_vals
        n₁_inds = sys.n₁_inds
        n₂_vals = sys.n₂_vals
        n₂_inds = sys.n₂_inds
        F_vals = sys.F_vals
        F_inds = sys.F_inds

        n = vcat(n₁_vals[n₁_inds] + 1.0im*n₂_vals[n₂_inds])
        ε = n.^2
        F = vcat(F_vals[F_inds])
        F[iszero.(F)] = TwoLevelSystemDefaults()
        sys.n = n
        sys.ε = ε
        sys.F = F

        sys.regions = which_region(sys, dis)
        ɛ_sm, F_sm = sub_pixel_smoothing(sys, dis)

        sys.ε_sm = ε_sm
        sys.F_sm = F_sm

        temp = vcat(ones(dN_PML[1],1)*transpose(ε_sm[1,:]), ε_sm, ones(dN_PML[2],1)*transpose(ε_sm[end,:]))
        ε_PML = hcat(temp[:,1]*ones(1,dN_PML[3]), temp, temp[:,end]*ones(1,dN_PML[4]))
        sys.ε_PML = ε_PML
        temp = vcat(ones(dN_PML[1],1)*transpose(F_sm[1,:]), F_sm, ones(dN_PML[2],1)*transpose(F_sm[end,:]))
        F_PML = hcat(temp[:,1]*ones(1,dN_PML[3]), temp, temp[:,end]*ones(1,dN_PML[4]))
        sys.F_PML = F_PML

        temp = zeros(Complex128, N[1], N[2])
        sct.ε₀ = zeros(Complex128, N[1], N[2])
        for wg in 1:length(wgs.dir)
            if wgs.dir[wg]=="x"
                temp += repeat(transpose(wgs.ε[wg]);inner=(N[1],1))
            elseif wgs.dir[wg]=="y"
                temp += wgs.ε[wg]*ones(1,N[2])
            end
        end
        sct.ε₀ = temp .+ 1
        # temp = vcat(ones(dN_PML[1],1)*transpose(ε_sm[1,:]), sct.ε₀, ones(dN_PML[2],1)*transpose(ε_sm[end,:]))
        # sct.ε₀_PML = hcat(temp[:,1]*ones(1,dN_PML[3]), temp, temp[:,end]*ones(1,dN_PML[4]))
        ε_temp = vcat(repmat(transpose(sct.ε₀[1,:]),dN_PML[1],1), sct.ε₀, repmat(transpose(sct.ε₀[1,:]),dN_PML[2],1))
        sct.ε₀_PML = hcat(repmat(ε_temp[:,1],1,dN_PML[3]), ε_temp, repmat(ε_temp[:,1],1,dN_PML[4]))

        return new(sys, bnd, dis, sct, tls)
    end
    function InputStruct(input::InputStruct)::InputStruct
        input = InputStruct(input.sys,input.bnd,input.dis,input.sct,input.tls)
        return input
    end
end # end of struct InputStruct

# ################################################################################
# ##### SYSTEM DEFINITION AUXILLIARIES
# ################################################################################
"""
r = which_region(XY, system, waveguides, discretization)
"""
function which_region(x::Array{Float64,1}, sys::SystemStruct,
    dis::DiscretizationStruct)::Array{Int,2}

    geometry = sys.geometry
    geoParams = sys.geoParams

    regions = zeros(Int,length(x))

    for i in 1:length(x)
        regions[i] = geometry(x[i], geoParams)
    end
    return regions
end
"""
r = which_region(system, waveguides, discretization)
"""
function which_region(sys::SystemStruct, dis::DiscretizationStruct)::Array{Int,2}
    regions = which_region(dis.x, sys, dis)
    return regions
end


"""
ε_sm, F_sm = sub_pixel_smoothing(system)
"""
function sub_pixel_smoothing(sys::SystemStruct, dis::DiscretizationStruct)::Tuple{Array{Complex128,1},Array{Float64,1}}

    x = dis.x
    sub_pixel_num = dis.sub_pixel_num

    r = sys.regions
    ε = sys.ε
    F = sys.F

    ɛ_sm = ɛ[r]
    F_sm = F[r]

    X = ones(Float64,sub_pixel_num)
    for i in 2:(length(x)-1)
        nearestNeighborFlag = (r[i]!==r[i+1]) | (r[i]!==r[i-1])
        if nearestNeighborFlag
            x_min = (x[i]+x[i-1])/2

            x_max = (x[i]+x[i+1])/2

            X = Array(linspace(x_min,x_max,sub_pixel_num))

            sub_regions = which_region(X, sys, dis)
            ɛ_sm[i] = mean(ɛ[sub_regions])
            F_sm[i] = mean(F[sub_regions])
        end
    end

    ε_sm[1] = ε_sm[2]
    ε_sm[end] = ε_sm[end-1]

    return ɛ_sm, F_sm
end

# ################################################################################
# ##### SYSTEM MODIFICATION
# ################################################################################
export update_input

"""
update_input(input, field, value)
update_input(input, fields, values)

    If changes were made to the ∂R, N, k₀, k, F, ɛ, Γ, bc, a, b, run update_input to
    propagate these changes through the system.
"""
function update_input(input::InputStruct, fields::Array{Symbol,1}, values::Any)::InputStruct

    f = fieldnames(input)
    for i in 1:length(f)
        g = fieldnames(getfield(input,f[i]))
        for m in 1:length(fields)
            if fields[m] in g
                setfield!(getfield(input,f[i]),fields[m],values[m])
            end
        end
    end

    if fields == :ω₀
        input.tls.k₀ = input.tls.ω₀
    elseif fields == :k₀
        input.tls.ω₀ = input.tls.k₀
    end

    flag = zeros(Bool,length(fields))
    for i in 1:length(fields)
        flag[i] = fields[i] in [:∂R, :N, :bc, :F_inds, :F_vals, :n₁_vals, :n₂_vals, :n₁_inds,
                        :n₂_inds, :sub_pixel_num, :geoParams]
    end

    if any(flag)
        sys = SystemStruct(input.sys)
        bnd = BoundaryStruct(input.bnd)
        dis = DiscretizationStruct(input.dis)
        sct = ScatteringStruct(input.sct)
        tls = TwoLevelSystemStruct(input.tls)
        input = InputStruct(sys, bnd, dis, sct, tls)
    end

    return input
end # end of function update_input
function update_input(input::InputStruct, field::Symbol, value::Any)::InputStruct
    input = update_input(input,[field],[value])
    return input
end # end of function update_input
function update_input(input::InputStruct, field::Symbol, ind::Int, value::Any)::InputStruct

    f = fieldnames(input)
    for i in 1:length(f)
        g = fieldnames(getfield(input,f[i]))
        if field in g
            temp = getfield(getfield(input,f[i]),field)
            temp[ind] = value
            input = update_input(input,field,temp)
        end
    end

    return input
end # end of function update_input


# ################################################################################
# ##### SYSTEM STANDARDIZATION
# ################################################################################
"""
fix_bc!(bc)
"""
function fix_bc!(bc::Array{String,1})::Array{String,1}
    for i in 1:4
        if bc[i] in ["pml_out", "PML_out"]
            bc[i] = "O"
        elseif bc[i] in ["pml_in", "PML_in"]
            bc[i] = "I"
        elseif bc[i] in ["d", "dirichlet", "Dirichlet", "hard", "h"]
            bc[i] = "d"
        elseif bc[i] in ["n", "neumann",   "Neumann",   "soft", "s"]
            bc[i] = "n"
        elseif bc[i] in ["o", "open", "Open"]
            bc[i] = "o"
        elseif bc[i] in ["p", "periodic", "Periodic"]
            bc[i] = "p"
        end
    end
    return bc
end






















# module Core
#
# export InputStruct, laplacian, whichRegion, trapz, dirac_δ, heaviside_Θ, subpixelSmoothing, processInputs, updateInputs!
#
# #############################################################
#
#
# """
# inputs = Inputs(ω, ω₀, k, k₀, N, ℓ, dx, x_ext, x_inds, ∂, ɛ, F, N_ext, ℓ_ext, ∂_ext, ɛ_ext, F_ext, x, xᵨ₋, xᵨ₊, γ⟂, D₀, a, b, Γ, dN, bc, subPixelNum, r_ext, ɛ_sm, F_sm)
#
# """
# mutable struct InputStruct
#     ω::Array{Complex128,1}
#     ω₀::Complex128
#     k::Array{Complex128,1}
#     k₀::Complex128
#     N::Int
#     ℓ::Float64
#     dx::Float64
#     x_ext::Array{Float64,1}
#     x_inds::Array{Int,1}
#     ∂::Array{Float64,1}
#     ɛ::Array{Complex128,1}
#     F::Array{Float64,1}
#     N_ext::Int
#     ℓ_ext::Float64
#     ∂_ext::Array{Float64,1}
#     ɛ_ext::Array{Complex128,1}
#     F_ext::Array{Float64,1}
#     x::Array{Float64,1}
#     xᵨ₋::Float64
#     xᵨ₊::Float64
#     γ⟂::Float64
#     D₀::Float64
#     a::Array{Complex128,1}
#     b::Array{Complex128,1}
#     Γ::Array{Complex128,1}
#     Γ_ext::Array{Complex128,1}
#     dN::Array{Int,1}
#     bc::Array{String,1}
#     subPixelNum::Int
#     r_ext::Array{Int,1}
#     ɛ_sm::Array{Complex128,1}
#     F_sm::Array{Float64,1}
# end
#
#
#
# """
# ∇ =  grad(N, dx)
#
# Gradient with N points, lattice spacing dx. It's the forward gradient (I think).
#
# sparse ∇[N,N+1]
#
# """
# function grad(N::Int, dx::Float64)
#
#     I₁ = Array(1:N)
#     J₁ = Array(1:N)
#     V₁ = fill(Complex(-1/dx), N)
#
#     I₂ = Array(1:N)
#     J₂ = Array(2:(N+1))
#     V₂ = fill(Complex(+1/dx), N)
#
#     ∇ = sparse( vcat(I₁,I₂), vcat(J₁,J₂), vcat(V₁,V₂), N, N+1, +)
#
# end # end of function grad
#
#
#
#
# """
# ∇² =  laplacian(k, inputs)
#
# Laplacian evaluated at k∈ℂ. Lattice size & spacing, bc's determined in inputs.
#
# sparse ∇²[# lattice pts, # lattice pts]
#
# """
# function laplacian(k::Complex128, inputs::InputStruct)
#
#     # definitions block#
#     N = inputs.N_ext
#     bc = inputs.bc
#     a = inputs.a
#     b = inputs.b
#     dx = inputs.dx
#     dx² = dx^2
#     ##
#
#     # generate robin parameter
#     λ = zeros(Complex128,2)
#     for i in 1:2
#         λ[i] = -(-1)^i*1im*k*(b[i]-a[i])/(b[i]+a[i])
#     end
#
#     ∇ = grad(N-1,dx)
#
#     # create PML layer
#     Σ = 1.+1im*σ(inputs)./real(k)
#     s₁ = sparse( 1:N-1, 1:N-1, 2./( Σ[1:end-1] .+ Σ[2:end] ), N-1, N-1)
#     s₂ = sparse( 1:N, 1:N, 1./Σ, N, N)
#     ∇² = -s₂*transpose(∇)*s₁*∇
#
#     # truncate
#     ind = [1, N, 1]
#     for i in 1:2
#         if bc[i] in ["pml_out" "pml_in"]
#             ∇²[ind[i],ind[i]]   += -2/dx²
#         elseif bc[i] in ["d" "dirichlet" "Dirichlet" "hard"]
#             ∇²[ind[i],ind[i]]   += -2/dx²
#         elseif bc[i] in ["n" "neumann" "Neumann" "soft"]
#             ∇²[ind[i],ind[i]]   += 0
#         elseif bc[i] in ["p" "periodic"]
#             ∇²[ind[i],ind[i]]   += -1/dx²
#             ∇²[ind[i],ind[i+1]] += +1/dx²
#         elseif bc[i] in ["o" "out" "outgoing"]
#             ∇²[ind[i],ind[i]]   += +(1im*k/dx)/(1-1im*dx*k/2)
#         elseif bc[i] in ["i" "in" "incoming" "incident"]
#             ∇²[ind[i],ind[i]]   += -(1im*k/dx)/(1+1im*dx*k/2)
#         elseif bc[i] in ["r" "robin"]
#             ∇²[ind[i],ind[i]]   += -(-1)^i*(1/dx)*(λ[i]/(1+(-1)^i*λ[i]*dx/2))
#         else
#             println("error in bc specification, not valid")
#             return
#         end
#     end
#
#     return ∇²
#
# end # end of function laplacian
#
#
#
# """
# Σ =  σ(inputs)
#
# Conductivity for PML layer.
#
# """
# function σ(inputs::InputStruct)
#
#     ####################
#     PML_extinction = 1e6
#     PML_ρ = 1/6
#     PML_power_law = 3
#     α_imag = -.15
#     ####################
#
#     x = inputs.x_ext
#     dx = inputs.dx
#     N = inputs.N_ext
#     ∂ = inputs.∂_ext
#
#     α = zeros(Complex128,2)
#     for i in 1:2
#         if inputs.bc[i] == "pml_out"
#             α[i] = +(1+α_imag*1im)*( (PML_ρ/dx)^(PML_power_law+1) )/ ( (PML_power_law+1)*log(PML_extinction) )^PML_power_law
#         elseif inputs.bc[i] == "pml_in"
#             α[i] = -(1+α_imag*1im)*( (PML_ρ/dx)^(PML_power_law+1) )/ ( (PML_power_law+1)*log(PML_extinction) )^PML_power_law
#         else
#             α[i] = 0
#         end
#     end
#
#     s = zeros(Complex128,N)
#
#     for i in 1:N
#         if x[i] ≤ ∂[2]
#             s[i] = α[1]*abs(x[i]-∂[2])^PML_power_law
#         elseif x[i] ≥ ∂[end-1]
#             s[i] = α[2]*abs(x[i]-∂[end-1])^PML_power_law
#         else
#             s[i] = 0
#         end
#     end
#
#     return s
#
# end # end of function σ
#
#
#
# """
# region = whichRegion(x, ∂)
#
# Takes vector x, gives vector of regions as defined in input file.
#
# """
# function whichRegion(x::Array{Float64,1}, ∂::Array{Float64,1})
#
#     region = similar(x,Int);
#
#     for i in 1:length(region), k in 1:length(∂)-1
#
#         if k+1 == length(∂)
#             if  ∂[k] ≤ x[i]
#                 region[i] = k
#             end
#         elseif k == 1
#             if  x[i] ≤ ∂[k+1]
#                 region[i] = k
#             end
#         else
#             if  ∂[k] ≤ x[i] ≤ ∂[k+1]
#                 region[i] = k
#             end
#         end
#
#     end
#
#     return region
#
# end # end of function whichRegion
#
#
#
# """
# ∫z_dx = trapz(z, dx)
#
# Trapezoidal sum of z.
#
# """
# function trapz(z::Array{Complex128,1},dx::Float64)::Complex128
#
#     ∫z_dx = dx*sum(z)-dx*(z[1]+z[end])/2
#
#     return ∫z_dx
#
# end # end of function trapz
#
#
#
# """
# δ, ∇δ = dirac_δ(x, x₀)
#
# δ sparse, dirac distribution weighted to be centered at x₀.
# ∇δ sparse, derivative of δ.
#
# """
# function dirac_δ(x::Array{Float64,1},x₀::Float64)
#
#     ind₁ = findmin(abs2.(x-x₀))[2]
#     ind₂ = ind₁ + (2*mod(findmin([findmin(abs2.(x[ind₁+1]-x₀))[1] findmin(abs2.(x[ind₁-1]-x₀))[1]])[2],2)[1] -1)
#     min_ind = min(ind₁,ind₂)
#     max_ind = max(ind₁,ind₂)
#
#     x₁ = x[min_ind]
#     x₂ = x[max_ind]
#     dx = abs(x₂-x₁)
#     dx² = dx^2
#
#     w₁ = abs(x₂-x₀)./dx²
#     w₂ = abs(x₀-x₁)./dx²
#
#     δ = sparsevec([min_ind,max_ind],[w₁,w₂],length(x),+)
#
#     δ1 = sparsevec([min_ind,max_ind]  ,[w₁,w₂],length(x)+1,+)
#     δ2 = sparsevec([min_ind,max_ind]+1,[w₁,w₂],length(x)+1,+)
#     ∇ = grad(length(x),dx)
#
#     ∇δ = ∇*(δ1.+δ2)/2
#
#     return δ,∇δ
#
# end # end of function dirac_δ
#
#
#
# """
# Θ,∇Θ,∇²Θ = heaviside_Θ(x, x₀)
#
# Θ,∇Θ,∇²Θ not sparse, weighted to be centered at x₀.
# """
# function heaviside_Θ(x::Array{Float64,1},x₀::Float64)
#
#     ind₁ = findmin(abs2.(x-x₀))[2]
#     ind₂ = ind₁ + (2*mod(findmin([findmin(abs2.(x[ind₁+1]-x₀))[1] findmin(abs2.(x[ind₁-1]-x₀))[1]])[2],2)[1] -1)
#     min_ind = min(ind₁,ind₂)
#     max_ind = max(ind₁,ind₂)
#
#     x₁ = x[min_ind]
#     x₂ = x[max_ind]
#     dx = x₂-x₁
#
#     Θ = zeros(length(x),1)
#     Θ[x .≥ x₀,1] = 1.
#     w₁ = (x₂-x₀)./dx
#     w₂ = (x₀-x₁)./dx
#     Θ[min_ind,1] = w₁
#     Θ[max_ind,1] = w₂
#
#     ∇Θ  = zeros(length(x),1)
#     ∇Θ[1,1]    = (Θ[2,1]-Θ[1,1])/dx
#     ∇Θ[end,1]  = (Θ[end,1]-Θ[end-1,1])/dx
#     ∇Θ[2:end-1,1]  = (Θ[3:end,1]-Θ[1:end-2,1])/2dx
#
#     ∇²Θ = zeros(length(x),1)
#     ∇²Θ[2:end-1,1] = (Θ[1:end-2,1] - 2Θ[2:end-1,1] + Θ[3:end,1])/dx^2
#
#     return (Θ,∇Θ,∇²Θ)
#
# end # end of function heaviside_Θ
#
#
#
# """
# ɛ_smoothed, F_smoothed = subpixelSmoothing(inputs; truncate = false, r = [])
#
# Sub-pixel smoothed ɛ & F.
#
# if truncate=true, then output is truncated to bulk region (sans PML)
#
# r is the output of whichRegion. If computing r elsewhere, can save time by using that output here.
#
# """
# function subpixelSmoothing(inputs::InputStruct; truncate::Bool = false, r::Array{Int,1} = Int[])
#     # for now it's for ɛ and F, could feasibly be extended eventually...
#
#     X = inputs.x
#     X_ext = inputs.x_ext
#     ∂ = inputs.∂_ext
#     ɛ = inputs.ɛ_ext
#     F = inputs.F_ext
#     subPixelNum = inputs.subPixelNum
#
#     ɛ_smoothed, F_smoothed = subpixelSmoothing_core(X, X_ext, ∂, ɛ, F, subPixelNum, truncate, r)
#
#     return ɛ_smoothed, F_smoothed
#
# end # end of function subpixelSmoothing
#
#
#
#
# """
# ɛ_smoothed, F_smoothed = subpixelSmoothing_core(inputs; truncate = false, r = [])
#
# Sub-pixel smoothed ɛ & F.
#
# heart of subpixelSmoothing routine, makes it usable in processInputs
# """
# function subpixelSmoothing_core(X::Array{Float64,1}, X_ext::Array{Float64,1}, ∂::Array{Float64,1}, ɛ::Array{Complex128,1}, F::Array{Float64,1}, subPixelNum::Int, truncate::Bool, r::Array{Int,1})
#     # for now it's for ɛ and F, could feasibly be extended eventually...
#
#     if truncate
#         x = X
#     else
#         x = X_ext
#     end
#
#     if isempty(r)
#         r = whichRegion(x, ∂)
#     end
#
#     ɛ_smoothed = deepcopy(ɛ[r])
#     F_smoothed = deepcopy(F[r])
#
#     for i in 2:(length(x)-1)
#
#         nearestNeighborFlag = (r[i]!==r[i+1]) | (r[i]!==r[i-1])
#
#         if nearestNeighborFlag
#
#             x_min = (x[i]+x[i-1])/2
#             x_max = (x[i]+x[i+1])/2
#
#             X = Array(linspace(x_min, x_max, subPixelNum))
#
#             subRegion = whichRegion(X, ∂)
#             ɛ_smoothed[i] = mean(ɛ[subRegion])
#             F_smoothed[i] = mean(F[subRegion])
#
#         end
#
#     end
#
#     return ɛ_smoothed, F_smoothed
#
# end # end of function subpixelSmoothing_core
#
#
#
#
#
# """
# processInputs(fileName = "./SALT_1d_Inputs.jl")
#
#
# """
# function processInputs(fileName::String = "./SALT_1d_Inputs.jl")::InputStruct
#
#     #################### See also definition block in sigma function
#     F_min = 1e-16
#     PML_extinction = 1e6
#     PML_ρ = 1/6
#     PML_power_law = 3
#     ####################
#
#     (N::Int, k₀::Complex128, k::Array{Complex128,1}, bc::Array{String,1}, ∂::Array{Float64,1}, Γ::Array{Complex128,1}, F::Array{Float64,1}, ɛ::Array{Complex128,1}, γ⟂::Float64, D₀::Float64, a::Array{Complex128,1}, b::Array{Complex128,1}, subPixelNum::Int) = evalfile(fileName)
#
#     ω₀ = k₀
#     ω  = k
#     ℓ = ∂[end] - ∂[1]
#
#     xᵨ₊ = (∂[1]+∂[2])/2
#     xᵨ₋ = (∂[end-1]+∂[end])/2
#
#     ##########################
#
#     dx = ℓ/(N-1);
#
#     dN = Int[0, 0]
#     for i in 1:2
#         if bc[i] in ["o", "out", "outgoing"]
#             dN[i] = 0
#         elseif bc[i] in ["i", "in", "incoming", "incident"]
#             dN[i] = 0
#         elseif bc[i] in ["d", "dirichlet", "Dirichlet", "hard"]
#             dN[i] = 0
#         elseif bc[i] in ["n", "neumann",   "Neumann",   "soft"]
#             dN[i] = 0
#         elseif bc[i] in ["p", "periodic"]
#             dN[i] = 0
#         elseif bc[i] in ["pml_out", "pml_in"]
#             dN[i] = ceil(Int,(PML_power_law+1)*log(PML_extinction)/PML_ρ)
#         elseif bc[i] in ["r", "robin"]
#             dN[i] = 0
#         else
#             println("error in bc specification, not valid")
#             return
#         end
#     end
#
#     N_ext = N + sum(dN)
#     ℓ_ext = dx*(N_ext-1)
#
#     x_ext = vcat(-[dN[1]:-1:1;]*dx+∂[1],  linspace(∂[1],∂[end],N), [1:dN[2];]*dx+∂[end])
#     x_inds = dN[1]+collect(1:N)
#     x = x_ext[x_inds]
#
#     ∂_ext = vcat([x_ext[1]-dx/2], ∂, [x_ext[end]+dx/2])
#
#     F[iszero.(F)] = F_min
#     F_ext = vcat([F_min], F, [F_min])
#
#     ɛ_ext = vcat([1], ɛ, [1])
#     Γ_ext = vcat([Inf], Γ, [Inf])
#
#     r_ext = whichRegion(x_ext, ∂_ext)
#     ɛ_sm, F_sm = subpixelSmoothing_core(x, x_ext, ∂_ext, ɛ_ext, F_ext, subPixelNum, false, r_ext)
#
#     inputs = InputStruct(ω, ω₀, k, k₀, N, ℓ, dx, x_ext, x_inds, ∂, ɛ, F, N_ext, ℓ_ext, ∂_ext, ɛ_ext, F_ext, x, xᵨ₋, xᵨ₊, γ⟂, D₀, a, b, Γ, Γ_ext, dN, bc, subPixelNum, r_ext, ɛ_sm, F_sm)
#
#
# end  # end of function processInputs
#
#
#
#
# """
# inputs = updateInputs(inputs)
#
# If changes were made to the ∂, N, k₀, k, F, ɛ, Γ, bc, a, b, run updateInputs to propagate these changes through the system.
#
# """
# function updateInputs!(inputs::InputStruct, field::Symbol, value::Any)
#
#     #################### See also definition block in sigma function
#     F_min = 1e-16
#     PML_extinction = 1e6
#     PML_ρ = 1/6
#     PML_power_law = 3
#     α_imag = -.15
#     ####################
#
#     fields = fieldnames(inputs)
#     setfield!(inputs,field,value)
#
#     if field == :ω₀
#         inputs.k₀ = inputs.ω₀
#     elseif field == :k₀
#         inputs.ω₀ = inputs.k₀
#     elseif field == :ω
#         inputs.k = inputs.ω
#     elseif field == :k
#         inputs.ω = inputs.k
#     end
#
#     if field in [:∂, :N, :bc]
#         ∂ = inputs.∂
#         N = inputs.N
#         bc = inputs.bc
#
#         xᵨ₊ = (∂[1]+∂[2])/2
#         xᵨ₋ = (∂[end-1]+∂[end])/2
#
#         ℓ = ∂[end] - ∂[1]
#
#         ##########################
#
#         dx = ℓ/(N-1);
#
#         dN = Int[0, 0]
#         for i in 1:2
#             if bc[i] in ["o" "out" "outgoing"]
#                 dN[i] = 0
#             elseif bc[i] in ["i" "in" "incoming" "incident"]
#                 dN[i] = 0
#             elseif bc[i] in ["d" "dirichlet" "Dirichlet" "hard"]
#                 dN[i] = 0
#             elseif bc[i] in ["n" "neumann"   "Neumann"   "soft"]
#                 dN[i] = 0
#             elseif bc[i] in ["p" "periodic"]
#                 dN[i] = 0
#             elseif bc[i] in ["pml_out" "pml_in"]
#                 dN[i] = ceil(Int,(PML_power_law+1)*log(PML_extinction)/PML_ρ)
#             elseif bc[i] in ["r" "robin"]
#                 dN[i] = 0
#             elseif bc[i] in ["CW" "cw" "clockwise" "Clockwise" "clock" "Clock"]
#                 dN[i] = 0
#             elseif bc[i] in ["CCW" "ccw" "counterclockwise" "CounterClockwise" "counter" "Counter"]
#                 dN[i] = 0
#             else
#                 println("error in bc specification, not valid")
#                 return
#             end
#         end
#
#         N_ext = N + sum(dN)
#         ℓ_ext = dx*(N_ext-1)
#
#         x_ext = vcat(-[dN[1]:-1:1;]*dx+∂[1],  linspace(∂[1],∂[end],N), [1:dN[2];]*dx+∂[end])
#         x_inds = dN[1]+collect(1:N)
#         x = x_ext[x_inds]
#
#         ∂_ext = vcat([x_ext[1]-dx/2], ∂, [x_ext[end]+dx/2])
#
#         inputs.xᵨ₊ = xᵨ₊
#         inputs.xᵨ₋ = xᵨ₋
#         inputs.ℓ = ℓ
#         inputs.dx = dx
#         inputs.dN = dN
#         inputs.N_ext = N_ext
#         inputs.ℓ_ext = ℓ_ext
#         inputs.x_ext = x_ext
#         inputs.x_inds = x_inds
#         inputs.x = x
#         inputs.∂_ext = ∂_ext
#
#     end
#
#     if  field in [:∂, :N, :bc, :F, :ɛ, :Γ, :subPixelNum]
#
#         F = inputs.F
#         ɛ = inputs.ɛ
#         Γ = inputs.Γ
#         subPixelNum = inputs.subPixelNum
#
#         F[iszero.(F)] = F_min
#         F_ext = vcat([F_min], F, [F_min])
#
#         ɛ_ext = vcat([1], ɛ, [1])
#         Γ_ext = vcat([Inf], Γ, [Inf])
#
#         r_ext = whichRegion(inputs.x_ext, inputs.∂_ext)
#         ɛ_sm, F_sm = subpixelSmoothing_core(inputs.x, inputs.x_ext, inputs.∂_ext, ɛ_ext, F_ext, inputs.subPixelNum, false, r_ext)
#
#         inputs.F_ext = F_ext
#         inputs.Γ_ext = Γ_ext
#         inputs.r_ext = r_ext
#         inputs.ɛ_sm = ɛ_sm
#         inputs.F_sm = F_sm
#
#     end
#
#     return inputs
#
# end # end of function updateInputs
#
# #######################################################################
#
# end # end of Module Core
