################################################################################
##### SYSTEM DEFINITION
################################################################################
(export SystemStruct, BoundaryStruct, DiscretizationStruct, ScatteringStruct,
            TwoLevelSystemStruct, ChannelStruct, InputStruct)

"""
sys = SystemStruct(geometry, geoParams, n₁_vals, n₁_inds, n₂_vals, n₂_inds,
    F_vals, F_inds, n, ε, F, ε_sm, F_sm, regions)
sys = SystemStruct(geometry, geoParams, n₁_vals, n₁_inds, n₂_vals, n₂_inds, F_vals, F_inds)
"""
mutable struct SystemStruct
    ∂C::Function
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
    ε_sm::Array{Complex128,1}
    F_sm::Array{Float64,1}
    ε_PML::Array{Complex128,1}
    F_PML::Array{Float64,1}
    regions::Array{Int,1}

    function SystemStruct(∂C::Function, geometry::Function, geoParams::Array{Float64,1},
        n₁_vals::Array{Float64,1}, n₁_inds::Array{Int,1}, n₂_vals::Array{Float64,1},
        n₂_inds::Array{Int,1}, F_vals::Array{Float64,1}, F_inds::Array{Int,1})::SystemStruct
        n = n₁_vals[n₁_inds]+1.0im*n₂_vals[n₂_inds]
        ε = n.^2
        F = F_vals[F_inds]
        ε_sm = Complex128[NaN]
        F_sm = Float64[NaN]
        ε_PML = Complex128[NaN]
        F_PML = Float64[NaN]
        regions = ones(Int,0)
        return new(∂C, geometry, geoParams, n₁_vals, n₁_inds, n₂_vals, n₂_inds,
                    F_vals, F_inds, n, ε, F, ε_sm, F_sm, ε_PML, F_PML, regions)
    end
    function SystemStruct(sys::SystemStruct)::SystemStruct
        return SystemStruct(sys.∂C, sys.geometry, sys.geoParams, sys.n₁_vals, sys.n₁_inds,
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
    incident_modes::Array{Int,1}

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
    dN_PML::Array{Int,1}

    dx::Float64

    x::Array{Float64,1}
    x_PML::Array{Float64,1}

    x_inds::Array{Int,1}

    sub_pixel_num::Int

    function DiscretizationStruct(N::Int, sub_pixel_num::Int)
        return new(N, N, [0, 0], NaN, [NaN], [NaN], [0], sub_pixel_num)
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
    ε₀::Array{Complex128,1}
    ε₀_PML::Array{Complex128,1}

    function ScatteringStruct(∂S::Array{Float64,1}, channels::Array{ChannelStruct,1})::ScattteringStruct
        return new(∂S,channels, Complex128[NaN],  Complex128[NaN])
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
        X = x

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
        X_PML = x_PML
        #xy_inds = reshape(kron(dN_PML[1]+(1:N[1]), ones(dN_PML[3]+(1:N[2]))), N[1],N[2])[:]
        # @time xy_inds = ( (A,B) -> sub2ind((N_PML...),A,B)).( repeat(dN_PML[1]+(1:N[1]); inner=(1,N[2])) , repeat(transpose(dN_PML[3]+(1:N[2])); outer=(N[1],1)))[:]
        x_inds = zeros(Int,N[1])
        for i in 1:N[1]
            x_inds[i] = sub2ind((N_PML...), dN_PML[1]+(1:N[1])[i])
        end

        bnd.∂R_PML = ∂R_PML
        bnd.ℓ = ℓ
        bnd.ℓ_PML = ℓ_PML

        dis.N_PML = N_PML
        dis.dN_PML = dN_PML
        dis.x = x
        dis.dx = dx
        dis.x_PML = x_PML
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

        ε_PML = vcat(ones(dN_PML[1])*ε_sm[1], ε_sm, ones(dN_PML[2])*ε_sm[end])

        sys.ε_PML = ε_PML
        F_PML = vcat(ones(dN_PML[1])*F_sm[1], F_sm, ones(dN_PML[2])*F_sm[end])
        sys.F_PML = F_PML

        temp = zeros(Complex128, N)
        sct.ε₀ = ones(Complex128, N)
        sct.ε₀_PML = ones(Complex128, N_PML)

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
r = which_region(X, system, waveguides, discretization)
"""
function which_region(x::Array{Float64,1}, sys::SystemStruct,
    dis::DiscretizationStruct)::Array{Int,1}

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
function which_region(sys::SystemStruct, dis::DiscretizationStruct)::Array{Int,1}
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
    for i in 1:2
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
