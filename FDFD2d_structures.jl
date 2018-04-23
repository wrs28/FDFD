(export SystemStruct, BoundaryStruct, DiscretizationStruct, ScatteringStruct,
            TwoLevelSystemStruct, ChannelStruct, WaveguideStruct, InputStruct)

################################################################################
##### SYSTEM DEFINITION
################################################################################
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
    ℓ::Array{Float64,1}
    ℓ_PML::Array{Float64,1}
    bc::Array{String,1}
    bc_sig::String
    bk::Array{Complex{Float64},1}
    incident_modes::Array{Array{Int,1},1}

    function BoundaryStruct(∂R::Array{Float64,1}, bc::Array{String,1})::BoundaryStruct
        ℓ = [∂R[2]-∂R[1], ∂R[4]-∂R[3]]
        ∂R_PML = ∂R
        ℓ_PML = [∂R_PML[2]-∂R_PML[1], ∂R_PML[4]-∂R_PML[3]]
        fix_bc!(bc)
        bc_sig = prod(bc)
        bk = [complex(float(0)),complex(float(0))]
        incident_modes = [Int[],Int[]]
        return new(∂R, ∂R_PML, ℓ, ℓ_PML, bc, bc_sig, bk, incident_modes)
    end
    function BoundaryStruct(∂R::Array{Float64,1}, bc::Array{String,1}, bk::Array{Complex128,1},
        incident_modes::Array{Array{Int,1},1})::BoundaryStruct
        x = BoundaryStruct(∂R, bc)
        x.bk = bk
        x.incident_modes = incident_modes
        return x
    end
    function BoundaryStruct(bnd::BoundaryStruct)::BoundaryStruct
        x = BoundaryStruct(bnd.∂R)
        return x
    end
end # end of struct BoundaryStruct

"""
dis = DiscretizationStruct(N, sub_pixel_num)
"""
mutable struct DiscretizationStruct
    N::Array{Int,1}
    N_PML::Array{Int,1}
    dN_PML::Array{Int,1}

    dx::Tuple{Float64,Float64}

    xy::Tuple{Array{Float64,1},Array{Float64,1}}
    XY::Tuple{Array{Float64,2},Array{Float64,2}}

    xy_PML::Tuple{Array{Float64,1},Array{Float64,1}}
    XY_PML::Tuple{Array{Float64,2},Array{Float64,2}}

    xy_inds::Array{Int,1}

    sub_pixel_num::Int

    function DiscretizationStruct(N::Array{Int,1}, sub_pixel_num::Int)
        return new(N, N, [0,0], (NaN,NaN), ([NaN],[NaN]), ([NaN NaN],[NaN NaN]), ([NaN],[NaN]),([NaN NaN],[NaN NaN]), [0], sub_pixel_num)
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
    ε₀::Array{Complex128,2}
    ε₀_PML::Array{Complex128,2}

    function ScatteringStruct(∂S::Array{Float64,1}, channels::Array{ChannelStruct,1})::ScattteringStruct
        return new(∂S,channels,Complex128[NaN NaN],Complex128[NaN NaN])
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
wgs = WaveguideStruct()
wgs = WaveguideStruct(wgd,wgp,wgt,wgn)
"""
mutable struct WaveguideStruct
    dir::String
    pos::Float64
    wdt::Float64
    ind::Float64
    ε::Array{Complex128,1}
    ε_PML::Array{Complex128,1}

    function WaveguideStruct()::WaveguideStruct
        return new("", 0., 0., 1.,complex([NaN]),complex([NaN]))
    end
    function WaveguideStruct(dir::String, pos::Float64, wdt::Float64, ind::Float64)::WaveguideStruct
        return new(dir, pos, wdt, ind, complex([NaN]), complex([NaN]))
    end
    function WaveguideStruct(wgs::WaveguideStruct)::WaveguideStruct
        return WaveguideStruct(wgs.dir, wgs.pos, wgs.wdt, wgs.ind)
    end
end # end of struct WaveguideStruct

"""
input = InputStruct(system, boundary, discretization, scattering, tls, channels, waveguides)
"""
mutable struct InputStruct
    sys::SystemStruct
    bnd::BoundaryStruct
    dis::DiscretizationStruct
    sct::ScatteringStruct
    tls::TwoLevelSystemStruct
    wgs::Array{WaveguideStruct,1}

    function InputStruct(sys, bnd, dis, sct, tls, wgs)

        # define discretization parameters and PML parameters
        N = dis.N
        ∂R = bnd.∂R
        bc = bnd.bc
        ℓ = bnd.ℓ

        dx = (ℓ[1]/N[1], ℓ[2]/N[2])
        xy = (∂R[1] + dx[1]/2 + dx[1]*(0:N[1]-1), ∂R[3] + dx[2]/2 + dx[2]*(0:N[2]-1))
        XY = ( repeat(xy[1];outer=(1,N[2])) , repeat(transpose(xy[2]);outer=(N[1],1)) )

        extinction, change_per_site, power_law, α_imag = PML_params()
        dN_PML = Int[0, 0, 0, 0]
        for i in 1:4
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
        N_PML = [N[1] + dN_PML[1] + dN_PML[2], N[2] + dN_PML[3] + dN_PML[4]]
        ∂R_PML = [∂R[1] - dx[1]*(dN_PML[1]) - dx[1]/2,
                ∂R[2] + dx[1]*(dN_PML[2]) + dx[1]/2,
                ∂R[3] - dx[2]*(dN_PML[3]) - dx[2]/2,
                ∂R[4] + dx[2]*(dN_PML[4]) + dx[2]/2]
        ℓ_PML = [∂R_PML[2]-∂R_PML[1], ∂R_PML[4]-∂R_PML[3]]
        xy_PML = (∂R_PML[1] + dx[1]/2 + dx[1]*(0:N_PML[1]-1), ∂R_PML[3] + dx[2]/2 + dx[2]*(0:N_PML[2]-1))
        XY_PML = ( repeat(xy_PML[1]; outer=(1,N_PML[2])) , repeat(transpose(xy_PML[2]); outer=(N_PML[1],1)) )
        #xy_inds = reshape(kron(dN_PML[1]+(1:N[1]), ones(dN_PML[3]+(1:N[2]))), N[1],N[2])[:]
        # @time xy_inds = ( (A,B) -> sub2ind((N_PML...),A,B)).( repeat(dN_PML[1]+(1:N[1]); inner=(1,N[2])) , repeat(transpose(dN_PML[3]+(1:N[2])); outer=(N[1],1)))[:]
        xy_inds = zeros(Int,N[1],N[2])
        for i in 1:N[1], j in 1:N[2]
            xy_inds[i,j] = sub2ind((N_PML...), dN_PML[1]+(1:N[1])[i],dN_PML[3]+(1:N[2])[j])
        end

        bnd.∂R_PML = ∂R_PML
        bnd.ℓ = ℓ
        bnd.ℓ_PML = ℓ_PML

        dis.N_PML = N_PML
        dis.dN_PML = dN_PML
        dis.xy = xy
        dis.XY = XY
        dis.dx = dx
        dis.xy_PML = xy_PML
        dis.XY_PML = XY_PML
        dis.xy_inds = xy_inds[:]

        # define system parameters
        n₁_vals = sys.n₁_vals
        n₁_inds = sys.n₁_inds
        n₂_vals = sys.n₂_vals
        n₂_inds = sys.n₂_inds
        F_vals = sys.F_vals
        F_inds = sys.F_inds

        n = vcat([wgs[m].ind for m in 1:length(wgs)], n₁_vals[n₁_inds] + 1.0im*n₂_vals[n₂_inds])
        ε = n.^2
        F = vcat(zeros(Float64,length(wgs)), F_vals[F_inds])
        F[iszero.(F)] = TwoLevelSystemDefaults()
        sys.n = n
        sys.ε = ε
        sys.F = F

        sys.regions = which_region(sys, wgs, dis)
        ɛ_sm, F_sm = sub_pixel_smoothing(sys, wgs, dis)

        sys.ε_sm = ε_sm
        sys.F_sm = F_sm

        temp = vcat(ones(dN_PML[1],1)*transpose(ε_sm[1,:]), ε_sm, ones(dN_PML[2],1)*transpose(ε_sm[end,:]))
        ε_PML = hcat(temp[:,1]*ones(1,dN_PML[3]), temp, temp[:,end]*ones(1,dN_PML[4]))
        sys.ε_PML = ε_PML
        temp = vcat(ones(dN_PML[1],1)*transpose(F_sm[1,:]), F_sm, ones(dN_PML[2],1)*transpose(F_sm[end,:]))
        F_PML = hcat(temp[:,1]*ones(1,dN_PML[3]), temp, temp[:,end]*ones(1,dN_PML[4]))
        sys.F_PML = F_PML

        for w in 1:length(wgs)
            wgs[w].ε = sub_pixel_smoothing(wgs[w], dis)
            if wgs[w].dir == "x"
                wgs[w].ε_PML = vcat(wgs[w].ε[1]*ones(Float64,dN_PML[3]),wgs[w].ε,wgs[w].ε[end]*ones(Float64,dN_PML[4]))
            elseif wgs[w].dir == "y"
                wgs[w].ε_PML = vcat(wgs[w].ε[1]*ones(Float64,dN_PML[1]),wgs[w].ε,wgs[w].ε[end]*ones(Float64,dN_PML[2]))
            end
        end

        temp = zeros(Complex128, N[1], N[2])
        sct.ε₀ = zeros(Complex128, N[1], N[2])
        for w in 1:length(wgs)
            if wgs[w].dir=="x"
                temp += repeat(transpose(wgs[w].ε);inner=(N[1],1))
            elseif waveguides[w].dir=="y"
                temp += wgs[w].ε*ones(1,N[2])
            end
        end
        sct.ε₀ = temp .- length(wgs) .+ 1
        # temp = vcat(ones(dN_PML[1],1)*transpose(ε_sm[1,:]), sct.ε₀, ones(dN_PML[2],1)*transpose(ε_sm[end,:]))
        # sct.ε₀_PML = hcat(temp[:,1]*ones(1,dN_PML[3]), temp, temp[:,end]*ones(1,dN_PML[4]))
        ε_temp = vcat(repmat(transpose(sct.ε₀[1,:]),dN_PML[1],1), sct.ε₀, repmat(transpose(sct.ε₀[1,:]),dN_PML[2],1))
        sct.ε₀_PML = hcat(repmat(ε_temp[:,1],1,dN_PML[3]), ε_temp, repmat(ε_temp[:,1],1,dN_PML[4]))

        return new(sys, bnd, dis, sct, tls, wgs)
    end
end # end of struct InputStruct

# ################################################################################
# ##### SYSTEM DEFINITION AUXILLIARIES
# ################################################################################
"""
r = which_region(XY, system, waveguides, discretization)
"""
function which_region(xy::Tuple{Array{Float64,1},Array{Float64,1}}, sys::SystemStruct,
    wgs::Array{WaveguideStruct,1}, dis::DiscretizationStruct)::Array{Int,2}

    x = xy[1]
    y = xy[2]
    geometry = sys.geometry
    geoParams = sys.geoParams

    regions = zeros(Int,length(x),length(y))

    for i in 1:length(x), j in 1:length(y)
        regions[i,j] = length(wgs) + geometry(x[i], y[j], geoParams)
        for w in 1:length(wgs)
            if wgs[w].dir in ["x", "X"]
                p = y[j]
            elseif wgs[w].dir in ["y", "Y"]
                p = x[i]
            else
                error("Invalid waveguide direction.")
            end
            if wgs[w].pos-wgs[w].wdt/2 < p < wgs[w].pos+wgs[w].wdt/2
                regions[i,j] = w
            end
        end
    end
    return regions
end
"""
r = which_region(system, waveguides, discretization)
"""
function which_region(sys::SystemStruct, wgs::Array{WaveguideStruct,1}, dis::DiscretizationStruct)::Array{Int,2}
    xy = (dis.xy[1], dis.xy[2])
    regions = which_region(xy, sys, wgs, dis)
    return regions
end
"""
r = which_region(X, waveguide, discretization)
"""
function which_region(x::Array{Float64,1}, wg::WaveguideStruct, dis::DiscretizationStruct)::Array{Int,1}
    regions = 2*ones(Int,length(x))
    wg_bool = wg.pos-wg.wdt/2 .< x .< wg.pos+wg.wdt/2
    regions[wg_bool] = 1
    return regions
end
"""
r = which_region(waveguide, discretization)
"""
function which_region(wg::WaveguideStruct, dis::DiscretizationStruct)::Array{Int,1}
    if wg.dir == "x"
        x = dis.xy[2]
    elseif wg.dir == "y"
        x = dis.xy[1]
    end
    regions = which_region(x, wg, dis)
    return regions
end # end of function which_region


"""
ε_sm, F_sm = sub_pixel_smoothing(system)
"""
function sub_pixel_smoothing(sys::SystemStruct, wgs::Array{WaveguideStruct,1},
    dis::DiscretizationStruct)::Tuple{Array{Complex128,2},Array{Float64,2}}

    x = dis.xy[1]
    y = dis.xy[2]
    sub_pixel_num = dis.sub_pixel_num

    r = sys.regions
    ε = sys.ε
    F = sys.F

    ɛ_sm = ɛ[r]
    F_sm = F[r]

    X = ones(Float64,sub_pixel_num)
    Y = ones(Float64,sub_pixel_num)
    for i in 2:(length(x)-1), j in 2:(length(y)-1)
        ( nearestNeighborFlag = (r[i,j]!==r[i,j+1]) | (r[i,j]!==r[i,j-1]) |
        (r[i,j]!==r[i+1,j]) | (r[i,j]!==r[i-1,j]) )

        ( nextNearestNeighborFlag = (r[i,j]!==r[i+1,j+1]) |
        (r[i,j]!==r[i-1,j-1]) | (r[i,j]!==r[i+1,j-1]) | (r[i,j]!==r[i-1,j+1]) )
        if nearestNeighborFlag | nextNearestNeighborFlag
            x_min = (x[i]+x[i-1])/2
            y_min = (y[j]+y[j-1])/2

            x_max = (x[i]+x[i+1])/2
            y_max = (y[j]+y[j+1])/2

            X = Array(linspace(x_min,x_max,sub_pixel_num))
            Y = Array(linspace(y_min,y_max,sub_pixel_num))

            sub_regions = which_region((X,Y), sys, wgs, dis)
            ɛ_sm[i,j] = mean(ɛ[sub_regions])
            F_sm[i,j] = mean(F[sub_regions])
        end
    end

    ε_sm[1,:] = ε_sm[2,:]
    ε_sm[end,:] = ε_sm[end-1,:]
    ε_sm[:,1] = ε_sm[:,2]
    ε_sm[:,end] = ε_sm[:,end-1]

    return ɛ_sm, F_sm
end
function sub_pixel_smoothing(wg::WaveguideStruct, dis::DiscretizationStruct)::Array{Complex128,1}

    if wg.dir == "x"
        x = dis.xy[2]
    elseif wg.dir == "y"
        x = dis.xy[1]
    end
    sub_pixel_num = dis.sub_pixel_num

    r = which_region(wg, dis)
    ε = complex([wg.ind^2,1.])

    ɛ_sm = ɛ[r]

    X = zeros(Float64,sub_pixel_num)
    for i in 2:(length(x)-1)
        ( nearestNeighborFlag = (r[i]!==r[i+1]) | (r[i]!==r[i-1]) |
        (r[i]!==r[i+1]) | (r[i]!==r[i-1]) )

        if nearestNeighborFlag
            x_min = (x[i]+x[i-1])/2
            x_max = (x[i]+x[i+1])/2

            X = Array(linspace(x_min,x_max,sub_pixel_num))

            sub_regions = which_region(X, wg, dis)
            ɛ_sm[i] = mean(ɛ[sub_regions])
        end
    end
    return ɛ_sm
end

# ################################################################################
# ##### SYSTEM MODIFICATION
# ################################################################################
(export update_input!)

"""
update_input!(input, field, value)
update_input!(input, fields, values)

    If changes were made to the ∂R, N, k₀, k, F, ɛ, Γ, bc, a, b, run update_input to
    propagate these changes through the system.
"""
function update_input!(input::InputStruct, fields::Array{Symbol,1}, values::Array{Any,1})::InputStruct

    for m in 1:length(fields)
        f = fieldnames(input)
        for i in 1:length(f)
            g = fieldnames(getfield(input,f[i]))
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
                        :n₂_inds, :ind, :wdt, :pos, :dir, :sub_pixel_num, :geoParams]
    end
    if any(flag)
        sys = SystemStruct(input.sys)
        bnd = BoundaryStruct(input.bnd)
        dis = DiscretizationStruct(input.dis)
        sct = ScatteringStruct(input.sct)
        tls = TwoLevelSystemStruct(input.tls)
        wgs = WaveguideStruct(input.wgs)
        input = InputStruct(sys, bnd, dis, sct, tls, wgs)
    end

    return input
end # end of function update_input
function update_input!(input::InputStruct, field::Symbol, value::Any)::InputStruct
    update_input!(input,[field],[value])
    return input
end # end of function update_input

#
#
#
# function set_bc(input::InputStruct, k_type::String, direction::Array{Int,1})::InputStruct
#
#     if k_type ∈ ["Pole","pole","P","p"]
#         input1 = open_to_pml_out(input)
#     elseif k_type ∈ ["Zero","zero","Z","z"]
#         input1 = open_to_pml_in(input)
#     elseif k_type ∈ ["UZR","uzr","U","u"]
#         input1 = deepcopy(input)
#         if (direction[1]==+1) && (input1.bc[1:2] !== ["I", "O"])
#             update_input!(input1, :bc, ["pml_in", "pml_out", input1.bc[3], input1.bc[4]])
#         elseif (direction[1]==-1) && (input1.bc[1:2] !== ["O", "I"])
#             update_input!(input1, :bc, ["pml_out", "pml_in", input1.bc[3], input1.bc[4]])
#         end
#         if (direction[2]==+1) && (input1.bc[3:4] !== ["I", "O"])
#             update_input!(input1, :bc, [input1.bc[1], input1.bc[2], "pml_in", "pml_out"])
#         elseif (direction[2]==-1) && (input1.bc[3:4] !== ["O", "I"])
#             update_input!(input1, :bc, [input1.bc[1], input1.bc[2], "pml_out", "pml_in"])
#         end
#     end
#     return input1
# end
#
#
# """
# input = open_to_pml_out(input)
#
#     converts "open" or "pml_in" to "pml_out", and creates a copy of input if it does so.
# """
# function open_to_pml_out(input1::InputStruct)::InputStruct
#     if any(input1.bc .== "o") || any(input1.bc .== "I")
#         input = deepcopy(input1)
#         for i in 1:4
#             if input.bc[i] in ["o", "I"]
#                 input.bc[i] = "pml_out"
#             end
#         update_input!(input, :bc, input.bc)
#         end
#     else
#         input = input1
#     end
#     return input
# end
# function open_to_pml_out(input1::InputStruct, flag::Bool)::InputStruct
#     input = deepcopy(input1)
#     for i in 1:4
#         if input.bc[i] in ["o", "I"]
#             input.bc[i] = "pml_out"
#         end
#         update_input!(input, :bc, input.bc)
#     end
#     return input
# end
#
#
# """
# input = open_to_pml_in(input)
#
#     converts "open" or "pml_out" to "pml_in", and creates a copy of input if it does so.
# """
# function open_to_pml_in(input1::InputStruct)::InputStruct
#     if any(input1.bc .== "o") || any(input1.bc .== "O")
#         input = deepcopy(input1)
#         for i in 1:4
#             if input.bc[i] in ["o", "O"]
#                 input.bc[i] = "pml_in"
#                 update_input!(input, :bc, input.bc)
#             end
#         end
#     else
#         input = input1
#     end
#     return input
# end
#
#
#
#
#
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
