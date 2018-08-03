@everywhere function geometry(x::Float64, y::Float64, geoParams::Array{Float64,1})::Int

    local region::Int

    r = sqrt(x^2 + y^2)

    R = geoParams[1]
    R_in = geoParams[2]
    d₁ = geoParams[3]
    r₁ = geoParams[4]
    d₂ = geoParams[5]
    r₂ = geoParams[6]
    β = geoParams[7]

    x₁ = cos(0.)*(R + d₁ + r₁)
    y₁ = sin(0.)*(R + d₁ + r₁)

    x₂ = cos(β)*(R + d₂ + r₂)
    y₂ = sin(β)*(R + d₂ + r₂)

    if (R_in<r≤R)
        region = 3
    elseif r ≤ R
        region = 2
    elseif (x-x₁)^2 + (y-y₁)^2 ≤ r₁^2
        region = 4
    elseif (x-x₂)^2 + (y-y₂)^2 ≤ r₂^2
        region = 4
    else
        region = 1
    end

    return region
end

    R = 1.
    R_in = R-.4
    d₁ = 0.04
    r₁ = 0.05
    d₂ = 0.0454
    r₂ = 0.05
    β = 156.3*π/180
geoParams = [R, R_in, d₁, r₁, d₂, r₂, β]





n₁_vals = [1.0, 1.5]
n₁_inds = [1, 2, 2, 2]

n₂_vals = [0.00, 0.002]
n₂_inds = [1, 2, 2, 1]

F_vals = [0.0, 0.0]
F_inds = [1, 2, 2, 1]

sys = SystemStruct(geometry, geoParams, n₁_vals, n₁_inds, n₂_vals, n₂_inds, F_vals, F_inds)

### BOUNDARY DEFINITION
∂R = [-1.4, 1.4, -1.4, 1.4]
bc = ["pml_in", "pml_out", "d", "d"]

bnd = BoundaryStruct(∂R, bc)

### DISCRETIZATION ###
N = [451,451]
sub_pixel_num = 50

dis = DiscretizationStruct(N, sub_pixel_num)

### SCATTERING ###
∂S = [1.35] # 1.35*[-1, 1, -1, 1]
L_max = 30
channels = [ChannelStruct(m,m,"") for m in -L_max:L_max]

sct = ScatteringStruct(∂S, channels)

### TWO LEVEL SYSTEM ###
k₀ = 20.
γ⟂ = 1e8
D₀ = 0.0

tls = TwoLevelSystemStruct(D₀, k₀, γ⟂)

### WAVEGUIDES

dir = String[] #["x", "x"]
pos = Float64[]#[+R+0.15, -R-0.15]
wdt = Float64[]#[ 0.04, 0.04]
ind = Float64[]#[n, n]

wgs = WaveguideStruct(dir, pos, wdt, ind)
