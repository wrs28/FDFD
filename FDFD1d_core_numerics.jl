################################################################################
####### DERIVATIVES
################################################################################
"""
∇ =  grad(N, dx)

    1-dim Gradient with N points, lattice spacing dx. It's the forward gradient (I think).

    sparse ∇[N,N+1]
"""
function grad(N::Int, dx::Float64)::SparseMatrixCSC{Complex128,Int}

    I₁ = Array(1:N)
    J₁ = Array(1:N)
    V₁ = fill(Complex(-1/dx), N)

    I₂ = Array(1:N)
    J₂ = Array(2:(N+1))
    V₂ = fill(Complex(+1/dx), N)

    ∇ = sparse(vcat(I₁,I₂), vcat(J₁,J₂), vcat(V₁,V₂), N, N+1, +)
end #end of function grad_1d

"""
∇² = laplacian(input, k)

    Computes 2-dim laplacian based on parameters and boundary conditions given in
    INPUTS, evaluated at (complex) frequency K.
"""
function laplacian(input::InputStruct, k::Complex128)::SparseMatrixCSC{Complex128,Int}

    # definitions block#
    ∂R = input.bnd.∂R_PML
    bc = input.bnd.bc

    k₁ = input.bnd.bk

    ℓ₁ = input.bnd.ℓ_PML

    N₁ = input.dis.N_PML

    dx₁ = input.dis.dx

    Σ₁ = σ(input)

    ∇ = grad(N₁-1,dx₁)

    s₁₁ = sparse(1:N₁-1,1:N₁-1,1./(1+.5im*(Σ₁[1:end-1] + Σ₁[2:end])/real(k)),N₁-1,N₁-1)
    s₁₂ = sparse(1:N₁,1:N₁,1./(1+1im*(Σ₁)/real(k)),N₁,N₁)

    ∇² = -(s₁₂*transpose(∇)*s₁₁*∇)
    ind = [1, N₁]

    for i in 1:2
        if bc[i] in ["O", "I"]
            ∇²[ind[i],ind[i]]   += -2/dx₁^2
        elseif bc[i] == "d"
            ∇²[ind[i],ind[i]]   += -2/dx₁^2
        elseif bc[i] == "n"
            ∇²[ind[i],ind[i]]   += 0
        elseif bc[i] == "p"
            ∇²[ind[i],ind[i]]   += -1/dx₁^2
            ∇²[ind[i],ind[i+1]] += +exp((-1)^(i+1)*1im*ℓ₁*k₁)/dx₁^2
        end
    end

    return ∇²
end # end of function laplacians

################################################################################
######## PML
################################################################################

"""
s = σ(input)

    Computes conductivity for PML layer in dimensions 1 and 2.
"""
function σ(input::InputStruct)::Array{Complex128,1}

    extinction, change_per_site, power_law, α_imag = PML_params()

    x₁ = input.dis.x_PML
    ∂R = input.bnd.∂R

    dx₁ = input.dis.dx

    α = zeros(Complex128,4)
    for i in 1:2
        if input.bnd.bc[i] == "O"
            ( α[i] = +(1+α_imag*1im)*( (change_per_site/dx₁)^(power_law+1) )/
                ( (power_law+1)*log(extinction) )^power_law )
        elseif input.bnd.bc[i] == "I"
            ( α[i] = -(1+α_imag*1im)*( (change_per_site/dx₁)^(power_law+1) )/
                ( (power_law+1)*log(extinction) )^power_law )
        else
            α[i] = 0
        end
    end

    s₁ = zeros(Complex128,length(x₁))

    for i in 1:length(x₁)
        if x₁[i] ≤ ∂R[1]
            s₁[i] = α[1]*abs(x₁[i]-∂R[1])^power_law
        elseif ∂R[2] ≤ x₁[i]
            s₁[i] = α[2]*abs(x₁[i]-∂R[2])^power_law
        end
    end

    return s₁
end # end of function σ

################################################################################
###### AUXILLIARIES
################################################################################
"""
γ(input, k) is the lorentzian gain curve
"""
function γ(input::InputStruct, k::Complex128)::Complex128
    return input.tls.γ⟂./(k-input.tls.k₀+1im*input.tls.γ⟂)
end

"""
∫z_dx = trapz(z, dx)
"""
function trapz(z::Array{Complex128,1}, dx::Float64)::Complex128

    ∫z_dx = prod(dx)*sum(z) # may have to address boundary terms later

    return ∫z_dx
end # end of function trapz
