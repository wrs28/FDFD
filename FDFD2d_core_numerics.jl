################################################################################
####### DERIVATIVES
################################################################################
"""
∇ =  grad_1d(N, dx)

    1-dim Gradient with N points, lattice spacing dx. It's the forward gradient (I think).

    sparse ∇[N,N+1]
"""
function grad_1d(N::Int, dx::Float64)::SparseMatrixCSC{Complex128,Int}

    I₁ = Array(1:N)
    J₁ = Array(1:N)
    V₁ = fill(Complex(-1/dx), N)

    I₂ = Array(1:N)
    J₂ = Array(2:(N+1))
    V₂ = fill(Complex(+1/dx), N)

    ∇ = sparse(vcat(I₁,I₂), vcat(J₁,J₂), vcat(V₁,V₂), N, N+1, +)
end #end of function grad_1d

"""
∇₁,∇₂ =  grad(N, dx)

    2-dim gradients with N[1],N[2] points, lattice spacing dx̄[1], dx̄[2].
    It's the forward gradient (I think).

    sparse (∇₁[N,N],∇₂[N,N])
"""
function grad(N::Array{Int,1}, dx::Array{Float64,1})::
    Tuple{SparseMatrixCSC{Complex128,Int},SparseMatrixCSC{Complex128,Int}}

    N₁ = N[1]
    dx₁ = dx[1]

    N₂ = N[2]
    dx₂ = dx[2]

    ∇₁ = grad_1d(N₁-1,dx₁)
    ∇₂ = grad_1d(N₂-1,dx₂)

    ∇₁ = kron(speye(N₂,N₂),∇₁)
    ∇₂ = kron(∇₂,speye(N₁,N₁))

    return ∇₁,∇₂
end # end of function grad

"""
∇₁², ∇₂² = laplacians(input, k)

    Computes 2-dim laplacian based on parameters and boundary conditions given in
    INPUTS, evaluated at (complex) frequency K.
"""
function laplacians(input::InputStruct, k::Complex128)::
    Tuple{SparseMatrixCSC{Complex128,Int},SparseMatrixCSC{Complex128,Int}}

    # definitions block#
    ∂R = input.bnd.∂R_PML
    bc = input.bnd.bc

    k₁ = input.bnd.bk[1]
    k₂ = input.bnd.bk[2]

    ℓ₁ = input.bnd.ℓ_PML[1]
    ℓ₂ = input.bnd.ℓ_PML[2]

    N₁ = input.dis.N_PML[1]
    N₂ = input.dis.N_PML[2]

    dx₁ = input.dis.dx[1]
    dx₂ = input.dis.dx[2]

    Σ₁,Σ₂ = σ(input)

    ∇₁ = grad_1d(N₁-1,dx₁)
    ∇₂ = grad_1d(N₂-1,dx₂)

    s₁₁ = sparse(1:N₁-1,1:N₁-1,1./(1+.5im*(Σ₁[1:end-1] + Σ₁[2:end])/real(k)),N₁-1,N₁-1)
    s₁₂ = sparse(1:N₁,1:N₁,1./(1+1im*(Σ₁)/real(k)),N₁,N₁)

    s₂₁ = sparse(1:N₂-1,1:N₂-1,1./(1+.5im*(Σ₂[1:end-1] + Σ₂[2:end])/real(k)),N₂-1,N₂-1)
    s₂₂ = sparse(1:N₂,1:N₂,1./(1+1im*(Σ₂)/real(k)),N₂,N₂)

    ∇₁² = -(s₁₂*transpose(∇₁)*s₁₁*∇₁)
    ∇₂² = -(s₂₂*transpose(∇₂)*s₂₁*∇₂)
    ind = [1, N₁, 1, N₂, 1]

    for i in 1:4
        if i ≤ 2
            if bc[i] in ["O", "I"]
                ∇₁²[ind[i],ind[i]]   += -2/dx₁^2
            elseif bc[i] == "d"
                ∇₁²[ind[i],ind[i]]   += -2/dx₁^2
            elseif bc[i] == "n"
                ∇₁²[ind[i],ind[i]]   += 0
            elseif bc[i] == "p"
                ∇₁²[ind[i],ind[i]]   += -1/dx₁^2
                ∇₁²[ind[i],ind[i+1]] += +exp((-1)^(i+1)*1im*ℓ₁*k₁)/dx₁^2
            end
        else
            if bc[i] in ["O", "I"]
                ∇₂²[ind[i],ind[i]]   += -2/dx₂^2
            elseif bc[i] == "d"
                ∇₂²[ind[i],ind[i]]   += -2/dx₂^2
            elseif bc[i] == "n"
                ∇₂²[ind[i],ind[i]]   += 0
            elseif bc[i] == "p"
                ∇₂²[ind[i],ind[i]]   += -1/dx₂^2
                ∇₂²[ind[i],ind[i+1]] += +exp((-1)^(i+1)*1im*ℓ₂*k₂)/dx₂^2
            end
        end
    end

    return ∇₁², ∇₂²
end # end of function laplacians

"""
∇² = laplacian(input, k)

    Computes 2-dim laplacian based on parameters and boundary conditions given in
    INPUTS, evaluated at (complex) frequency K.
"""
function laplacian(input::InputStruct, k::Complex128)::SparseMatrixCSC{Complex128,Int}

    # definitions block#
    ∂R = input.bnd.∂R_PML
    bc = input.bnd.bc

    k₁ = input.bnd.bk[1]
    k₂ = input.bnd.bk[2]

    ℓ₁ = input.bnd.ℓ_PML[1]
    ℓ₂ = input.bnd.ℓ_PML[2]

    N₁ = input.dis.N_PML[1]
    N₂ = input.dis.N_PML[2]

    dx₁ = input.dis.dx[1]
    dx₂ = input.dis.dx[2]

    ∇₁², ∇₂² = laplacians(input,k)

    if any(bc .== "o")

        if any(bc[1:2] .== "o") && !(bc[3:4]⊆["d", "n"])
            error("Inconsistent boundary conditions. Cannot have an open side without Dirichlet/Neumann top and bottom.")
        elseif any(bc[3:4] .== "o") && !(bc[1:2]⊆["d", "n"])
            error("Inconsistent boundary conditions. Cannot have open top or bottom without Dirichlet sides.")
        end

        if any(input.bc[1:2] .== "o")
            N⟂ = N₂
            Φ = zeros(Complex128,N⟂,N⟂,2)
            x⟂ = x₂ - ∂R[3]
            ℓ⟂ = ℓ₂
        else
            N⟂ = N₁
            Φ = zeros(Complex128,N⟂,N⟂,2)
            x⟂ = x₁ - ∂R[1]
            ℓ⟂ = ℓ₁
        end

        m_cutoff = 2*floor(Int,ℓ⟂*sqrt(real(k^2))/π)
    end

    for i in 1:4
        if i ≤ 2
            if bc[i] == "o"
                for m in 1:m_cutoff
                    if m in input.input_modes[i]
                        M = +1
                    else
                        M = -1
                    end
                    k⟂ = M*sqrt(k^2 - (m*π/ℓ⟂)^2)
                    Φ[:,:,i] += (1im*dx₁*dx₂/ℓ⟂)*k⟂*sin.(m*π*x⟂/ℓ⟂)*sin.(m*π*transpose(x⟂)/ℓ⟂)
                end
                Φ[:,:,i] = -(eye(Complex128,N⟂,N⟂)+Φ[:,:,i])\(Φ[:,:,i]*2/dx₁^2)
            end
        else
            if bc[i] == "o"
                for m in 1:m_cutoff
                    if m in input.input_modes[i]
                        M = +1
                    else
                        M = -1
                    end
                    k⟂ = M*sqrt(k^2 - (m*π/ℓ⟂)^2)
                    Φ[:,:,i-2] += (1im*dx₁*dx₂/ℓ⟂)*k⟂*sin.(m*π*x⟂/ℓ⟂)*sin.(m*π*transpose(x⟂)/ℓ⟂)
                end
                Φ[:,:,i-2] = -(eye(Complex128,N⟂,N⟂)+Φ[:,:,i-2])\(Φ[:,:,i-2]*2/dx₂^2)
            end
        end
    end

    if !any(bc .== "o")
        ∇² = kron(speye(N₂,N₂),∇₁²) + kron(∇₂²,speye(N₁,N₁))
    elseif any(bc[1:2] .== "o")
        ( ∇² = kron(speye(N₂,N₂),∇₁²) + kron(∇₂²,speye(N₁,N₁)) + kron(Φ[:,:,1],
            sparse([1],[1],[1],N₁,N₁)) + kron(Φ[:,:,2],sparse([N₁],[N₁],[1],N₁,N₁)) )
    elseif any(bc[3:4] .== "o")
        ( ∇² = kron(speye(N₂,N₂),∇₁²) + kron(∇₂²,speye(N₁,N₁)) +
            kron(sparse([1],[1],[1],N₂,N₂),Φ[:,:,1]) +
            kron(sparse([N₂],[N₂],[1],N₂,N₂),Φ[:,:,2]) )
    end

    return ∇²
end # end of function laplacian

################################################################################
######## PML
################################################################################

"""
s₁, s₂ = σ(input)

    Computes conductivity for PML layer in dimensions 1 and 2.
"""
function σ(input::InputStruct)::Tuple{Array{Complex128,1},Array{Complex128,1}}

    extinction, change_per_site, power_law, α_imag = PML_params()

    x₁ = input.dis.xy_PML[1]
    x₂ = input.dis.xy_PML[2]
    ∂R = input.bnd.∂R

    dx₁ = input.dis.dx[1]
    dx₂ = input.dis.dx[2]

    α = zeros(Complex128,4)
    for i in 1:4
        if i in 1:2
            dx = dx₁
        else
            dx = dx₂
        end
        if input.bnd.bc[i] == "O"
            ( α[i] = +(1+α_imag*1im)*( (change_per_site/dx)^(power_law+1) )/
                ( (power_law+1)*log(extinction) )^power_law )
        elseif input.bnd.bc[i] == "I"
            ( α[i] = -(1+α_imag*1im)*( (change_per_site/dx)^(power_law+1) )/
                ( (power_law+1)*log(extinction) )^power_law )
        else
            α[i] = 0
        end
    end

    s₁ = zeros(Complex128,length(x₁))
    s₂ = zeros(Complex128,length(x₂))

    for i in 1:length(x₁)
        if x₁[i] ≤ ∂R[1]
            s₁[i] = α[1]*abs(x₁[i]-∂R[1])^power_law
        elseif ∂R[2] ≤ x₁[i]
            s₁[i] = α[2]*abs(x₁[i]-∂R[2])^power_law
        end
    end

    for j in 1:length(x₂)
        if x₂[j] ≤ ∂R[3]
            s₂[j] = α[3]*abs(x₂[j]-∂R[3])^power_law
        elseif ∂R[4] ≤ x₂[j]
            s₂[j] = α[4]*abs(x₂[j]-∂R[4])^power_law
        end
    end

    return s₁,s₂
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
function trapz(z::Array{Complex128,1}, dx::Tuple{Float64,Float64})::Complex128

    ∫z_dx = prod(dx)*sum(z) # may have to address boundary terms later

    return ∫z_dx
end # end of function trapz
