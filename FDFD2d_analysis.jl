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
