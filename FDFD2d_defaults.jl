"""
PML_params()
"""
function PML_params()::Tuple{Float64,Float64,Float64,Float64}

    extinction = 1e3
    change_per_site = 1/3
    power_law = 2
    α_imag = -.25

    return extinction, change_per_site, power_law, α_imag
end # end function PML_params


function TwoLevelSystemDefaults()::Float64
    F_min = 1e-15
end

"""
wg_transverse_y_params(q)
    q is the transverse quantum number
    nev is the number of eigenvalues computed for transverse problem
    k_scale*k² defines the center of the region in which transverse eigensolver
        works
"""
function wg_transverse_y_params(q::Int)::Tuple{Int,Int}

    nev = 4 + 2*q
    k_scale = 3

    return nev, k_scale
end # end function wg_transverse_y_params


function analyze_into_angular_momentum_defaults()::Int

    nθ = Int(7e3)+1

end
