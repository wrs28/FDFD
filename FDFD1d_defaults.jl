"""
PML_params()
"""
function PML_params()::Tuple{Float64,Float64,Float64,Float64}

    extinction = 5e2
    change_per_site = 2/3
    power_law = 4
    α_imag = -.25

    return extinction, change_per_site, power_law, α_imag
end # end function PML_params


function TwoLevelSystemDefaults()::Float64
    F_min = 1e-16
end
