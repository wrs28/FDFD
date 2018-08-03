"""
PML_params()
"""
function PML_params()::Tuple{Float64,Float64,Float64,Float64}

    extinction = 1e7
    change_per_site = 1/50
    power_law = 1
    α_imag = -.15

    return extinction, change_per_site, power_law, α_imag
end # end function PML_params


function TwoLevelSystemDefaults()::Float64
    F_min = 1e-16
end
