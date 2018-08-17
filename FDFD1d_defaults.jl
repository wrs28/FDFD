"""
PML_params()
"""
function PML_params()::Tuple{Float64,Float64,Float64,Float64}

    extinction = 1e6
    change_per_site = 1/4
    power_law = 2
    α_imag = -.1

    return extinction, change_per_site, power_law, α_imag
end # end function PML_params


function TwoLevelSystemDefaults()::Float64
    F_min = 1e-16
end
