module FDFD2d_plus

using Formatting
using Interpolations
using SpecialFunctions
using PyPlot

################################################################################
include("FDFD2d_defaults.jl")
include("FDFD2d_structures.jl")
include("FDFD2d_core_numerics.jl")
include("FDFD2d_scattering.jl")
include("FDFD2d_synthesis.jl")
include("FDFD2d_smatrix.jl")
include("FDFD2d_analysis.jl")
include("FDFD2d_eigensolvers.jl")
include("FDFD2d_parallel.jl")

################################################################################
##### WAVE PLOT
################################################################################
export wave_plot
"""
wavePlot(ψ, inputs)
"""
function wave_plot(input::InputStruct, ψ::Union{Array{Complex{Float64},1}, Array{Complex{Float64},2}};
    array::Bool = true, how::Function = abs2)::
    Tuple{Array{Float64,1},Array{Float64,1},Array{Complex128,3}}

    N = size(ψ,2)

    if size(ψ,1) == prod(input.dis.N)
        M = input.dis.N
        x₁ = input.dis.xy[1]
        x₂ = input.dis.xy[2]
        ∂R = input.bnd.∂R
        ε_sm = input.sys.ε_sm
        F_sm = input.sys.F_sm
    else
        M = input.dis.N_PML
        x₁ = input.dis.xy_PML[1]
        x₂ = input.dis.xy_PML[2]
        ∂R = input.bnd.∂R_PML
        ε_sm = input.sys.ε_PML
        F_sm = input.sys.F_PML
    end

    ψ_plot = NaN*zeros(Complex128, M[1], M[2], N)

    for i in 1:N
        ψ_plot[:,:,i] = reshape(ψ[:,i], (M[1],M[2]) )
    end

    if array

        figure(figsize=4.8*[3,(N+1)*(∂R[4]-∂R[3])/(∂R[2]-∂R[1])])
        subplots_adjust(hspace=0.0)
        subplots_adjust(wspace=0.0)

        subplot(N+1,3,1); axt = gca()
        pcolormesh(x₁, x₂, transpose(real(sqrt.(ɛ_sm-1im*input.tls.D₀*F_sm)))) #transpose(r))
        xlim( [ ∂R[1],∂R[2] ] )
        ylim( [ ∂R[3],∂R[4] ] )
        axis("tight")
        xlabel("Real")
        axt[:xaxis][:tick_top]()
        axt[:xaxis][:set_label_position]("top")

        subplot(N+1,3,2); axt = gca()
        pcolormesh(x₁, x₂, transpose(imag(sqrt.(ɛ_sm-1im*input.tls.D₀*F_sm))), cmap="bwr")
        xlim( [ ∂R[1],∂R[2] ] )
        ylim( [ ∂R[3],∂R[4] ] )
        clim([-1,1]*findmax(abs.(imag(sqrt.(ɛ_sm-1im*input.tls.D₀*F_sm))))[1])
        axis("tight")
        xlabel("Imag")
        setp(axt[:get_xticklabels](),visible=false)
        setp(axt[:get_yticklabels](),visible=false)
        axt[:xaxis][:tick_top]()
        axt[:xaxis][:set_label_position]("top")

        subplot(N+1,3,3); axt = gca()
        pcolormesh(x₁, x₂, transpose(F_sm))
        xlim( [ ∂R[1],∂R[2] ] )
        ylim( [ ∂R[3],∂R[4] ] )
        axis("tight")
        xlabel("F/abs2")
        setp(axt[:get_xticklabels](),visible=false)
        setp(axt[:get_yticklabels](),visible=false)
        axt[:xaxis][:tick_top]()
        axt[:xaxis][:set_label_position]("top")

        for i in 1:N

            subplot(N+1,3,3+3*(i-1)+1); axt = gca()
            pcolormesh(x₁, x₂,transpose(real(ψ_plot[:,:,i])),cmap = "bwr")
            xlim( [ ∂R[1],∂R[2] ] )
            ylim( [ ∂R[3],∂R[4] ] )
            axis("tight")
            clim([-1,1])
            setp(axt[:get_xticklabels](),visible=false)
            setp(axt[:get_yticklabels](),visible=false)

            subplot(N+1,3,3+3*(i-1)+2); axt = gca()
            pcolormesh(x₁, x₂,transpose(imag(ψ_plot[:,:,i])),cmap = "bwr")
            xlim( [ ∂R[1],∂R[2] ] )
            ylim( [ ∂R[3],∂R[4] ] )
            axis("tight")
            clim([-1,1])
            setp(axt[:get_xticklabels](),visible=false)
            setp(axt[:get_yticklabels](),visible=false)

            subplot(N+1,3,3+3*(i-1)+3); axt = gca()
            pcolormesh(x₁, x₂,transpose(abs2.(ψ_plot[:,:,i])),cmap = "gray_r")
            xlim( [ ∂R[1],∂R[2] ] )
            ylim( [ ∂R[3],∂R[4] ] )
            axis("tight")
            clim([0,1])
            setp(axt[:get_xticklabels](),visible=false)
            setp(axt[:get_yticklabels](),visible=false)

        end

    else

        for i in 1:N
            figure(i, figsize=8*[1,(∂R[4]-∂R[3])/(∂R[2]-∂R[1])])
            if how==abs2
                cmap = "gray_r"
            else
                cmap = "bwr"
            end
            pcolormesh(x₁, x₂,transpose(how.(ψ_plot[:,:,i])),cmap=cmap)
            xlim( [ ∂R[1],∂R[2] ] )
            ylim( [ ∂R[3],∂R[4] ] )
            axis("equal")
        end

    end

    return (x₁,x₂,ψ_plot)
end # end of function wavePlot

end #end of module SALT_2d_plus
