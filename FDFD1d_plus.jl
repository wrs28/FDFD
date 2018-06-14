module FDFD1d_plus

using Formatting
using Interpolations
using PyPlot

################################################################################
include("FDFD1d_defaults.jl")
include("FDFD1d_structures.jl")
include("FDFD1d_core_numerics.jl")
include("FDFD1d_scattering.jl")
include("FDFD1d_synthesis.jl")
include("FDFD1d_smatrix.jl")
include("FDFD1d_analysis.jl")
include("FDFD1d_eigensolvers.jl")
include("FDFD1d_parallel.jl")

################################################################################
##### WAVE PLOT
################################################################################
export wave_plot
"""
wavePlot(input, ψ)
"""
function wave_plot(input::InputStruct, ψ::Union{Array{Complex{Float64},1}, Array{Complex{Float64},2}};
    array::Bool = true, how::Function = abs2)::Tuple{Array{Float64,1},Array{Complex128,2}}

    N = size(ψ,2)

    if size(ψ,1) == input.dis.N
        M = input.dis.N
        x₁ = input.dis.x
        ∂R = input.bnd.∂R
        ε_sm = input.sys.ε_sm
        F_sm = input.sys.F_sm
    else
        M = input.dis.N_PML
        x₁ = input.dis.x_PML
        ∂R = input.bnd.∂R_PML
        ε_sm = input.sys.ε_PML
        F_sm = input.sys.F_PML
    end

    ψ_plot = NaN*zeros(Complex128, M, N)

    for i in 1:N
        ψ_plot[:,i] = ψ[:,i]
    end

    if array

        figure(figsize=4.8*[3,(N+1)])
        subplots_adjust(hspace=0.0)
        subplots_adjust(wspace=0.0)

        subplot(N+1,3,1);
        plot(x₁, real(sqrt.(ɛ_sm-1im*input.tls.D₀*F_sm)))
        axt = gca()
        xlim( [ ∂R[1],∂R[2] ] )
        # ylim( [ ∂R[3],∂R[4] ] )
        axis("tight")
        xlabel("Real")
        axt[:xaxis][:tick_top]()
        axt[:xaxis][:set_label_position]("top")

        subplot(N+1,3,2); axt = gca()
        plot(x₁, imag(sqrt.(ɛ_sm-1im*input.tls.D₀*F_sm)))
        xlim( [ ∂R[1],∂R[2] ] )
        # ylim( [ ∂R[3],∂R[4] ] )
        axis("tight")
        xlabel("Imag")
        setp(axt[:get_xticklabels](),visible=false)
        setp(axt[:get_yticklabels](),visible=false)
        axt[:xaxis][:tick_top]()
        axt[:xaxis][:set_label_position]("top")

        subplot(N+1,3,3); axt = gca()
        plot(x₁, F_sm)
        xlim( [ ∂R[1],∂R[2] ] )
        # ylim( [ ∂R[3],∂R[4] ] )
        axis("tight")
        xlabel("F/abs2")
        setp(axt[:get_xticklabels](),visible=false)
        setp(axt[:get_yticklabels](),visible=false)
        axt[:xaxis][:tick_top]()
        axt[:xaxis][:set_label_position]("top")

        for i in 1:N

            subplot(N+1,3,3+3*(i-1)+1); axt = gca()
            plot(x₁, real(ψ_plot[:,i]), "b")
            xlim( [ ∂R[1],∂R[2] ] )
            # ylim( [ ∂R[3],∂R[4] ] )
            axis("tight")
            setp(axt[:get_xticklabels](),visible=false)
            setp(axt[:get_yticklabels](),visible=false)

            subplot(N+1,3,3+3*(i-1)+2); axt = gca()
            plot(x₁, imag(ψ_plot[:,i]), "r")
            xlim( [ ∂R[1],∂R[2] ] )
            # ylim( [ ∂R[3],∂R[4] ] )
            axis("tight")
            setp(axt[:get_xticklabels](),visible=false)
            setp(axt[:get_yticklabels](),visible=false)

            subplot(N+1,3,3+3*(i-1)+3); axt = gca()
            plot(x₁, abs2.(ψ_plot[:,i]), "k")
            xlim( [ ∂R[1],∂R[2] ] )
            # ylim( [ ∂R[3],∂R[4] ] )
            axis("tight")
            setp(axt[:get_xticklabels](),visible=false)
            setp(axt[:get_yticklabels](),visible=false)

        end

    else

        for i in 1:N
            figure(i, figsize=8*[1,(∂R[4]-∂R[3])/(∂R[2]-∂R[1])])
            if how==abs2
                cmap = "k"
            elseif how==imag
                cmap = "r"
            else
                cmap = "b"
            end
            plot(x₁, how.(ψ_plot[:,i]),color=cmap)
            xlim( [ ∂R[1],∂R[2] ] )
            # ylim( [ ∂R[3],∂R[4] ] )
            axis("equal")
        end

    end

    return (x₁, ψ_plot)
end # end of function wave_plot

end #end of module FDFD1d_plus
