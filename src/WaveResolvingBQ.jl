module WaveResolvingBQ

import Statistics
import Base.@kwdef

include("Grid/Grid.jl")

export Grid

include("Parameters/Parameters.jl")

export SimulationParameters
export check_parameters

include("TimeStepping/TimeStepping.jl")

export TimeStepper
export timestepping

include("Boundaries/Boundaries.jl")

export Boundaries

include("Discretization/Discretization.jl")

export SpatialScheme
export discretize

include("Model/Model.jl")

export Model
export modelize 

include("Utils/Utils.jl")

export window

function run(M::Model,T::TimeStepper,S::SpatialScheme,P::SimulationParameters,B::Boundaries,G::Grid,η)
    TV = zeros(P.Nt)

    for step in 1:P.Nt

        #c=-cos.(2*π/1*P.Δt*step)*cmax
        #function deriv(η)
        #    return modelize(M,S,B,η,P)
        #end
        deriv = modelize(M,S,G,B,η)
        
        η_new = timestepping(T,P,η,deriv)

        # Check for NaN or Inf
        if any(isnan, η_new) || any(isinf, η_new)
            println("Numerical instability detected at step $step, time $(step * P.Δt)")
            break
        end
        
        # Update previous and current states
        η_old = copy(η)
        η = copy(η_new)
        TV[step] = Statistics.var(η)
        TV[step] = Statistics.mean(η)
        
        # (Optional) Visualization or data storage
        if step % 1000 == 0
            println("Time: $(step * P.Δt)")
            # You can plot η here using Plots.jl or other packages
        end
    end
    return η,TV
end

export run

end # module WaveResolvingBQ