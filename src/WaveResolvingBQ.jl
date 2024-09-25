module WaveResolvingBQ

import Statistics
import Base.@kwdef

include("Parameters/Parameters.jl")

export SimulationParameters
export check_parameters

include("TimeStepping/TimeStepping.jl")

export TimeStepper
export timestepping

include("Discretization/Discretization.jl")

export SpatialScheme
export discretize

include("Model/Model.jl")

export Model
export modelize 

include("Utils/Utils.jl")

export window

function run(M::Model,T::TimeStepper,S::SpatialScheme,P::SimulationParameters,η)
    TV = zeros(P.Nt)
    a=0.0
    b=0.0
    for step in 1:P.Nt

        function deriv(η, Δx, c, a, b)
            return modelize(M,S,η, Δx, c, a, b)
        end
        η_new = timestepping(T,η, P.Δt, P.Δx, P.c, a, b,deriv)

        # Check for NaN or Inf
        if any(isnan, η_new) || any(isinf, η_new)
            println("Numerical instability detected at step $step, time $(step * Δt)")
            break
        end
        
        # Update previous and current states
        η_old = copy(η)
        η = copy(η_new)
        TV[step] = Statistics.var(η)
        
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