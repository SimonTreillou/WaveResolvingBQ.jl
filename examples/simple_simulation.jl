using WaveResolvingBQ
using Plots

# Simulation parameters
L = 1000.0     # Domain length (m)
Nx = 1000      # Number of spatial points
h = 10.0       # Undisturbed water depth (m)
g = 9.81       # Gravity acceleration (m/s²)
T = 10.0       # Total simulation time (s)
dt = 0.01      # Time step size (s)

# Run simulation
state = WaveResolvingBQ.run_simulation(L, Nx, h, g, T, dt)

# Plot final surface elevation
plot(state.x, state.η, title="Surface Elevation at t=$(state.t) s",
     xlabel="x (m)", ylabel="η (m)", legend=false)
savefig("surface_elevation_final.png")