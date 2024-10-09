using Plots
using WaveResolvingBQ
using Statistics
using BenchmarkTools
using LaTeXStrings

G=WaveResolvingBQ.Grid(Lx=10.0,Δx=0.01)
S=WaveResolvingBQ.WENO5()
S=WaveResolvingBQ.UP1()
B=WaveResolvingBQ.Periodic1D(4)


############### SINUS ##################
x=G.x
y=sin.(x)
dydx=cos.(x)
dy2dx2=-sin.(x)

############### POLYNOM ################
x=G.x
y=x.^3 
dydx=3 .* x.^2
dy2dx2=6 .* x


###########################################
############### PLOTTING ##################

plot(layout=(2,1))
plot!(x,dydx,subplot=1,label=L"$\partial y/\partial x$ théorique",line=(4,:red))
dyn =WaveResolvingBQ.discretize(S,y,G,1,velocity=y*0 .-1)
#WaveResolvingBQ.apply_boundaries!(B,G,dyn)
plot!(x,dyn,label=L"$\partial y/\partial x$ numérique",line=(2,:black))
plot!(x,dy2dx2,subplot=2,label=L"$\partial^2 y/\partial x^2$ théorique",line=(4,:red))
dyn2 =WaveResolvingBQ.discretize(S,y,G,2,velocity=y*0 .-1)
#WaveResolvingBQ.apply_boundaries!(B,G,dyn)
plot!(x,dyn2,subplot=2,label=L"$\partial^2 y/\partial x^2$ numérique",line=(2,:black))
xlabel!("x (m)",subplot=2)
ylabel!("y (m)")