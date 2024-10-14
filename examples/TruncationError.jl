using Plots
using WaveResolvingBQ
using Statistics
using BenchmarkTools
using LaTeXStrings

S=WaveResolvingBQ.UP1(2)
B=WaveResolvingBQ.Periodic1D(4)


############### FUNCTION ##################
# See https://oceanrep.geomar.de/id/eprint/52876/1/K%C3%A4mpf2009_Book_OceanModellingForBeginners.pdf (p. 67)
#
A=0.5
λ=5.0
x=G.x
function f(x)
    return A .* sin.(x .* 2*π/λ)
end
function dfdx(x)
    return 2*π*A/λ .* cos.(x .* 2*π/λ)
end


############### TRUNCATION ERROR ##################
#
function ϵ(Δx)
    return 1 .- sin(2*π*Δx/λ) / (2*π*Δx/λ)
end

############### COMPUTING ################
ΔX=[i*λ for i=0.1:0.01:0.5]
ϵ_num=zero(ΔX)
ϵ_ana=zero(ΔX)
n=1
for dx=ΔX
    G=WaveResolvingBQ.Grid(Lx=10.0,Δx=dx)
    df_ana=dfdx(G.x)
    df_num=0.5*(WaveResolvingBQ.discretize(S,f(G.x),G,1,-1.0)+WaveResolvingBQ.discretize(S,f(G.x),G,1,1.0)) # centered scheme
    #df_num[2:end-1]=(f(G.x[3:end]) - f(G.x[1:end-2]))/(2*G.Δx)
    tmp=(df_ana .-df_num)./df_ana
    ϵ_num[n]=mean(tmp[2:end-2])
    ϵ_ana[n]=ϵ(G.Δx) 
    n+=1
end

############### PLOTTING ##################
#
plot(ΔX/λ,ϵ_ana*100,title="Truncation error test case",ylabel="Relative error (%)",xlabel="Δx/λ",line=(4,:black),label="Theoretical ϵ")
plot!(ΔX/λ,ϵ_num*100,label="Numerical ϵ",line=(4,:red,:dash))