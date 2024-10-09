
using NetCDF

lonatts = Dict("longname" => "x-coordinates",
               "units" => "m")
latatts = Dict("longname" => "y-coordinates",
               "units" => "m")
timeatts = Dict("longname" => "Time",
                "units" => "s")
varatts = Dict("longname" => "zeta",
                "units" => "m")


y=[1,2,3]
zeta=zeros(length(P.x),1,length(P.t)) .+ rand()
zeta = rand(length(P.x),length(P.t))
fn = joinpath("test.nc")
isfile(fn) && rm(fn)
nccreate(fn, "zeta","x",P.x,lonatts,"time",P.t,timeatts,atts=varatts)
ncwrite(zeta,fn,"zeta")
ncclose(fn)
ncinfo(fn)

tt=ncread(fn,"x")
