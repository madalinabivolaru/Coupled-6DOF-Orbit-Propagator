using LinearAlgebra
using OrdinaryDiffEq
using DifferentialEquations
using StaticArrays
using Plots
using BenchmarkTools
include("./PropagatorConstants.jl")
include("Sun_position.jl")
include("Moon_position.jl")
using .PropagatorConstants
#include("./JulDat.jl")
include("asrp.jl")
#include("RVtoCOE.jl")
using DataFrames
using CSV
using CubicSplines
using SatelliteToolbox
#include("Geogrid.jl")
using Interpolations

#eop_IAU1980 = get_iers_eop_iau_2000A(force_download=false)
pos = 1.;
jd0 = 2456727.500;
p =   @SVector [10.0, 7991.0, PropagatorConstants.J2, PropagatorConstants.J3, PropagatorConstants.J4, 
     PropagatorConstants.μ_e, PropagatorConstants.μ_m, PropagatorConstants.μ_s, 
     PropagatorConstants.R_e, PropagatorConstants.ω_e, PropagatorConstants.C_d, pos, jd0]

function KepCowell!(dy::Array{Float64, 1}, y::Array{Float64, 1}, p::SVector{13, Float64}, t::Float64)#::Vector{Float64},y::Vector{Float64},p,t)
    
    A, m, J2,J3,J4, μe, μm, μs, Re, ωe, Cd, pos, jd0 = p
    jd = jd0 + t/3600/24


    #Reference densities for LEO
    rho0 = @SVector [3.725*10^-12, 1.585*10^-12,6.967*10^-13, 1.454*10^-13,
    3.614*10^-14,1.17*10^-14, 5.245*10^-15,3.019*10^-15,
    9.518*10^-12, 2.418*10^-11, 7.248*10^-11, 2.789*10^-10, 5.464*10^-10 ]

    #Access data directly from y, dy vector with view macro
    r = @view y[1:3]
    v = @view y[4:6]
    #print(convert(Vector{Float64}, Geogrid(r, jd, eop_IAU1980)))
    rsu = @view y[7:9]
    vsu = @view y[10:12]
    rmoo = @view y[13:15]
    vmoo = @view y[16:18]

    dr = @view dy[1:3] 
    dv = @view dy[4:6] 
    drs = @view dy[7:9]
    dvs = @view dy[10:12] 
    drm = @view dy[13:15]
    dvm = @view dy[16:18] 

    r_norm = norm(r)
    #gg = Geogrid(r, jd, eop_IAU1980)
    #print(gg)
    if (pos == 1)
        rsn, rmn = norm(rsu), norm(rmoo)
        rmMB, rsMB = rmoo, rsu
    else
        rmMB = rmoon(jd)
        rsMB = rsun(jd)
        rsn, rmn = norm(rsMB), norm(rmMB)
    end

    #r3ms = rsu - rmoo
    rsat3m = rmMB - r
    rsat3s = rsMB -r

    #Third body accelerations
    #a3m = PropagatorConstants.μ_m * (rsat3m/norm(rsat3m)^3 - rmoon/rmn^3)
    #a3s = PropagatorConstants.μ_s * (rsat3s/norm(rsat3s)^3 - rsun/rsn^3)
    #a3m = μm * (rsat3m/norm(rsat3m)^3 - rmMB/norm(rmMB)^3) #CK
    #a3s = μs * (rsat3s/norm(rsat3s)^3 - rsMB/norm(rsMB)^3) #CK
    a3m = μm * (rsat3m/norm(rsat3m)^3 - rmMB/rmn^3) #CK
    a3s = μs * (rsat3s/norm(rsat3s)^3 - rsMB/rsn^3) 
    #a3ms = μs * (r3ms/norm(r3ms)^3 - rsMB/rsn^3)
    #Calculate J2 influence

    a_J2 = 3/2*J2*μe*Re^2/r_norm^4/r_norm*[r[1]*(5*r[3]^2/r_norm^2-1), r[2]*(5*r[3]^2/r_norm^2-1), r[3]*(5*r[3]^2/r_norm^2-3)]
    #Calculate J3 influence
    d = 1/(2*r_norm^7);
    dd = -5*J3*μe*Re^3;
    ddd = 3*r[3]-7*r[3]^3/r_norm^2;

    a_J3 = d*dd*[r[1]*ddd,r[2]*ddd,6*r[3]^2 - 7*r[3]^4/r_norm^2 - 3*r_norm^2/5]; #checked

    #Calculate J4 influence
    g = 1/(8*r_norm^7);
    gg = 15*J4*μe*Re^4;
    ggg = 1 - 14*r[3]^2/r_norm^2 + 21*r[3]^4/r_norm^4;

    a_J4 = g*gg*[r[1]*ggg,r[2]*ggg,r[3]*(5 - 70*r[3]^2/(3*r_norm^2) + 21*r[3]^4/r_norm^4)];
    
    #Calculate J22 influence
    #aJ22 = J22(r, jd, eop_IAU1980)=#
    altitude = r_norm - Re

    #Exponential atmospheric model as per Vallado
    @inbounds rho =  
    altitude ≤ 200 ? rho0[12]*exp(-(altitude-180)/29.740) :
    altitude ≤ 250 ? rho0[11]*exp(-(altitude-200)/37.105) :
    altitude ≤ 300 ? rho0[10]*exp(-(altitude-250)/45.546) :
    altitude ≤ 350 ? rho0[9]*exp(-(altitude-300)/53.628) :
    altitude ≤ 400 ? rho0[8]*exp(-(altitude-350)/53.298) :
    altitude ≤ 450 ? rho0[1]*exp(-(altitude-400)/58.515) :
    altitude ≤ 500 ? rho0[2]*exp(-(altitude-450)/60.828) :
    altitude ≤ 600 ? rho0[3]*exp(-(altitude-500)/63.822) :
    altitude ≤ 700 ? rho0[4]*exp(-(altitude-600)/71.835) :
    altitude ≤ 800 ? rho0[5]*exp(-(altitude-700)/88.667) :
    altitude ≤ 900 ? rho0[6]*exp(-(altitude-800)/124.64) :
    altitude ≤ 1000 ? rho0[7]*exp(-(altitude-900)/181.05) : rho0[8]*exp(-(altitude-1000)/268) #CK
 
    #Obtain relative velocity to Earth's atmosphere
    vr  = [v[1] + ωe*r[2], v[2]- ωe*r[1], v[3]] #CK
    v_r = norm(vr)

    #Calculate drag
    a_d = -10^3*0.5*rho*Cd*A/m * v_r .* vr#v_r^2 * vr/v_r;
    #gg = Geogrid(r, jd)
    #g = gg
    #print(gg)
    #Broadcast to existing arrays
    @. dr = v
    #print(dr)
    @. dv = -r.*μe./r_norm^3 .+ a_J2 .+ a_d .+ a3m .+ a3s .+ asrp(rsu, r) .+ a_J3 .+ a_J4 #gg# + gg ##g## #a_d+ a3m + a3s+ asrp(rsMB, r)# #a3m + a3s+ a_J2 #+ aJ22 #@+ aJ22 #+ asrp(rsun, r) #+ a_J4# ++ asrp(rsun, r)  a_J3+ asrp(rsun, r)  #A+ a_d+dd aJ22 when ECEF to ECI solved
    @. drs = vsu
    @. dvs = -rsu.*μs./rsn^3
    @. drm = vmoo
    @. dvm = -rmoo.*μe./rmn^3 
    
end

#Starting conditions
#=
y0 = [1315.0823093944996, -6765.511976126983, 0.0,-0.9692300284946672, -0.1883992325643531, 7.540501269182395, 
1.3219e8,-0.6101e8, -0.2645e8, 13.8614, 24.5178 ,  10.6295, 
2.6183e5, 2.6184e5, 0.757e5, -0.7733, 0.6312, 0.2803] =#


y0 = [-2.402864124972275e+03,8.538336244770163e+02,-6.682344154832135e+03,-4.951327157415951,5.020437431604642,2.422461551425946,
1.463957088922231e+08, -2.335131546485938e+07, -1.012403398015358e+07, 5.576319298565035, 27.026588503272446, 11.716155914748752,
-1.170981644822336e+05, 3.683752523717154e+05, 1.208032774797189e+05, -0.927440482124745, -0.241519315834115, -0.125577393037498]

tspan = (0.0, 24*3600*455)
prob = ODEProblem(KepCowell!,y0,tspan, p)
sol = solve(prob, VCABM(),  reltol=1e-13, abstol=1e-13, maxiters = 1e8)#, save_everystep = false)#, save_everystep = false)#save_at = 30*60) #save_everystep = false


