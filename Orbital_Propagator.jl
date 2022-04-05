using LinearAlgebra
using OrdinaryDiffEq
using DifferentialEquations
using StaticArrays
#using Plots
using BenchmarkTools
include("./PropagatorConstants.jl")
using .PropagatorConstants
#include("./JulDat.jl")
#include("asrp.jl")
#include("RVtoCOE.jl")
using DataFrames
using CSV


function KepCowell!(dy::Vector{Float64},y::Vector{Float64},p,t)

    A = 0.7;
    m = 400;

    #Reference densities for LEO
    rho0 = [3.735*10^-12, 1.585*10^-12,6.967*10^-13, 1.454*10^-13,
    3.614*10^-14,1.17*10^-14, 5.245*10^-15,3.019*10^-15,
    9.518*10^-12, 2.418*10^-11, 7.248*10^-11, 2.789*10^-10, 5.464*10^-10 ]

    #Access data directly from y, dy vector with view macro
    r = @view y[1:3]
    v = @view y[4:6]
    rsun = @view y[7:9]
    vsun = @view y[10:12]
    rmoon = @view y[13:15]
    vmoon = @view y[16:18]

    dr = @view dy[1:3] 
    dv = @view dy[4:6] 
    drs = @view dy[7:9]
    dvs = @view dy[10:12] 
    drm = @view dy[13:15]
    dvm = @view dy[16:18] 


    r_norm, rsn, rmn = norm(r), norm(rsun), norm(rmoon)

    #Obtain relative vectors 
    rsat3m = rmoon - r 
    rsat3s = rsun - r 
    #rm3s = rsun - rmoon 

    #Third body accelerations
    a3m = PropagatorConstants.μ_m * (rsat3m/norm(rsat3m)^3 - rmoon/rmn^3)
    a3s = PropagatorConstants.μ_s * (rsat3s/norm(rsat3s)^3 - rsun/rsn^3)

    #Calculate J2 influence
    c = 1/(2*r_norm^5);
    cc = 5*r[3]^2/r_norm^2 ;
    ccc= -3*PropagatorConstants.J2*PropagatorConstants.μ_e*PropagatorConstants.R_e^2;

    a_J2 = ccc*c*[r[1] * (1-cc), r[2] * (1-cc), r[3]* (3-cc)]; 

    #Calculate J3 influence
    d = 1/(2*r_norm^7);
    dd = -5*PropagatorConstants.J3*PropagatorConstants.μ_e*PropagatorConstants.R_e^3;
    ddd = 3*r[3]-7*r[3]^3/r_norm^2;

    a_J3 = d*dd*[r[1]*ddd,r[2]*ddd,6*r[3]^2 - 7*r[3]^4/r_norm^2 - 3*r_norm^2/5]; #checked

    #Calculate J4 influence
    g = 1/(8*r_norm^7);
    gg = 15*PropagatorConstants.J4*PropagatorConstants.μ_e*PropagatorConstants.R_e^4;
    ggg = 1 - 14*r[3]^2/r_norm^2 + 21*r[3]^4/r_norm^4;

    a_J4 = g*gg*[r[1]*ggg,r[2]*ggg,r[3]*(5 - 70*r[3]^2/(3*r_norm^2) + 21*r[3]^4/r_norm^4)];

    #Calculate J22 influence
    #aJ22 = J22(r)
    altitude = r_norm - PropagatorConstants.R_e

    #Exponential atmospheric model as per Vallado
    rho =  
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
    altitude ≤ 1000 ? rho0[7]*exp(-(altitude-900)/181.05) : rho0[8]*exp(-(altitude-1000)/268)
 
    #Obtain relative velocity to Earth's atmosphere
    vr  = [v[1]+PropagatorConstants.ω_e*r[2], v[2]-PropagatorConstants.ω_e*r[1], v[3]]
    v_r = norm(vr)

    #Calculate drag
    a_d = -10^3*0.5*rho*PropagatorConstants.C_d*A/m * v_r^2 * vr/v_r;

    #Broadcast to existing arrays
    @. dr = v
    @. dv = -r*PropagatorConstants.μ_e/r_norm^3 + a3m + a3s + a_d + a_J2 + a_J3 + a_J4 + asrp(rsun, r)  #Add aJ22 when ECEF to ECI solved
    @. drs = vsun
    @. dvs = -rsun*PropagatorConstants.μ_s/rsn^3
    @. drm = vmoon
    @. dvm = -rmoon*PropagatorConstants.μ_e/rmn^3 
    
end

#Starting conditions
y0 = [1315.0823093944996, -6765.511976126983, 0.0,-0.9692300284946672, -0.1883992325643531, 7.540501269182395, 
1.3219e8,-0.6101e8, -0.2645e8, 13.8614, 24.5178 ,  10.6295, 
2.6183e5, 2.6184e5, 0.757e5, -0.7733, 0.6312, 0.2803] 


tspan = (0.0, 24*3600*365)
prob = ODEProblem(KepCowell!,y0,tspan)
sol = solve(prob, VCABM(), reltol=1e-13, abstol=1e-13, save_everystep = false , maxiters = 1e8) #save_everystep = false

#plot(sol,vars=(1,2,3))

