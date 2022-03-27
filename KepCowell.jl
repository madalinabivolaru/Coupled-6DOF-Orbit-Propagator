using LinearAlgebra
using OrdinaryDiffEq
using DifferentialEquations
#using Plots
using BenchmarkTools
include("./PropagatorConstants.jl")
using .PropagatorConstants
include("./JulDat.jl")
include("asrp.jl")
include("RVtoCOE.jl")
using DataFrames
using CSV

function KepCowell!(dy::Vector{Float64},y::Vector{Float64},p,t)
    #A = 0.52^2;
    #m = 100;
    A = 0.7;
    m = 400;
    rho0 = [3.735*10^-12, 1.585*10^-12,6.967*10^-13, 1.454*10^-13,
    3.614*10^-14,1.17*10^-14, 5.245*10^-15,3.019*10^-15]
    r = y[1:3]
    v = y[4:6]
    rsun = y[7:9]
    vsun = y[10:12]
    rmoon = y[13:15]
    vmoon = y[16:18]
    v_norm, r_norm, rsn, rmn = norm(v),norm(r), norm(rsun), norm(rmoon)

    rsat3m = rmoon - r #correct orientation
    rsat3s = rsun - r #correct orientation
    rm3s = rsun - rmoon #correct orientation
    a3m = PropagatorConstants.μ_m * (rsat3m/norm(rsat3m)^3 - rmoon/rmn^3)
    a3s = PropagatorConstants.μ_s * (rsat3s/norm(rsat3s)^3 - rsun/rsn^3)
    am3s = PropagatorConstants.μ_s * (rm3s/norm(rm3s)^3 - rsun/rsn^3)

    c = 1/(2*r_norm^5);
    cc = 5*r[3]^2/r_norm^2 ;
    ccc= -3*PropagatorConstants.J2*PropagatorConstants.μ_e*PropagatorConstants.R_e^2;

    a_J2 = ccc*c*[r[1] * (1-cc),r[2] * (1-cc),r[3]* (3-cc)];

    d = 1/(2*r_norm^7);
    dd = -5*PropagatorConstants.J3*PropagatorConstants.μ_e*PropagatorConstants.R_e^3;
    ddd = 3*r[3]-7*r[3]^3/r_norm^2;


    a_J3 = d*dd*[r[1]*ddd,r[2]*ddd,6*r[3]^2-7*r[3]^4/r_norm^2-3*r_norm^2/5];

    g = 1/(8*r_norm^7);
    gg = 15*PropagatorConstants.J4*PropagatorConstants.μ_e*PropagatorConstants.R_e^4;
    ggg = 1 - 14*r[3]^2/r_norm^2 + 21*r[3]^4/r_norm^4;


    a_J4 = g*gg*[r[1]*ggg,r[2]*ggg,r[3]*(5 - 70*r[3]^2/(3*r_norm^2) + 21*r[3]^4/r_norm^4)];

    x22 = -3*PropagatorConstants.R_e^2*PropagatorConstants.μ_e/2/r_norm^7*
    (PropagatorConstants.C22*(r_norm^2*r[1]+5*r[1]^3-15*r[1]*r[2]^2-5*r[1]*r[3]^2)+
    PropagatorConstants.S22*(r_norm^2*r[2]+15*r[1]^2*r[2]-5*r[2]^3-5*r[2]*r[3]^2));
    y22 = 3*PropagatorConstants.R_e^2*PropagatorConstants.μ_e/2/r_norm^7*
    (PropagatorConstants.C22*(r_norm^2*r[2]-15*r[1]^2*r[2]+5*r[2]^3-5*r[2]*r[3]^2)+
    PropagatorConstants.S22*(-r_norm^2*r[1]+5*r[1]^3-15*r[1]*r[2]^2+5*r[1]*r[3]^2));
    z22 = -15*PropagatorConstants.R_e^2*PropagatorConstants.μ_e*r[3]/r_norm^7*
    (PropagatorConstants.C22*(r[1]^2-r[2]^2)+PropagatorConstants.S22*r[1]*r[2]);
    aJ22 = [x22,y22,z22];
    rr = r_norm - PropagatorConstants.R_e

    rho =  rr ≤ 450 ? rho0[1]*exp(-(rr-400)/58.515) :
    #400 < rr ≤ 450 ? rho0[1]*exp(-(rr-400)/58.515) :
    450 < rr ≤ 500 ? rho0[2]*exp(-(rr-450)/60.828) :
    500 < rr ≤ 600 ? rho0[3]*exp(-(rr-500)/63.822) :
    600 < rr ≤ 700 ? rho0[4]*exp(-(rr-600)/71.835) :
    700 < rr ≤ 800 ? rho0[5]*exp(-(rr-700)/88.667) :
    800 < rr ≤ 900 ? rho0[6]*exp(-(rr-800)/124.64) :
    900 < rr ≤ 1000 ? rho0[7]*exp(-(rr-900)/181.05) : rho0[8]*exp(-(rr-1000)/268)

    vr = v - cross([0,0,PropagatorConstants.ω_e],r)
    v_r = norm(vr)
    Drag = -10^3*0.5*rho*PropagatorConstants.C_d*A/m * v_r^2 * vr/v_r;
    #Drag = -1000* 0.5*rho*PropagatorConstants.C_d*A/m*v_norm^2*v/v_norm;
    #println(Drag)
    dy[1:3] = y[4:6] 
    dy[4:6] = -r*PropagatorConstants.μ_e/r_norm^3 + a3m + a3s + a_J2  + Drag + aJ22+ a_J3 + a_J4 #+ asrp(rsun,r,false)
    dy[7:9] = y[10:12]
    dy[10:12] = -rsun*PropagatorConstants.μ_s/rsn^3
    dy[13:15] = y[10:12]
    dy[16:18] = -rmoon*PropagatorConstants.μ_e/rmn^3 + am3s
    
end
y0 = [1315.0823093944996, -6765.511976126983, 0.0,-0.9692300284946672, -0.1883992325643531, 7.540501269182395, 1.3219e8,-0.6101e8, -0.2645e8, 13.8614, 24.5178 ,  10.6295,     0.1583e5,    3.6069e5,   1.6382e5,  -0.9845,   -0.0440,    0.0730]
#y0 = [ 1131.340, - 2282.343, 6672.423,- 5.64305, 4.30333, 2.42879,-7.926e7,-1.1448e8,-4.963e7, 25.6188,-14.5466, -6.3048,1.757e5,3.3491e5,1.473e5,-0.8700, 0.3558, 0.2520]
tspan = (0.0,24*3600*365)
prob = ODEProblem(KepCowell!,y0,tspan)
@benchmark sol = solve(prob, VCABM(), reltol=1e-14, abstol=1e-14,  save_everystep=false)#saveat = 36000)

#f = DataFrame(sol)
#CSV.write("all9.csv",f[:,2:7])
#=
function Anal(sol)
    len = 1000
    a = ones(1,len)
    e = ones(1,len)
    inc = ones(1,len)
    Om = ones(1,len)
    om = ones(1,len)
    for i in 1:2:4
        a[i] = RVtoCOE(sol(i)[1:3],sol(i)[4:6])[1]
        e[i] = RVtoCOE(sol(i)[1:3],sol(i)[4:6])[2]
        inc[i] RVtoCOE(sol(i)[1:3],sol(i)[4:6])[3]
        Om[i] = RVtoCOE(sol(i)[1:3],sol(i)[4:6])[4]
        om[i] = RVtoCOE(sol(i)[1:3],sol(i)[4:6])[5] 
    #println(RVtoCOE(sol(i)[1:3],sol(i)[4:6]))
    end
end
=#