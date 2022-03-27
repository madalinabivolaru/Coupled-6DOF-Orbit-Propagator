using LinearAlgebra
using DifferentialEquations
using Plots
using BenchmarkTools

const R_e = 6378.137; #km
const μ = 398600.4418;#km^3/s^2

function Cowell!(dy,y,p,t)
    #Propagator with perturbation fue to J2 influence
    μ = 398600.4418; 
    R_e = 6378.137;
    J2 =  0.0010826269;
    r = y[1:3]
    v = y[4:6]

    v_norm = norm(v)
    r_norm = norm(r) 
   
    c = 2*r_norm^5
    cc = 5*r[3]^2/r_norm^2 
    ccc= -3*J2*μ*R_e^2

    ax = ccc*r[1] / (c) * (1-cc)
    ay = ccc*r[2] / (c) * (1-cc)
    az = ccc*r[3] / (c) * (3-cc)
    a_J2 = [ax,ay,az];

    dy[1:3] = y[4:6]
    dy[4:6] = -r*μ/r_norm^3 + a_J2
end



y0 = [ 1131.340, - 2282.343, 6672.423,- 5.64305, 4.30333, 2.42879]
tspan = (0.0,2480942.0)

prob = ODEProblem(Cowell!,y0,tspan)
sol = solve(prob, DP8(), reltol=1e-10, abstol=1e-10, save_everystep = false) 

 