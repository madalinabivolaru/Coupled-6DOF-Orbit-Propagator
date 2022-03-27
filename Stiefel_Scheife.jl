using LinearAlgebra
using DifferentialEquations
using OrdinaryDiffEq
using Plots
using BenchmarkTools
#using StaticArrays


function Stiefel(dy::Float64,y::Float64,p::Float64,t::Float64)
    μ = 398601.0; 
    R_e = 6371.22;
    J2 =  0.00108265;
    r = y[1:3]
    v = y[4:6]
    v_norm, r_norm= norm(v), norm(r) 

    c = 2*r_norm^5
    cc = 5*r[3]^2/r_norm^2 
    ccc= -3*J2*μ*R_e^2
    
    #J2 acceleration
    ax = ccc*r[1] / (c) * (1-cc)
    ay = ccc*r[2] / (c) * (1-cc)
    az = ccc*r[3] / (c) * (3-cc)
    a_J2 = [ax, ay, az];

    #Introducing Lunar influence by Stiefel & Scheifele
    rm = 384400.0;
    Ωm = 2.665315780887 * 10^-6;
    rmoon = rm.*[sin(Ωm*t), -0.5*cos(Ωm*t)*sqrt(3),-0.5*cos(Ωm*t)];
   

    μm = 4902.66;
    rs3m = rmoon - r;
    a3m = μm * (rs3m/norm(rs3m)^3 - rmoon/norm(rmoon)^3);

    dy[1:3] = y[4:6];
    dy[4:6]= -r*μ/r_norm^3 + a_J2 + a3m ;
    
end

#Start conditions for the test
y0 = [0.0, -5888.9727, -3400.0, 10.691338, 0.0, 0.0];
tspan = (0.00,24894232.3650240);

#Define and solve ODE
prob = ODEProblem(Stiefel, y0, tspan);
toler = 10^-13
sol = solve(prob, DP8(), reltol = toler, abstol = toler)

