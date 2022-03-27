using LinearAlgebra
using DifferentialEquations
using Plots

 R_e = 6378.137; #km
 μ = 398600.4418;#km^3/s^2

function Cowell!(dy,y,p,t)
    #Propagator with accelerations due to J2, J3 and J4
    μ = 398600.4418; 
    R_e = 6378.137;
    J2 = 0.00108262617;
    J3 = -2.53241052e-06;
    J4 = -1.6198976e-06;
    

    r = [y[1],y[2],y[3]]
    v = [y[4],y[5],y[6]]
    v_norm = norm(v)
    r_norm = norm(r)

    ax_J2 = -3*J2*μ*R_e^2*r[1] / (2*r_norm^5) * (1-5*r[3]^2/r_norm^2)
    ay_J2 = -3*J2*μ*R_e^2*r[2] / (2*r_norm^5) * (1-5*r[3]^2/r_norm^2)
    az_J2 = -3*J2*μ*R_e^2*r[3] / (2*r_norm^5) * (3-5*r[3]^2/r_norm^2)
    a_J2 = [ax_J2,ay_J2,az_J2];

    ax_J3 = -5*J3*μ*R_e^3*r[1] / (2*r_norm^7) * (3*r[3]-7*r[3]^3/r_norm^2)
    ay_J3 = -5*J3*μ*R_e^3*r[2] / (2*r_norm^7) * (3*r[3]-7*r[3]^3/r_norm^2)
    az_J3 = -5*J3*μ*R_e^3 / (2*r_norm^7) * (6*r[3]^2-7*r[3]^4/r_norm^2-3*r_norm^2/5)
    a_J3 = [ax_J3,ay_J3,az_J3];

    ax_J4 = 15*J4*μ*R_e^4*r[1] / (8*r_norm^7) *(1 - 14*r[3]^2/r_norm^2 + 21*r[3]^4/r_norm^4)
    ay_J4 = 15*J4*μ*R_e^4*r[2] / (8*r_norm^7) *(1 - 14*r[3]^2/r_norm^2 + 21*r[3]^4/r_norm^4)
    az_J4 = 15*J4*μ*R_e^4*r[3] / (8*r_norm^7) *(5 - 70*r[3]^2/(3*r_norm^2) + 21*r[3]^4/r_norm^4)
    a_J4 = [ax_J4,ay_J4,az_J4]

    dy[1:3] = y[4:6]
    dy[4:6] = -r*μ/r_norm^3 + a_J2 + a_J3 + a_J4
end

y0 = [ 1131.340, - 2282.343, 6672.423,- 5.64305, 4.30333, 2.42879]
tspan = (0.0,604800.0)

prob = ODEProblem(Cowell!,y0,tspan)
sol = solve(prob, Tsit5(),reltol=1e-14, abstol=1e-14) 
