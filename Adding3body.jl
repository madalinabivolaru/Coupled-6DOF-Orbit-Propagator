using LinearAlgebra
using DifferentialEquations
using Plots

const R_e = 6378.137; #km
const μ = 398600.4418;#km^3/s^2

function JulDat(yr,mo,d,h,min,s)
    JD = 367*yr - floor(7/4*(yr+floor((mo+9)/12)))+
    floor(275*mo/9)+d+1721013.5+((s/60 + min)/60 + h)/24;
    return JD
end
function rmoon(JD)
    Ttdb = (JD-2451545.0)/36525
    s(x) = sind(x%360)
    c(x) = cosd(x%360)
    λecl = (218.32 + 481267.8813 * Ttdb + 6.29 * s(134.9 + 477198.85*Ttdb)- 
    1.27*s(259.2 - 413335.38*Ttdb) + 0.66*s(235.7 + 890534.23* Ttdb)+ 
    0.21*s(269.9 + 954397.70*Ttdb) - 0.19*s(357.5 + 35999.05 * Ttdb)-
    0.11*s(186.6 + 966404.05*Ttdb)) % 360 + 360

    Φecl = 5.13*s(93.3 + 483202.03 * Ttdb) + 0.28*s(228.2 + 960400.87 * Ttdb)-
    0.28*s(318.3 + 6003.18 * Ttdb) - 0.17*s(217.6 - 407332.20 *Ttdb)
    P = 0.9508 + 0.0518 * c(134.9 + 477198.85 * Ttdb) +
    0.0095 * c(259.2 - 413335.38 * Ttdb) + 0.0078 * c(235.7 + 890534.23 * Ttdb) + 
    0.0028 * c(269.9 + 954397.70 * Ttdb)
    Ε = 23.439291 - 0.0130042 * Ttdb
    r_moon = 6378.137/sind(P) 
    R_moon = r_moon * [cosd(Φecl)*cosd(λecl); 
    cosd(Ε)*cosd(Φecl)*sind(λecl) - sind(Ε)*sind(Φecl);
    sind(Ε)*cosd(Φecl)*sind(λecl) + cosd(Ε)*sind(Φecl)]
    return transpose(R_moon)[:]
end
function rsun(JD)
    Tut1 = (JD-2451545.0)/36525
    AU = 149597870.7
    λM = (280.460 + 36000.771*Tut1)%360
    Ttdb = Tut1
    M = (357.5291092 + 35999.05034*Ttdb)%360
    λecliptic = (λM +1.914666471*sind(M) + 0.019994643 * sind(2*M))%360
    r = 1.000140612 - 0.016708617*cosd(M) - 0.0001395892 * cosd(2*M)
    ϵ = 23.439291 - 0.0130042*Ttdb
    R =  AU * [r*cosd(λecliptic); r*cosd(ϵ)*sind(λecliptic);r*sind(ϵ)*sind(λecliptic)]
    return transpose(R)[:]
end
function Rdbody(rsat,rmoon,rsun)
    μ = 398600.4418;
    μm = 4904.8695;
    μs = 1.32712440018*10^11
    rsat3m = rmoon - rsat
    rsat3s = rsun - rsat
    a3m = μm * (rsat3m/norm(rsat3m)^3 - rmoon/norm(rmoon)^3)
    a3s = μs * (rsat3s/norm(rsat3s)^3 - rsun/norm(rsun)^3)
    return a3m + a3s
end
function asrp(JD, rsat)
    p_srp = 4.57 * 10^-6;   #next you should determine Daphelion
    Cr = 1.2;
    As = 0.5*0.5;
    m = 100;
    rsats = rsun(JD) - rsat
    asrp = - (p_srp*Cr*As/m)*rsats/norm(rsats)
    return asrp
end

function Cowell!(dy,y,p,t)
    μ = 398600.4418; 
    R_e = 6378.137;
    J2 =  0.0010826269;
    r = [y[1],y[2],y[3]]
    v = [y[4],y[5],y[6]]
    v_norm = norm(v)
    r_norm = norm(r) 
    JD = JulDat(2021,11,20,0,0,t) 

    ax = -3*J2*μ*R_e^2*r[1] / (2*r_norm^5) * (1-5*r[3]^2/r_norm^2)
    ay = -3*J2*μ*R_e^2*r[2] / (2*r_norm^5) * (1-5*r[3]^2/r_norm^2)
    az = -3*J2*μ*R_e^2*r[3] / (2*r_norm^5) * (3-5*r[3]^2/r_norm^2)
    a_J2 = [ax,ay,az];

    dy[1], dy[2], dy[3] = y[4], y[5], y[6]
    dy[4], dy[5], dy[6] = -r*μ/r_norm^3 + a_J2 + Rdbody(r,rmoon(JD),rsun(JD)) + asrp(JD,r)
end

y0 = [ 1131.340, - 2282.343, 6672.423,- 5.64305, 4.30333, 2.42879]
tspan = (0.0,2500.0)
prob = ODEProblem(Cowell!,y0,tspan)
sol = solve(prob, Tsit5(),reltol=1e-8, abstol=1e-8) 