#First attempt at a coupled orbit and attitude PropagatorConstants
#Third body was eliminated on account of complexity
using LinearAlgebra
using OrdinaryDiffEq
using DifferentialEquations
using Plots
using BenchmarkTools
include("./PropagatorConstants.jl")
using .PropagatorConstants
include("./JulDat.jl")
include("asrp.jl")
include("RVtoCOE.jl")
using DataFrames
using CSV

function A2(r, v)
    #Function which defines rotatiuon matrix A2 as ECI = A2*NTW
    T = v/(norm(v))
    W = cross(r,v)/norm(cross(r,v))
    N = cross(T,W)
    #Get att matrix
    A2 = [N, T, W]

    return inv(reduce(vcat,transpose.(A2)))
end

function Beta2Cos(q)
    #function which returns A1 for quaternions propagated; Here A1 is defines as NTW = A1*Body
    Cfin = [q[1]^2+q[2]^2-q[3]^2-q[4]^2 2*(q[2]*q[3]+q[1]*q[4]) 2*(q[2]*q[4]-q[1]*q[3]);
            2*(q[2]*q[3]-q[1]*q[4]) q[1]^2-q[2]^2+q[3]^2-q[4]^2 2*(q[3]*q[4]+q[1]*q[2]);
            2*(q[2]*q[4]+q[1]*q[3]) 2*(q[3]*q[4]-q[1]*q[2]) q[1]^2-q[2]^2-q[3]^2+q[4]^2]; #checked
    #if abs(Cfin[3,3] - 1)<0.0001
        #Cfin[3,3] = 1
    #end
end

function KepCowell!(dy::Vector{Float64},y::Vector{Float64},p,t)
    #Define constants, body vectors and area of panels
    Lsp = 1.32
    A = [1.11,  1,  1,  1  ];
    nx = [1,0,0]
    nz = [0,0,1]
    ny = [0,1,0]

    rpanel = Lsp*ny
    ry = 0.5*ny
    rx = 0.5*nx
    rz = 0.5*nz
    m = 100;
    #I = diagm([42.24,104.93, 105.82])
    I = diagm([104.93, 42.24, 105.82])
    rho0 = [3.735*10^-12, 1.585*10^-12,6.967*10^-13, 1.454*10^-13,
    3.614*10^-14,1.17*10^-14, 5.245*10^-15,3.019*10^-15,
    9.518*10^-12, 2.418*10^-11, 7.248*10^-11, 2.789*10^-10, 5.464*10^-10 ]
    SF = 1367
    p_srp = SF/3*10^-8;
    Cr = 1.2;

    r = y[1:3]
    v = y[4:6]
    rsun = y[7:9]
    vsun = y[10:12]
    rmoon = y[13:15]
    vmoon = y[16:18]    
    q = y[19:22]
    ω = y[23:25]

    dr = dy[1:3] 
    dv = dy[4:6] 
    drs = dy[7:9]
    dvs = dy[10:12] 
    drm = dy[13:15]
    dvm = dy[16:18] 
    dq = dy[19:22]
    dω = dy[23:25]

    #Get total Attitude matrix ECI in Body
    At2 = A2(r,v)
    At1 = Beta2Cos(q)


    #v_norm, r_norm, rsn, rmn = norm(v), norm(r), norm(rsun), norm(rmoon)
    r_norm, rsn, rmn = norm(r), norm(rsun), norm(rmoon)

    rsat3m = rmoon - r #correct orientation
    rsat3s = rsun - r #correct orientation
    rm3s = rsun - rmoon #correct orientation
    #Qm = (r_norm^2+2*dot(r,rsat3m))*(rmn^2+rmn*norm(rsat3m)+norm(rsat3m)^2)/(rmn^3*norm(rsat3m)^3)./(rmoon+rsat3m)
    #a3m =(PropagatorConstants.μ_m*(rsat3m.*Qm-r/rmn^3))

    #Qs = (r_norm^2+2*dot(r,rsat3s))*(rsn^2+rsn*norm(rsat3s)+norm(rsat3s)^2)/(rsn^3*norm(rsat3s)^3)./(rsun+rsat3s)
    #a3s = (PropagatorConstants.μ_s*(rsat3s.*Qs-r/rsn^3))
    a3m = PropagatorConstants.μ_m * (rsat3m/norm(rsat3m)^3 - rmoon/rmn^3)
    a3s = PropagatorConstants.μ_s * (rsat3s/norm(rsat3s)^3 - rsun/rsn^3)
    #Q = (rmn^2+2*dot(rmoon,rm3s))*(rsn^2+rsn*norm(rm3s)+norm(rm3s)^2)/(rsn^3*norm(rm3s)^3)./(rsun+rm3s)
    #am3s= rm3s.*Q-rmoon/rmn^3
    #am3s = PropagatorConstants.μ_s * (rm3s/norm(rm3s)^3 - rsun/rsn^3) previous version

    c = 1/(2*r_norm^5);
    cc = 5*r[3]^2/r_norm^2 ;
    ccc= -3*PropagatorConstants.J2*PropagatorConstants.μ_e*PropagatorConstants.R_e^2;

    a_J2 = ccc*c*[r[1] * (1-cc), r[2] * (1-cc), r[3]* (3-cc)]; #checked

    d = 1/(2*r_norm^7);
    dd = -5*PropagatorConstants.J3*PropagatorConstants.μ_e*PropagatorConstants.R_e^3;
    ddd = 3*r[3]-7*r[3]^3/r_norm^2;


    a_J3 = d*dd*[r[1]*ddd,r[2]*ddd,6*r[3]^2 - 7*r[3]^4/r_norm^2 - 3*r_norm^2/5]; #checked

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

    rho =  
    rr ≤ 200 ? rho0[12]*exp(-(rr-180)/29.740) :
    rr ≤ 250 ? rho0[11]*exp(-(rr-200)/37.105) :
    rr ≤ 300 ? rho0[10]*exp(-(rr-250)/45.546) :
    rr ≤ 350 ? rho0[9]*exp(-(rr-300)/53.628) :
    rr ≤ 400 ? rho0[8]*exp(-(rr-350)/53.298) :
    rr ≤ 450 ? rho0[1]*exp(-(rr-400)/58.515) :
    rr ≤ 500 ? rho0[2]*exp(-(rr-450)/60.828) :
    rr ≤ 600 ? rho0[3]*exp(-(rr-500)/63.822) :
    rr ≤ 700 ? rho0[4]*exp(-(rr-600)/71.835) :
    rr ≤ 800 ? rho0[5]*exp(-(rr-700)/88.667) :
    rr ≤ 900 ? rho0[6]*exp(-(rr-800)/124.64) :
    rr ≤ 1000 ? rho0[7]*exp(-(rr-900)/181.05) : rho0[8]*exp(-(rr-1000)/268)
    #rho = expatmosphere(rr*1000)
    #vr = v - cross([0,0,PropagatorConstants.ω_e],r)
    vr  = [v[1]+PropagatorConstants.ω_e*r[2], v[2]-PropagatorConstants.ω_e*r[1], v[3]]
    v_r = norm(vr)

    #Perform rotation from ECI to Body to get Drag torque
    ECItoBody = inv(At2*At1)
    vr_B = ECItoBody*vr
    vs_B = ECItoBody*(-rsat3s)/norm(rsat3s)
    cosSx = dot(nx,vs_B)
    cosSy = dot(ny,vs_B)
    cosSz = dot(nz,vs_B)


    cosTx = dot(nx,vr_B)/v_r
    cosTy = dot(ny,vr_B)/v_r
    cosTz = dot(nz,vr_B)/v_r

    #print(cosTx)
    #Drag coupled with attitude
    #Change As of asrp
    #asrp = -10^(-3)*(p_srp*Cr*As/m)*rsun/norm(rsun)^3*PropagatorConstants.AUtoKM^2
    Drag = -10^3*0.5*rho*PropagatorConstants.C_d*(2*A[1]*max(cosTx,0)+2*A[1]*max(-cosTx,0) + A[3]*max(cosTx,0) + A[3]*max(-cosTx,0) +
            A[2]*max(cosTy,0) + A[2]*max(-cosTy,0) + A[4]*max(cosTz,0) + A[4]*max(-cosTz,0))/m * v_r^2 * vr/v_r;

    #Fsrp = -10^(-3)*p_srp*Cr*vs_B/norm(vs_B)*(4*max(cosSx,0)+2*max(cosSy,0)+2*max(cosSx,0))
    #asrp = inv(ECItoBody)*(Fsrp/m)

    rsats = r-rsun
    asrp = 10^(-3)*(p_srp*Cr*(2*A[1]*max(cosSx,0) + 2*A[1]*max(-cosSx,0) + A[3]*max(cosSx,0) + A[3]*max(-cosSx,0) +
    A[2]*max(cosSy,0) + A[2]*max(-cosSy,0) + A[4]*max(cosSz,0) + A[4]*max(-cosSz,0))/m)*rsats/norm(rsats)
    Fs = 10^(-3)*[(-p_srp*Cr*vs_B)*A[1]*max(cosSx,0), (-p_srp*Cr*vs_B)*A[1]*max(-cosSx,0), (-p_srp*Cr*vs_B)*A[3]*max(cosSx,0),
           (-p_srp*Cr*vs_B)*A[3]*max(-cosSx,0), (-p_srp*Cr*vs_B)*A[2]*max(cosSy,0), (-p_srp*Cr*vs_B)*A[2]*max(-cosSy,0),
           (-p_srp*Cr*vs_B)*A[4]*max(cosSz,0), (-p_srp*Cr*vs_B)*A[4]*max(-cosSz,0)]
    
           #Drag components for torques
    Fd = [-10^6*0.5*rho*PropagatorConstants.C_d*v_r*vr_B*A[1]*max(cosTx,0), -10^6*0.5*rho*PropagatorConstants.C_d*v_r*vr_B*A[1]*max(-cosTx,0),
           -10^6*0.5*rho*PropagatorConstants.C_d*v_r*vr_B*A[3]*max(cosTx,0), -10^6*0.5*rho*PropagatorConstants.C_d*v_r*vr_B*A[3]*max(-cosTx,0),
           -10^6*0.5*rho*PropagatorConstants.C_d*v_r*vr_B*A[2]*max(cosTy,0), -10^6*0.5*rho*PropagatorConstants.C_d*v_r*vr_B*A[2]*max(-cosTy,0),
           -10^6*0.5*rho*PropagatorConstants.C_d*v_r*vr_B*A[4]*max(cosTz,0), -10^6*0.5*rho*PropagatorConstants.C_d*v_r*vr_B*A[4]*max(-cosTz,0)]

    Ld = 1/m*(cross(rpanel, Fd[1])+cross(-rpanel, Fd[2])+cross(ry, Fd[3]) + cross(-ry, Fd[4]) + cross(rx, Fd[5]) + cross(-rx, Fd[6]) + cross(rz, Fd[7]) + cross(-rz, Fd[8]))
    #print(Ld)
    Lsrp = cross(rpanel, Fs[1])+cross(-rpanel, Fs[2])+cross(ry, Fs[3]) + cross(-ry, Fs[4]) + cross(rx, Fs[5]) + cross(-rx, Fs[6]) + cross(rz, Fs[7]) + cross(-rz, Fs[8])
    #print(Fs)
    #Attitude 
    ωx, ωy, ωz = ω[1], ω[2], ω[3]

    Ω = [0.0 -ωx -ωy -ωz;
        ωx 0 ωz -ωy;
        ωy -ωz 0 ωx;
        ωz ωy -ωx 0.0];
    #println((1.11*2*max(cosSx,0)+ max(cosSx,0) + (1.11*2+1)*max(cosSxx,0)+max(cosSy,0)+max(cosSyy,0)+max(cosSz,0)+max(cosSzz,0))/m)
    dy[1:3] = v
    dy[4:6] = -r*PropagatorConstants.μ_e/r_norm^3+ a3m+ a3s + Drag+ a_J2 + a_J3 + a_J4+ asrp   # #    #+ asrp(rsun,r)  + aJ22
    dy[7:9] = vsun
    dy[10:12] = -rsun*PropagatorConstants.μ_s/rsn^3
    dy[13:15] = vmoon
    dy[16:18] = -rmoon*PropagatorConstants.μ_e/rmn^3 #+ am3s
    #print(length(transpose(0.5 .* Ω* transpose(q)[:])[:]))
    #print(length(transpose(I^-1*(transpose(Ld)[:]-cross(transpose(ω)[:],I*transpose(ω)[:])))[:]))
    dy[19:22] = transpose(0.5 .* Ω* transpose(q)[:])[:]
    dy[23:25] = transpose(I^-1*(transpose(Ld+Lsrp)[:]-cross(transpose(ω)[:],I*transpose(ω)[:])))[:]
    
end
y0 = [1315.0823093944996, -6765.511976126983, 0.0,-0.9692300284946672, -0.1883992325643531, 7.540501269182395, 
1.3219e8,-0.6101e8, -0.2645e8, 13.8614, 24.5178 ,  10.6295, 
2.6183e5, 2.6184e5, 0.757e5, -0.7733, 0.6312, 0.2803,
1,0,0,0,
0.0,0.,0.] 
tspan = (0.0,100000)
prob = ODEProblem(KepCowell!,y0,tspan)
sol = solve(prob, VCABM(), reltol=1e-8, abstol=1e-8, saveat = 10)#,  save_everystep=false)#saveat = 36000)

tim = sol.t[:]
R = transpose(sol[1:3,:])[:,:]
V = transpose(sol[4:6,:])[:,:]
a = ones(length(tim),1)
e_mod = ones(length(tim),1)
i = ones(length(tim),1)
ν = ones(length(tim),1)
Ω  =ones(length(tim),1)
ω = ones(length(tim),1)

for j in 1:1:length(tim)
    a[j],e_mod[j],i[j],Ω[j],ω[j],ν[j] = RVtoCOE(R[j,:],V[j,:])
    #print(a[j],e_mod[j],i[j],Ω[j],ω[j],ν[j])
end 
plot(tim, transpose(sol[23:25,:])[:,:], label = ["ωx" "ωy" "ωz"])
#plot(tim, transpose(sol[19:22,:])[:,:], label = ["q0" "q1" "q2" "q3"])
#=
p1 = plot(tim,a)
p2 = plot(tim,e_mod)
p3 = plot(tim,i)
p4 = plot(tim,Ω)
p5 = plot(tim, ω)
p6 = plot(tim,ν)
plot(p1,p2,p3,p4,p5,p6, layout = (3,2), label = ["a" "e" "i" "Ω" "ω" "ν"])
=#