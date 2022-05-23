#Propagator test on ENVISAT Geometry as proposed by Tsinghua Paper

using LinearAlgebra
using OrdinaryDiffEq
using DifferentialEquations
using BenchmarkTools
include("./PropagatorConstants.jl")
using .PropagatorConstants
include("./JulDat.jl")
include("asrp.jl")
include("RVtoCOE.jl")
using StaticArrays


function A2(r, v)
    #Function which defines rotation matrix A2 as ECI = A2*NTW
    T = v./(norm(v))
    W = cross(r,v)/norm(cross(r,v))
    N = cross(T,W)

    return [N T W]
end

function CosTo123(Cfin)
    #Gives the 321 Euler angles of the DCM
    ψ, θ, ϕ = atan(-Cfin[3,2],Cfin[3,3]), asin(Cfin[3,1]), atan(-Cfin[2,1],Cfin[1,1]) #checked
end

function CosTo321(Cfin)
    #Gives the 321 Euler angles of the DCM
    ψ, θ, ϕ = atan(Cfin[1,2],Cfin[1,1]), asin(-Cfin[1,3]), atan(Cfin[2,3],Cfin[3,3]) #checked
end
function Beta2Cos(q::Vector{Float64})
    #function which returns A1 for quaternions propagated; Here A1 is defines as NTW = A1*Body
    q = q./norm(q)
    @inbounds Cfin = [q[1]^2+q[2]^2-q[3]^2-q[4]^2 2*(q[2]*q[3]+q[1]*q[4]) 2*(q[2]*q[4]-q[1]*q[3]);
            2*(q[2]*q[3]-q[1]*q[4]) q[1]^2-q[2]^2+q[3]^2-q[4]^2 2*(q[3]*q[4]+q[1]*q[2]);
            2*(q[2]*q[4]+q[1]*q[3]) 2*(q[3]*q[4]-q[1]*q[2]) q[1]^2-q[2]^2-q[3]^2+q[4]^2]; #checked
    #if abs(Cfin[3,3] - 1)<0.0001
        #Cfin[3,3] = 1
    #end
end

function KepCowell!(dy::Vector{Float64},y::Vector{Float64},p::SVector{13, Float64},t::Float64)
    #Define constants, body vectors and area of panels
    Lsp = 1.32
    A, m, J2,J3,J4, μe, μm, μs, Re, ωe, Cd, pos, jd0 = p
    jd = jd0 + t/3600/24
    A = @SVector [71.12, 15.64, 22.92, 38.26 ];

    #Normal in satellite reference frame
    nx = [1,0,0]
    nz = [0,0,1]
    ny = [0,1,0]
    npanel = [-sind(22), cosd(22), 0]

    rpanel = -6*ny #solar panel offset
    ry = 1.7*ny 
    rx = 1.0*nx
    rz = 1.0*nz
   
    #I = diagm([42.24,104.93, 105.82])
    #I = diagm([129180.25,16979.74, 124801.21]) #km m^2
    I = @SVector [129180.25,16979.74, 124801.21]
    rho0 = @SVector [3.735*10^-12, 1.585*10^-12,6.967*10^-13, 1.454*10^-13,
    3.614*10^-14,1.17*10^-14, 5.245*10^-15,3.019*10^-15,
    9.518*10^-12, 2.418*10^-11, 7.248*10^-11, 2.789*10^-10, 5.464*10^-10 ]
    SF =1367
    p_srp = SF/3*10^-8;
    Cr = 1.2;

    r = @view y[1:3]
    v = @view y[4:6]
    rsun = @view y[7:9]
    vsun = @view y[10:12]
    rmoon = @view y[13:15]
    vmoon = @view y[16:18]    
    q = @view y[19:22]
    ω = @view y[23:25]

    dr = @view dy[1:3] 
    dv = @view dy[4:6] 
    drs = @view dy[7:9]
    dvs = @view dy[10:12] 
    drm = @view dy[13:15]
    dvm = @view dy[16:18] 
    dq = @view dy[19:22]
    dω = @view dy[23:25]

    #Normalise quaternion vector for long term propagations
    q = q./norm(q)
    #Get total Attitude matrix ECI in Body
    At2 = A2(r,v)
    At1 = Beta2Cos(q)


    r_norm = norm(r)
    #gg = Geogrid(r, jd, eop_IAU1980)
    #print(gg)
    if (pos == 1)
        rsn, rmn = norm(rsun), norm(rmoon)
        rmMB, rsMB = rmoon, rsun
    else
        rmMB = rmoon(jd)
        rsMB = rsun(jd)
        rsn, rmn = norm(rsMB), norm(rmMB)
    end

    #r3ms = rsu - rmoo
    rsat3m = rmMB .- r
    rsat3s = rsMB .-r
    
    #Third body accelerations
    a3m = μm .* (rsat3m./norm(rsat3m)^3 .- rmMB./rmn^3) #CK
    a3s = μs .* (rsat3s./norm(rsat3s)^3 .- rsMB./rsn^3) 

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

    ECItoBody = Array{Float64}(undef, 3,3)
    #Perform rotation from ECI to Body to get Drag torque
    mul!(ECItoBody, At2,At1)
    ECItoBody = inv(ECItoBody)
    vr_B = Vector{Float64}(undef, 3)
    mul!(vr_B,ECItoBody, vr./v_r) 
    vs_B = Vector{Float64}(undef, 3)
    mul!(vs_B, ECItoBody, (-rsat3s)./norm(rsat3s))
    
    cosSx = dot(nx,vs_B)
    cosSy = dot(ny,vs_B)
    cosSz = dot(nz,vs_B)
    cosSp = dot(npanel,vs_B)


    cosTx = dot(nx,vr_B)
    cosTy = dot(ny,vr_B)
    cosTz = dot(nz,vr_B)
    cosTp = dot(npanel,vr_B)

    #print(cosTx)
    #Drag coupled with attitude
    #Change As of asrp
    #asrp = -10^(-3)*(p_srp*Cr*As/m)*rsun/norm(rsun)^3*PropagatorConstants.AUtoKM^2
    @inbounds Drag = -10^3*0.5*rho*Cd*(A[1]*max(cosTp,0)+A[1]*max(-cosTp,0) + A[3]*max(cosTx,0) + A[3]*max(-cosTx,0) +
            A[2]*max(cosTy,0) + A[2]*max(-cosTy,0) + A[4]*max(cosTz,0) + A[4]*max(-cosTz,0))/m * v_r^2 .* vr/v_r;

    #Fsrp = -10^(-3)*p_srp*Cr*vs_B/norm(vs_B)*(4*max(cosSx,0)+2*max(cosSy,0)+2*max(cosSx,0))
    #asrp = inv(ECItoBody)*(Fsr@p/m)
    #As = 1
    rsats = r-rsun
    @inbounds asrp = 10^(-3).*(p_srp*Cr*(A[1]*max(cosSp,0)+A[1]*max(-cosSp,0) + A[3]*max(cosSx,0) + A[3]*max(-cosSx,0) +
            A[2]*max(cosSy,0) + A[2]*max(-cosSy,0) + A[4]*max(cosSz,0) + A[4]*max(-cosSz,0))/m).*rsats./norm(rsats)

    @inbounds Fs = [(-p_srp*Cr*vs_B)*A[1]*max(cosSp,0), (-p_srp*Cr*vs_B)*A[1]*max(-cosSp,0), (-p_srp*Cr*vs_B)*A[3]*max(cosSx,0),
           (-p_srp*Cr*vs_B)*A[3]*max(-cosSx,0), (-p_srp*Cr*vs_B)*A[2]*max(cosSy,0), (-p_srp*Cr*vs_B)*A[2]*max(-cosSy,0),
           (-p_srp*Cr*vs_B)*A[4]*max(cosSz,0), (-p_srp*Cr*vs_B)*A[4]*max(-cosSz,0)]
    
           #Drag components for torques
    @inbounds Fd = [-10^6*0.5*rho*Cd*v_r.*vr_B*A[1]*max(cosTp,0), -10^6*0.5*rho*Cd*v_r.*vr_B*A[1]*max(-cosTp,0),
            -10^6*0.5*rho*Cd*v_r.*vr_B*A[3]*max(cosTx,0), -10^6*0.5*rho*Cd*v_r.*vr_B*A[3]*max(-cosTx,0),
            -10^6*0.5*rho*Cd*v_r.*vr_B*A[2]*max(cosTy,0), -10^6*0.5*rho*Cd*v_r.*vr_B*A[2]*max(-cosTy,0),
            -10^6*0.5*rho*Cd*v_r.*vr_B*A[4]*max(cosTz,0), -10^6*0.5*rho*Cd*v_r.*vr_B*A[4]*max(-cosTz,0)]

    @inbounds Ld = 1/m.*(cross(rpanel, Fd[1]) .+ cross(rpanel, Fd[2]) .+ cross(ry, Fd[5]) .+ cross(-ry, Fd[6]) .+ cross(rx, Fd[3]) .+ cross(-rx, Fd[4]) .+
            cross(rz, Fd[7]) .+ cross(rz, Fd[8]))
    #print(Ld)
    @inbounds Lsrp = cross(rpanel, Fs[1]) .+ cross(rpanel, Fs[2]) .+ cross(ry, Fs[5]) .+ cross(-ry, Fs[6]) .+ cross(rx, Fs[3]) .+ cross(-rx, Fs[4]) .+
             cross(rz, Fs[7]) .+ cross(rz, Fs[8])
    #print(Lsrp)
    #Attitude 

    #println((1.11*2*max(cosSx,0)+ max(cosSx,0) + (1.11*2+1)*max(cosSxx,0)+max(cosSy,0)+max(cosSyy,0)+max(cosSz,0)+max(cosSzz,0))/m)
    @. dr = v
    @. dv = -r*μe/r_norm^3 .+ a3m .+ a3s .+ a_J2 .+ a_J3 .+ a_J4 .+ asrp .+ Drag
    @. drs = vsun
    @. dvs = -rsun*μs./rsn^3
    @. drm = vmoon
    @. dvm = -rmoon*μe./rmn^3 
    #print(length(transpose(0.5 .* Ω* transpose(q)[:])[:]))
    #print(length(transpose(I^-1*(transpose(Ld)[:]-cross(transpose(ω)[:],I*transpose(ω)[:])))[:]))
    #@. dq = 0.5.*Ω*q
    #@. dω = I^-1*(Ld + Lsrp-cross(ω,I*ω))
    
    @. dq = 0.5.*[-q[2]*ω[1]-q[3]*ω[2]-q[4]*ω[3], q[1]*ω[1]-q[4]*ω[2]+q[3]*ω[3], q[4]*ω[1]+q[1]*ω[2]-q[2]*ω[3], q[2]*ω[2]-q[3]*ω[1]+q[1]*ω[3]]
    @. dω = (Ld + Lsrp - [(I[3]-I[2])*ω[2]*ω[3], (I[1]-I[3])*ω[1]*ω[3], (I[2]-I[1])*ω[2]*ω[1]] )./I #[nn[1]/I[1], nn[2]/I[2],nn[3]/I[3]]
    
end


y0 = [-2.402809058043280E+03,  8.537777890343411E+02, -6.682371096072878E+03, -4.951358974727093E+00,  5.020450518834991E+00,  2.422381912322207E+00,
1.463957088922231e+08, -2.335131546485938e+07, -1.012403398015358e+07, 5.576319298565035, 27.026588503272446, 11.716155914748752,
-1.170981644822336e+05, 3.683752523717154e+05, 1.208032774797189e+05, -0.927440482124745, -0.241519315834115, -0.125577393037498,
#1.4679e8,-0.2135e8, -0.0926e8, 5.1388, 27.0915 , 11.7445, 
#-1.8342e5, 3.4449e5, 1.0956e5, -0.8613, -0.4026, -0.1776,
1,0,0,0,
0.0,0.0,0]
#0.906, 0, 0.423, 0,
#0.,0.,0.1745] 


tolx = 1e-11#=[1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,
        1e-10,1e-10,1e-10,1e-9,1e-9,1e-9,#1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,#
        1e-10,1e-10,1e-10,1e-9,1e-9,1e-9,#1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,
        1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,1e-11]#1e-10,1e-10,1e-10,1e-10,1e-8,1e-8,1e-8,1e-8]=#

pos = 1.;
jd0 = 2456727.500;

p =   @SVector [12.0, 7991.0, PropagatorConstants.J2, PropagatorConstants.J3, PropagatorConstants.J4, 
        PropagatorConstants.μ_e, PropagatorConstants.μ_m, PropagatorConstants.μ_s, 
        PropagatorConstants.R_e, PropagatorConstants.ω_e, PropagatorConstants.C_d, pos, jd0]

tspan = (0.0,3600*24*36.5*20)

prob = ODEProblem(KepCowell!,y0,tspan,p)
@time sol = solve(prob, VCABM(), reltol=tolx, abstol=tolx, maxiters = 1e8, save_everystep=false)#,  maxiters = 1e8,)#, save_everystep = false)#,  save_everystep=false)#saveat = 36000)
#=
for i in 1:10
    global y0
    local tspan = (0.,3600*24*36.5*i)
    local prob = ODEProblem(KepCowell!,y0,tspan,p)
    local sol = solve(prob, VCABM(), reltol=tolx, abstol=tolx, maxiters = 1e8, save_everystep=false)
    y0 = sol(3600*24*36.5*i)
end=#