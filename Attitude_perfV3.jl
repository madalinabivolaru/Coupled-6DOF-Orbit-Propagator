using LinearAlgebra
using OrdinaryDiffEq
using DifferentialEquations
using StaticArrays
#using Plots
using BenchmarkTools
#using LaTeXStrings
#using DataFrames
#using CSV
#using RecursiveArrayTools


function Beta2Cos(q)
    #Gives the DCM of the quaternion vector q
    Cfin = [q[1]^2+q[2]^2-q[3]^2-q[4]^2 2*(q[2]*q[3]+q[1]*q[4]) 2*(q[2]*q[4]-q[1]*q[3]);
            2*(q[2]*q[3]-q[1]*q[4]) q[1]^2-q[2]^2+q[3]^2-q[4]^2 2*(q[3]*q[4]+q[1]*q[2]);
            2*(q[2]*q[4]+q[1]*q[3]) 2*(q[3]*q[4]-q[1]*q[2]) q[1]^2-q[2]^2-q[3]^2+q[4]^2]; #checked
end

function CosTo321(Cfin)
    #Gives the 321 Euler angles of the DCM
    ψ, θ, ϕ = atan(Cfin[1,2],Cfin[1,1]), asin(-Cfin[1,3]), atan(Cfin[2,3],Cfin[3,3]) #checked
end

function CosTo313(Cfin)
    #Gives the 313 Euler angles of the DCM
    if abs(Cfin[3,3] - 1)<0.0001
        Cfin[3,3] = 1
    end
    ϕ, θ, ψ = atan(Cfin[1,3],Cfin[2,3]), acos(Cfin[3,3]), atan(Cfin[3,1],-Cfin[3,2])
end

function ()
    
end

c(x) = cosd(x)
s(x) = sind(x)
#=
C = [c(θ)*c(ψ) c(θ)*s(ψ) -s(θ); 
    s(ϕ)*s(θ)*c(ψ)-c(ϕ)*s(ψ) s(ϕ)*s(θ)*s(ψ)+c(ϕ)*c(ψ) s(ϕ)*c(θ);
    c(ϕ)*s(θ)*c(ψ)+s(ϕ)*s(ψ) c(ϕ)*s(θ)*s(ψ)-s(ϕ)*c(ψ) c(ϕ)*c(θ)]; #Cosine matrix for 321 euler angles
=#

#q = [c(ϕ)*c(θ)*c(ψ) - s(ϕ)*c(θ)*s(ψ), c(ϕ)*s(θ)*c(ψ) + s(ϕ)*s(θ)*s(ψ), c(ϕ)*s(θ)*s(ψ) s(ϕ)*s(θ)*c(ψ), c(ϕ)*c(θ)*s(ψ) + s(ϕ)*c(θ)*c(ψ)]
#=
C = [c(ϕ)*c(ψ) - s(ϕ)*c(θ)*s(ψ) c(ϕ)*s(ψ) + s(ϕ)*c(θ)*c(ψ) s(ϕ)*s(θ);
    -s(ϕ)*c(ψ) - c(ϕ)*c(θ)*s(ψ) -s(ϕ)*s(ψ) + c(ϕ)*c(θ)*c(ψ) c(ϕ)*s(θ);
    s(θ)*s(ψ) -s(θ)*c(ψ) c(θ)]; #Cosine matrix for 313 Euler angles =#

#quaternion definition from C is checked against Schaub example
#q = [0.5*sqrt(C[1,1]+C[2,2]+C[3,3]+1) 1/4/(0.5*sqrt(C[1,1]+C[2,2]+C[3,3]+1))*(C[2,3]-C[3,2]) 1/4/(0.5*sqrt(C[1,1]+C[2,2]+C[3,3]+1))*(C[3,1]-C[1,3]) 1/4/(0.5*sqrt(C[1,1]+C[2,2]+C[3,3]+1))*(C[1,2]-C[2,1])]

Ωq = ones(1,4)
II = ones(1,4)

function Attitude!(dy::Vector{Float64},y::Vector{Float64},p,t)
    Ωq = ones(4,1)
    I = diagm([42.24,104.93, 105.82])
    #I = diagm([400,750,850])

    #Normalise quaternion 
    q =  @view y[1:4]
    ω =  @view y[5:7]
    dq = @view dy[1:4]
    dω = @view dy[5:7]
    #print(q[1])
    #Normalise quaternion    
    q = q./norm(q)
    

    Ω = Float64[0.0 -ω[1] -ω[2] -ω[3];
        ω[1] 0 ω[3] -ω[2];
        ω[2] -ω[3] 0 ω[1];
        ω[3] ω[2] -ω[1] 0.0];

    n = norm([1.3219e8,0.0, 0.2645e8])
    Tsp = 1367/3/10^8*1.6*(1.32*1.11*10^3*[1.3219e8,0, 0.2645e8] + (1.32*1.11*10^3*[1.3219e8, 0 , 0.2645e8]) + 0.01*1.11*10^3*[1.3219e8,0, 0.2645e8])
    Tsp = 1000*1367/3/10^8*1.6*(cross([0.1,0.01,0.01],1.11*[1.3219e8,0, 0.2645e8]/n)+cross([0.1,0.01,0.01],1.11*[1.3219e8,0, 0.2645e8]/n)+cross([0,0,0],10^3*[1.3219e8,0, 0.2645e8]/n))
    #=
    if t>= 0 && t<0.001
        T = [1000,1000,0]
    else
        T = [0,0,0]
    end=#
    Ta = 0.5*0.1*(3.019*10^-15)*2.2*1.2*7500^2
    v = 1000*[-0.9692300284946672, -0.1883992325643531, 7.540501269182395]
    r = 1000*[1315.0823093944996, -6765.511976126983, 0.0]

    vr  = [v[1]+0.7292*10^-4*r[2], v[2]-0.7292*10^-4*r[1], v[3]]
    v_r = norm(vr)
    Drag = 1000*[0.1,0,0].*(0.5*3.019*10^-15*2.2*1/100 * v_r^2 * vr/v_r);
    Drag = [0,0,0]
    TSP = [0,0,0]
    #nn = [0,0,0]-[(I[3]-I[2])*ω[2]*ω[3], (I[1]-I[3])*ω[1]*ω[3], (I[2]-I[1])*ω[2]*ω[1]] 
    #n = 0.5.*[-q[2]*ωx-q[3]*ωy-q[4]*ωz, q[1]*ωx-q[4]*ωy+q[3]*ωz, q[4]*ωx+q[1]*ωy-q[2]*ωz, q[2]*ωy-q[3]*ωx+q[1]*ωz]
    #print(n)
    
    @. dq = 0.5*mul!([0,0,0,0], Ω, q)
    @. dω = I^-1*([0,0,0]-cross(ω,I*ω))#mul!(II,I^-1,([0,0,0]-cross(ω,I*ω))) #[nn[1]/I[1], nn[2]/I[2],nn[3]/I[3]]
    #@SVector [dy[1:4], dy[5:7]]
 
end

#q0 = [1, 0, 0, 0]
q0 = [1, 0, 0, 0]
ω0 = [0, 0.5, 1]
#Attitude0 = @SVector [q0,ω0]
#Attitude = ArrayPartition(x::AbstractArray...)
#Attitude0 = reduce(vcat, [q0,ω0])
Attitude0 = [1, 0, 0, 0, 0, 0.05, 0.175]

tspan = (0.0, 3600*24*10)

Xtol = 1e-13#[1e-10,1e-10,1e-10,1e-10,1e-8,1e-8,1e-8] #
probA = ODEProblem(Attitude!,Attitude0,tspan)
sol2 = solve(probA, VCABM(),reltol = Xtol, abstol = Xtol, save_everystep=false, maxiters=1e8)#reltol=[1e-10,1e-10,1e-10,1e-10,1e-12,1e-12,1e-12], abstol=[1e-10,1e-10,1e-10,1e-10,1e-12,1e-12,1e-12], maxiters = 1e8


#plot(tim, Q)
#=
plot(sol2.t[:], transpose(sol1[1:4,:])[:,:]-transpose(sol2[1:4,:])[:,:], labels = false)
ylabel!("Quaternion Difference, (-)")
xlabel!("Time, (s)") =#
#=
tim = sol1.t[:]
Q = transpose(sol1[1:4,:])[:,:]
W = transpose(sol1[5:7,:])[:,:]=#
#=
tim = sol1.t[:]
Q = transpose(sol1[1:4,:])[:,:]
W = transpose(sol1[5:7,:])[:,:]

qq = zeros(length(tim),3)
Ww = zeros(length(tim),3)=#
#=
for i in 1:1:length(Q[:,1])
    qq[i,:] .= broadcast(rad2deg,CosTo313(Beta2Cos(Q[i,:]./norm(Q[i,:]))))
    #println(CosTo313(Beta2Cos(Q[i,:]./norm(Q[i,:]))))
    #println(broadcast(rad2deg,CosTo313(Beta2Cos(Q[i,:]))))
end 

for i in 1:1:length(W[:,1])
    Ww[i,:] .= W[i,:]
    #println(broadcast(rad2deg,W[i,:]))
end =#
#=
for i in 1:1:length(W[:,1])
    #println(transpose(Ww[i,:])[:,:]*diagm([0.234375, 0.46875, 0.234375])*Ww[i,:])
    println( dot(-Ww[i,:],cross(Ww[i,:],diagm([0.234375, 0.46875, 0.234375])*Ww[i,:])))
end=#

#plot(tim, W, label =[L"$q_0$" L"$q_1$" L"$q_2$" L"$q_3$"])
#plot(tim, W, label =[L"$\omega_x$" L"$\omega_y$" L"$\omega_z$"])
#ylabel!("Quaternion, (-)")
#ylabel!("Angular velocity, (rad/s)")
#xlabel!("Time, (s)")
#plot(sol.t[:],sol[1,:])
#println((Beta2Cos(Q[26,:])))
#println((CosTo321(transpose(Beta2Cos(Q[2,:]))[:,:])))
#println(Beta2Cos(Q[3,:]))
#println(norm(Q[11,:]))





