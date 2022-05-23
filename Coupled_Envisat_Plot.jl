#=
Angles1 = zeros(length(tim),3)
Angles2 = zeros(length(tim),3)
tim= sol.t[:]
for j in 1:1:length(tim)
    Angles1[j,:] .= rad2deg.(CosTo321(Beta2Cos(sol[19:22,j])))
end
for j in 1:1:length(tim)
    Angles2[j,:] .= rad2deg.(CosTo123(Beta2Cos(sol[19:22,j])))
end



p1 = plot(tim, Angles1, label = ["yaw" "pitch" "roll"])
p2 = plot(tim, transpose(sol[19:22,:])[:,:], label = ["q0" "q1" "q2" "q3"])
p3 = plot(tim, transpose(rad2deg.(sol[23:25,:]))[:,:],["ωx" "ωy" "ωz"])# label = [L"ω\psi" L"ω\phi"  L"ω\theta"])
p4 = plot(tim, Angles2, label = [L"\phi" L"\theta" L"\psi"])
plot(p2,p3,p4, layout = @layout [a;b;c])=#
#plot(tim, transpose(sol[19:22,:])[:,:], label = ["q0" "q1" "q2" "q3"])
#=
AMATO = CSV.read("AMATO.csv", header=false, DataFrame)[1:182627,:]
H3 =  CSV.read("horizons3.csv", header=false, DataFrame)[:,:]=#
#=
#Th = LinRange(0, 8928*300, 8929)
#=
Xsp = CubicSpline(AMATO[:,8],AMATO[:,2])
Ysp = CubicSpline(AMATO[:,8],AMATO[:,3])
Zsp = CubicSpline(AMATO[:,8],AMATO[:,4])
vxsp = CubicSpline(AMATO[:,8],AMATO[:,5])
vysp = CubicSpline(AMATO[:,8],AMATO[:,6])
vzsp = CubicSpline(AMATO[:,8],AMATO[:,7])=#
Xsp = CubicSpline(AMATO[:,8],AMATO[:,2])
Ysp = CubicSpline(AMATO[:,8],AMATO[:,3])
Zsp = CubicSpline(AMATO[:,8],AMATO[:,4])
vxsp = CubicSpline(AMATO[:,8],AMATO[:,5])
vysp = CubicSpline(AMATO[:,8],AMATO[:,6])
vzsp = CubicSpline(AMATO[:,8],AMATO[:,7])

AMA = Xsp[Th]
#plot(AMATO[:,8],AMATO[:,2])
#l = @layout [a b; c d; e f]
l = @layout[a ;b;c]
p1=plot(Th[:], abs.(sol(Th)[1,:] - Xsp[Th]), label = "rx", xlabel = "time, s", ylabel = "Position, km")
p2=plot(Th[:], abs.( sol(Th)[2,:]- Ysp[Th]), label = "ry",xlabel = "time, s", ylabel = "Position, km")
p3=plot(Th[:], abs.(sol(Th)[3,:]- Zsp[Th]), label = "rz",xlabel = "time, s", ylabel = "Position, km")
p4=plot(Th[:], abs.(sol(Th)[4,:] - vxsp[Th]), label = "vx",xlabel = "time, s", ylabel = "velocity, km/s")
p5=plot(Th[:], abs.(sol(Th)[5,:] - vysp[Th]), label = "vy",xlabel = "time, s", ylabel = "velocity, km/s")
p6=plot(Th[:], abs.(sol(Th)[6,:] - vzsp[Th]), label = "vz",xlabel = "time, s", ylabel = "velocity, km/s")
#plot(p1,p4,p2,p5,p3,p6,layout = l)
#plot(p4,p5,p6,layout = l)
plot(p1,p2,p3,layout = l)
=#


#=
tim = sol.t[:]
R = transpose(sol[1:3,:])[:,:]
V = transpose(sol[4:6,:])[:,:]
a = ones(length(tim),1)
e_mod = ones(length(tim),1)
i = ones(length(tim),1)
ν = ones(length(tim),1)
Ω  =ones(length(tim),1)
ω = ones(length(tim),1)=#
#=
for j in 1:1:length(tim)
    a[j],e_mod[j],i[j],Ω[j],ω[j],ν[j] = RVtoCOE(R[j,:],V[j,:])
    print(a[j],e_mod[j],i[j],Ω[j],ω[j],ν[j])
end 
=#
#plot(tim, transpose(sol[23:25,:])[:,:], label = ["ωx" "ωy" "ωz"])
#plot(tim, transpose(sol[1:3,:])[:,:], label = ["x" "y" "z"])
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