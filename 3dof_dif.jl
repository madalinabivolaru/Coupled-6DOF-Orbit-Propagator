#Th = sol.t[:]
#Th = LinRange(8)
using CSV
AMATO = CSV.read("AMATO.csv", DataFrame)[1:182627,:]
AM = CSV.read("cartJ3J4.csv", header=false, DataFrame)[1:182626,2:7]
AM =  CSV.read("horizons3.csv", header=false, DataFrame)[:,:]
#Th = LinRange(0, 31550169.600, 182583)
Th = LinRange(0, 0.002*24*3600*182625,182626)
Th = LinRange(0, 5*60*8928,8929)
#=
Xsp = CubicSpline(AMATO[:,8],AMATO[:,2])
Ysp = CubicSpline(AMATO[:,8],AMATO[:,3])
Zsp = CubicSpline(AMATO[:,8],AMATO[:,4])
vxsp = CubicSpline(AMATO[:,8],AMATO[:,5])
vysp = CubicSpline(AMATO[:,8],AMATO[:,6])
vzsp = CubicSpline(AMATO[:,8],AMATO[:,7])=#
#AMA = Xsp[Th]
#plot(AMATO[:,8],AMATO[:,2])
#l = @layout [a b; c d; e f]
#Th =convert.(Float16, AM[:,8])
p1=plot(Th[:], abs.(sol(Th)[1,:]- AM[:,1]), label = "rx", xlabel = "time, s", ylabel = "δr, km")
p2=plot(Th[:], abs.( sol(Th)[2,:]- AM[:,2]), label = "ry",xlabel = "time, s", ylabel = "δr, km")
p3=plot(Th[:], abs.(sol(Th)[3,:]- AM[:,3]), label = "rz",xlabel = "time, s", ylabel = "δr, km")
p4=plot(Th[:], abs.(sol(Th)[4,:] - AM[:,4]), label = "vx",xlabel = "time, s", ylabel = "δv, km/s")
p5=plot(Th[:], abs.(sol(Th)[5,:] - AM[:,5]), label = "vy",xlabel = "time, s", ylabel = "δv, km/s")
p6=plot(Th[:], abs.(sol(Th)[6,:] - AM[:,6]), label = "vz",xlabel = "time, s", ylabel = "δv, km/s")


l = @layout[a ;b;c]#=
p1=plot(Th[:], abs.(sol(Th)[1,:] - Xsp[Th]), label = "rx", xlabel = "time, s", ylabel = "Position, km")
p2=plot(Th[:], abs.( sol(Th)[2,:]- Ysp[Th]), label = "ry",xlabel = "time, s", ylabel = "Position, km")
p3=plot(Th[:], abs.(sol(Th)[3,:]- Zsp[Th]), label = "rz",xlabel = "time, s", ylabel = "Position, km")
p4=plot(Th[:], abs.(sol(Th)[4,:] - vxsp[Th]), label = "vx",xlabel = "time, s", ylabel = "velocity, km/s")
p5=plot(Th[:], abs.(sol(Th)[5,:] - vysp[Th]), label = "vy",xlabel = "time, s", ylabel = "velocity, km/s")
p6=plot(Th[:], abs.(sol(Th)[6,:] - vzsp[Th]), label = "vz",xlabel = "time, s", ylabel = "velocity, km/s")
#plot(p1,p4,p2,p5,p3,p6,layout = l)=#
plot(p1,p2,p3,layout = l)