using LinearAlgebra, Plots
H3 =  CSV.read("horizons3.csv", header=false, DataFrame)[:,:]
Th = LinRange(0, 8299*300, 8300)
p1=plot(Th[:], (H3[1:8300,1] - sol(Th)[1,:]), label = "rx", xlabel = "time, s", ylabel = "Position, km")
p2=plot(Th[:], (H3[1:8300,2] - sol(Th)[2,:]), label = "ry", xlabel = "time, s", ylabel = "Position, km")
p3=plot(Th[:], (H3[1:8300,3] - sol(Th)[3,:]), label = "rz", xlabel = "time, s", ylabel = "Position, km")
p4=plot(Th[:], (H3[1:8300,4] - sol(Th)[4,:]), label = "vx", xlabel = "time, s", ylabel = "v, km/s")
p5=plot(Th[:], (H3[1:8300,5] - sol(Th)[5,:]), label = "vy", xlabel = "time, s", ylabel = "v, km/s")
p6=plot(Th[:], (H3[1:8300,6] - sol(Th)[6,:]), label = "vz", xlabel = "time, s", ylabel = "v, km/s")
l = @layout[a ;b;c]
plot(p1,p2,p3,layout = l)