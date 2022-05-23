#This is an attempt at linearising Attitude dynamics
#define ODE to be used:
#=
function LinearAttitude!(dy::Vector{Float64},y::Vector{Float64},p,t)

    I = [42.24, 104.93, 105.82]
#I = diagm([400,750,850])

#Normalise quaternion 
    X = y[1:7]
    Ẋ = dy[1:7]

    ω = X[5:7]
    q = X[1:4]
#print(q[1])
#Normalise quaternion    
    q = q ./norm(q)

    #def torques
    Torques = [0., 0., 0.]

    A = [0. -ω[1]/2 -ω[2]/2 -ω[3]/2 -q[2]/2 -q[3]/2 -q[4]/2;
        ω[1]/2 0. ω[3]/2 -ω[2]/2 q[1]/2 -q[4]/2 q[3]/2;
        ω[2]/2 -ω[3]/2 0. ω[1]/2 q[4]/2 q[1]/2 -q[2]/2;
        ω[3]/2 ω[2]/2 -ω[1]/2 0. -q[3]/2 q[2]/2 q[1]/2;
        0. 0. 0. 0. 0. (Torques[1]-(I[3]-I[2])*ω[3])/I[1]  (Torques[1]-(I[3]-I[2])*ω[2])/I[1];
        0. 0. 0. 0. (Torques[2]-(I[1]-I[3])*ω[3])/I[2] 0. (Torques[2]-(I[1]-I[3])*ω[1])/I[2];
        0. 0. 0. 0. (Torques[3]-(I[2]-I[1])*ω[2])/I[2] (Torques[2]-(I[2]-I[1])*ω[1])/I[2] 0.]

    #expA = exp(A)
    dy[1:7] = exp(A)*X
end=#
using IterativeSolvers
function fifi(X)
    X = zeros(3601,7)
    X[1,1:7] = [1, 0, 0, 0, 0.9, 0, 0.]
    tim = 3600
    I = [42.24, 104.93, 105.82]
    for j in 1:1:tim
        Xx = X[j, :]
        ω = Xx[5:7]
        q = Xx[1:4]
   
        q = q ./norm(q)

    #def torques
        Torques = [0., 0., 0.]
        A = [0. -ω[1]/2 -ω[2]/2 -ω[3]/2 -q[2]/2 -q[3]/2 -q[4]/2;
        ω[1]/2 0. ω[3]/2 -ω[2]/2 q[1]/2 -q[4]/2 q[3]/2;
         ω[2]/2 -ω[3]/2 0. ω[1]/2 q[4]/2 q[1]/2 -q[2]/2;
        ω[3]/2 ω[2]/2 -ω[1]/2 0. -q[3]/2 q[2]/2 q[1]/2;
         0. 0. 0. 0. 0. (Torques[1]-(I[3]-I[2])*ω[3])/I[1]  (Torques[1]-(I[3]-I[2])*ω[2])/I[1];
        0. 0. 0. 0. (Torques[2]-(I[1]-I[3])*ω[3])/I[2] 0. (Torques[2]-(I[1]-I[3])*ω[1])/I[2];
        0. 0. 0. 0. (Torques[3]-(I[2]-I[1])*ω[2])/I[2] (Torques[2]-(I[2]-I[1])*ω[1])/I[2] 0.]

        X[j+1,1:7] = exp(A)*X[j,1:7]
    end
    return X
end

fifi([1, 0, 0, 0, 0.9, 0, 0.])

#=
q0 = [1 0 0 0]
ω0 = [0 0 0.175]
X0 = [1, 0, 0, 0, 0.9, 0, 0.]
tspan = (0.0, 360)


Xtol = 1e-10#[1e-10,1e-10,1e-10,1e-10,1e-8,1e-8,1e-8] #
prob_lin = ODEProblem(LinearAttitude!,X0,tspan)
sol_linear = solve(prob_lin, RK4(),dt =1, reltol = Xtol, abstol = Xtol) #save_everystep=false, maxiters=1e8)=#