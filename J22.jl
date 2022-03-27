using LinearAlgebra

function J22(r::Vector{Float64})
    #Calculate acceleration due to J22 influence as per Montenbruck p66 and described in Cunningham (1970)
    #Using recurrence relations for the evaluation of Legendre polynomials
    #Earth-fixed coordinate system as a function of the Earth-fixed position vector r = (x, y, z).
    #Conversion to inertial frame system is required
    include("./PropagatorConstants.jl")

    #Still need to transform r in ECI to ECEF
    r_norm = norm(r)
    println(r_norm)
    x22 = -3*PropagatorConstants.R_e^2*PropagatorConstants.μ_e/2/r_norm^7*
            (PropagatorConstants.C22*(r_norm^2*r[1]+5*r[1]^3-15*r[1]*r[2]^2-5*r[1]*r[3]^2)+
            PropagatorConstants.S22*(r_norm^2*r[2]+15*r[1]^2*r[2]-5*r[2]^3-5*r[2]*r[3]^2));

    y22 = 3*PropagatorConstants.R_e^2*PropagatorConstants.μ_e/2/r_norm^7*
            (PropagatorConstants.C22*(r_norm^2*r[2]-15*r[1]^2*r[2]+5*r[2]^3-5*r[2]*r[3]^2)+
            PropagatorConstants.S22*(-r_norm^2*r[1]+5*r[1]^3-15*r[1]*r[2]^2+5*r[1]*r[3]^2));

    z22 = -15*PropagatorConstants.R_e^2*PropagatorConstants.μ_e*r[3]/r_norm^7*
            (PropagatorConstants.C22*(r[1]^2-r[2]^2)+PropagatorConstants.S22*r[1]*r[2]);
    
    println(norm([x22,y22,z22]))
    aJ22 = [x22,y22,z22]
    #perform ECEF to ECI tranformation to get aJ22
end

J22([6672.423,1131.340,- 2282.343 ])