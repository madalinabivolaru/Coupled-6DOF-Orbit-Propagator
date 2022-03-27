
function rsun(JD::Float64)
    #Function determining the sun position vector from the earth to the sun following the Astronomical Almanac 2006 expressed in T_{ut1}
    Tut1 = (JD-2451545.0)/36525
    λM = (280.460 + 36000.771*Tut1)%360
    Ttdb = Tut1
    M = (357.5291092 + 35999.05034*Ttdb)%360
    λecliptic = (λM +1.914666471*sind(M) + 0.019994643 * sind(2*M))%360
    r = 1.000140612 - 0.016708617*cosd(M) - 0.0001395892 * cosd(2*M)
    ϵ = 23.439291 - 0.0130042*Ttdb
    R =  PropagatorConstants.AUtoKM * [r*cosd(λecliptic); r*cosd(ϵ)*sind(λecliptic);r*sind(ϵ)*sind(λecliptic)]
    return transpose(R)[:]
end




