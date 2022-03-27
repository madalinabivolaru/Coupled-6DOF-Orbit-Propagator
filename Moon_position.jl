
function rmoon(JD::Float64)
    #Moon position vector following the Astronomical Almanac 2006 expressed in T_{ut1}
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

    r_moon = PropagatorConstants.R_e/sind(P) 
    R_moon = r_moon * [cosd(Φecl)*cosd(λecl); 

    cosd(Ε)*cosd(Φecl)*sind(λecl) - sind(Ε)*sind(Φecl);
    sind(Ε)*cosd(Φecl)*sind(λecl) + cosd(Ε)*sind(Φecl)]

    return transpose(R_moon)[:]
end

