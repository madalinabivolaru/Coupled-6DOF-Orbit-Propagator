function COEtoRV(p,e,i,Ω,ω,ν)
    μ = 398600.4418; 
    Rpqw = [p*cosd(ν)/(1+e*cosd(ν)), p*sind(ν)/(1+e*cosd(ν)), 0 ]
    Vpqw = [-sqrt(μ/p)*sind(ν), sqrt(μ/p)*(e+cosd(ν)), 0 ]
    c(α) = cosd(α)
    s(α) = sind(α)

    PQWtoIJK = [c(Ω)*c(ω) - s(Ω)*s(ω)*c(i) -c(Ω)*s(ω) - s(Ω)*c(ω)*c(i) s(Ω)*s(i)
                s(Ω)*c(ω) + c(Ω)*s(ω)*c(i) -s(Ω)*s(ω) + c(Ω)*c(ω)*c(i) -c(Ω)*s(i)
                s(ω)*s(i) c(ω)*s(i) c(i)]
    Rijk, Vijk =  PQWtoIJK*Rpqw, PQWtoIJK*Vpqw
end

COEtoRV(11067.790,0.83285,87.87,227.89,53.38,92.335)
                  