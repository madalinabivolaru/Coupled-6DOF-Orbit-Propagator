
#using SatelliteToolbox
#using LinearAlgebra
function Geogrid(r, jd)
    #LOD = eop.LOD[]
    mu =  3.986004414498200E+05
    Re = 6378.136460
    ω = 0.7292115*10^-4#*(1-LOD/86400)
    #ω = 86400/(2*pi/ω)
    #-r*mu/r_norm^3
    gmst = jd_to_gmst(jd)
    #D_eci_ecef = r_eci_to_ecef(GCRF(),ITRF(),  jd, eop_IAU1980)
    #r =  D_eci_ecef*r
    r = [r[1]*cos(gmst)+r[2]*sin(gmst),
         r[2]*cos(gmst)-r[1]*sin(gmst),
         r[3]]
    #v =  D_eci_ecef*(v - cross([0,0,ω],r))

    r_norm = norm(r)
    #r = r/r_norm
    #print(r)
    #C20 = -4.84165115570149964E-04
    S20 = 0
    #=
    C21 = -1.62865090839200010E-10
    S21 = 1.32488349710740023E-09
    C22 = 2.43931947200109987E-06
    S22 = -1.40027653068339989E-06=#
    C20 = -1.08E-03
    C21 = -2.66739475237484E-10
    S21 = 1.78727064852404E-09
    C22 = 1.57461532572292E-06
    S22 = -9.03872789196567E-07
     #=
    C21 = -2.06615509074176E-10
    S21 = 1.38441389137979E-09
    C22 = 2.43938357328313E-06
    S22 = -1.40027653068339989E-06=#
    #C21 = -2.66739475237484E-10
    #S21 = 1.78727064852404E-09
    #C20 = -4.84165115570149964E-04
   # C20 = -1.08E-03
    #S20 = 0
    #C21 = 0
    #S21 = 0
    #C22 = 1.57461532572292E-06
    #S22 = -9.03872789196567E-07


    V = zeros(3,3)
    W = V
    V00 = Re/r_norm
    W00 = 0
    m = 1
    #print(r)
    @inbounds V[1,1] = (2*m-1)*(Re*r[1]/r_norm^2*V00- Re*r[2]/r_norm^2*W00)
    @inbounds W[1,1] = (2*m-1)*(Re*r[1]/r_norm^2*W00 +  Re*r[2]/r_norm^2*V00)
    m = 2
    @inbounds V[2,2] = (2*m-1)*(Re*r[1]/r_norm^2*V[m-1,m-1]- Re*r[2]/r_norm^2*W[m-1,m-1])
    @inbounds W[2,2] = (2*m-1)*(Re*r[1]/r_norm^2*W[m-1,m-1] +  Re*r[2]/r_norm^2*V[m-1,m-1])
    m = 3
    @inbounds V[3,3] = (2*m-1)*(Re*r[1]/r_norm^2*V[m-1,m-1]- Re*r[2]/r_norm^2*W[m-1,m-1])
    @inbounds W[3,3] = (2*m-1)*(Re*r[1]/r_norm^2*W[m-1,m-1] +  Re*r[2]/r_norm^2*V[m-1,m-1])

    n =1
    @inbounds V10 = (2*n-1)/n*r[3]*Re/r_norm^2*V00 - (n-1)/n*Re^2/r_norm^2
    @inbounds W10 = 0#(2*n-1)/n*r[3]*Re/r_norm^2*W00 - (n-1)/n*Re^2/r_norm^2

    n =2 
    @inbounds V20 = (2*n-1)/n*r[3]*Re/r_norm^2*V10 - (n-1)/n*Re^2/r_norm^2*V00
    @inbounds W20 = 0#(2*n-1)/n*r[3]*Re/r_norm^2*W10 - (n-1)/n*Re^2/r_norm^2*W00

    n=3
    @inbounds V30 = (2*n-1)/n*r[3]*Re/r_norm^2*V20 - (n-1)/n*Re^2/r_norm^2*V10
    @inbounds W30 = 0#(2*n-1)/n*r[3]*Re/r_norm^2*W20 - (n-1)/n*Re^2/r_norm^2*W10

    n=2
    m=1
    @inbounds V[2,1] = (2*n-1)/(n-m)*r[3]*Re/r_norm^2*V[1,1] 
    @inbounds W[2,1] = (2*n-1)/(n-m)*r[3]*Re/r_norm^2*W[1,1] 

    n=3
    m=1
    @inbounds V[3,1] = (2*n-1)/(n-m)*r[3]*Re/r_norm^2*V[2,1] - (n+m-1)/(n-m)*Re^2/r_norm^2*V[1,1]
    @inbounds W[3,1] = (2*n-1)/(n-m)*r[3]*Re/r_norm^2*W[2,1] - (n+m-1)/(n-m)*Re^2/r_norm^2*W[1,1]

    n=3
    m=2
    @inbounds V[3,2] = (2*n-1)/(n-m)*r[3]*Re/r_norm^2*V[2,2] 
    @inbounds W[3,2] = (2*n-1)/(n-m)*r[3]*Re/r_norm^2*W[2,2] 

    #create accelerations from newly found parameters

    @inbounds  a20 = mu/Re^2*[-C20*V[3,1],
         -C20*W[3,1],
         3*(-C20*V30-S20*W30)]

    @inbounds a21 = mu/Re^2*[0.5*(-C21*V[3,2] - S21*W[3,2] + factorial(2-1+2)/factorial(2-1) * (C21*V30 + S21*W30)),
         0.5*(-C21*W[3,2] + S21*V[3,2] + factorial(2-1+2)/factorial(2-1) * (-C21*W30 + S21*V30)),
         (2-1+1)*(-C21*V[3,1] - S21*W[3,1])]

    @inbounds a22 = mu/Re^2*[0.5*(-C22*V[3,3] - S22*W[3,3] + factorial(2-2+2)/factorial(2-2) * (C22*V[3,1] + S22*W[3,1])),
         0.5*(-C22*W[3,3] + S22*V[3,3] + factorial(2-2+2)/factorial(2-2) * (-C22*W[3,1] + S22*V[3,1])),
         (2-2+1)*(-C22*V[3,2] - S22*W[3,2])]
    
    #D_ecef_eci = r_ecef_to_eci(ITRF(), GCRF(), jd, eop_IAU1980)
   
    #a =  D_ecef_eci*((a20) +  [0;0;ω] × r + 2*[0;0;ω] × v)  
    #v =  D_ecef_eci*(v + cross([0,0,ω],r))
    #r =  D_ecef_eci*r
    #print(a)
    #print(a20) 
    #a = transpose(D_eci_ecef)[:,:]*(a20)  FIn(1) = F(1)*cosGMST - F(2)*sinGMST
  #FIn(2) = F(1)*sinGMST + F(2)*cosGMST
    a20=a22
    a = [a20[1]*cos(gmst)-a20[2]*sin(gmst), a20[1]*sin(gmst)+a20[2]*cos(gmst),a20[3]]
    #print(a)
    return a
end
#=
eop_IAU1980 = get_iers_eop_iau_2000A()
jd = 2456727.500 
r = [-2.402809058043280E+03,  18.537777890343411E+02, -6.682371096072878E+03]
v = [-4.951358974727093E+00,  5.020450518834991E+00,  2.422381912322207E+00]
r_norm = norm(r)
Re = 6378.136460
μe = 3.986004414498200E+05
Geogrid(r, jd, eop_IAU1980)
a_J2 = 3/2*J2*μe*Re^2/r_norm^4/r_norm*[r[1]*(5*r[3]^2/r_norm^2-1), r[2]*(5*r[3]^2/r_norm^2-1), r[3]*(5*r[3]^2/r_norm^2-3)]
=#
#×