function asrp(rsun, rsat)
    #Computing solar pressure with spherical body assumption
    #=if Wertz == true
        D_aphe = (JD - JulDat(2021.0,7.0,4.0,0.0,0.0,0.0))/365
        SF =  1358/(1.004+0.0334*cos(2*pi*D_aphe))
    else 
        SF =1367
    end=#
    
    SF = 1367
    p_srp = SF/3*10^-8;
    Cr = 1.2;
  
    As = 0.7;
    m = 400;
    rsats = rsat - rsun
    #asrp = 10^(-3)*(p_srp*Cr*As/m)*rsats/norm(rsats)   #positive as rsunsat is used 

    asrp = -10^(-3)*(p_srp*Cr*As/m)*rsun/norm(rsun)^3*PropagatorConstants.AUtoKM^2

    #Check umbra penumbra instances!
    α_umb = 0.264121687
    α_pen = 0.269007205
    
    cos_ξ = dot(rsun, rsat)/norm(rsun)/norm(rsat)
    sat_h =  norm(rsat)*cos_ξ
    sat_v = norm(rsat)*sqrt(abs(1-cos_ξ^2))
    x = PropagatorConstants.R_e/sind(α_pen)
    pen_v = tand(α_pen)*(x + sat_h)
   #= if sat_v<=pen_v
        y = PropagatorConstants.R_e/sind(α_umb)
        umb_v = tand(α_umb)*(y-sat_h)
        if sat_v<=umb_v
            asrp = [0.0,0.0,0.0]
       # else
            #asrp = asrp.*((sat_v-umb_v)/(pen_v-umb_v))
        end
    end =#
    return asrp
end
include("./PropagatorConstants.jl")
asrp([ 1.3219e8,-0.6101e8, -0.2645e8],[1315.0823093944996, -6765.511976126983, 4000.0])