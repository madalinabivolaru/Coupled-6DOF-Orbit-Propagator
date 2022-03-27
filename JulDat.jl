
function JulDat(yr::Float64,mo::Float64,d::Float64,h::Float64,min::Float64,s::Float64)
    #Julian date calculator defined as per Vallado for years 1900 to 2100 

    JD = 367.0*yr - floor(7/4*(yr+floor((mo+9)/12)))+
    floor(275.0*mo/9)+d+1721013.5+((s/60 + min)/60 + h)/24;
    
    return JD
end

    
