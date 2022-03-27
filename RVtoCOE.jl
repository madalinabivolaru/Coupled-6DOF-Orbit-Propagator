function RVtoCOE(r,v)
    #Function changing state vector [r, v] to osculating orbital elements
    μ = 398600.4418;  
    h = cross(r,v) 
    h_mod = norm(h)
    r_mod = norm(r) 
    v_mod = norm(v)
    n = cross([0.0,0.0,1.0], h)
    n_mod = norm(n)
    e = ((v_mod^2-μ/r_mod)*r-(r ⋅ v)*v)/μ
    e_mod = norm(e)
    ξ = 0.5*v_mod^2 - μ/r_mod

    if e != 1
        a = -0.5*μ/ξ
        p = a*(1 - e_mod^2)
    else
        p = h^2/μ 
        a = "infinity"
    end

    i = acosd(h[3]/h_mod)
    ω = acosd(n ⋅ e/(n_mod*e_mod))

    if e[3]<0
        ω = 360-ω
    end

    Ω = acosd(n[1]/ n_mod)

    if n[2] < 0 
        Ω = 360 - Ω
    end

    ν = acosd(e ⋅ r/(e_mod*r_mod))

    if r ⋅ v < 0
        ν = 360 - ν
    end
        
    return a,e_mod,i,Ω,ω,ν
end






    
