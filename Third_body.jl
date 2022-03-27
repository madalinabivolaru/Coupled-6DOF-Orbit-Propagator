#Computing Third body acceleration
function Rdbody(rsat::Vector{Float64},rmoon::Vector{Float64},rsun::Vector{Float64})
    rsat3m = rmoon - rsat
    rsat3s = rsun - rsat
    a3m = PropagatorConstants.μ_m * (rsat3m/norm(rsat3m)^3 - rmoon/norm(rmoon)^3)
    a3s = PropagatorConstants.μ_s * (rsat3s/norm(rsat3s)^3 - rsun/norm(rsun)^3)
    return a3m + a3s
end
