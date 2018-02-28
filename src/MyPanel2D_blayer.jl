#=##############################################################################
# DESCRIPTION
    Methods for calculating boundary layer and viscous effects.
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Feb 2018
  * License   : MIT License
=###############################################################################

"""
  `calc_blayer(body, mu, nu; theta0="blunt")`

Calculates properties of the boundary layer at each control point using
Thwaites' method for the laminar region and Head's method for the turbulent
region, identifying transition through Michel's method (incompressibility is
assumed throughout).

  **Arguments**
  * `body::Body`                            : Paneled body.
  * `mu::Float64`                           : Dynamic viscosity.
  * `nu::Float64`                           : Kinematic viscosity.

  **Optional Arguments**
  * `theta0::Any`       : Indicates initial condition at leading edge. Give it
                              "flat" for a flat plate, "blunt" for an airfoil,
                              or a number for manually setting the initial
                              condition.

Returns `(theta, H, cf)`, with `theta` momentum thickness, `H` shape factor,
`cf` friction coefficient at the wall.
"""
function calc_blayer(body::Body, mu::Real, nu::Real; theta0::Any="blunt",
                                  more_outputs=nothing, t::Real=0, debug=false)
  if !body.solved
    error("Body's panels must be solved before calculating boundary layer."*
            "Call `solve()` before calling this function.")
  end

  # Control points and velocities at each control point
  CPs = [ get_panel(body, i)[4] for i in 1:body.n ]
  vels = [ norm(Vind(body, CP) + body.Vinf(CP, t)) for CP in CPs ]

  # Stagnation point
  _, stg_i = findmin(vels)

  # Splits the surface according to stagnation point
  low_points = CPs[1:stg_i]; reverse!(low_points);
  low_vels = vels[1:stg_i]; reverse!(low_vels);
  up_points = CPs[stg_i:end]
  up_vels = vels[stg_i:end]

  if more_outputs!=nothing
    push!(more_outputs, stg_i)
    push!(more_outputs, low_points)
    push!(more_outputs, up_points)
  end

  # Solves lower and upper surfaces separately
  low_theta, low_H, low_cf = calc_blayer(low_points, low_vels, mu, nu;
                          theta0=theta0, debug=debug, more_outputs=more_outputs)
  up_theta, up_H, up_cf = calc_blayer(up_points, up_vels, mu, nu;
                          theta0=theta0, debug=debug, more_outputs=more_outputs)

  # Put them back together
  theta = vcat(reverse(low_theta), up_theta[2:end])
  H = vcat(reverse(low_H), up_H[2:end])
  cf = vcat(reverse(low_cf), up_cf[2:end])


  return theta, H, cf
end


"""
  `calc_blayer(points, vels, mu, nu; theta0="blunt")`

Given a collection of points from leading to trailing edge and their respective
velocities at the edge of the boundary layer, it calculates properties of the
boundary layer at each point using Thwaites' method for the laminar region and
Head's method for the turbulent region, identifying transition through Michel's
method (incompressibility is assumed throughout).

  **Arguments**
  * `points::Array{Array{Float64, 1}, 1}`   : Collection of points.
  * `vels::Array{Float64, 1}`               : Magnitude of boundary layer edge
                                              velocity at each point.
  * `mu::Float64`                           : Dynamic viscosity.
  * `nu::Float64`                           : Kinematic viscosity.

  **Optional Arguments**
  * `theta0::Any`       : Indicates initial condition at leading edge. Give it
                              "flat" for a flat plate, "blunt" for an airfoil,
                              or a number for manually setting the initial
                              condition.

Returns `(theta, H, cf)`, with `theta` momentum thickness, `H` shape factor,
`cf` friction coefficient at the wall.
"""
function calc_blayer(points::Array{T, 1} where {T<:AbstractArray},
                      vels::Array{T, 1} where {T<:Real}, mu::Real, nu::Real;
                      theta0::Any="blunt", more_outputs=nothing, debug=false)

  npoints = size(points, 1)
  rho = mu/nu

  # Solves it as fully laminar (Thwaites method)
  th_lambda, th_theta, th_H, th_tau_w = thwaites(points, vels, mu, nu; theta0=theta0)
  th_cf = th_tau_w./(0.5*rho*vels.^2)

  # Checks for flow separation
  sep_i = nothing                       # Index of separation point
  if NaNMath.minimum(th_lambda) <= -0.09

    # Finds separation point
    for i in 1:npoints
      if th_lambda[i] <= -0.09
        sep_i = i
        break
      end
    end

    warn("Flow separation detected at X=$(points[sep_i]) (lambda=$(th_lambda[sep_i]))."*
          "No flow separation method is currently implement, calculations will"*
          "proceed ignoring separation.")
  end

  # Checks for transition (Michel's method)
  tran_i = nothing                      # Index of transition point
                                        # Distance between points
  dist_x = vcat(0,[norm(points[i]-points[i-1]) for i in 2:npoints])
  pos_x = [sum(dist_x[1:i]) for i in 1:npoints]   # Position along surface
  Rex = vels.*pos_x/nu                  # Position-based Reynolds number
  Retheta = vels.*th_theta/nu           # Momentum-thickness-based Reynolds numb
  tran_crit = 1.174*(1+22400./Rex).*(Rex.^0.46)   # Transition criteria
  tran_crit = Retheta-tran_crit

  if NaNMath.maximum(tran_crit)>=0
    # Finds transition point
    for i in 1:npoints
      if tran_crit[i]>=0
        tran_i = i
        break
      end
    end
  end

  if tran_i!=nothing
    hd_theta, hd_H, hd_cf = head(points[tran_i:end], vels[tran_i:end], nu,
                                th_theta[tran_i]; H0=th_H[tran_i], debug=debug)

    theta = vcat(th_theta[1:tran_i-1], hd_theta)
    H = vcat(th_H[1:tran_i-1], hd_H)
    cf = vcat(th_cf[1:tran_i-1], hd_cf)

  else
    theta = th_theta
    H = th_H
    cf = th_cf

  end

  if more_outputs!=nothing
    push!(more_outputs, sep_i)
    push!(more_outputs, tran_i)
    push!(more_outputs, th_lambda)
    push!(more_outputs, tran_crit)
  end

  return theta, H, cf
end

"""
  `thwaites(points, vels, mu, nu; theta0="blunt")`

Given a collection of points from leading to trailing edge and their respective
velocities at the edge of the boundary layer, it calculates properties of the
boundary layer at each point using Thwaites method (assumes laminar,
incompressible flow).

  **Arguments**
  * `points::Array{Array{Float64, 1}, 1}`   : Collection of points.
  * `vels::Array{Float64, 1}`               : Magnitude of boundary layer edge
                                              velocity at each point.
  * `mu::Float64`                           : Dynamic viscosity.
  * `nu::Float64`                           : Kinematic viscosity.

  **Optional Arguments**
  * `theta0::Any`       : Indicates initial condition at leading edge. Give it
                              "flat" for a flat plate, "blunt" for an airfoil,
                              or a number for manually setting the initial
                              condition.

Returns `(lambda, theta, H, tau_w)`, with `lambda` Thwaites' non-dimensional
parameter, `theta` momentum thickness, `H` shape factor, `tau_w` shear stress
at the wall.
"""
function thwaites(points::Array{T, 1} where {T<:AbstractArray},
                  vels::Array{T, 1} where {T<:Real}, mu::Real, nu::Real;
                  theta0::Any="blunt")

  # Initial condition (leading edge)
  if typeof(theta0)==String
    if theta0=="flat"
      _theta0 = 0
    elseif theta0=="blunt"
      dVdx0 = (vels[2]-vels[1]) / norm(points[2]-points[1])
      _theta0 = sqrt(0.075*nu/dVdx0)
      if dVdx0==0
        error("dVdx0 must be different than zero for blunt mode."*
                " Try flat plate mode instead.")
      end
    else
      error("Invalid initial condition $theta0.")
    end
  else
    _theta0 = theta0
  end

  npoints = size(points, 1)
  if size(vels,1) != npoints
    error("Number of points and velocities must match!")
  end

  dx = [norm(points[i]-points[i-1]) for i in 2:npoints]   # Steps

  # Calculates the theta factor at each point
  theta2V6 = zeros(Real, npoints)
  theta2V6[1] = _theta0^2*vels[1]^6
  for i in 2:npoints

    dtheta2V6 = dx[i-1]*0.45*nu*vels[i]^5
    theta2V6[i] = theta2V6[i-1] + dtheta2V6

  end

  # Converts it into theta
  theta = sqrt.(theta2V6./(vels.^6))

  # Converts momentum thickness theta into lambda
  dVdx = [(vels[i]-vels[i-1])/dx[i-1] for i in 2:npoints] # Velocity change
  dVdx = vcat(dVdx[1], dVdx)                              # Duplicates it at LE
  lambda = (theta.^2).*dVdx/nu

  # Shape factor
  z = 0.25 - lambda
  H = 2.0 + 4.14*z - 83.5*z.^2 + 854*z.^3 - 3337*z.^4 + 4576*z.^5

  # Shear stress at wall
  tau_w = (mu*vels./theta).*(lambda+0.09+0*im).^0.62

  return lambda, theta, H, tau_w
end



"""
  `head(points, vels, nu, theta0; H0=1.4)`

Given a collection of points from leading to trailing edge and their respective
velocities at the edge of the boundary layer, it calculates properties of the
boundary layer at each point using Head method (assumes turbulent,
incompressible flow).

  **Arguments**
  * `points::Array{Array{Float64, 1}, 1}`   : Collection of points.
  * `vels::Array{Float64, 1}`               : Magnitude of boundary layer edge
                                              velocity at each point.
  * `nu::Float64`                           : Kinematic viscosity.
  * `theta0::Float64`                       : Momentum thickness at first point.

  **Optional Arguments**
  * `H0::Float64`             : Indicates shape factor at first point. The default
                              is the value of a flat plate's leading edge, which
                              should be good enough for an airfoil.

Returns `(theta, H, cf)`, with `theta` momentum thickness, `H` shape factor,
`cf` friction coefficient at the wall.
"""
function head(points::Array{T, 1} where {T<:AbstractArray},
                  vels::Array{T, 1} where {T<:Real}, nu::Real, theta0::Real;
                  H0::Real=1.4, debug=false)

  npoints = size(points, 1)
  if size(vels,1) != npoints
    error("Number of points and velocities must match!")
  end

  dx = [norm(points[i]-points[i-1]) for i in 2:npoints]   # Steps
  dVdx = [(vels[i]-vels[i-1])/dx[i-1] for i in 2:npoints] # Velocity change

  theta = zeros(Real, npoints)        # Momentum thickness
  H = zeros(Real, npoints)            # Shape factor
  cf = zeros(Real, npoints)           # Friction coefficient

  calc_cf(this_H, Retheta) = 0.246 * 10^(-0.678*this_H) * Retheta^(-0.268)

  # Initial conditions
  theta[1] = theta0
  H[1] = H0
  cf[1] = calc_cf(H[1], vels[1]*theta[1]/nu)

  if debug
    println("theta0=$(theta[1])")
    println("H0=$(H[1])")
    println("cf0=$(cf[1])")
  end

  for i in 2:npoints

    dthetadx = cf[i-1]/2 - theta[i-1]*(H[i-1]+2)*dVdx[i-1]/vels[i-1]
    H1 = _H2H1(H[i-1])
    fH1 = 0.0306*(H1-3)^(-0.6169)
    dH1dx = fH1/theta[i-1] - H1/vels[i-1]*dVdx[i-1] - H1/theta[i-1]*dthetadx

    if debug
      println("i=$i")
      println("\tdthetadx=$dthetadx")
      println("\t\tcf[i-1]/2=$(cf[i-1]/2)")
      println("\t\ttheta[i-1]*(H[i-1]+2)*dVdx[i-1]/vels[i-1]=$(theta[i-1]*(H[i-1]+2)*dVdx[i-1]/vels[i-1])")
      println("\tprev_H1=$H1")
      println("\tfH1=$fH1")
      println("\tdH1dx=$dH1dx")
      println("\t\tfH1/theta[i-1]=$(fH1/theta[i-1])")
      println("\t\tH1/vels[i-1]*dVdx[i-1]=$(H1/vels[i-1]*dVdx[i-1])")
      println("\t\tH1/theta[i-1]*dthetadx=$(H1/theta[i-1]*dthetadx)")
      println("\tdx=$(dx[i-1])")
      println("\tnew_H1=$(H1 + dH1dx*dx[i-1])")
    end

    theta[i] = theta[i-1] + dthetadx*dx[i-1]
    H[i] = real(_H12H(H1 + dH1dx*dx[i-1]))
    cf[i] = calc_cf(H[i], vels[i]*theta[i]/nu)

  end

  return theta, H, cf
end



# ------ INTERNAL FUNCTIONS ----------------------------------------------------
function _H12H(H1)
  if H1>=5.3
    return 1.1 + 0.86*(H1-3.3+0*im)^(-0.777)
  else
    return 0.6778 + 1.1538*(H1-3.3+0*im)^(-0.326)
  end
end

function _H2H1(H)
  if H<=1.6
    return 0.8234*(H-1.1)^(-1.287) + 3.3
  else
    return 1.5501*(H-0.6778)^(-3.064) + 3.3
  end
end












#
