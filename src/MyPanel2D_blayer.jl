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
  * `Vinf::Array{Float64, 1}`               : Freestream for nondimensionalizing
                                               and direction of drag.

  **Optional Arguments**
  * `theta0::Any`       : Indicates initial condition at leading edge. Give it
                              "flat" for a flat plate, "blunt" for an airfoil,
                              or a number for manually setting the initial
                              condition.

Returns `(theta, H, cf, CDs)`, with `theta` momentum thickness, `H` shape factor,
`cf` friction coefficient at the wall, and `CDs=[CD, CDf, CDp]` total viscous,
friction, and pressure drag.
"""
function calc_blayer(body::Body, mu::Real, nu::Real,
                      Vinf::Array{T,1} where T<:Real;
                      theta0::Any="blunt", more_outputs=nothing, t::Real=0,
                      show_plot::Bool=false, plot_title::String="Boundary layer",
                      debug=false, verbose=false)
  if !body.solved
    error("Body's panels must be solved before calculating boundary layer."*
            "Call `solve()` before calling this function.")
  end

  # --------- PREPARATION ------------------------------------------------------
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

  # More outputs and parameters for plot
  if more_outputs!=nothing || show_plot
    _more_outputs=[]
  else
    _more_outputs=nothing
  end
  if _more_outputs!=nothing
    push!(_more_outputs, stg_i)
    push!(_more_outputs, low_points)
    push!(_more_outputs, up_points)
  end

  # --------- SOLVES BOUNDARY LAYER --------------------------------------------
  # Solves lower and upper surfaces separately
  low_theta, low_H, low_cf = calc_blayer(low_points, low_vels, mu, nu;
                          theta0=theta0, debug=debug, more_outputs=_more_outputs)
  up_theta, up_H, up_cf = calc_blayer(up_points, up_vels, mu, nu;
                          theta0=theta0, debug=debug, more_outputs=_more_outputs)

  # Put them back together
  theta = vcat(reverse(low_theta), up_theta[2:end])
  H = vcat(reverse(low_H), up_H[2:end])
  cf = vcat(reverse(low_cf), up_cf[2:end])


  # Checks for flow separation
  low_sep_i = _more_outputs[4]        # Index of lower point of separation
  up_sep_i = _more_outputs[8]         # Index of upper point of separation
  low_sep_i = low_sep_i==nothing ? size(low_points,1) : low_sep_i-1^(low_sep_i!=1)
  up_sep_i = up_sep_i==nothing ? size(up_points,1) : up_sep_i-1^(up_sep_i!=1)

  # Calculates total viscous drag using Squire-Young's formula
  magVinf = norm(Vinf)                # Freestream velocity
  dirVinf = Vinf/magVinf              # Freestream direction
  c = get_c(body)                     # Chord length
  CD = 2*up_theta[up_sep_i]/c * (
                            up_vels[up_sep_i]/magVinf)^( (up_H[up_sep_i]+5)/2 )
  CD += 2*low_theta[low_sep_i]/c * (
                        low_vels[low_sep_i]/magVinf)^( (low_H[low_sep_i]+5)/2 )
  CD = CD*dirVinf

  if debug
    println("up_theta=$(up_theta[up_sep_i])")
    println("low_theta=$(low_theta[low_sep_i])")
    println("up_H=$(up_H[up_sep_i])")
    println("low_H=$(low_H[low_sep_i])")
    println("up_vel=$(up_vels[up_sep_i])")
    println("low_vel=$(low_vels[low_sep_i])")
    println("magVinf=$magVinf")
  end

  # Calculates friction drag
  CDf = sum([(up_points[i]-up_points[i-1])*(up_cf[i]+up_cf[i-1])/2
                                                for i in 2:up_sep_i])
  CDf += sum([(low_points[i]-low_points[i-1])*(low_cf[i]+low_cf[i-1])/2
                                                for i in 2:low_sep_i])
  # (sanity check and gets rid of imaginary part)
  check_crit = norm(imag.(CDf))/norm(real.(CDf))
  if check_crit>1e-1
    error("Sanity check on CDf failed! (Imaginary part is too large"*
            ", crit=$check_crit) \n CDf=$CDf")
  end
  CDf = real.(CDf)

  # Calculates pressure drag
  CDp = CD - dot(CDf,dirVinf)*dirVinf

  # --------- EXTRA STUFF ------------------------------------------------------
  # Plot
  if show_plot
    _plot_blayer(body, theta, H, _more_outputs;
                              verbose=verbose, str_title=plot_title)
  end

  # More outputs
  if more_outputs!=nothing
    for elem in _more_outputs
      push!(more_outputs, elem)
    end
  end

  return theta, H, cf, [CD, CDf, CDp]
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

  # Puts together laminar and turbulent solutions
  if tran_i!=nothing  # Case of turbulent transition
    # H0 = th_H[tran_i]
    H0 = 1.4  # Resets the shape factor to flat plate initial condition
    hd_theta, hd_H, hd_cf = head(points[tran_i:end], vels[tran_i:end], nu,
                                th_theta[tran_i]; H0=H0, debug=debug)

    theta = vcat(th_theta[1:tran_i-1], hd_theta)
    H = vcat(th_H[1:tran_i-1], hd_H)
    cf = vcat(th_cf[1:tran_i-1], hd_cf)

  else # Case of fully laminar
    theta = th_theta
    H = th_H
    cf = th_cf

  end


  sep_i = nothing                       # Index of separation point
  sep_crit = 3.518333702400003          # Separation criterium, lambda=-0.09
  H_fl = [Float64(Hi) for Hi in H]      # Converts Real array to Float array
                                        #  in order to use NaNMath package
  # Checks for flow separation
  if NaNMath.maximum(H_fl) >= sep_crit

    # Finds separation point
    for i in 1:npoints
      if H[i] >= sep_crit && pos_x[i]/pos_x[end]>0.20
        sep_i = i

        warn("Flow separation detected at X=$(points[sep_i]) (H=$(H[sep_i]))."*
              "No flow separation method is currently implement, calculations will"*
              "proceed ignoring separation.")
        break
      end
    end
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

function _plot_blayer(body::Body, theta, H, more_outputs;
                          verbose=false, str_title::String="Boundary layer")

  stg_i, low_points, up_points = more_outputs[1:3]
  low_sep_i, low_tran_i, low_th_lambda, low_tran_crit = more_outputs[4:7]
  up_sep_i, up_tran_i, up_th_lambda, up_tran_crit = more_outputs[8:11]


  # Plots geometry and boundary layer properties
  CPs = [ get_panel(body, i)[4] for i in 1:body.n ] # Control points
  theta_points = []                                 # Momentum thickness
  deltastar_points = []                             # Displacement thickness
  for i in 1:body.n
    deltastar = H[i]*theta[i]
    t, n = get_tn(body, i)
    push!(theta_points, CPs[i] + 20*theta[i]*n)
    push!(deltastar_points, CPs[i] + 20*deltastar*n)
  end

  fig = figure("blayer")
  y_lims = norm(CPs[stg_i]-CPs[1])/2
  xlim([-y_lims/4, y_lims*2*9/8])
  ylim([-y_lims, y_lims])
  plot([p[1] for p in CPs], [p[2] for p in CPs], "-k", label="Body")
  plot([p[1] for p in deltastar_points], [p[2] for p in deltastar_points],
  "-.b", label=L"Displacement thickness $20\times\delta^*$")
  plot([p[1] for p in theta_points], [p[2] for p in theta_points],
  "--g", label=L"Momentum thickness $20\times\theta$")
  plot([up_points[1][1]], [up_points[1][2]], "*c", label="Stagnation point")


  if verbose; println("\tStagnation point X=$(round.(up_points[1],2))"); end;
  if low_tran_i!=nothing
    X = low_points[low_tran_i]
    plot([X[1]], [X[2]], "oy")
    if verbose
      println("\tLower surface turbulent transition at"*
      " X=$(round.(X,2)) (tran_crit=$(round(low_tran_crit[low_tran_i],1)))")
    end
  end
  if up_tran_i!=nothing
    X = up_points[up_tran_i]
    plot([X[1]], [X[2]], "oy", label="Transition point")
    if verbose
      println("\tUpper surface turbulent transition at"*
      " X=$(round.(X,2))  (tran_crit=$(round(up_tran_crit[up_tran_i],1)))")
    end
  end
  if low_sep_i!=nothing
    X = low_points[low_sep_i]
    plot([X[1]], [X[2]], "dr")
    if verbose
      println("\tLower surface flow separation at"*
      " X=$(round.(X,2)) (H=$(round(H[length(low_th_lambda)-low_sep_i+1],1)))")
    end
  end
  if up_sep_i!=nothing
    X = up_points[up_sep_i]
    plot([X[1]], [X[2]], "dr", label="Separation point")
    if verbose
      println("\tUpper surface flow separation at"*
      " X=$(round.(X,2)) (H=$(round(H[length(low_th_lambda)+up_sep_i-1],1)))")
    end
  end

  xlabel("x")
  ylabel("y")
  legend(loc="best")
  title(str_title)
end











#
