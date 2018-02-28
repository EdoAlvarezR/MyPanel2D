"Test boundary layer laminar-turbulence transition on NACA2412 airfoil"
function test_blayer_naca2412(; verbose=true, show_plots=false,
                              npoints::Signed=150, debug=false)

  Re = 1e6        # Reynolds number based on chord length
  c = 1           # Chord length
  mu = 1          # Dynamic viscosity
  rho = 1         # Density
  nu = mu/rho     # Kinematic viscosity
  U = Re*nu/c     # Freestream velocity

  # Freestream
  AOA = 5.0                       # (deg) angle of attack
  magVinf = U                     # (m/s) freestream velocity
  Vinf(X, t) = magVinf*[cos(AOA*pi/180), sin(AOA*pi/180)]

  # Airfoil geometry
  if verbose; println("Defining geometry..."); end;
  A, B, CC = 2, 4, 12
  x, y = p2d.generate_naca4(A, B, CC, npoints; sharp=true)

  # Creates the paneled body
  points = hcat(x, y)
  types = ["source", "vortex"]
  body = p2d.Body(points, types)

  # Solves the paneled body
  if verbose; println("Solving paneled body..."); end;
  p2d.solve(body, Vinf)

  # Calculates boundary layer
  if verbose; println("Calculating boundary layer..."); end;
  more_outputs = []
  theta, H, cf = p2d.calc_blayer(body, mu, nu; debug=debug,
                                                  more_outputs=more_outputs)
  stg_i, low_points, up_points = more_outputs[1:3]
  low_sep_i, low_tran_i, low_th_lambda, low_tran_crit = more_outputs[4:7]
  up_sep_i, up_tran_i, up_th_lambda, up_tran_crit = more_outputs[8:11]


  # Plots geometry and boundary layer properties
  if verbose; println("Plotting and generating good looking stuff..."); end;
  CPs = [ p2d.get_panel(body, i)[4] for i in 1:body.n ] # Control points
  theta_points = []                                     # Momentum thickness
  deltastar_points = []                                 # Displacement thickness
  for i in 1:body.n
    deltastar = H[i]*theta[i]
    t, n = p2d.get_tn(body, i)
    push!(theta_points, CPs[i] + 10*theta[i]*n)
    push!(deltastar_points, CPs[i] + 10*deltastar*n)
  end

  plot([p[1] for p in CPs], [p[2] for p in CPs], "-k", label="Body")
  plot([p[1] for p in deltastar_points], [p[2] for p in deltastar_points],
                            "-.b", label=L"Displacement thickness $10\times\delta^*$")
  plot([p[1] for p in theta_points], [p[2] for p in theta_points],
                                "--g", label=L"Momentum thickness $10\times\theta$")
  plot([up_points[1][1]], [up_points[1][2]], "*r", label="Stagnation point")


  if verbose; println("\tStagnation point X=$(up_points[1])"); end;
  if low_sep_i!=nothing
    X = low_points[low_sep_i]
    plot([X[1]], [X[2]], "^y")
    if verbose
      println("\tLower surface flow separation at"*
          " X=$(X) (lambda=$(low_th_lambda[low_sep_i]))")
    end
  end
  if up_sep_i!=nothing
    X = up_points[up_sep_i]
    plot([X[1]], [X[2]], "^y", label="Separation point")
    if verbose
      println("\tUpper surface flow separation at"*
          " X=$(X) (lambda=$(up_th_lambda[up_sep_i]))")
    end
  end
  if low_tran_i!=nothing
    X = low_points[low_tran_i]
    plot([X[1]], [X[2]], "oy")
    if verbose
      println("\tLower surface turbulent transition at"*
      " X=$(X) (tran_crit=$(low_tran_crit[low_tran_i]))")
    end
  end
  if up_tran_i!=nothing
    X = up_points[up_tran_i]
    plot([X[1]], [X[2]], "oy", label="Transition point")
    if verbose
      println("\tUpper surface turbulent transition at"*
        " X=$(X)  (tran_crit=$(up_tran_crit[up_tran_i]))")
    end
  end

  xlabel("x")
  ylabel("y")
  y_lims = maximum([p[1] for p in CPs])/2
  ylim([-y_lims, y_lims])
  legend(loc="best")
  title("Boundary layer visualization on NACA2412")
end
