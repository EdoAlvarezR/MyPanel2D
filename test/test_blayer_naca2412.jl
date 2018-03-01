"Test boundary layer laminar-turbulence transition on NACA2412 airfoil"
function test_blayer_naca2412(; verbose=true, show_plots=false,
                              npoints::Signed=150, debug=false, AOA=5.0)

  plot_title = "Boundary layer visualization on NACA2412, AOA=$AOA"

  Re = 1e6        # Reynolds number based on chord length
  c = 1           # Chord length
  mu = 1          # Dynamic viscosity
  rho = 1         # Density
  nu = mu/rho     # Kinematic viscosity
  U = Re*nu/c     # Freestream velocity

  # Freestream
  # AOA = 5.0                       # (deg) angle of attack
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
  theta, H, cf = p2d.calc_blayer(body, mu, nu; debug=debug, show_plot=true,
                                        plot_title=plot_title, verbose=verbose)

  if verbose; println("Done."); end;
end
