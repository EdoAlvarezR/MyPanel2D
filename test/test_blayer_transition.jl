"Tests implementation of Michael's boundary-layer transition method"
function test_blayer_transition(; verbose=true, show_plots=false, npoints=500)

  # Re = 1e6       # Reynolds number based on length of plate
  Re = 5e6       # Reynolds number based on length of plate
  L = 1           # Length of plate
  mu = 1          # Dynamic viscosity
  rho = 1         # Density
  nu = mu/rho     # Kinematic viscosity
  U = Re*nu/L

  points = [[xi, 0] for xi in linspace(0,L,npoints)]
  xs = [p[1] for p in points]
  Us = U*ones(npoints)

  aux1 = []
  theta, H, cf = p2d.calc_blayer(points, Us, mu, nu; theta0="flat",
                                                      more_outputs=aux1)
  sep_i, tran_i, lambda, tran_crit = aux1

  # println("\tsep_i=$sep_i")
  # println("\ttran_i=$tran_i")
  # println("\tlambda=$lambda")
  # println("\ttran_crit=$tran_crit")

  # Creates the boundary layer geometry
  theta_points = []               # Momentum thickness
  deltastar_points = []           # Displacement thickness
  prev_t, prev_n = p2d.get_tn(points[1], points[2])
  for i in 1:npoints
    deltastar = H[i]*theta[i]

    if i!=npoints
      t, n = p2d.get_tn(points[i], points[i+1])
    else
      t, n = prev_t, prev_n
    end

    push!(theta_points, points[i] + theta[i]*(prev_n+n)/2)
    push!(deltastar_points, points[i] + deltastar*(prev_n+n)/2)

    prev_t, prev_n = t, n
  end

  # Plots geometry and boundary layer properties
  plot([p[1] for p in points], [p[2] for p in points], "-k", label="Body")
  plot([p[1] for p in deltastar_points], [p[2] for p in deltastar_points],
                            "-.b", label=L"Displacement thickness $\delta^*$")
  plot([p[1] for p in theta_points], [p[2] for p in theta_points],
                                "--g", label=L"Momentum thickness $\theta$")
  xlabel("x")
  ylabel("y")
  legend(loc="best")
  title("Boundary layer visualization on flat plate")

end
