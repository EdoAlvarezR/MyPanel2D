"Compares implementation of Head's method of boundary layer with Schlichting
empiricial relation on a flat plate"
function test_blayer_head(; verbose=true, show_plots=false,
                              npoints::Signed=500, npoints_plot::Signed=25,
                              debug=false, more_outputs=nothing)

  Re = 1e6       # Reynolds number based on length of plate
  L = 1           # Length of plate
  mu = 1          # Dynamic viscosity
  rho = 1         # Density
  nu = mu/rho     # Kinematic viscosity
  U = Re*nu/L

  # Schlichting relation
  sc_npoints = 1000
  sc_xs = [xi for xi in linspace(0,L,sc_npoints) ]
  Rex = U*sc_xs/nu
  sc_deltastar = 0.046*sc_xs./(Rex.^0.2)
  sc_theta = 0.036*sc_xs./(Rex.^0.2)
  sc_H = sc_deltastar./sc_theta
  sc_cf = 0.0592./(Rex.^0.2)
  # sc_Cf = sum( sc_cf[2:end].*[sc_xs[i]-sc_xs[i-1] for i in 2:sc_npoints] )
  sc_Cf = 0.074/Rex[end]^0.2

  # Head's method
  i0 = 10
  if debug; println("sc_cf0 = $(sc_cf[i0])"); end;
  x0 = sc_xs[i0]
  theta0 = sc_theta[i0]
  H0 = sc_H[i0]
  points = [[xi, 0] for xi in linspace(x0,L,npoints)]
  xs = [p[1] for p in points]
  Us = U*ones(npoints)
  theta, H, cf = p2d.head(points, Us, nu, theta0; H0=H0, debug=debug)
  deltastar = H.*theta
  Cf = sum( cf[2:end].*[xs[i]-xs[i-1] for i in 2:npoints] )

  # Plots
  if show_plots
    is = [i+1 for i in 0:npoints-1 if i%(ceil(npoints/npoints_plot))==0]
    xs_plot = [xs[i] for i in is]
    fig = figure("thwaites_flat", figsize=(7*3,5*1))
    subplot(131)
    plot(sc_xs, sc_deltastar, "--k", label=L"$\delta^*$ Schlichting")
    plot(xs_plot, [deltastar[i] for i in is], "or", label=L"$\delta^*$ Head")
    plot(sc_xs, sc_theta, "-.k", label=L"$\theta$ Schlichting")
    plot(xs_plot, [theta[i] for i in is], "ob", label=L"$\theta$ Head")
    xlabel(L"x")
    ylabel(L"y")
    legend(loc="best")
    grid(true, color="0.8", linestyle="--")
    title("Boundary layer growth")


    subplot(132)
    plot(sc_xs, sc_H, "-k", label=L"$H$ Schlichting")
    plot(xs_plot, [H[i] for i in is], "or", label=L"$H$ Head")
    ylim([mean(H)*(1-0.5), mean(H)*(1+0.5)])
    xlabel(L"x")
    ylabel(L"H")
    grid(true, color="0.8", linestyle="--")
    legend(loc="best")
    title("Shape factor")

    subplot(133)
    plot(sc_xs, sc_cf, "--k", label=L"$c_f$ Schlichting")
    plot(xs_plot, [cf[i] for i in is], "or", label=L"$c_f$ Head")
    xlabel(L"x")
    ylabel(L"c_f")
    legend(loc="best")
    grid(true, color="0.8", linestyle="--")
    title("Skin friction")
  end

  # Test criteria
  err = (sc_Cf-Cf)/sc_Cf
  flag = abs(err)<2e-1

  if verbose
    println("Schilchting total skin friction Cf:\t$sc_Cf")
    println("Head total skin friction Cf:\t\t$Cf")
    println("Error:\t\t\t\t$(err*100)\%")
    if flag
      println("Passed!")
    else
      println("*****FAILED!****")
    end
  end

  if more_outputs!=nothing
    push!(more_outputs, err)
  end

  return flag
end
