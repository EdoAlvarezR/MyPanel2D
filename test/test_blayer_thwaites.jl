"Compares implementation of Thwaites' method of boundary layer with Blasius
solution on a flat plate"
function test_blayer_thwiates(; verbose=true, show_plots=false, npoints=500)

  Re = 1e6       # Reynolds number based on length of plate
  L = 1           # Length of plate
  mu = 1          # Dynamic viscosity
  rho = 1         # Density
  nu = mu/rho     # Kinematic viscosity
  U = Re*nu/L

  # npoints = 250
  points = [[xi, 0] for xi in linspace(0,L,npoints)]
  xs = [p[1] for p in points]
  Us = U*ones(npoints)

  # Thwaites method
  lambda, theta, H, tau_w = p2d.thwaites(points, Us, mu, nu; theta0="flat")
  deltastar = H.*theta
  cf = tau_w / (0.5*rho*U^2)
  Cf = sum( cf[2:end].*[xs[i]-xs[i-1] for i in 2:npoints] )

  # Blasius solution
  bl_npoints = 1000
  bl_xs = [xi for xi in linspace(0,L,bl_npoints) ]
  Rex = U*bl_xs/nu
  bl_deltastar = 1.7208 * bl_xs ./ sqrt.(Rex)
  bl_theta = 0.664 * bl_xs ./ sqrt.(Rex)
  bl_H = bl_deltastar./bl_theta
  bl_cf = 0.664 ./ sqrt.(Rex)
  bl_Cf = sum( bl_cf[2:end].*[bl_xs[i]-bl_xs[i-1] for i in 2:bl_npoints] )

  # Plots
  if show_plots
    fig = figure("thwaites_flat", figsize=(7*3,5*1))
    subplot(131)
    plot(bl_xs, bl_deltastar, "--k", label=L"$\delta^*$ Blasius")
    plot(xs, deltastar, "or", label=L"$\delta^*$ Thwaites")
    plot(bl_xs, bl_theta, "-.k", label=L"$\theta$ Blasius")
    plot(xs, theta, "ob", label=L"$\theta$ Thwaites")
    xlabel(L"x")
    ylabel(L"y")
    grid(true, color="0.8", linestyle="--")
    legend(loc="best")
    title("Boundary layer growth")

    subplot(132)
    plot(bl_xs, bl_H, "-k", label=L"$H$ Blasius")
    plot(xs, H, "or", label=L"$H$ Thwaites")
    ylim([mean(H)*(1-0.1), mean(H)*(1+0.1)])
    xlabel(L"x")
    ylabel(L"H")
    grid(true, color="0.8", linestyle="--")
    legend(loc="best")
    title("Shape factor")

    subplot(133)
    plot(bl_xs, bl_cf, "--k", label=L"$c_f$ Blasius")
    plot(xs, cf, "or", label=L"$c_f$ Thwaites")
    xlabel(L"x")
    ylabel(L"c_f")
    legend(loc="best")
    grid(true, color="0.8", linestyle="--")
    title("Skin friction")
  end

  # Test criteria
  err = (bl_Cf-Cf)/bl_Cf
  flag = abs(err)<1e-3

  if verbose
    println("Blasius total skin friction Cf:\t\t$bl_Cf")
    println("Thwaites total skin friction Cf:\t$Cf")
    println("Error:\t\t\t$(err*100)\%")
    if flag
      println("Passed!")
    else
      println("*****FAILED!****")
    end
  end

  return flag
end
