"Test boundary layer laminar-turbulence transition on NACA2412 airfoil"
function test_blayer_naca2412(; verbose=true, show_plots=false,
                              npoints::Signed=150, debug=false, AOA=5.0,
                              plot_blayer=false)

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
  more_outputs = []
  theta, H, cf, CDs = p2d.calc_blayer(body, mu, nu, Vinf(0,0); debug=debug,
                        show_plot=true, plot_title=plot_title, verbose=verbose,
                        more_outputs=more_outputs)
  cf = real.(cf)
  CD, CDf, CDp = CDs
  if verbose
    println("\tTotal viscous drag coeff (CD):\t$(round(norm(CD),4))")
    println("\tFriction drag coeff (CDf):\t$(round.(CDf,4))")
    println("\tPressure drag coeff (CDp):\t$(round(norm(CDp),4))")
  end

  # Plots boundary layer properties
  if plot_blayer
    stg_i, low_points, up_points = more_outputs[1:3]
    low_sep_i, low_tran_i, low_th_lambda, low_tran_crit = more_outputs[4:7]
    up_sep_i, up_tran_i, up_th_lambda, up_tran_crit = more_outputs[8:11]

    up_tran_i_f = up_tran_i!=nothing ? stg_i+up_tran_i-1 : nothing
    low_tran_i_f = low_tran_i!=nothing ? stg_i-(low_tran_i-1) : nothing
    up_sep_i_f = up_sep_i!=nothing ? stg_i+up_sep_i-1 : nothing
    low_sep_i_f = low_sep_i!=nothing ? stg_i-(low_sep_i-1) : nothing

                                    # Distances between points
    low_ds = [norm(low_points[i]-low_points[i-1]) for i in 2:size(low_points,1)]
    up_ds = [norm(up_points[i]-up_points[i-1]) for i in 2:size(up_points,1)]
                                    # Surface position of each point
    low_s = vcat(0, [sum(low_ds[1:i]) for i in 1:size(low_ds,1)])
    up_s = vcat(0, [sum(up_ds[1:i]) for i in 1:size(up_ds,1)])
    rev_low_s = reverse(low_s)

    deltastar = H.*theta

    fig = figure("blayer_props", figsize=(7*3, 5*1))
    subplot(131)
    plot(up_s, deltastar[stg_i:end], "-b",
                                  label=L"$\delta^*$ upper surface")
    plot(up_s, theta[stg_i:end], "-g",
                                  label=L"$\theta$ upper surface")
    plot(rev_low_s, deltastar[1:stg_i], "--b",
                                  label=L"$\delta^*$ lower surface")
    plot(rev_low_s, theta[1:stg_i], "--g",
                                  label=L"$\theta$ lower surface")
    if up_tran_i!=nothing
      plot(up_s[up_tran_i]*ones(1), deltastar[up_tran_i_f], "oy")
    end
    if low_tran_i!=nothing
      plot(low_s[low_tran_i]*ones(1), deltastar[low_tran_i_f], "oy")
    end
    if up_sep_i!=nothing
      plot(up_s[up_sep_i]*ones(1), deltastar[up_sep_i_f], "dr")
    end
    if low_sep_i!=nothing
      plot(low_s[low_sep_i]*ones(1), deltastar[low_sep_i_f], "dr")
    end
    xlim([0, max(up_s[end], low_s[end])])
    ylim([0, 0.015])
    xlabel(L"s")
    ylabel(L"y")
    grid(true, color="0.8", linestyle="--")
    legend(loc="best")
    title("Boundary layer growth")

    subplot(132)
    aux2 = [minimum(H), maximum(H)]
    plot(up_s, H[stg_i:end], "-k", label=L"$H$ upper surface")
    plot(rev_low_s, H[1:stg_i], "--k", label=L"$H$ lower surface")
    if up_tran_i!=nothing
      plot(up_s[up_tran_i]*ones(1), H[up_tran_i_f], "oy")
    end
    if low_tran_i!=nothing
      plot(low_s[low_tran_i]*ones(1), H[low_tran_i_f], "oy")
    end
    if up_sep_i!=nothing
      plot(up_s[up_sep_i]*ones(1), H[up_sep_i_f], "dr")
    end
    if low_sep_i!=nothing
      plot(low_s[low_sep_i]*ones(1), H[low_sep_i_f], "dr")
    end
    xlim([0, max(up_s[end], low_s[end])])
    ylim([0, maximum(H)])
    xlabel(L"s")
    ylabel(L"H")
    grid(true, color="0.8", linestyle="--")
    legend(loc="best")
    title("Shape factor")

    subplot(133)
    plot(up_s, cf[stg_i:end], "-k", label=L"$c_f$ upper surface")
    plot(rev_low_s, cf[1:stg_i], "--k", label=L"$c_f$ lower surface")
    if up_tran_i!=nothing
      plot(up_s[up_tran_i]*ones(1), cf[up_tran_i_f], "oy")
    end
    if low_tran_i!=nothing
      plot(low_s[low_tran_i]*ones(1), cf[low_tran_i_f], "oy")
    end
    if up_sep_i!=nothing
      plot(up_s[up_sep_i]*ones(1), cf[up_sep_i_f], "dr")
    end
    if low_sep_i!=nothing
      plot(low_s[low_sep_i]*ones(1), cf[low_sep_i_f], "dr")
    end
    xlim([0, max(up_s[end], low_s[end])])
    ylim([0, 0.05])
    xlabel(L"s")
    ylabel(L"c_f")
    grid(true, color="0.8", linestyle="--")
    legend(loc="best")
    title("Skin friction")
  end


  if verbose; println("Done."); end;
end
