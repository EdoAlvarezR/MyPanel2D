"Test generating the Cl, Cm, and Cd polar curves of a NACA2412 airfoil"
function test_cl_naca2412(; verbose=true, show_plots=false,
                              npoints::Signed=150, debug=false, Re=1e6,
                              dalpha=1.0, show_ite_plot=true, polar=true,
                              check_stall=true, save_path=nothing, prompt=true)

  c = 1           # Chord length
  mu = 1          # Dynamic viscosity
  rho = 1         # Density
  nu = mu/rho     # Kinematic viscosity
  U = Re*nu/c     # Freestream velocity
  magVinf = U                     # (m/s) freestream velocity


  # Airfoil geometry
  if verbose; println("Defining geometry..."); end;
  A, B, CC = 2, 4, 12
  x, y = p2d.generate_naca4(A, B, CC, npoints; sharp=true)

  # Creates the paneled body
  points = hcat(x, y)
  types = ["source", "vortex"]
  body = p2d.Body(points, types)

  if save_path!=nothing
    p2d.vtk.create_path(save_path, prompt)
  end

  if verbose; println("Generating curves..."); end;
  alphas = []
  cls = []
  cds = []
  cms = []
  stall = zeros(Int64, 2)
  fig_num = 0
  for (j,ite_alphas) in enumerate([-dalpha:-dalpha:-20, 0:dalpha:20])
    for (i,alpha) in enumerate(ite_alphas)
      fig_num += 1

      if verbose
        println("\n\t*******************************************")
        println("\talpha=$alpha")
        println("\t*******************************************")
      end
      push!(alphas, alpha)

      # Freestream
      Vinf(X, t) = magVinf*[cos(alpha*pi/180), sin(alpha*pi/180)]

      # Solves the paneled body
      if verbose; println("\t\tSolving paneled body..."); end;
      p2d.solve(body, Vinf)

      if verbose; println("\t\tCalculating lift..."); end;
      push!(cls, p2d.calc_liftKJ(body, magVinf; c=c))

      if verbose; println("\t\tCalculating drag..."); end;
      if show_ite_plot
        plot_title = "Boundary layer on NACA2412, AOA=$alpha"
        fig = figure("blayer")
        clf()
      end
      more_outputs = []
      _, _, _, CDs = p2d.calc_blayer(body, mu, nu, Vinf(0,0); debug=debug,
                                          show_plot=show_ite_plot, verbose=verbose,
                                          more_outputs=more_outputs,
                                          plot_title=plot_title)
      push!(cds, CDs)

      if verbose; println("\t\tCalculating moment..."); end;
      Cm = p2d.calc_moment(body, magVinf)
      push!(cms, Cm)

      if save_path!=nothing
        num_str = (fig_num<10 ? "000" : fig_num<100 ? "00" : fig_num<1000 ? "0" : "")
        num_str *= "$fig_num"
        savefig(joinpath(save_path, "test_cl_naca2412."*num_str*".png"))
      end


      if verbose; print("\t\tChecking stall: "); end;
      low_points, up_points = more_outputs[2:3]
      low_sep_i = more_outputs[4]
      up_sep_i = more_outputs[8]

      # Stall condition: Separation before 80% of chord
      if up_sep_i!=nothing && up_points[up_sep_i][1]<0.77 && check_stall
        if verbose
          println("Upper surface stall at X=$(round.(up_points[up_sep_i],1))")
        end
        stall[j] = i
        break
      elseif low_sep_i!=nothing && low_points[low_sep_i][1]<0.77 && check_stall
        if verbose
          println("Lower surface stall at X=$(round.(low_points[low_sep_i],1))")
        end
        stall[j] = i
        break
      else
        if verbose; println("nothing"); end;
      end
    end
  end

  # println("alphas=$alphas")
  # println("cls=$cls")
  # println("cds=$cds")
  # println("stall=$stall")

  # XFOIL data
  data = CSV.read(joinpath(data_path, "xf-naca2412-il-1000000.csv"))

  # Plots
  alphas = vcat(reverse(alphas[1:stall[1]]), alphas[stall[1]+1:end])
  cls = vcat(reverse(cls[1:stall[1]]), cls[stall[1]+1:end])
  cds = vcat(reverse(cds[1:stall[1]]), cds[stall[1]+1:end])
  # cms = vcat(reverse(cms[1:stall[1]]), cms[stall[1]+1:end])
  if verbose; println("\t\tPlotting curves..."); end;
  fig = figure("curves", figsize=(7*3, 5*1))

  subplot(131)
  title("Lift Curve")
  plot(data[1], data[2], "-k", label="XFOIL")
  plot(alphas, cls, "--or", label="MyPanel2D")
  xlabel(L"Angle of attack $\alpha$ ($^\circ$)")
  ylabel(L"C_l")
  legend(loc="best")
  grid(true, color="0.8", linestyle="--")

  subplot(132)
  if polar
    title("Drag Polar")
    plot(data[3], data[2], "-k", label="XFOIL")
    plot([norm(cd[1]) for cd in cds], cls, "--or", label="MyPanel2D")
    ylabel(L"C_l")
    xlabel(L"C_d")
  else
    title("Drag Curve")
    plot(data[1], data[3], "-k", label="XFOIL")
    plot(alphas, [norm(cd[1]) for cd in cds], "--or", label="MyPanel2D")
    xlabel(L"Angle of attack $\alpha$ ($^\circ$)")
    ylabel(L"C_d")
  end
  legend(loc="best")
  grid(true, color="0.8", linestyle="--")

  subplot(133)
  title("Moment Curve")
  plot(data[1], data[5], "-k", label="XFOIL")
  plot(alphas, cms, "--or", label="MyPanel2D")
  xlabel(L"Angle of attack $\alpha$ ($^\circ$)")
  ylabel(L"C_m")
  ylim([-1, 1])
  legend(loc="best")
  grid(true, color="0.8", linestyle="--")



  if verbose; println("Done."); end;
end
