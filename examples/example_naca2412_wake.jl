
vpm_path = "/home/user/Dropbox/FLOWResearch/MyCodes/MyVPM"
include(joinpath(vpm_path, "src/MyVPM.jl"))
vpm = MyVPM

p2d_path = "/home/user/Dropbox/FLOWResearch/MyCodes/MyPanel2D"
include(joinpath(p2d_path, "src/MyPanel2D.jl"))
p2d = MyPanel2D

function example_VPMwake_naca2412(; AOA=2.0, verbose=true, show_plots=true,
                            npoints::Signed=75, debug=false,
                            save_path="temps/naca2412wake15", prompt=true)


  # ---------- SIMULATION PARAMETERS -------------------------------------------
  run_name = "naca2412wake"
  plot_title = "Boundary layer visualization on NACA2412, AOA=$AOA"

  # Physical parameters
  Re = 1e6                    # Reynolds number based on chord length
  c = 1                       # Chord length
  rho = 1.225                 # (kg/m^3) air density
  mu = 1.789e-5               # (kg/ms) dynamic viscosity
  nu = mu/rho                 # Kinematic viscosity
  U = Re*nu/c                 # Freestream velocity
  K = 0.6                     # Circulation reduction factor
  # AOA = 5.0                 # (deg) angle of attack
  magVinf = U                 # (m/s) freestream velocity
  Vinf(X, t) = magVinf*[cos(AOA*pi/180), sin(AOA*pi/180)]

  # VPM parameters
  solver = "ExaFMM"           # VPM solver
  nsteps = 250               # Number of time steps
  dt = 0.025*c/magVinf        # Time step
  n_steps_shedding = 1        # Time steps in between sheddings
  sigmap = 3*magVinf*dt       # Particle smoothing radius
  nsteps_save = 1             # Time steps in between VPM vtk outputs
  max_particles = 2*nsteps    # Max particles
  wake_max = 5.0              # Radii from origin to keep wake



  # ---------- PANEL SETUP -----------------------------------------------------
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
                        show_plot=show_plots, plot_title=plot_title,
                        verbose=verbose, more_outputs=more_outputs)

  # Identify separation point on lower and upper surface
  stg_i, low_points, up_points = more_outputs[1:3]
  low_sep_i = more_outputs[4]
  up_sep_i = more_outputs[8]

  # # Transforms the indices to body indices
  # blow_sep_i = low_sep_i==nothing? 1 : stg_i-((low_sep_i-1)-1)
  # bup_sep_i = up_sep_i==nothing? body.n : stg_i+(up_sep_i-1)-1 -5

  new_low_sep_i = low_sep_i==nothing ? size(low_points,1) : low_sep_i - 1
  new_up_sep_i = up_sep_i==nothing ? size(up_points,1) : up_sep_i - 1 -(up_sep_i>15? 15 : 0)

  # Transforms the indices to body indices
  blow_sep_i = stg_i-(new_low_sep_i-1)
  bup_sep_i = stg_i+new_up_sep_i-1


  # Shedding positions
  CP_up = p2d.get_panel(body, bup_sep_i)[4]
  X_up = CP_up + 0.5*p2d.get_tn(body, bup_sep_i)[2]*theta[bup_sep_i]*H[bup_sep_i]
  CP_low = p2d.get_panel(body, blow_sep_i)[4]
  X_low = CP_low + 0.5*p2d.get_tn(body, blow_sep_i)[2]*theta[blow_sep_i]*H[blow_sep_i]
  X_up3D = vcat(X_up, 0)
  X_low3D = vcat(X_low, 0)

  # Velocity at separation points
  # V_up = norm(Vinf(X_up, 0.0) + p2d.Vind(body, X_up))
  # V_low = norm(Vinf(X_low, 0.0) + p2d.Vind(body, X_low))
  V_up = norm(Vinf(X_up, 0.0) + p2d.Vind(body, CP_up))
  V_low = norm(Vinf(X_low, 0.0) + p2d.Vind(body, CP_low))

  println(V_up)
  println(V_low)



  # ---------- VPM SETUP -------------------------------------------------------
  if verbose; println("Setting up VPM..."); end;

  # vpm_Vinf(X, t) = vcat(Vinf(X[1:2], t) + p2d.Vind(body, X[1:2]), 0)

  # The panel method is way to slow, so let's define a radius where we will
  # actually calculate induced velocity
  function vpm_Vinf(X, t)
    this_Vinf = Vinf(X[1:2], t)
    this_Vind = norm(X)<1.1 && abs(X[2])<0.15 ? p2d.Vind(body, X[1:2]) : zeros(2)
    return vcat( this_Vinf+this_Vind, 0)
  end

  pfield = vpm.ParticleField(max_particles, vpm_Vinf, nothing, solver)




  # --------------- SETS UP RUNTIME ROUTINE ------------------------------------
  # Delete function
  function delete_particles(PFIELD, del_r)
    for i in 1:PFIELD.np
      if i>PFIELD.np
        break
      end
      this_X = vpm.get_x(PFIELD, i)
      if norm(this_X)>=del_r || abs(this_X[3])>0.05
        vpm.delparticle(PFIELD, i)
      end
    end
  end
  # Runtime function
  function runtime_function(PFIELD, T, DT)
    # Adds particles
    if PFIELD.nt%n_steps_shedding==0
      # Vectorial circulations
      Gamma_up = -[0, 0, K/2*V_up^2*DT]
      Gamma_low = [0, 0, K/2*V_low^2*DT]
      # Adds
      vpm.addparticle(pfield, vcat(X_up3D, Gamma_up, sigmap))
      vpm.addparticle(pfield, vcat(X_low3D, Gamma_low, sigmap))
    end
    # Deletes particles
    if PFIELD.nt%n_steps_shedding==0
      delete_particles(PFIELD, wake_max)
    end
    # Outputs panel
    if PFIELD.nt==0
      p2d.save(body, run_name; path=save_path)
    end

    return false
  end


  # --------------- RUN THE VPM ------------------------------------------------
  vpm.run_vpm!(pfield, dt, nsteps; save_path=save_path, run_name=run_name,
                    runtime_function=runtime_function, solver_method=solver,
                    nsteps_save=nsteps_save,
                    save_code="examples/example_naca2412_wake.jl",
                    prompt=prompt)


  return pfield, body
end
