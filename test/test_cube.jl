
"Testing that calculation of inviscid flow (no dettachement) around a cube
predicts zero circulation"
function test_cube(; save_path=nothing, run_name="test_cube", delete=false,
                      show_plots=false, verbose=true)

  # Generates the nodes around each face of the cube
  lcube = sqrt(2)             # Face length
  Ocube = [0, 0]              # Cube centroid
  npanels = 10                # Panel per face
  x,y = zeros(4*npanels+1), zeros(4*npanels+1)
  i = 1
  for face_i in 1:4

    if face_i==1
      t = [-1, -1]
      A = lcube/sqrt(2)*[1,0]
    elseif face_i==2
      t = [-1, 1]
      A = lcube/sqrt(2)*[0,-1]
    elseif face_i==3
      t = [1, 1]
      A = lcube/sqrt(2)*[-1,0]
    else
      t = [1, -1]
      A = lcube/sqrt(2)*[0,1]
    end
    t = t/norm(t)

    for s in linspace(0, 1, npanels+1)[1:end-1]
      p = Ocube + A + (s*lcube)*t
      x[i] = p[1]
      y[i] = p[2]
      i += 1
    end
  end
  x[end] = x[1]
  y[end] = y[1]

  # Creates the paneled body
  points = hcat(x, y)
  types = ["source", "vortex"]
  body = p2d.Body(points, types)

  # Freestream
  AOA = 0.0                      # (deg) angle of attack
  magVinf = 10.0                 # (m/s) freestream velocity
  Vinf(X, t) = magVinf*[cos(AOA*pi/180), sin(AOA*pi/180)]

  # Solves the paneled body
  p2d.solve(body, Vinf)

  if show_plots
    p2d.plot_Vtan(body)
    p2d.plot_P(body)
  end

  # Run tests
  flag = true

  # Verification of flow tangency
  if verbose
    println("\n---------- FLOW TANGENCY TEST ------------------------")
    println("Panel\tt\thatV\tTangent?")
  end
  for p_i in 1:body.n
    A, B, _, CP = p2d.get_panel(body, p_i)
    this_Vinf = Vinf(CP, 0.0)
    this_Vind = p2d.Vind(body, CP)
    this_V = this_Vinf+this_Vind
    t = (B-A)/norm(B-A)
    hatV = this_V / norm(this_V)
    tang = norm(cross(vcat(t,0), vcat(hatV,0))) < 1e-4

    flag *= tang
    if verbose; println("$p_i\t$t\t$hatV\t$tang"); end;
  end


  # Verification of Kutta condition
  A, B, strengths, CP = p2d.get_panel(body, 1)
  aux1 = []
  p2d.sourcevortex2D_V(A, B, 1.0, 1.0, CP;  more_outputs=aux1)
  t, n, _, _ = aux1
  this_Vinf = dot(t, Vinf(CP, 0.0))
  this_Vind = dot(t, p2d.Vind(body, CP))
  magVt1 = this_Vinf + this_Vind
  t1 = t
  CP1 = CP
  A, B, strengths, CP = p2d.get_panel(body, body.n)
  aux1 = []
  p2d.sourcevortex2D_V(A, B, 1.0, 1.0, CP;  more_outputs=aux1)
  t, n, _, _ = aux1
  this_Vinf = dot(t, Vinf(CP, 0.0))
  this_Vind = dot(t, p2d.Vind(body, CP))
  magVtN = this_Vinf + this_Vind
  tN = t
  CPN = CP
  kutta_test = abs(magVtN+magVt1)/abs(magVtN)<1e-6
  if verbose
    println("\n---------- KUTTA CONDITION TEST ------------------------")
    println("\tCP1=$CP1\tCPN=$CPN")
    println("\tVt1: $magVt1\tt1=$t1")
    println("\tVtN: $magVtN\ttN=$tN")
    if kutta_test
      println("Passed!")
    else
      println("*****FAILED!****")
    end
  end
  flag *= kutta_test

  # Verification of Kutta condition on pressure
  _,_,_,CP1 = p2d.get_panel(body, 1)
  _,_,_,CPend = p2d.get_panel(body, body.n)
  p1 = p2d.p_coeff(body, CP1)
  pend = p2d.p_coeff(body, CPend)
  kutta_test2 = abs(p1-pend)/abs(p1)<1e-6
  if verbose
    println("\n---------- KUTTA CONDITION TEST2 ------------------------")
    println("\tCP1: $(p1)")
    println("\tCPN: $(pend)")
    if kutta_test2
      println("Passed!")
    else
      println("*****FAILED!****")
    end
  end
  flag *= kutta_test2

  # Verification of circulation
  gamma_exp = 0.0
  gamma_calc = p2d.calc_Gamma(body)
  gamma_norm = gamma_exp!=0 ? gamma_exp : 1
  gamma_test = abs((gamma_exp-gamma_calc)/gamma_norm)<1e-6
  if verbose
    println("\n---------- CIRCULATION TEST ------------------------")
    println("\tExpected Circulation: $gamma_exp")
    println("\tCalculated Circulation: $gamma_calc")
    if gamma_test
      println("Passed!")
    else
      println("*****FAILED!****")
    end
  end
  flag *= gamma_test

  # Visualization
  if save_path!=nothing
      p2d.vtk.create_path(save_path, false)
      p2d.save(body, run_name; path=save_path)

      strn = "$(run_name)_pnl.vtk;"
      run(`paraview --data="$save_path$strn"`)

      if delete
          run(`rm -rf save_path`)
      end
  end

  return flag
end
