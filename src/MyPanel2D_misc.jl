
"Generates a NACA 4-series airfoil of the form ABCC where A is the max camber
as percetange of chord, B is the location of max camber in tenths of chord,
and CC is the thickness of the airofil as percentage of chord."
function generate_naca4(A::Int64, B::Int64, CC::Int64, n::Int64;
                              sharp::Bool=false, r::Float64=100.0)
  eps = A/100
  p = B/10
  tau = CC/100

  function thickness(x)
    return 10*tau*( 0.2969*sqrt(x) - 0.1260*x -
                        (sharp ? 0.3537:0.3516)*x^2 + 0.2843*x^3 - 0.1015*x^4)
  end
  function ybar(x)
    if 0<=x<=p
      return eps/p^2 * (2*p*x-x^2)
    else
      return eps/(1-p)^2 * (1-2*p+2*p*x-x^2)
    end
  end
  yu(x) = ybar(x) + thickness(x)/2
  yl(x) = ybar(x) - thickness(x)/2

  f(x) = x
  xu = vtk.discretize(f, 0.0, 1.0, n, r; central=true)
  xu = Float64[xi for xi in xu]
  xl = reverse(xu[2:end])
  yl = [yl(x) for x in xl]
  yu = [yu(x) for x in xu]

  x = vcat(xl, xu)
  y = vcat(yl, yu)

  return x, y
end

"Calculates the total circulation of the body by integrating the vortex strength density
along the surface assuming it to be constant along each panel."
function calc_Gamma(body::Body)

  circ = 0

  if body.types in [["source", "vortex"], ["source", "vortex2"]]
    for i in 1:body.n
      A, B, strengths, CP = get_panel(body, i)
      circ += body.strengths[i, 2]*norm(B-A)
    end

  else
    error("Circulation can't be calculated for panel type $(body.types)")
  end

  return circ
end

"Calculates lift coefficient using Kutta-Joukowski theorem"
function calc_liftKJ(body::Body, magVinf::Real; c="automatic")
  _c = c=="automatic" ? get_c(body) : c
  Gamma = calc_Gamma(body)
  cl = 2*Gamma/(magVinf*_c)
  return cl
end

"Integrates the pressure along the surface to obtain a force"
function calc_pressureforce(body::Body, magVinf::Real; c="automatic", t::Real=0)
  # Integrates the pressure
  Cf = zeros(2)
  for i in body.n
    A, B, strengths, CP = get_panel(body, i)
    _, n = get_tn(body, i)
    p = p_coeff(body, CP, magVinf; t=t)
    Cf += -p*norm(B-A)*n
  end

  # Normalizes it by the chord length
  Cf = Cf / (c=="automatic" ? get_c(body) : c)

  return Cf
end

"Calculates the moment coefficient about the quarter chord by integrating the
pressure distribution along the surface of the body"
function calc_moment(body::Body, magVinf::Real; c="automatic", pos::Real=1/4,
                                                                    t::Real=0)
  # Defines the moment center of reference
  LE = get_LE(body)
  X0 = LE + pos*(get_TE(body)-LE)

  # Integrates the pressure
  Cm = 0
  for i in 1:body.n
    A, B, strengths, CP = get_panel(body, i)
    _, n = get_tn(body, i)
    p = p_coeff(body, CP, magVinf; t=t)   # Pressure
    Cf = p*norm(B-A)*n                    # Force
    l = CP-X0                             # Moment arm
    Cm += cross(vcat(l, 0), vcat(Cf, 0))[3]     # Moment
  end

  Cm = Cm / (c=="automatic" ? get_c(body) : c)^2    # Normalizes by chord square

  return Cm
end
