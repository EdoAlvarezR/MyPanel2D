#=##############################################################################
# DESCRIPTION
    Solvers for a variety of panel types.
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Feb 2018
  * License   : MIT License
=###############################################################################

# ------- INTERNAL FUNCTION ----------------------------------------------------
"Solves the body of source-vortex panels with uniform vortex strength density,
and imposing the Kutta condition of equal tangent velocity between first and
last panel (TE)"
function _solve_sourcevortex(body::Body, Vinf; time::Real=0)
  N = body.n                      # Number of panels
  G = zeros(Real, N+1, N+1)       # Geometric matrix
  Q = zeros(Real, N+1)            # Source strengths and vortex
  V = zeros(Real, N+1)            # Freestream normal

  for i in 1:N # Iterates over target panels

    A, B, _, CP = get_panel(body, i)  # Target panel

    # Gets tangent, normal, and self-induced velocities of this panel
    aux1 = []
    sourcevortex2D_V(A, B, 1.0, 1.0, CP; more_outputs=aux1)
    t, n, Gs, Gv = aux1

    # Adds self induced terms
    G[i, i] = dot(n, Gs)          # Source
    G[i, N+1] = dot(n, Gv)        # Vortex

    # Kutta condition on self induced terms
    if i==1
      G[N+1, i] = -dot(t, Gs)
      G[N+1, N+1] = -dot(t, Gv)
    elseif i==N
      G[N+1, i] -= dot(t, Gs)
      G[N+1, N+1] -= dot(t, Gv)
    end

    for j in 1:N  # Velocity induced by all other panels on target
      if i!=j
        A, B, _, _ = get_panel(body, j)
        aux1 = []
        sourcevortex2D_V(A, B, 1.0, 1.0, CP; more_outputs=aux1)
        _, _, Gs, Gv = aux1

        G[i, j] = dot(n, Gs)
        G[i, N+1] += dot(n, Gv)

        # Kutta condition
        if i==1
          G[N+1, j] -= dot(t, Gs)
          G[N+1, N+1] -= dot(t, Gv)
        elseif i==N
          G[N+1, j] -= dot(t, Gs)
          G[N+1, N+1] -= dot(t, Gv)
        end

      end

    end

    # Freestream normal component
    V[i] = -dot(n, Vinf(CP, time))
    if i==1
      V[N+1] = dot(t, Vinf(CP, time))
    elseif i==N
      V[N+1] += dot(t, Vinf(CP, time))
    end

  end

  # Calculates strengths
  Q = G\V

  return Q
end

"Solves the body of source-vortex panels with uniform vortex strength density,
and imposing the Kutta condition of stagnation at the first node (TE)"
function _solve_sourcevortex2(body::Body, Vinf;
                                  time::Real=0, kt_i::Int64=1)

  # Checks that the body has a closed trailing edge
  if get_panel(body, 1)[1]!=get_panel(body, body.n)[2]
    error("First and last nodes must be coincident!")
  end

  N = body.n                      # Number of panels
  G = zeros(Real, N+1, N+1)       # Geometric matrix
  Q = zeros(Real, N+1)            # Source strengths and vortex
  V = zeros(Real, N+1)            # Freestream normal

  for i in 1:N # Iterates over target panels

    A, B, _, CP = get_panel(body, i)  # Target panel

    # Gets tangent, normal, and self-induced velocities of this panel
    aux1 = []
    sourcevortex2D_V(A, B, 1.0, 1.0, CP; more_outputs=aux1)
    t, n, Gs, Gv = aux1

    # Adds self induced terms
    G[i, i] = dot(n, Gs)          # Source
    G[i, N+1] = dot(n, Gv)        # Vortex

    for j in 1:N  # Velocity induced by all other panels on target
      if i!=j
        A, B, _, _ = get_panel(body, j)
        aux1 = []
        sourcevortex2D_V(A, B, 1.0, 1.0, CP; more_outputs=aux1)
        _, _, Gs, Gv = aux1

        G[i, j] = dot(n, Gs)
        G[i, N+1] += dot(n, Gv)
      end
    end

    # Freestream normal component
    V[i] = -dot(n, Vinf(CP, time))

  end

  # Kutta condition (stagnation at trailing edge node)
  TE = get_panel(body, 1)[1]
  Vinf_TE = Vinf(TE, time)
  G[N+1, :] = 0
  G[kt_i, :] = 0
  V[N+1] = -Vinf_TE[1]
  V[kt_i] = -Vinf_TE[2]

  for i in 1:N # Iterates over inducing panels
      A, B, _, _ = get_panel(body, i)

      # Gets induced-velocity factors of this panel
      aux1 = []
      sourcevortex2D_V(A, B, 1.0, 1.0, TE; more_outputs=aux1)
      _, _, Gs, Gv = aux1

      G[N+1, i] = Gs[1]     # x-component of source
      G[N+1, N+1] += Gv[1]  # x-component of vortex

      # NOTE: These line overwrites one of the no-flow equations
      G[kt_i, i] = Gs[2]     # y-component of source
      G[kt_i, N+1] += Gv[2]  # y-component of vortex
  end

  # Calculates strengths
  Q = G\V

  return Q
end

"Solves the body of source panels (constant strengths)"
function _solve_source(body::Body, Vinf; time::Float64=0.0)
  N = body.n                      # Number of panels
  G = zeros(Real, N, N)           # Geometric matrix
  Q = zeros(Real, N)              # Source strengths
  V = zeros(Real, N)              # Freestream normal

  for i in 1:N # Iterates over target panels

    A, B, _, CP = get_panel(body, i)  # Target panel

    # Gets tangent, normal, and self-induced velocities of this panel
    aux1 = []
    sourcevortex2D_V(A, B, 1.0, 0.0, CP; more_outputs=aux1)
    t, n, Gs, _ = aux1

    # Adds self induced terms
    G[i, i] = dot(n, Gs)

    for j in 1:N  # Velocity induced by all other panels on target
      if i!=j
        A, B, _, _ = get_panel(body, j)
        aux1 = []
        sourcevortex2D_V(A, B, 1.0, 0.0, CP; more_outputs=aux1)
        _, _, Gs, _ = aux1

        G[i, j] = dot(n, Gs)
      end
    end

    # Freestream normal component
    V[i] = -dot(n, Vinf(CP, time))
  end

  # Calculates strengths
  Q = G\V

  return Q
end
