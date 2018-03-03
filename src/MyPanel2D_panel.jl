#=##############################################################################
# DESCRIPTION
    Definition of a paneled body and solver interface.
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Feb 2018
  * License   : MIT License
=###############################################################################

"""
  `body(points, types)`
Receives a collection of points defining a panelized surface and creates a 2D
body.

  **Arguments**
  * `points::Array{Float64, 2}`   : Nodes paneling the body.
  * `types::Array{String, 1}`     : Type of panels (see global variable `TYPES`)


NOTE: Following the convention of an airfoil, the points must be given in
clockwise order (starting at TE, going from bottom to top surface and back to
TE). First and last may or may not need to be coincident (depending on the panel
types), but either way will apply the Kutta condition between first and last
panels.
"""
type Body
  # USER INPUTS
  nodes::Array{T, 2} where {T<:Real}          # Nodes on the body
  types::Array{String, 1}                     # Type of panels

  # Properties
  n::Signed                                   # Number of panels
  strengths::Array{T, 2} where {T<:Real}      # Stengths of each panel
  Vinf::Any                                   # Freestream used for solution
  solved::Bool                                # Flag whether solution exists

  Body(nodes, types=_check_types(types),
        n=size(nodes,1)-1,
        strengths=zeros(size(nodes,1)-1, length(types)),
        Vinf=nothing, solved=false
      ) = new(nodes, types,
        n,
        strengths,
        Vinf, solved)
end

"Sets a freestream velocity field `Vinf(X, t)`"
function setVinf(self::Body, Vinf)
  _reset(self)
  self.Vinf = Vinf
end

"Solves for the strengths of a paneled body with a given freestream Vinf(X,t)"
function solve(self::Body, Vinf; t::Float64=0.0)
  setVinf(self, Vinf)

  if self.types==["source", "vortex"]
    Q = _solve_sourcevortex(self, self.Vinf; time=t)
    self.strengths[:, 1] = Q[1:end-1]
    self.strengths[:, 2] = Q[end]

  elseif self.types==["source", "vortex2"]
    Q = _solve_sourcevortex2(self, self.Vinf; time=t, kt_i=Int64(ceil(self.n/3)))
    self.strengths[:, 1] = Q[1:end-1]
    self.strengths[:, 2] = Q[end]

  elseif self.types==["source"]
    Q = _solve_source(self, self.Vinf; time=t)
    self.strengths[:, 1] = Q[1:end]

  elseif !(self.types in keys(TYPES))
    error("Solution for panel type $(self.types) not implemented yet."*
            " Types available: $(keys(TYPES))")

  else
    error("Logic error!")
  end

  self.solved = true
end

"Returns nodes A and B, strength densities, and control point of the i-th panel"
function get_panel(self::Body, i::Int64)
  if i>self.n || i<=0
    error("Invalid panel $i ($i>$(self.n) || $i<=0).")
  end

  A = self.nodes[i, :]
  B = self.nodes[i+1, :]
  strengths = self.strengths[i,:]
  CP = (A+B)/2

  return A, B, strengths, CP
end

"Returns tangent and normal vectors of the i-th panel"
function get_tn(self::Body, i::Int64)
  A, B, _, _ = get_panel(self, i)
  return get_tn(A, B)
end

"Returns point of trailing edge"
function get_TE(self::Body)
  return self.nodes[1,:]
end

"Returns point of leading edge"
function get_LE(self::Body)
  TE = get_TE(self)   # Trailing edge
  cmax = -Inf         # Maximum chord
  cmax_i = nothing    # Point index of maximum chord

  # Iterates over all points searching the point of maximum chord
  for i in 1:self.n
    this_c = norm(self.nodes[i, :] - TE)
    if this_c >= cmax
      cmax = this_c
      cmax_i = i
    end
  end

  if cmax_i==nothing
    error("Logic error: No chord found!")
  end

  return self.nodes[cmax_i,:]
end

"Returns chord length"
function get_c(self::Body)
  xmax = maximum(self.nodes[:,1])
  ymin = maximum(self.nodes[:,1])
end

"Returns the velocity induced at `X` by this paneled body.
NOTE: This velocity doesn't include the freestream, only induced"
function Vind(self::Body, X; debug=false)
  if !self.solved
    error("Induced velocity requested on an unsolved body."*
            " Call `solve()` before calling this function.")
  end

  V = zeros(2)

  if self.types in [["source", "vortex"], ["source", "vortex2"], ["source"]]

    for i in 1:self.n
      A, B, strengths, _ = get_panel(self, i)
      if self.types==["source"]
        q, gamma = vcat(strengths, 0)
      else
        q, gamma = strengths
      end

      V += sourcevortex2D_V(A, B, q, gamma, X)
    end

  else
    error("Velocity calculation for panel type $(self.types) not"*
                                                            "implemented yet!")
  end

  return V
end

"Returns the pressure coefficient at point X"
function p_coeff(self::Body, X, magVinf::Real; t::Real=0)
  this_Vinf = self.Vinf(X, t)
  return 1 - ( norm( this_Vinf+Vind(self, X) ) / norm(magVinf) )^2
end

function save(self::Body, filename::String;
                              num=nothing, time=nothing, path="", comments="",
                              t::Real=0, magVinf="automatic")

  # Formats the body as a VTK
  points = [ vcat(self.nodes[i,:], 0) for i in 1:self.n+1 ]
  lines = [[i-1 for i in 1:self.n+1]]

  if self.solved

    # Formats the strength densities as a data field
    data = [Dict(
                  "field_name" => self.types[i],
                  "field_type" => "scalar",
                  "field_data" => vcat(self.strengths[:,i],
                                                    self.strengths[end,i])
                  )
              for i in 1:length(self.types)
            ]

    # Creates a pressure distribution data field
    _magVinf = magVinf=="automatic" ? norm(self.Vinf(zeros(2), 0)) : magVinf
    arr_p_field = [ p_coeff(self, CP, _magVinf; t=t) for CP in
                                [get_panel(self, i)[4] for i in 1:self.n]]
    arr_p_field = vcat(arr_p_field, p_coeff(self, get_panel(self, self.n)[2],
                                                                _magVinf; t=t))
    p_field = Dict(
        "field_name" => "Cp",
        "field_type" => "scalar",
        "field_data" =>  arr_p_field
        )
    push!(data, p_field)

    # Creates velocity distribution data fields
    arr_Vind_field = Array{Real, 1}[]
    arr_Vinf_field = Array{Real, 1}[]
    arr_Vtot_field = Array{Real, 1}[]
    for i in 1:self.n+1
      if i==(self.n+1)
        CP = get_panel(self, i-1)[2]
      else
        CP = get_panel(self, i)[4]
      end
      this_Vind = vcat(Vind(self, CP), 0)
      this_Vinf = vcat(self.Vinf(CP, t), 0)
      push!(arr_Vind_field, this_Vind)
      push!(arr_Vinf_field, this_Vinf)
      push!(arr_Vtot_field, this_Vind+this_Vinf)
    end
    for (field_name, field_data) in [("Vind", arr_Vind_field),
                                      ("Vinf", arr_Vinf_field),
                                      ("Vtot", arr_Vtot_field)]
      push!(data, Dict(
          "field_name" => field_name,
          "field_type" => "vector",
          "field_data" => field_data))
    end

  else
    data = []

  end


  # Creates vtk
  vtk.generateVTK(filename*"_pnl", points; lines=lines, point_data=data,
                          num=num, time=time, path=path, comments=comments)
end

"Plots the pressure distribution along the body"
function plot_P(self::Body; magVinf="automatic", t::Real=0)
  xs = []
  ys = []
  _magVinf = magVinf=="automatic" ? norm(self.Vinf(zeros(0), 0)) : magVinf
  for i in 1:self.n
    _, _, _, CP = get_panel(self, i)
    push!(xs, CP[1])
    push!(ys, p_coeff(self, CP, _magVinf; t=t))
  end

  fig = plt.figure("Cp")
  this_plot = plt.plot(xs, ys, "--.k")
  plt.xlabel("x-position")
  plt.ylabel(L"C_p")
  plt.title("Pressure distribution")
  plt.grid(true, color="0.8", linestyle="--")
  this_plot[1]["axes"]["invert_yaxis"]()
end


"Plots the tangential velocity along the body"
function plot_Vtan(self::Body; t::Real=0)
  xs = []
  ys = []
  for i in 1:self.n
    _, _, _, CP = get_panel(self, i)
    push!(xs, CP[1])
    push!(ys, norm(self.Vinf(CP,t) + Vind(self, CP)))
  end

  fig = plt.figure("Vtan")
  plt.plot(xs, ys, "--.k")
  plt.xlabel("x-position")
  plt.ylabel(L"\hat{t}\cdot \vec{V}")
  plt.title("Tangential velocity")
  plt.grid(true, color="0.8", linestyle="--")
end

# ------- INTERNAL FUNCTION ----------------------------------------------------
function _check_types(types::Array{String, 1})
  if !(types in keys(TYPES))
    error("Invalid type $types. Please use one of the following: $(keys(TYPES))")
  end
  return types
end

"Removes previous solutions"
function _reset(self::Body)
  self.strengths = zeros(size(self.nodes,1)-1, length(self.types))
  self.Vinf = nothing
  self.solved = false
end
#
