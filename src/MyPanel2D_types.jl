#=##############################################################################
# DESCRIPTION
    Calculation of induced velocity of a variety of panel types.
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Feb 2018
  * License   : MIT License
=###############################################################################

"""
Velocity induced on point `C` by a dual source/vortex panel of ends `A`,`B`,
constant source strength density `q`, and uniform vortex stength density
`gamma`. The panel panel is expected to point from A to B such as to define the
tangent vector of that panel in that direction and the normal as the cross product
with an implicit z-axis coming out of the page (i.e., if the tanget points from
left to right, the normal will point upwards)
<!-- NOTE: Defining  ẑ*  as coming out of the page implies that a clockwise
circulation will have to be negative -->
"""
function sourcevortex2D_V(A::Array{Float64, 1}, B::Array{Float64, 1},
                              q::Float64, gamma::Float64, C::Array{Float64, 1};
                              more_outputs=nothing, reg_radius=1e-13, smth=1e-6)
  t = (B-A); t = t/norm(t);   # Tangent vector
  n = [-t[2], t[1]]           # Normal vector

  # Distances
  r_j = (C-A)
  magr_j = norm(r_j)
  r_jp1 = (C-B)
  magr_jp1 = norm(r_jp1)

  # Angles
  a_j = atan2( cross(vcat(t,0), vcat(r_j,0))[3] ,  dot(t, r_j) )
  a_jp1 = atan2( cross(vcat(t,0), vcat(r_jp1,0))[3] ,  dot(t, r_jp1) )

  # NOTE: I was having issues with numerically detecting collinearity, so
  #         here I am forcing colinearity (and its sign) when C is close enough
  #         to be colinear with AB, alomg with forcing aux2=1/2 on self induced
  a_j = (abs(a_j)<smth ? 0.0 : (abs(abs(a_j)-pi)<smth ? pi : a_j))
  a_jp1 = (abs(a_jp1)<smth ? 0.0 : (abs(abs(a_jp1)-pi)<smth ? pi : a_jp1))

  beta = a_jp1-a_j

  # # This line is necessary to force aux2 = 1/2 numerically on self induced velocity
  # if beta<0 && abs(beta+pi)<=1e-2
  #   beta = pi # If beta=-pi, it is forced to be +pi
  # end

  # # NOTE: This formulation puts singularities at both A and B
  # aux1 = -1/(2*pi) * log(magr_jp1 / magr_j)
  # aux2 = beta/(2*pi)

  # NOTE: This formulation avoid the singularities at A and B
  if magr_jp1<reg_radius
    magr_jp1 = reg_radius
  end
  if magr_j<reg_radius
    magr_j = reg_radius
  end

  aux1 = -1/(2*pi) * log(magr_jp1 / magr_j)
  aux2 = beta/(2*pi)

  # Source contribution
  sourceu = aux1
  sourcev = aux2
  sourceV = q*(sourceu*t + sourcev*n)

  # Vortex contribution
  vortexu = aux2
  vortexv = -aux1
  vortexV = gamma*(vortexu*t + vortexv*n)

  if typeof(more_outputs)==Array{Any, 1}
    push!(more_outputs, t)
    push!(more_outputs, n)
    push!(more_outputs, sourceV )
    push!(more_outputs, vortexV )
  end

  return sourceV+vortexV
end
