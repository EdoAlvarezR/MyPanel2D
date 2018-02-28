"""
  Two dimensional panel methods for aerodynamics calculations.

  # AUTHORSHIP
    * Author    : Eduardo J Alvarez
    * Email     : Edo.AlvarezR@gmail.com
    * Created   : Feb 2018
    * License   : MIT License

  # DEVELOPMENT HISTORY
    * 20180224  : Define architecture and write source+vortex method.
"""
module MyPanel2D

# ------------ GENERIC MODULES -------------------------------------------------
using PyPlot
plt = PyPlot

import NaNMath

# ------------ FLOW MODULES ----------------------------------------------------
# VTK geometry processing tools https://github.com/byuflowlab/VTKtools.jl
vtk_path = "/home/user/Dropbox/FLOWResearch/FLOWCodes/VTKtools/"
include(joinpath(vtk_path, "src/VTKtools.jl"))
vtk = VTKtools

# ------------ HEADERS ---------------------------------------------------------
for header_name in ["types", "panel", "solver", "misc", "blayer"]
  include("MyPanel2D_"*header_name*".jl")
end


# ------------ GLOBAL VARIABLES ------------------------------------------------
global module_path; module_path,_ = splitdir(@__FILE__);   # Path to this module
global data_path = module_path*"/../data/"       # Path to data folder

# Types of panels currently implemented
global TYPES = Dict( # [type1, type2, ...] => []

  # Source panels of constant strengths
  ["source"]              => [],

  # Dual source+vortex panels of constant source strength and uniform vortex
  # strength with tangential velocity TE Kutta condition
  ["source", "vortex"]    => [],

  # Dual source+vortex panels of constant source strength and uniform vortex
  # strength with stagnation TE Kutta condition.
  ["source", "vortex2"]   => []
)


end # END OF MODULE
