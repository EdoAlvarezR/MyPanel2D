include("../src/MyPanel2D.jl")
p2d = MyPanel2D

using PyPlot
import CSV

module_path,_ = splitdir(@__FILE__);   # Path to this module
data_path = joinpath(module_path, "data")

for test_name in ["cube", "naca0012", "blayer_thwaites", "blayer_head",
                  "blayer_transition", "blayer_naca2412",
                  "cl_naca2412"]
  include("test_"*test_name*".jl")
end
