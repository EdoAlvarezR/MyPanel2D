include("../src/MyPanel2D.jl")
p2d = MyPanel2D

using PyPlot

for test_name in ["cube", "naca0012", "blayer_thwaites", "blayer_head",
                  "blayer_transition", "blayer_naca2412"]
  include("test_"*test_name*".jl")
end
