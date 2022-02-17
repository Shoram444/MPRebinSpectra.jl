module MPRebinSpectra

using DataFrames, DataFramesMeta, Random, Distributions
# Write your package code here.

include("Line.jl")
export  get_line_params, 
        get_line_point

include("Volume.jl")
export  get_area,
        get_row_volume,
        get_total_volume

include("Rebin2D.jl")
export  rebin,
        rebin2D,
        normalize2D!

end
