module MPRebinSpectra

using DataFrames, DataFramesMeta, Random, Distributions, StatsPlots, ColorSchemes
# Write your package code here.

include("Line.jl")
export  get_line_params, 
        get_line_point,
        get_integral_linear,
        #solvequadratic,
        plot_lines

include("Volume.jl")
export  get_area,
        get_row_volume,
        get_total_volume,
        get_segment_volume,
        get_volume_matrix


include("Rebin2D.jl")
export  rebin,
        rebin2D,
        normalize2D!,
        get_cdf,
        get_cdf!


end
