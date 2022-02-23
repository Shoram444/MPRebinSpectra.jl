module MPRebinSpectra

using DataFrames, DataFramesMeta, Random, Distributions
# Write your package code here.

include("Line.jl")
export  get_line_params, 
        get_line_point,
        get_integral_linear,
        solvequadratic

include("Volume.jl")
export  get_area,
        get_row_volume,
        get_total_volume,
        get_segment_volume,
        get_volume_matrix


include("Rebin2D.jl")
export  rebin,
        rebin2D,
        normalize2D!

include("Sampling.jl")
export  sample_rejection,
        sample_discrete_CDF,
        sample_energies_old,
        get_cdf,
        get_cdf!,
        sample_energies

end
