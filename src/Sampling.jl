function sample_rejection(df::DataFrame; 
    minX = minimum(df.minE), maxX = maximum(df.maxE), 
    minY = minimum(df.minG), maxY = maximum(df.maxG))

    ξ   = rand(Uniform(minX,maxX),1)[1]
    η   = rand(Uniform(minY,maxY),1)[1]
    a,b = get_line_params(ξ, df)

    while η > get_line_point(ξ, a, b)
    ξ   = rand(Uniform(minX,maxX),1)[1]
    η   = rand(Uniform(minY,maxY),1)[1]
    a,b = get_line_params(ξ, df)
    end

    return ξ
end

function sample_discrete_CDF(gamma::Float64, v::Vector{Float64}, cdf::Vector{Float64})
    
    gamma < cdf[1] && return 1 # if gamma falls in the first bin move on right away
    
    for i in 1:length(v)-1
        if( cdf[i] = gamma <= cdf[i+1] )
            return i+1
            
        end
    end
end    

function sample_energies(df::DataFrame, volumesE1::Vector{Float64})
    cdf = cumsum(v)
    
    gamma = rand(Uniform())
    E1 = unique(df.E1)[sample_discrete_CDF(gamma, volumesE1, cdf)]
    E2 = sample_rejection(df[df.E1 .== E1, :])
    
    return E1, E2
end
