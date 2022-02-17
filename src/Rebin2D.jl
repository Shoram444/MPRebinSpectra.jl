function rebin2D(df::DataFrame, _prec::Float64 = 0.001)
    df_rebinned = DataFrame(E1 = Float64[],
                            minE = Float64[], maxE = Float64[],            
                            minG = Float64[], maxG = Float64[],
                            a = Float64[], b = Float64[])

    for e in unique(df.E1)
        df_rebinned = vcat(df_rebinned, rebin(df[df.E1 .== e, 2], df[df.E1 .== e, 3], e, _prec))
    end
    
    return df_rebinned
end

function rebin(_x::Vector{Float64}, _y::Vector{Float64}, _E1::Float64, _prec::Float64 = 0.001)
    e0=1 #initialize index of the first point in fit 
    ef=2 #initialize index of the second point in fit
    
    df_rebinned = DataFrame(E1 = Float64[],
                            minE = Float64[], maxE = Float64[],            
                            minG = Float64[], maxG = Float64[],
                            a = Float64[], b = Float64[] )
    
    if length(_x) == 1
        push!(df_rebinned, [_E1,
                                _x[e0], _x[ef-1],    #refit line with previous point 
                                _y[e0], _y[ef-1],    #unless it's last point of the series
                                0     , _x[e0]  ] )
        return df_rebinned
    end
    
    a,b = get_line_params(_x[e0], _y[e0], _x[ef], _y[ef])  #first fit points 1,2
    Δy  = abs( (get_line_point(_x[ef], a, b) - _y[ef]) / _y[ef] )  #Δy of points 1,2
    
    
                    
    for _ in 1:length(_x)

        while Δy <= _prec && ef < length(_x)
            Δy = abs((get_line_point(_x[ef], a, b) - _y[ef])/_y[ef])
            ef += 1
        end
        
        if ef != length(_x)
            a, b = get_line_params(_x[e0], _y[e0], _x[ef-1], _y[ef-1]) 
                                                                    
            push!(df_rebinned, [_E1, 
                                _x[e0], _x[ef-1],    #refit line with previous point 
                                _y[e0], _y[ef-1],    #unless it's last point of the series
                                a     , b       ] )
        else
            a, b = get_line_params(_x[e0], _y[e0], _x[ef], _y[ef]) 
            push!(df_rebinned, [_E1, 
                                _x[e0], _x[ef], 
                                _y[e0], _y[ef],        #for last point line stops at the endpoint
                                a     , b        ] )
        end


        if ef >= length(_x) # if the next point is more than the size of data, check ends
            break
        end

        e0  = ef -1 # last point from fit becomes new point in next fit
        a,b = get_line_params(_x[e0], _y[e0], _x[ef], _y[ef])  #find new line from the next set of points
        Δy  = abs((get_line_point(_x[ef], a, b) - _y[ef])/_y[ef]) #find new Δy

    end
    
    return df_rebinned
end


function normalize2D!(df_rebinned::DataFrame; thickness = 0.0001)
    volume_rebinned = get_total_volume(df_rebinned; thickness)
    
    df_normed = select(df_rebinned,
                    :E1            => :E1,
                    :minE          => :minE,
                    :maxE          => :maxE,
                    :minG=> ByRow(g->g/volume_rebinned) => :minG,
                    :maxG=> ByRow(g->g/volume_rebinned) => :maxG,
                    :a   => ByRow(a->a/volume_rebinned) => :a,
                    :b   => ByRow(b->b/volume_rebinned) => :b,
                    )
    return df_normed
end