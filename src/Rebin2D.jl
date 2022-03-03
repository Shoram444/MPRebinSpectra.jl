function rebin2D(df::DataFrame, _prec::Float64 = 0.001)
    df_rebinned = DataFrame(E1 = Float64[],   
                            minE = Float64[], maxE = Float64[],            
                            minG = Float64[], maxG = Float64[],
                            a = Float64[], b = Float64[])

    for e in unique(df.E1)
        df_rebinned = vcat(df_rebinned, rebin1D(df[df.E1 .== e, 2], df[df.E1 .== e, 3], e, _prec))
    end
    
    return df_rebinned
end

function rebin1D(_x::Vector{Float64}, _y::Vector{Float64}, _E1::Float64, _prec::Float64 = 0.001)
    e0=1 #initialize index of the first point in fit 
    ef=2 #initialize index of the second point in fit

    a,b = get_line_params(0.0, 0.0, _x[e0], _y[e0])                 # first approximation is made from origin to first point in spectrum
    
    df_rebinned = DataFrame(E1 = _E1,                               # initialized DataFrame with first approximation
                            minE = 0.0, maxE = _x[e0],            
                            minG = 0.0, maxG = _y[e0],
                            a = a, b = b )
    
    if length(_x) == 1                                              # if there is only one point in the vector, first point is approximated from (0,0) to (x1,y1)
        a,b = get_line_params(_x[e0], _y[e0], 2*_x[e0], 0.0)        # second line is approximated from (x1,y1) to (x1+dx, 0) creating a triangle with area equal to 
                                                                    # rectangle with sides dx and y1 (dx is the step in x-direction...we assume x-array is evenly spaced)
        push!(df_rebinned, [_E1,
                            _x[e0], 2*_x[e0],    
                            _y[e0], 0.0     ,     
                            a     , b       ] )
        return df_rebinned
    end
    
    a,b = get_line_params(_x[e0], _y[e0], _x[ef], _y[ef])           # first fit points 1,2
    Δy  = abs( (get_line_point(_x[ef], a, b) - _y[ef]) / _y[ef] )   # Δy of points 1,2
    
    for _ in 1:length(_x)

        while Δy <= _prec && ef < length(_x)
            Δy = abs((get_line_point(_x[ef], a, b) - _y[ef])/_y[ef])
            ef += 1
        end
        
        if ef != length(_x)
            a, b = get_line_params(_x[e0], _y[e0], _x[ef-1], _y[ef-1]) 
                                                                    
            push!(df_rebinned, [_E1, 
                                _x[e0], _x[ef-1],                   # refit line with previous point 
                                _y[e0], _y[ef-1],                   # unless it's last point of the series
                                a     , b       ] )
        else
            a, b = get_line_params(_x[e0], _y[e0], _x[ef], _y[ef]) 
            push!(df_rebinned, [_E1, 
                                _x[e0], _x[ef], 
                                _y[e0], _y[ef],                     # for last point line stops at the endpoint
                                a     , b       ] )
        end

 

        if ef >= length(_x)                                         # if the next point is more than the size of data, check ends
            break
        end

        e0  = ef -1                                                 # last point from fit becomes new point in next fit
        a,b = get_line_params(_x[e0], _y[e0], _x[ef], _y[ef])       # find new line from the next set of points
        Δy  = abs((get_line_point(_x[ef], a, b) - _y[ef])/_y[ef])   # find new Δy

    end

    a, b = get_line_params(_x[ef], _y[ef], _x[ef]+_x[1], 0.0)       # finally to add very last approximation by going from (xf, yf) to (xf+dx, 0.0)
    push!(df_rebinned, [_E1, 
                        _x[ef], _x[ef]+_x[1], 
                        _y[ef], 0.0         ,                       # for last point line stops at the endpoint
                        a     , b           ] )       
    
    return df_rebinned
end


function normalize2D!(df_rebinned::DataFrame, thickness = 0.001)
    volume_rebinned = get_total_volume(df_rebinned, thickness)
    
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


"""
#### function ```get_cdf(df::DataFrame, thickness::Real = 0.001)```
<br>

    Description of ```get_cdf```
    ------------------------------
Calculate CDF for pdf = a + bx. CDF is cummulative density function defined as integral from -inf to x of pdf dx. 

"""
function get_cdf(df::DataFrame, thickness::Real = 0.001)
    n = nrow(df)
    cdf = Vector{Float64}(undef, n)
    cdf[1] = get_integral_linear(df[1,2], df[1,3], df[1,6], df[1,7])*thickness
    
    for i in 2:n
        cdf[i] = cdf[i-1] + get_integral_linear(df[i,2], df[i,3], df[i,6], df[i,7])*thickness
    end
    return cdf
end

function get_cdf!(df::DataFrame, thickness::Real = 0.001)
    n = nrow(df)
    cdf = Vector{Float64}(undef, n)
    cdf[1] = get_integral_linear(df[1,2], df[1,3], df[1,6], df[1,7])*thickness
    
    for i in 2:n
        cdf[i] = cdf[i-1] + get_integral_linear(df[i,2], df[i,3], df[i,6], df[i,7])*thickness
    end
    
    @transform! df :cdf = cdf
end
