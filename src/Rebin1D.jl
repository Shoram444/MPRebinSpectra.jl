"""
function rebin1D(_x::Vector{Float64}, _y::Vector{Float64}, _prec::Float64 = 0.001)

Description of ```rebin1D```
------------------------------
Reshapes the 1D spectrum based on linear approximations. 

Arguments in are:
    + _x::Vector{Float64} : vector of x-values
    + _y::Vector{Float64} : vector of y-values
    + _prec::Float64      : floating point precision (default is 1 promile, _prec= 0.001)
------------------------------

Returns a DataFrame which has the following columns:
    + firstFitPoint	: the first point of the fit within the given dataframe row
    + lastFitPoint	: the last point of the fit within the given dataframe row
    + E0	        : energy of the first point
    + Ef	        : energy of the last point
    + G0	        : y-value of the first point (used for finding the minimum of the "box" in rejection method sampling)
    + Gf	        : y-value of the last point (used for finding the maximum of the "box" in rejection method sampling)
    + a             : a parameter of the fit (y = b +ax)
    + b             : b parameter of the fit (y = b +ax)

"""
function rebin1D(_x::Vector{Float64}, _y::Vector{Float64}, _prec::Float64 = 0.001)
    e0=1 #initialize index of the first point in fit 
    ef=2 #initialize index of the second point in fit
    
    a,b = get_line_params(_x[e0], _y[e0], _x[ef], _y[ef])  #first fit points 1,2
    Δy  = abs( (get_line_point(_x[ef], a, b) - _y[ef]) / _y[ef] )  #Δy of points 1,2
    
    df_rebinned = DataFrame(firstFitPoint = e0,lastFitPoint=ef,  # initialize DataFrame to be returned             
                            E0 = _x[e0], Ef = _x[ef],            # first row is a dummy and will be deleted 
                            G0 = _y[e0], Gf = _y[ef],
                            a = a, b = b, Δy = Δy)
                    
    for _ in 1:length(_x)

        while Δy <= _prec && ef < length(_x)
            Δy = abs((get_line_point(_x[ef], a, b) - _y[ef])/_y[ef])
            ef += 1
        end
        
        if ef != length(_x)
            a, b = get_line_params(_x[e0], _y[e0], _x[ef-1], _y[ef-1]) 
                                                                    
            push!(df_rebinned, [e0    ,     ef-1, 
                                _x[e0], _x[ef-1],    #refit line with previous point 
                                _y[e0], _y[ef-1],    #unless it's last point of the series
                                a, b, Δy])           # add row into df 
        else
            a, b = get_line_params(_x[e0], _y[e0], _x[ef], _y[ef]) 
            push!(df_rebinned, [e0    ,     ef,
                                _x[e0], _x[ef], 
                                _y[e0], _y[ef],        #for last point line stops at the endpoint
                                a     ,      b, Δy])   # add row into df 
        end


        if ef >= length(_x) # if the next point is more than the size of data, check ends
            break
        end

        e0  = ef -1 # last point from fit becomes new point in next fit
        a,b = get_line_params(_x[e0], _y[e0], _x[ef], _y[ef])  #find new line from the next set of points
        Δy  = abs((get_line_point(_x[ef], a, b) - _y[ef])/_y[ef]) #find new Δy

    end
    
    delete!(df_rebinned, [1]) 
    return df_rebinned
end