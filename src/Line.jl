"""
#### function ```get_line_params(xi::Float64, df::DataFrame)``` 
<br>

    Description of ```get_line_params```
    -----------------------------------------
Returns the parameters a,b of linear approximation y = ax +b. For the segment where the random number ξ falls into. 
``
"""
function get_line_params(xi::Float64, df::DataFrame)
    row = df[(df.minE .< xi) .& (df.maxE .> xi),:]
    return (row.a[1], row.b[1])
end


"""
#### function ```get_line_params(x1::Float64, y1::Float64, x2::Float64, y2::Float64)``` 
<br>

    Description of ```get_line_params```
    -----------------------------------------
Returns the parameters a,b of linear approximation y = ax +b. 
``
"""
function get_line_params(x1::Float64, y1::Float64, x2::Float64, y2::Float64) 
    a = (y1 - y2) / (x1 - x2)
    b = y1 - a*x1
    return a,b #returns the parameters of line given by y = b + ax
end

get_line_point(x::Float64,a::Float64,b::Float64) = b + a*x

"""
#### function ```get_integral_linear(minE::Real, maxE::Real, a::Real, b::Real)``` 
<br>

    Description of ```get_integral_linear```
    -----------------------------------------
Return the definitive integral of ``
    \\int_{minE}^{maxE} (a_iE + b_i) dE= \\frac{a_i}{2}(maxE^2 - minE^2) + b(maxE - minE)
``
"""
function get_integral_linear(minE::Real, maxE::Real, a::Real, b::Real)    
    return 0.5*a * (maxE^2 - minE^2) + b * ( maxE - minE )
end


function solvequadratic(a, b, c)
    d  = sqrt(b^2 - 4*a*c)
    return (-b - d) / (2*a), (-b + d) / (2*a)
end

