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
    a = (x1-x2) == 0.0 ? (0.0) : (y1 - y2) / (x1 - x2)
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
    if (a == 0) 
        @warn "a = ($a); not a quadratic formula!"
        return 0.0
    end

    d  = sqrt(b^2 - 4*a*c)
    return (-b - d) / (2*a), (-b + d) / (2*a)
end



function plot_lines(step::Real, df::DataFrame, c::ColorPalette)
    p::Plots.Plot = plot()
    cp::Float64 = 1.0
    for e1 in 1:step:length(unique(df.E1))
        for row in eachrow(df[df.E1 .== unique(df.E1)[e1],:])
            xs = row.minE : 1e-3 : row.maxE
            line2(x) = get_line_point(x, row.a, row.b)
            plot!(xs, line2.(xs), lw = 2, alpha = 0.3, c = c[ceil(Int,cp)], legend = :false,
                    xlabel = "E2 [MeV]", ylabel ="dGdE", ylims = (minimum(df.minG), 1.1*maximum(df.maxG)),
                    title = "projection of the linear approximations, \n every $step")
        end
        if cp < 250.0
            cp += step/10.0
        end

    end
    return p
end