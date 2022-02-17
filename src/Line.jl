function get_line_params(xi::Float64, df_normed::DataFrame)
    row = df_normed[(df_normed.E0 .< xi) .& (df_normed.Ef .> xi),:]
    return (row.a_normed[1], row.b_normed[1])
end

function get_line_params(x1::Float64, y1::Float64, x2::Float64, y2::Float64) 
    a = (y1 - y2) / (x1 - x2)
    b = y1 - a*x1
    return a,b #returns the parameters of line given by y = b + ax
end

get_line_point(x::Float64,a::Float64,b::Float64) = b + a*x



