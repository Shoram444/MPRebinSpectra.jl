function get_area(x1::Float64, x2::Float64, a::Float64, b::Float64)
    # area under the trapezoid is A = Δx * ⟨y⟩
    # where Δx is the width of the base of trapezoid
    # and ⟨y⟩ is mean height
    
    Δx = abs(x2-x1) 
    mean_y = abs( get_line_point(x2,a,b) + get_line_point(x1,a,b)) / 2
    return Δx * mean_y
end

function get_row_volume(df_rebinned::DataFrame, thickness = 0.01)
    totalArea = sum(get_area.(df_rebinned.minE, df_rebinned.maxE, df_rebinned.a, df_rebinned.b))
    volume = totalArea*thickness
    return volume
end

function get_total_volume(df_rebinned::DataFrame, thickness = 0.01)
    volume = 0.0
    for e in unique(df_rebinned.E1)
        volume += get_row_volume(df_rebinned[df_rebinned.E1 .== e, :], thickness)
    end
    return volume
end