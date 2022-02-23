"""
#### function ```get_area(x1::Float64, x2::Float64, a::Float64, b::Float64)```
<br>

    Description of ```get_area```
    ----------------------------------------
Returns the area under trapezoid, where "roof" is defined by a linear function y = ax +b. 

"""
function get_area(x1::Float64, x2::Float64, a::Float64, b::Float64)
    # area under the trapezoid is A = Δx * ⟨y⟩
    # where Δx is the width of the base of trapezoid
    # and ⟨y⟩ is mean height
    
    Δx = abs(x2-x1) 
    mean_y = abs( get_line_point(x2,a,b) + get_line_point(x1,a,b)) / 2
    return Δx * mean_y
end

"""
#### function ```get_row_volume(df_row::DataFrame, thickness = 0.001)```
<br>

    Description of ```get_row_volume```
    ----------------------------------------
Returns the volume of trapezoid, where "roof" is defined by a linear function y = ax +b and thickness is the step between rows - ΔE. 
```df_row``` is the df subset where df.E1 .== e1, e1 is the constant energy for which the linear approximation was made. 

"""
function get_row_volume(df_row::DataFrame, thickness = 0.001)
    totalArea = sum(get_area.(df_row.minE, df_row.maxE, df_row.a, df_row.b))
    volume = totalArea*thickness
    return volume
end

"""
#### function ```get_total_volume(df_rebinned::DataFrame, thickness = 0.001)```
<br>

    Description of ```get_total_volume```
    ----------------------------------------
Returns the total volume as a sum of all rows. Thickness is the step between energies, ΔE. 

"""
function get_total_volume(df_rebinned::DataFrame, thickness = 0.001)
    volume = 0.0
    for e in unique(df_rebinned.E1)
        volume += get_row_volume(df_rebinned[df_rebinned.E1 .== e, :], thickness)
    end
    return volume
end

"""
#### function ```get_segment_volume(df::DataFrame, X1::Float64, X2::Float64, thickness = 0.001)```
<br>

    Description of ```get_segment_volume```
    ----------------------------------------
Returns the volume of a segment enclosing x1, x2. 

"""
function get_segment_volume(df::DataFrame, X1::Float64, X2::Float64, thickness = 0.001)
    X1 < minimum(df.minE) && (X1 = minimum(df.minE))  #if x1 is less than the lowest minE value in df, set X1 to the lowest value 
    X2 > maximum(df.maxE) && (X2 = maximum(df.maxE))  #if x2 is more than the highest minE value in df, set X2 to the highest value 
    

    df_segment = subset(df, :minE => e -> e .<= X1, :maxE => e -> e .> X1, view = :true) # return df subset where minE < X1 < maxE

    if nrow(df_segment) == 0
        return 0.0
    end

    A = 0.0
     
    if(X1 > df_segment.minE[1]) && (X2 < df_segment.maxE[1] )                     # if the bounds are within one line region
        A += get_area(X1, X2, df_segment.a[1], df_segment.b[1])
        return A*thickness
    end
    
    while df_segment.maxE[1] - X2 < 0
        A += get_area(X1, df_segment.maxE[1], df_segment.a[1], df_segment.b[1])
        X1 = df_segment.maxE[1]
        df_segment = subset(df, :minE => e -> e .<= X1, :maxE => e -> e .> X1, view = :true)

    end
    A += get_area(df_segment.minE[1], X2, df_segment.a[1], df_segment.b[1])
    return A*thickness
end


"""
#### function ```get_volume_matrix(df::DataFrame, eMin, eMax, dE, thickness = 0.001)```
<br>

    Description of ```get_volume_matrix```
    ----------------------------------------
Returns a square matrix of volumes. Each cell defines volume of a segment defined by dE1*dE2 = dE^2.

"""
function get_volume_matrix(df::DataFrame, eMin, eMax, dE, thickness = 0.001) #returns square matrix with volumes
    dims = length(minimum(df.E1):dE:maximum(df.E1))  # dimension of the square matrix
        
    volumesTruncated = zeros(dims,dims) # final volume to output with dimensions n x n 
    temp_volumes = zeros(dims, length(unique(df.E1)))  # 30x2997
    
    for (idxE1, uniqueE1) in enumerate(unique(df.E1))  #1:2997
        df_e1 = subset(df, :E1 => e -> e .== uniqueE1)
        
        for (idxEnergy, energy) in enumerate(eMin:dE:eMax) # 1:30
            temp_volumes[idxEnergy, idxE1] = get_segment_volume(df_e1, energy, energy+dE, thickness)
        end
    end
    
    n = floor(Int,length(unique(df.E1))/ dims )  # how many cols to truncate
    for r in 1:length(temp_volumes[:,1])              # truncating n cols into 1
        volumesTruncated[r,:] = [sum(temp_volumes[r,i:i+n]) for i in 1:n:length(temp_volumes[r,:])-n]
    end
    return volumesTruncated
end