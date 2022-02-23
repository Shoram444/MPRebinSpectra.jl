function sample_rejection(df::DataFrame; 
    minX = minimum(df.minE), maxX = maximum(df.maxE), 
    minY = minimum(df.minG), maxY = maximum(df.maxG))

    ξ   = rand(Uniform(minX,maxX))
    η   = rand(Uniform(minY,maxY))
    a,b = get_line_params(ξ, df)

    while η > get_line_point(ξ, a, b)
        ξ   = rand(Uniform(minX,maxX))
        η   = rand(Uniform(minY,maxY))
        a,b = get_line_params(ξ, df)
    end

    return ξ
end

function sample_discrete_CDF(gamma::Float64, v::Vector{Float64}, cdf::Vector{Float64})
    # returns the index (i+1) where gamma surpasses cdf[i] see:
    # Oleg N. Vassiliev
    # Monte Carlo Methods for Radiation Transport
    # page 20

    gamma < cdf[1] && return 1 # if gamma falls in the first bin move on right away
    
    return findfirst(p -> p .>= gamma, cdf) + 1 

end    


"""
## DEPRECATED DO NOT USE! Use ```sample_energies(df::DataFrame, thickness = 0.001)```
#### function ```sample_energies_old(df::DataFrame, volumesE1::Vector{Float64})```
<br>

    Description of ```sample_energies_old```
    ----------------------------------------
Returns an array of cumulative probabilities, where pdf is given by a linear function: `` y = aᵢx + bᵢ `` 
where each i corresponds to the segment where linear approximation was used on the discrete data. (i.e. 2ββ spectrum of E1,E2 was approximated by lines).

"""
function sample_energies_old(df::DataFrame, volumesE1::Vector{Float64})
    cdf = cumsum(volumesE1)
    
    gamma = rand(Uniform())
    E1 = unique(df.E1)[sample_discrete_CDF(gamma, volumesE1, cdf)]
    E2 = sample_rejection(df[df.E1 .== E1, :])
    
    return E1, E2
end


"""
#### function ```get_cdf(df::DataFrame, thickness::Real = 0.001)```
<br>

    Description of ```get_cdf```
    ------------------------------
Returns a tuple of energies (E, E2).

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


"""
#### function ```sample_energies(df::DataFrame, thickness = 0.001)```
<br>

    Description of ```sample_energies```
    ------------------------------
Returns a tuple of energies (E1, E2) given by the distribution specified in dataframe ```df```. 

```df``` must have fields: ```:E1, :minE, :maxE, :a, :b, :cdf (given by get_cdf!(df))```

Sampling is done in the following way. 
E2 is sampled using inverse CDF method - where cdf is obtained from ```get_cdf!(df)```. 
E1 is dependent on E2 as it can only be sampled in the given row where E2 was chosen.
E1 is sampled as a random uniform number between (row of E1, next row).

<br>
It takes roughly 100s to sample 1e5 energies. 

"""
function sample_energies(df::DataFrame, thickness = 0.001)
    gamma = rand(Uniform()) # random uniform number from 0 to 1
    
    if gamma > maximum(df.cdf) || gamma < minimum(df.cdf) # sanity check that gamma isnt some crazy number
        @show "sanity check"
        sample_energies1(df)
        
    end
    
    r_id = findfirst(x -> x .>= gamma,df.cdf)
        
    r_id > 1 ? Pᵢ₋₁ = df[r_id-1, 8] : Pᵢ₋₁ = 0 
    
    df_Pᵢ= df[r_id, :]    
    
    a = 0.5*df_Pᵢ.a * thickness
    b = df_Pᵢ.b * thickness
    c = Pᵢ₋₁ - gamma - a*df_Pᵢ.minE^2 - b*df_Pᵢ.minE
    
    x1, x2 = solvequadratic(a, b, c)
    
    if (df_Pᵢ.minE < x1) && (df_Pᵢ.maxE) > x1
        return rand(Uniform(df_Pᵢ.E1, df_Pᵢ.E1+thickness)), x1  #E1 is sampled as uniform number between steps, 
                                                                #E2 is sampled by CDF
    else
        return rand(Uniform(df_Pᵢ.E1, df_Pᵢ.E1+thickness)), x2 
    end
    
end

function solvequadratic(a, b, c)
    d  = sqrt(b^2 - 4*a*c)
    return (-b - d) / (2*a), (-b + d) / (2*a)
end