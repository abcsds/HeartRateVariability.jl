module Preprocessing

import Interpolations
import Statistics

#=
This function turns any zero value into a NaN value.
:param n: the array that contains the NN-intervals
:return: the array without zeros
=#
function replace_zeros(n::Array{T,1}) where T<:Real
    return Float64[e==0 ? NaN : e for e in n]
end # replace_zeros

#=
This function turns any rr-interval outside biologically plausible limits, i.e., less than 300 ms or greater than 2000 ms, into a NaN value.
:param n: the array that contains the NN-intervals
:param min: the minimum value of the interval, default=300 (200 BPM)
:param max: the maximum value of the interval, default=2000 (30 BPM)
:return: the array without the intervals that are less than 300 ms or greater than 2000 ms
=#
function replace_bio_outliers(n::Array{T,1};min=300,max=2000) where T<:Real
    return Float64[e<min || e>max ? NaN : e for e in n]
end # replace_bio_outliers

#=
This function turns any rr-interval outside the 2.5th and 97.5th percentiles into a NaN value.
:param n: the array that contains the NN-intervals
:return: the array without the intervals outside the 2.5th and 97.5th percentiles
=#
function replace_statistical_outliers(n::Array{T,1};low::Float64=0.025,high::Float64=0.975) where T<:Real
    return Float64[e<Statistics.quantile(n,low) || e>Statistics.quantile(n,high) ? NaN : e for e in n]
end # replace_statistical_outliers

#=
This function interpolates nan values
:param n: the array that contains the NN-intervals
:param method: the interpolation method, default=:linear (options: :linear, :spline)
:return: the array with the interpolated values
=#
function interpolate(n::Array{Float64,1}; method::Symbol=:linear)
    x=1:length(n)
    # itp=Interpolations.LinearInterpolation(x,n)
    # return [isnan(e) ? itp(i) : e for (i,e) in enumerate(n)]
    # Check if the number of missing values is valid
    n_nan = count(isnan, n)
    @info "Number of NaN values: $n_nan"
    if n_nan > length(n) ÷ 2
        throw(ArgumentError("Too many missing values: $n_nan out of $(length(n))"))
    end

    # Identify indices of valid values and NaN values
    valid_indices = findall(!isnan, n)
    NaN_idx = findall(isnan, n)

    # Extract valid values
    valid_values = n[valid_indices]
    x = collect(1:length(n))
    println("Length of valid values: ", length(valid_values))
    println("Length of invalid_values: ", length(NaN_idx))
    # Select interpolation method
    if method == :constant
        itp = Interpolations.ConstantInterpolation(x[valid_indices], valid_values) # Nearest neighbor
    elseif method == :linear
        itp = Interpolations.LinearInterpolation(x[valid_indices], valid_values)
    elseif method == :quadratic
        # itp = Interpolations.interpolate([x[valid_indices] valid_values], Interpolations.BSpline(Interpolations.Quadratic(Interpolations.InPlace(Interpolations.OnCell()))))
        throw(ArgumentError("Quadratic interpolation is not implemented yet"))
    elseif method == :cubic
        # itp = Interpolations.interpolate(valid_values, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Free(Interpolations.OnCell()))))
        throw(ArgumentError("Cubic interpolation is not implemented yet"))
    else
        throw(ArgumentError("Unsupported interpolation method: $method"))
    end
    # Replace missing values with interpolated values
    println("Interpolated values: ", itp(NaN_idx))
    n[NaN_idx] .= itp(NaN_idx)
    return n
end # interpolate

end # module