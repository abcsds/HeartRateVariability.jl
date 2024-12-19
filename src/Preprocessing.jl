module Preprocessing

import DataInterpolations
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
    l = Statistics.quantile(n,low)
    h = Statistics.quantile(n,high)
    return Float64[e<l || e>h ? NaN : e for e in n]
end # replace_statistical_outliers

#=
This function strips any NaN values from the extremes of the array of NN-intervals.
:param n: the array that contains the NN-intervals
:return: the array without the NaN values at the extremes
=#
function strip_extremes(n::Array{Float64,1})
    return n[findfirst(!isnan, n):findlast(!isnan, n)]
end # strip_extremes

#=
This function interpolates nan values
:param n: the array that contains the NN-intervals
:param method: the interpolation method, default=:linear (options: :linear, :spline)
:return: the array with the interpolated values
=#
function interpolate(n::Array{Float64,1}; method::Symbol=:linear)
    n = strip_extremes(n)
    # Check if the number of missing values is valid
    n_nan = count(isnan, n)
    n_nan > length(n) / 2 && throw(ArgumentError("Too many missing values: $n_nan out of $(length(n))"))

    # Identify indices of valid values and NaN values
    valid_indices = findall(!isnan, n)
    NaN_idx = findall(isnan, n)

    # Extract valid values
    valid_values = n[valid_indices]
    x = Float64.(collect(1:length(n)))
    # Select interpolation method
    if method == :constant
        itp = DataInterpolations.ConstantInterpolation(x[valid_indices], valid_values)
    elseif method == :linear
        itp = DataInterpolations.LinearInterpolation(x[valid_indices], valid_values)
    elseif method == :quadratic
        itp = DataInterpolations.QuadraticInterpolation(x[valid_indices], valid_values)
    elseif method == :cubic
        itp = DataInterpolations.CubicSpline(x[valid_indices], valid_values)
    else
        throw(ArgumentError("Unsupported interpolation method: $method"))
    end
    # Replace missing values with interpolated values
    n[NaN_idx] .= itp(NaN_idx)
    return n
end # interpolate

#=
This function returns views of the array of NN-intervals in a sliding window. The sliding window can be defined in beats or milliseconds.
:param n: the array that contains the NN-intervals
:param window_size: the size of the window, default=60 (beats)
:param stride: the step size of the sliding window, default=1 (beat)
:param time: the time unit of the window, default=:beats (options: :beats, :ms)
:param f: reducing function to apply to the window, default=identity
:return: the views of the array in a sliding window, with the reducing function applied to each window
=#
function windowed(n::Array{Float64,1};window_size::Int=60, stride::Int=1, time::Symbol=:beats, f::Function=identity)
    if time == :beats
        return [f(view(n,i:i+window_size-1)) for i in 1:stride:length(n)-window_size+1]
    elseif time == :ms
        t = cumsum(n)
        max_t = t[end] # Record duration in ms
    else
        throw(ArgumentError("Unsupported time unit: $time"))
    end
    window_starts = 1:stride:max_t
    window_ends = window_starts .+ window_size
    views = Vector{Any}(undef, length(window_starts))
    for (i, (start, stop)) in enumerate(zip(window_starts, window_ends))
        idx = findall(x -> x >= start && x < stop, t)
        if !isempty(idx)
            views[i] = f(view(n, idx))
        end
    end
    return views
end # windowed

end # module