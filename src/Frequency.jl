module Frequency

import LombScargle
import Trapz
import DSP
import ..Preprocessing
import StatsBase

#=
This function calculates a lomb scargle transformation
:param n: is the array that contains the NN-intervals
:return: the result of the lomb scargle transformation
=#
function lomb_scargle(n)
    t=cumsum(n).-n[1]
    t=t./1000
    plan=LombScargle.plan(t,n,normalization=:psd,minimum_frequency=0.003,maximum_frequency=0.4)
    return LombScargle.lombscargle(plan)
end # lomb_scargle

#=
This function calculates the power of a frequency band between two given frequencies
:param freq: The frequency of a lomb scargle transformation
:param power: The power of a lomb scargle transformation
:param min: The minimum value of the frequency band to be calculated
:param max: The maximum value of the frequency band to be calculated
:return p: The power of the frequency band
=#
function get_power(freq,power,min,max)
    index=findall(x->x>=min && x<max,freq)
    isempty(index) && return NaN
    return Trapz.trapz(freq[index[1]:index[end]],power[index[1]:index[end]])
end # get_power

function welch(n::Array{Float64,1}; method::Symbol=:linear, fs::Int=4, kwargs...)
    method in [:constant, :linear, :quadratic, :cubic] || throw(ArgumentError("Unsupported interpolation method: $method"))
    s = Preprocessing.interpolate(n; method=method, fs=fs)
    s = s .- StatsBase.mean(s)
    # Default call:
    # welch_pgram(s::AbstractVector, n=div(length(s), 8), noverlap=div(n, 2); onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window)
    return DSP.welch_pgram(s, fs=fs, kwargs...)
end

function get_power(pgram::DSP.Periodograms.Periodogram, min, max)
    freqs = pgram.freq
    power = pgram.power
    index = findall(x->x>=min && x<max, freqs)
    return Trapz.trapz(freqs[index[1]:index[end]], power[index[1]:index[end]])
end # get_power

end # module
