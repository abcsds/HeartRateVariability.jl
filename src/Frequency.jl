module Frequency

import LombScargle
import Trapz

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
    p=Trapz.trapz(freq[index[1]:index[end]],power[index[1]:index[end]])
    return p
end # get_power

end # module
