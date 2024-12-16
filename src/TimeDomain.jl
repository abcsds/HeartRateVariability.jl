module TimeDomain

import Statistics

#=
This function calculates the standard deviation of the NN intervals
:param n: is the array that contains the NN-intervals
:return: the standard deviation
=#
sdnn(n) = Statistics.std(n)

#=
This function calculates the root mean square of successive differences
:param diff: is the array containing the differences between the NN intervals
:return: the rmssd
=#
rmssd(diff) = sqrt(Statistics.mean(diff.^2))

#=
This function calculates the standard deviation of successive differences
:param diff: is the array containing the differences between the NN intervals
:return: the sdsd
=#
sdsd(diff) = Statistics.std(diff)

#=
This function calculates the percentage of successive NN intervals,
with an interval smaller than x ms
:param diff: is the array containing the differences between the NN intervals
:param x: is the number of milliseconds the intervals may differ
:return: the percentage of successive intervals with a difference < x ms
=#
pnn(diff,x) = nn(diff,x)/(length(diff)+1)*100

#=
This function calculates the number of successive NN intervals,
with an interval smaller than x ms
:param diff: is the array containing the differences between the NN intervals
:param x: is the number of milliseconds the intervals may differ
:return: the number of successive intervals with a difference < x ms
=#
nn(diff,x) = sum(abs.(diff) .> x)

#=
This function calculates the mean of the NN intervals
:param n: is the array that contains the NN-intervals
:return: the mean value
=#
mean(n) = Statistics.mean(n)

#=
This function calculates the median of the NN intervals
:param n: is the array that contains the NN-intervals
:return: the median value
=#
median(n) = Statistics.median(n)

#=
This function calculates the range of the NN intervals
:param n: is the array that contains the NN-intervals
:return: the range
=#
range(n) = Statistics.maximum(n) - Statistics.minimum(n)

#=
This function calculates the relative RR
:param n: is the array that contains the NN-intervals
:return: the relative RR
=#
function rRR(n)
    rr = 2 .* diff(n) ./ (n[1:end-1] .+ n[2:end])
    m = Statistics.mean(rr)
    d = [sqrt((m-rr[i])^2+(m-rr[i+1])^2) for i in 1:length(rr)-1]
    return Statistics.median(d)*100
end # rRR

#=
This function calculates the coefficient of variation of the successive differences of the NN intervals.
:param n: is the array that contains the NN-intervals
:return: the CVSD
=#
cvsd(n) = sqrt(Statistics.mean(diff(n).^2))/Statistics.mean(n)

#=
This function calculates the mean heart rate from the NN intervals
:param n: is the array that contains the NN-intervals
:return: the mean heart rate
=#
mean_hr(n) = 60000 / Statistics.mean(n)

#=
This function calculates the standard deviation of the heart rate from the NN intervals
:param n: is the array that contains the NN-intervals
:return: the standard deviation of the heart rate
=#
sd_hr(n) = Statistics.std(60000 ./ n)

#=
This function calculates the maximum heart rate from the NN intervals
:param n: is the array that contains the NN-intervals
:return: the maximum heart rate
=#
max_hr(n) = 60000 / Statistics.minimum(n)

#=
This function calculates the minimum heart rate from the NN intervals
:param n: is the array that contains the NN-intervals
:return: the minimum heart rate
=#
min_hr(n) = 60000 / Statistics.maximum(n)

end # module
