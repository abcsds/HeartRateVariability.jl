module TimeDomain

import Statistics

#=
This function calculates the differences between the NN intervals
:param n: is the array that contains the NN-intervals
:return diff: is an array containing the differences
=#
function nn_diff(n)
    diff=[]
    for i in 1:length(n)-1
        push!(diff,n[i+1]-n[i])
    end
    return diff
end #nn_diff

#=
This function calculates the standard deviation of the NN intervals
:param n: is the array that contains the NN-intervals
:return: the standard deviation
=#
function sdnn(n)
    return Statistics.std(n)
end # sdnn

#=
This function calculates the root mean square of successive differences
:param diff: is the array containing the differences between the NN intervals
:return: the rmssd
=#
function rmssd(diff)
    return sqrt(Statistics.mean(diff.^2))
end # rmssd

#=
This function calculates the standard deviation of successive differences
:param diff: is the array containing the differences between the NN intervals
:return: the sdsd
=#
function sdsd(diff)
    return Statistics.std(diff)
end # sdsd

#=
This function calculates the percentage of successive NN intervals,
with an interval smaller than x ms
:param diff: is the array containing the differences between the NN intervals
:param x: is the number of milliseconds the intervals may differ
:return: the percentage of successive intervals with a difference < x ms
=#
function pnn(diff,x)
    return nn(diff,x)/(length(diff)+1)*100
end # pnn

#=
This function calculates the number of successive NN intervals,
with an interval smaller than x ms
:param diff: is the array containing the differences between the NN intervals
:param x: is the number of milliseconds the intervals may differ
:return: the number of successive intervals with a difference < x ms
=#
function nn(diff,x)
    count=0
    for d in diff
        if abs(d)>x
            count+=1
        end
    end
    return count
end # nn

#=
This function calculates the mean of the NN intervals
:param n: is the array that contains the NN-intervals
:return: the mean value
=#
function mean_nn(n)
    return Statistics.mean(n)
end # mean_nn

#=
This function calculates the median of the NN intervals
:param n: is the array that contains the NN-intervals
:return: the median value
=#
function median_nn(n)
    return Statistics.median(n)
end # median_nn

#=
This function calculates the range of the NN intervals
:param n: is the array that contains the NN-intervals
:return: the range
=#
function range_nn(n)
    return Statistics.maximum(n)-Statistics.minimum(n)
end # range_nn

#=
This function calculates the relative RR
:param n: is the array that contains the NN-intervals
:return: the relative RR
=#
function rRR(n)
    rr=[]
    for i in 2:length(n)
        r=(2*(n[i]-n[i-1])/(n[i]+n[i-1]))
        push!(rr,r)
    end
    m=sum(rr)/length(rr)
    d=[]
    for i in 1:length(rr)-1
        push!(d,sqrt((m-rr[i])^2+(m-rr[i+1])^2))
    end
    return Statistics.median(d)*100
end # rRR

#=
This function calculates the CVSD of the NN intervals
:param n: is the array that contains the NN-intervals
:return: the CVSD
=#
function cvsd(n)
    return sqrt(Statistics.mean(diff(n).^2))/Statistics.mean(n)
end # cvsd

#=
This function calculates the mean heart rate from the NN intervals
:param n: is the array that contains the NN-intervals
:return: the mean heart rate
=#
function mean_hr(n)
    return 60000 / Statistics.mean(n)
end # mean_hr

#=
This function calculates the standard deviation of the heart rate from the NN intervals
:param n: is the array that contains the NN-intervals
:return: the standard deviation of the heart rate
=#
function sd_hr(n)
    return Statistics.std(60000 ./ n)
end # sd_hr

#=
This function calculates the maximum heart rate from the NN intervals
:param n: is the array that contains the NN-intervals
:return: the maximum heart rate
=#
function max_hr(n)
    return 60000 / Statistics.minimum(n)
end # max_hr

#=
This function calculates the minimum heart rate from the NN intervals
:param n: is the array that contains the NN-intervals
:return: the minimum heart rate
=#
function min_hr(n)
    return 60000 / Statistics.maximum(n)
end # min_hr

end # module
