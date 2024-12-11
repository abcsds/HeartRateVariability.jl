module Geometric

import Plots
import Images
import Statistics

#=
This function creates a Poincaré plot
:param n: is the array that contains the NN-intervals
:return: a plot object
=#
function poincare(n)
    x=[]
    y=[]
    for i in 1:length(n)-1
        push!(x,n[i])
        push!(y,n[i+1])
    end
    p=Plots.scatter(x,y,xlabel="RRn",ylabel="RRn+1",legend=false);
    return p;
end # poincare

#=
This function creates a recurrence plot
:param n: is the array that contains the NN-intervals
:param e: the maximum distance between two intervals, default="mean" == the mean value of successive differences
:return: a plot object
=#
function recurrence(n,e)
    if e=="mean"
        diff=[]
        for i in 1:length(n)-1
            push!(diff,abs(n[i+1]-n[i]))
        end
        e=sum(diff)/length(diff)
    end
    x=zeros(length(n),length(n))
    for i in 1:length(n)
        for j in i:length(n)
            if sqrt((n[i]-n[j])^2)<=e
                x[i,j]=1
                x[j,i]=1
            end
        end
    end
    img=Images.Gray.(x)
    r=Plots.plot(img);
    return r;
end # recurrence

#=
This function calculates the standard deviation of projection of the 
Poincaré plot on the line perpendicular to the line of identity, SD1,
or the width of the Poincaré plot.
:param n: the array that contains the NN-intervals
=#
function sd1(n)
    # return sqrt(Statistics.std(diff(n))^2/2)
    x=[n[i] for i in 1:length(n)-1]
    y=[n[i+1] for i in 1:length(n)-1]
    sd1=sqrt(Statistics.var((x-y)/sqrt(2)))
end # sd1

#=
This function calculates the standard deviation of projection of the
Poincaré plot on the line along the line of identity, SD2, or the length
of the Poincaré plot.
:param n: the array that contains the NN-intervals
=#
function sd2(n)
    # return sqrt(2*Statistics.std(n)^2 - 0.5*Statistics.std(diff(n))^2)
    x=[n[i] for i in 1:length(n)-1]
    y=[n[i+1] for i in 1:length(n)-1]
    sd2=sqrt(Statistics.var((x+y)/sqrt(2)))
end # sd2

#=
This function calculates the ratio of SD2 and SD1.
:param n: the array that contains the NN-intervals
=#
function sd2_sd1(n)
    return sd2(n)/sd1(n)
end # sd2_sd1

#=
This function calculates the area covered by the Poincaré ellipse.
:param n: the array that contains the NN-intervals
=#
function sd1_sd2_area(n)
    return π*sd1(n)*sd2(n)
end # sd1_sd2_area

#=
This function calculates the Cardiac Sympathetic Index (CSI).
:param n: the array that contains the NN-intervals
=#
function csi(n)
    return sd2(n)/sd1(n)
end # csi

#=
This function calculates the Cardiac Vagal Index (CVI).
:param n: the array that contains the NN-intervals
=#
function cvi(n)
    return log10(sd1(n)*sd2(n)*16)
end # cvi

#=
This function calculates the corrected Cardiac Sympathetic Index (CSI).
:param n: the array that contains the NN-intervals
=#
function ccsi(n)
    return ((4*sd2(n))^2)/sd1(n)
end # ccsi

end # module
