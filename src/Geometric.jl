module Geometric

import Plots
import Images
import Statistics
import StatsBase

#=
This function creates a Poincaré plot
:param n: is the array that contains the NN-intervals
:return: a plot object
=#
function poincare(n)
    return Plots.scatter(n[1:end-1],n[2:end],xlabel="RRn",ylabel="RRn+1",legend=false);
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

#=
This function returns the histogram of the array x given a range of bins.
:param x: the array that contains the NN-intervals
:param bins: the range of bins
=#
function histogram(x, bins)
    h = StatsBase.fit(StatsBase.Histogram, x, bins)
    return h
end

#=
This function calculates the triangular index of the HRV: Integral of the density of the RR interval histogram
divided by its height at the mode.
:param n: the array that contains the NN-intervals
=#
function triangular_index(n)
    return length(n) / maximum(histogram(n, range(300, 2000, step=8)).weights)
end # triangular_index

#=
This function calculates the TINN index of the HRV: The width of the RR interval histogram.
[Read more](https://www.ahajournals.org/doi/full/10.1161/01.CIR.93.5.1043)[^1]
[^1]: Heart Rate Variability | Circulation. (n.d.). Retrieved December 17, 2024, from https://www.ahajournals.org/doi/full/10.1161/01.CIR.93.5.1043
:param n: the array that contains the NN-intervals
=#
function tinn(nn)
    h = histogram(nn, range(300, 2000, step=8))
    iX = argmax(h.weights)
    X = [h.edges[1];][iX]
    Y = maximum(h.weights)
    N_range = X<=300 ? [300] : range(300, X, step=8)
    M_range = X>=2000 ? [2000] : range(X, 2000, step=8)
    min_sse = 1e10
    vars = []
    for (i,n) in enumerate(N_range)
        for (j,m) in enumerate(M_range)
            D_edges = h.edges[1]
            D_weights = h.weights
            # @assert length(D_edges) == length(D_weights)
            iN = findfirst(x->x==n, D_edges)
            iM = findfirst(x->x==m, D_edges)
            @assert iM <= length(D_edges)

            Q_weights = zeros(length(D_weights))
            a = iX - iN
            b = iM - iX
            a == 0 && (redg = [])
            a == 1 && (Q_weights[iX] = Y)
            a >= 2 && (Q_weights[iN:iX] = [i for i in LinRange(0, Y, iX-iN+1)])

            b == 0 && (fedg = [])
            b == 1 && (Q_weights[iX] = Y)
            b >= 2 && (Q_weights[iX:iM-1] = [i for i in LinRange(Y, 0, iM-iX)])
            Q_edges = D_edges
            @assert length(D_weights) == length(Q_weights)
            # Catch Float overflow
            sse = Inf
            try
                sse = sum((Q_weights .- D_weights) .^ 2)
            catch
                @debug "Overflow, M: $m, N: $n"
            end
            if sse < min_sse
                min_sse = sse
                vars = [n, m, i, j]
            end
        end
    end
    N, M = vars[1], vars[2]
    return M - N
end # tinn

end # module
