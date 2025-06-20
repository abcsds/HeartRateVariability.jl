module Nonlinear

import StatsBase
import Statistics
import DFA
import EntropyHub

#=
This function calculates the approximate entropy
:param n: the array that contains the NN-intervals
:param m: the embedding dimension, default=2
:param r: the tolerance, default=6
:return: the approximate entropy
=#
function _apen(n,m,r)
    c1=get_apen_dist(n,m,r)
    c2=get_apen_dist(n,m+1,r)
    @info "Using custom approximate entropy calculation"
    @info "m: $m, r: $r"
    @info "n length: $(length(n))"
    @info "c1: $c1"
    @info "c2: $c2"
    return log(c1/c2)
end # apen

function apen(n,m,r)
    # apen1, _ = EntropyHub.ApEn(n, m=m, r=r)
    # apen2, _ = EntropyHub.ApEn(n, m=m+1, r=r)
    # return -log(sum(apen1)/sum(apen2))
    apens, _ = EntropyHub.ApEn(n, m=m+1, r=r)
    @info "Using EntropyHub for approximate entropy calculation"
    @info "m: $m, r: $r"
    @info "n length: $(length(n))"
    @info "apens: $apens"
    @info "apens: $apens"
    @info "apens length: $(length(apens))"
    @info "c1: $(apens[end-1])"
    @info "c2: $(apens[end])"
    @info "apens ratio: $(apens[end-1]/apens[end])"
    @info "apens log ratio: $(log(apens[end-1]/apens[end]))"
    @info apens
    return log(apens[end-1]/apens[end])
end # apen

#=
This function calculates the sample entropy
:param n: the array that contains the NN-intervals
:param m: the embedding dimension, default=2
:param r: the tolerance, default=6
:return: the sample entropy
=#
function _sampen(n,m,r)
    c1=get_sampen_dist(n,m,r,1)
    c2=get_sampen_dist(n,m+1,r,0)
    @info "Using custom sample entropy calculation"
    @info "m: $m, r: $r"
    @info "n length: $(length(n))"
    @info "c1: $c1"
    @info "c2: $c2"
    @info "c1/c2: $(c2/c1)"
    @info "log(c2/c1): $(log(c2/c1))"
    @info "sampen: $(-log(c2/c1))"
    return -log(c2/c1)
end # sampen

function sampen(n,m,r)
    sampen1, _ = EntropyHub.SampEn(n, m=m, r=r)
    sampen2, _ = EntropyHub.SampEn(n, m=m+1, r=r)
    @info "Using EntropyHub for sample entropy calculation"
    @info "m: $m, r: $r"
    @info "n length: $(length(n))"
    @info "sampen1: $sampen1"
    @info "sampen2: $sampen2"
    @info "sampen1 length: $(length(sampen1))"
    @info "sampen2 length: $(length(sampen2))"
    @info "sampen1 last: $(sampen1[end])"
    @info "sampen2 last: $(sampen2[end])"
    @info "sampen ratio: $(sampen1[end]/sampen2[end])"
    @info "sampen log ratio: $(log(sampen1[end]/sampen2[end]))"
    return log(sampen1[end]/sampen2[end])
end # sampen

#=
This function creates a template of a given array over an embedding dimension
:param n: the array that contains the NN-intervals
:param m: the embedding dimension, default=2
:return template: the created template
=#
function get_template(n,m)
    return [n[i:i+m-1] for i in 1:length(n)-m+1]
end # get_template

#=
This function calculates the distances for the approximate entropy
:param n: the array that contains the NN-intervals
:param m: the embedding dimension, default=2
:param r: the tolerance, default=6
:return: the distance for the approximate entropy
=#
function get_apen_dist(n,m,r)
    template=get_template(n,m)
    count=zeros(length(template))
    for i in 1:length(template)
        for j in i+1:length(template)
            if maximum(abs.(template[i].-template[j]))<=r
                count[i]+=1
                count[j]+=1
            end
        end
    end
    return sum(count./(length(n)-m+1))/(length(n)-m+1)
end # get_apen_dist

#=
This function calculates the distances for the sample entropy
:param n: the array that contains the NN-intervals
:param m: the embedding dimension, default=2
:param r: the tolerance, default=6
:param l: a value to limit the for-loops
:return: the distance for the sample entropy
=#
function get_sampen_dist(n,m,r,l)
    template=get_template(n,m)
    counts=[]
    count=0
    for i in 1:length(template)-l
        for j in 1:length(template)-l
            if maximum(abs.(template[i].-template[j]))>=r || i==j
                push!(counts,count)
                count=0
            else
                count+=1
            end
        end
    end
    return sum(counts)
end # get_sampen_dist

#=
This function calculates the renyi entropy of a given order
:param n: the array that contains the NN-intervals
:param a: the order of the renyi entropy
:return: the calculated renyi entropy
=#
renyi(n,a) = StatsBase.renyientropy(n,a)

#=
This function calculates the hurst coefficient
It was inspired by the python hurst package by Dmitry A. Mottl (https://github.com/Mottl/hurst)
:param n: the array that contains the NN-intervals
:return H: the hurst coefficient
=#
function hurst(n)
    ws=Array(range(log10(10),stop=log10(length(n)),step=0.25))
    window = [round(Int64,exp10(x),RoundDown) for x in ws]
    if !(length(n) in window)
        push!(window,length(n))
        push!(ws,log10(length(n)))
    end
    RS=[]
    for w in window
        rs=[]
        for start in (range(0,stop=length(n),step=w))
            if (start+w)>length(n)
                break
            end
            RS_part= get_rs(n[start+1:start+w])
            if RS_part != 0
                push!(rs,RS_part)
            end
        end
        if length(rs)>0
            push!(RS,Statistics.mean(rs))
        end
    end
    A=Array{Float64}([ws ones(length(RS))])
    RSlog=log10.(RS)
    B=Array{Float64}(RSlog)
    H,c=A\B
    c=exp10(c)
    return H
end # hurst

#=
This function calculates the rescaled range of a time series
It was inspired by the python hurst package by Dmitry A. Mottl (https://github.com/Mottl/hurst)
:param n: the array that contains the NN-intervals
:return: the rescaled range
=#
function get_rs(n)
    incs=n[2:end].-n[1:end-1]
    mean_inc=(n[end]-n[1])/length(incs)
    deviations=incs.-mean_inc
    Z=cumsum(deviations)
    R=maximum(Z)-minimum(Z)
    S=Statistics.std(incs)
    if R==0 || S==0
        return 0
    else
        return R/S
    end
end # get_rs

#=
This function calculates the detrended fluctuation analysis
:param n: the array that contains the NN-intervals
:param window_size: the size of the window, default=10
:return: the detrended fluctuation analysis
=#
function dfa(n)
    scales, fluc = DFA.dfa(n, boxmax=64, boxmin=4, boxratio=2, overlap=0.0)
    log_scales = log10.(scales)
    log_fluc = log10.(fluc)
    ntercept, α1 = DFA.polyfit(log_scales, log_fluc)
        
    scales, fluc = DFA.dfa(n, boxmax=16, boxmin=4, boxratio=2, overlap=0.0)
    log_scales = log10.(scales)
    log_fluc = log10.(fluc)
    intercept, α2 = DFA.polyfit(log_scales, log_fluc)
    return α1, α2
end # dfa

end # module
