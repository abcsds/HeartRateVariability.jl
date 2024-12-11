module HeartRateVariability

include("TimeDomain.jl")
include("Input.jl")
include("Frequency.jl")
include("Nonlinear.jl")
include("Geometric.jl")

"""
    geometric_plots(n,e="mean")

Arguments:
- n: the array that contains the NN-intervals
- e: the maximum distance between two intervals, default="mean" (the mean value of the successive differences), has to be "mean" or a number

Results:
- poincare: the Poincaré plot
- recurrence: the recurrence plot
"""
function geometric_plots(n::Array{Float64,1},e="mean")
    if (e!="mean" && !isa(e,Number))
           error("e has to be a numerical value or 'mean'")
       end
    return (poincare=Geometric.poincare(n),recurrence=Geometric.recurrence(n,e))
end # geometric_plots

"""
    geometric(n)

Arguments:
- n: the array that contains the NN-intervals

Results:
- sd1: the width of the Poincaré plot
- sd2: the length of the Poincaré plot
- sd2_sd1: the ratio of sd2 and sd1
- sd1_sd2_area: the area of the Poincaré plot
- csi: cardiac sympathetic index
- cvi: cardiac vagal index
- ccsi: corrected cardiac sympathetic index
"""
function geometric(n::Array{Float64,1})
    return (sd1=Geometric.sd1(n),
            sd2=Geometric.sd2(n),
            sd2_sd1=Geometric.sd2_sd1(n),
            sd1_sd2_area=Geometric.sd1_sd2_area(n),
            csi=Geometric.csi(n),
            cvi=Geometric.cvi(n),
            ccsi=Geometric.ccsi(n))
end # geometric

"""
    nonlinear(n,m=2,r=6)

Arguments:
- n: the array that contains the NN-intervals
- m: the embedding dimension, default=2
- r: the tolerance, default=6

Results:
- apen: the approximate entropy
- sampen: the sample entropy
- hurst: the hurst exponent (only valid if the length of n is >= 100)
- renyi0, renyi1, renyi2: the Rényi entropy of order 0,1 and 2
"""
function nonlinear(n::Array{Float64,1},m::Int64=2,r::Number=6)
    if length(n)<100
        @warn("To obtain a valid value for the hurst coefficient, the length of the data series must be greater than or equal to 100.")
    end
    return (apen=Nonlinear.apen(n,m,r), sampen=Nonlinear.sampen(n,m,r),
            hurst=Nonlinear.hurst(n), renyi0=Nonlinear.renyi(n,0),
            renyi1=Nonlinear.renyi(n,1), renyi2=Nonlinear.renyi(n,2))
end # nonlinear

"""
    frequency(n)

Arguments:
- n: the array that contains the NN-intervals

Results:
- vlf: the very low-frequency power
- lf: the low-frequency power
- hf: the high-frequency power
- lfhf_ratio: the lf/hf ratio
- tp: the total power
"""
function frequency(n::Array{Float64,1})
    ls=Frequency.lomb_scargle(n)
    vlf=Frequency.get_power(ls.freq,ls.power,0.003,0.04)
    lf=Frequency.get_power(ls.freq,ls.power,0.04,0.15)
    hf=Frequency.get_power(ls.freq,ls.power,0.15,0.4)
    tp=vlf+lf+hf
    return (vlf=vlf, lf=lf, hf=hf, lfhf_ratio=lf/hf, tp=tp)
end # frequency

"""
    time_domain(n)

Arguments:
- n: the array that contains the NN-intervals

Results:
- mean: the mean value
- median: the median value of NN intervals
- range: the range of NN intervals
- sdnn: the standard deviation
- rmssd: the root mean square of successive differences
- sdsd: the standard deviation of successive differences
- nn50: the number of successive NN intervals with an interval smaller than 50 ms
- pnn50: the percentage of successive NN intervals with an interval smaller than 50 ms
- nn20: the number of successive NN intervals with an interval smaller than 20 ms
- pnn20: the percentage of successive NN intervals with an interval smaller than 20 ms
- rRR: the percentage of relative RR intervals
- cvsd: the coefficient of variation of the successive differences
- mean_hr: the mean heart rate
- sd_hr: the standard deviation of the heart rate
- max_hr: the maximum heart rate
- min_hr: the minimum heart rate
"""
function time_domain(n::Array{Float64,1})
    diff=TimeDomain.nn_diff(n)
    return (mean=TimeDomain.mean_nn(n), median=TimeDomain.median_nn(n),
            range=TimeDomain.range_nn(n), sdnn=TimeDomain.sdnn(n),
            rmssd=TimeDomain.rmssd(diff), sdsd=TimeDomain.sdsd(diff),
            nn50=TimeDomain.nn(diff,50), pnn50=TimeDomain.pnn(diff,50),
            nn20=TimeDomain.nn(diff,20), pnn20=TimeDomain.pnn(diff,20),
            rRR=TimeDomain.rRR(n), cvsd=TimeDomain.cvsd(n),
            mean_hr=TimeDomain.mean_hr(n), sd_hr=TimeDomain.sd_hr(n),
            max_hr=TimeDomain.max_hr(n), min_hr=TimeDomain.min_hr(n))
end # time_domain

"""
    infile(file)

This function reads the data from a txt or csv file.

Arguments:
- file: is the path of the input file

"""
function infile(file::String)
    return Input.read_txt(file)
end # infile

"""
    infile(record,annotator)

This function reads the data from a wbdb file.

Arguments:
- record: is the name of the record
- annotator: is the annotator of the record

!!! note
    In order to use the `infile` function for wfdb files, the WFDB Software Package from
    Pysionet is required. See [Installation](installation) for more information.
"""
function infile(record::String,annotator::String)
    return Input.read_wfdb(record,annotator)
end # infile

export nonlinear, frequency, time_domain, infile, geometric
end # module
