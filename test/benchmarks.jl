using BenchmarkTools
import HeartRateVariability
import Statistics

url = "https://physionet.org/files/rr-interval-healthy-subjects/1.0.0/4016.txt"
hrv = parse.(Float64, filter!(e->e!="", split(String(HTTP.get(url).body), r"[^\d.]")))
n = hrv
N = length(n)
println(N)
println(sizeof(n)) # 1.8 MB

println("Preprocessing")
@btime HeartRateVariability.preprocess(n);
x = n .* (rand(N) .< 0.1);
@btime HeartRateVariability.Preprocessing.replace_zeros(x); # 103.044 μs (3 allocations: 1.73 MiB)
@btime HeartRateVariability.Preprocessing.replace_bio_outliers(n); # 101.456 μs (3 allocations: 1.73 MiB)
@btime HeartRateVariability.Preprocessing.replace_statistical_outliers(n); # 
x = n .* (rand([NaN,1,1,1,1,1,1,1],N));
@btime HeartRateVariability.Preprocessing.interpolate(x); # 

@btime HeartRateVariability.Preprocessing.windowed(n; window_size=1000, stride=1000, time=:ms, f=x->sum(x)/length(x)); # 4.433 s (603532 allocations: 2.29 GiB)
@btime HeartRateVariability.Preprocessing.windowed(n; window_size=1000, stride=1000, time=:ms, f=Statistics.mean(x));
@btime HeartRateVariability.Preprocessing.windowed(n; window_size=60, stride=1, time=:beats, f=x->sum(x)); # 3.087 ms (3 allocations: 1.73 MiB)


println("Time Domain")
@btime HeartRateVariability.time_domain(n);
@btime HeartRateVariability.TimeDomain.mean(n);
@btime HeartRateVariability.TimeDomain.median(n);
@btime HeartRateVariability.TimeDomain.range(n);
@btime HeartRateVariability.TimeDomain.sdnn(n);
@btime diff(n);
dn = diff(n)
@btime HeartRateVariability.TimeDomain.rmssd(dn);
@btime HeartRateVariability.TimeDomain.sdsd(dn);
@btime HeartRateVariability.TimeDomain.sdann(n);
@btime HeartRateVariability.TimeDomain.nn(dn,50);
@btime HeartRateVariability.TimeDomain.pnn(dn,50);
@btime HeartRateVariability.TimeDomain.nn(dn,20);
@btime HeartRateVariability.TimeDomain.pnn(dn,20);
@btime HeartRateVariability.TimeDomain.rRR(n);
@btime HeartRateVariability.TimeDomain.cvsd(n);
@btime HeartRateVariability.TimeDomain.mean_hr(n);
@btime HeartRateVariability.TimeDomain.sd_hr(n);
@btime HeartRateVariability.TimeDomain.max_hr(n);
@btime HeartRateVariability.TimeDomain.min_hr(n);

println("Frequency Domain")
@btime HeartRateVariability.frequency(n);
@btime HeartRateVariability.Frequency.lomb_scargle(n);
ls = HeartRateVariability.Frequency.lomb_scargle(n);
@btime HeartRateVariability.Frequency.get_power(ls.freq, ls.power, 0.04, 0.15);

println("Nonlinear")
@btime HeartRateVariability.nonlinear(n);
@btime HeartRateVariability.Nonlinear.apen(n,2,6); # Singled out the problem
@btime HeartRateVariability.Nonlinear.sampen(n,2,6); # Singled out the problem
@btime HeartRateVariability.Nonlinear.hurst(n);
# @btime HeartRateVariability.Nonlinear.dfa(n);
# @btime HeartRateVariability.Nonlinear.lyap(n);
@btime HeartRateVariability.Nonlinear.renyi(n,0);
@btime HeartRateVariability.Nonlinear.renyi(n,1);
@btime HeartRateVariability.Nonlinear.renyi(n,2);

println("Geometric")
@btime HeartRateVariability.geometric(n);
@btime HeartRateVariability.Geometric.sd1(n);
@btime HeartRateVariability.Geometric.sd2(n);
@btime HeartRateVariability.Geometric.sd2_sd1(n);
@btime HeartRateVariability.Geometric.sd1_sd2_area(n);
@btime HeartRateVariability.Geometric.csi(n);
@btime HeartRateVariability.Geometric.ccsi(n);
