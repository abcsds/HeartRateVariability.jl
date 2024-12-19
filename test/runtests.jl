using HeartRateVariability
using Test
using HTTP
using Statistics
using BenchmarkTools

n=HeartRateVariability.infile("e1304.txt")
td=HeartRateVariability.time_domain(n)
fd=HeartRateVariability.frequency(n)
nl=HeartRateVariability.nonlinear(n)
g=HeartRateVariability.geometric(n)
gp=HeartRateVariability.geometric_plots(n)

@testset verbose=true "HeartRateVariability.jl" begin

    @testset "HeartRateVariability.infile" begin
        @test HeartRateVariability.infile("e1304","atr")==n
    end

    @testset "HeartRateVariability.time_domain" begin
        @test td.mean≈917.24 atol=0.1
        @test td.median≈968.0 atol=0.1
        @test td.range≈808.0 atol=0.1
        @test td.sdnn≈137.19 atol=0.1
        @test td.rmssd≈27.85 atol=0.1
        @test td.sdsd≈27.85 atol=0.1
        @test td.nn50≈342 atol=1
        @test td.pnn50≈4.41 atol=0.1
        @test td.nn20≈2831 atol=1
        @test td.pnn20≈36.53 atol=0.1
        @test td.rRR≈2.67 atol=0.1
        @test td.cvsd≈0.0303 atol=0.1
        @test td.mean_hr≈65.4137 atol=0.1
        @test td.sd_hr≈14.0535 atol=0.1
        @test td.max_hr≈178.5714 atol=0.1
        @test td.min_hr≈52.4476 atol=0.1
    end

    @testset "HeartRateVariabiliy.time_domain.verify" begin
        @test td.mean≈mean(n) atol=0.1
        @test td.median≈median(n) atol=0.1
        @test td.range≈maximum(n) - minimum(n) atol=0.1
        @test td.sdnn≈std(n) atol=0.1
        @test td.rmssd≈sqrt(mean(diff(n).^2)) atol=0.1
        @test td.sdsd≈std(diff(n)) atol=0.1
        @test td.nn50≈sum(abs.(diff(n)) .> 50) atol=0.1
        @test td.pnn50≈sum(abs.(diff(n)) .> 50) / length(n) * 100 atol=0.1
        @test td.nn20≈sum(abs.(diff(n)) .> 20) atol=0.1
        @test td.pnn20≈sum(abs.(diff(n)) .> 20) / length(n) * 100 atol=0.1
        rr = [2*(n[i]-n[i-1])/(n[i]+n[i-1]) for i in 2:length(n)]
        d = [sqrt((Statistics.mean(rr)-rr[i])^2+(Statistics.mean(rr)-rr[i+1])^2) for i in 1:length(rr)-1]
        @test td.rRR≈Statistics.median(d)*100 atol=0.1
        @test td.cvsd≈sqrt(mean(diff(n).^2)) / mean(n) atol=0.1
        @test td.mean_hr≈60000/mean(n) atol=0.1
        @test td.sd_hr≈std(60000 ./ n) atol=0.1
        @test td.max_hr≈60000/minimum(n) atol=0.1
        @test td.min_hr≈60000/maximum(n) atol=0.1
    end

    @testset "HeartRateVariability.frequency" begin
        @test fd.vlf≈1317.96 atol=0.01*fd.vlf
        @test fd.lf≈90.36 atol=0.01*fd.lf
        @test fd.hf≈186.93 atol=0.01*fd.hf
        @test fd.lfhf_ratio≈0.48 atol=0.01*fd.lfhf_ratio
        @test fd.tp≈1584.35 atol=0.01*fd.tp
    end

    @testset "HeartRateVariability.nonlinear" begin
        @test nl.apen≈2.16 atol=0.1
        @test nl.sampen≈2.16 atol=0.1
        @test nl.hurst≈0.37 atol=0.1
        @test nl.renyi0≈-6.82 atol=0.1
        @test nl.renyi1≈-6.83 atol=0.1
        @test nl.renyi2≈-6.84 atol=0.1

        #testing if get_rs from module Nonlinear returns 0 when S or R is 0
        @test HeartRateVariability.Nonlinear.get_rs(ones(100))==0

        #testing if warning is thrown
        @test_logs (:warn,"To obtain a valid value for the hurst coefficient, the length of the data series must be greater than or equal to 100.") HeartRateVariability.nonlinear([1.0,2.0,1.0,2.0,1.0,2.0,1.0,2.0,1.0,2.0])
    end
    @testset "HeartRateVariability.geometric" begin
        @test g.sd1≈19.695 atol=0.1
        @test g.sd2≈193.017 atol=0.1
        @test g.sd2_sd1≈9.8005 atol=0.1
        @test g.sd1_sd2_area≈11943.4333 atol=0.1
        @test g.csi≈9.8 atol=0.1
        @test g.cvi≈4.7840 atol=0.1
        @test g.ccsi≈30268.0552 atol=0.1
        @test g.ti≈20.8868 atol=0.1
        @test g.tinn≈288 atol=0.1
    end
    @testset "HeartRateVariability.geometric_plots" begin
        @test gp.poincare!=nothing
        @test gp.recurrence!=nothing
        @test_throws ErrorException HeartRateVariability.geometric_plots(n,"error")
    end

    @testset "rr-interval-healthy-subjects" begin
        """
        """
        url = "https://physionet.org/files/rr-interval-healthy-subjects/1.0.0/000.txt"
        hrv = parse.(Float64, filter!(e->e!="", split(String(HTTP.get(url).body), r"[^\d.]")))
        n = hrv[1:1000]
        td = HeartRateVariability.time_domain(n)
        fd = HeartRateVariability.frequency(n)
        nl = HeartRateVariability.nonlinear(n)
        g = HeartRateVariability.geometric(n)
        gp = HeartRateVariability.geometric_plots(n)

        @testset "rr-interval-healthy-subjects.time_domain" begin
            @test td.mean≈mean(n) atol=0.1
            @test td.median≈median(n) atol=0.1
            @test td.sdnn≈std(n) atol=0.1
            @test td.rmssd≈sqrt(mean(diff(n).^2)) atol=0.1
            @test td.sdsd≈std(diff(n)) atol=0.1
            @test td.rmssd≈td.sdsd atol=5
            @test td.nn50≈sum(abs.(diff(n)) .> 50) atol=1
            @test td.pnn50≈sum(abs.(diff(n)) .> 50) / length(n) * 100 atol=0.1
            @test td.nn20≈sum(abs.(diff(n)) .> 20) atol=1
            @test td.pnn20≈sum(abs.(diff(n)) .> 20) / length(n) * 100 atol=0.1
            rr = [2*(n[i]-n[i-1])/(n[i]+n[i-1]) for i in 2:length(n)]
            d = [sqrt((Statistics.mean(rr)-rr[i])^2+(Statistics.mean(rr)-rr[i+1])^2) for i in 1:length(rr)-1]
            @test td.rRR≈Statistics.median(d)*100 atol=0.01
            @test td.cvsd≈sqrt(mean(diff(n).^2)) / mean(n) atol=0.1
            @test td.mean_hr≈60000/mean(n) atol=0.1
            @test td.sd_hr≈std(60000 ./ n) atol=0.1
            @test td.max_hr≈60000/minimum(n) atol=0.1
            @test td.min_hr≈60000/maximum(n) atol=0.1

            @test td.mean≈932.165 atol=0.1
            @test td.median≈922.0 atol=0.1
            @test td.sdnn≈65.3648 atol=0.1
            @test td.rmssd≈57.4586 atol=0.1        
            @test td.sdsd≈57.4873 atol=0.1
            @test td.nn50≈298 atol=1
            @test td.pnn50≈29.7999 atol=0.1
            @test td.nn20≈708 atol=1
            @test td.pnn20≈70.8 atol=0.1
            @test td.rRR≈5.9890 atol=0.1
            @test td.cvsd≈0.0616 atol=0.1
            @test td.mean_hr≈64.3663 atol=0.1
            @test td.sd_hr≈4.6177 atol=0.1
            @test td.max_hr≈116.2791 atol=0.1
            @test td.min_hr≈47.6947 atol=0.1
        end
        @testset "rr-interval-healthy-subjects.frequency" begin
            @test fd.vlf≈HeartRateVariability.Frequency.get_power(HeartRateVariability.Frequency.lomb_scargle(n).freq, HeartRateVariability.Frequency.lomb_scargle(n).power, 0.003, 0.04) atol=0.01*fd.vlf
            @test fd.lf≈HeartRateVariability.Frequency.get_power(HeartRateVariability.Frequency.lomb_scargle(n).freq, HeartRateVariability.Frequency.lomb_scargle(n).power, 0.04, 0.15) atol=0.01*fd.lf
            @test fd.hf≈HeartRateVariability.Frequency.get_power(HeartRateVariability.Frequency.lomb_scargle(n).freq, HeartRateVariability.Frequency.lomb_scargle(n).power, 0.15, 0.4) atol=0.01*fd.hf
            @test fd.lfhf_ratio≈fd.lf / fd.hf atol=0.01*fd.lfhf_ratio
            @test fd.tp≈fd.vlf + fd.lf + fd.hf atol=0.01*fd.tp

            @test fd.vlf≈713.8915 atol=0.01*fd.vlf
            @test fd.lf≈197.8093 atol=0.01*fd.lf
            @test fd.hf≈336.9385 atol=0.01*fd.hf
            @test fd.lfhf_ratio≈0.58707842 atol=0.01*fd.lfhf_ratio
            @test fd.tp≈1248.6395 atol=0.01*fd.tp
        end
        @testset "rr-interval-healthy-subjects.nonlinear" begin
            # @test nl.apen≈HeartRateVariability.Nonlinear.apen(n) atol=0.1
            # @test nl.sampen≈HeartRateVariability.Nonlinear.sampen(n) atol=0.1
            # @test nl.hurst≈HeartRateVariability.Nonlinear.hurst(n) atol=0.1
            # @test nl.renyi0≈HeartRateVariability.Nonlinear.renyi(n, 0) atol=0.1
            # @test nl.renyi1≈HeartRateVariability.Nonlinear.renyi(n, 1) atol=0.1
            # @test nl.renyi2≈HeartRateVariability.Nonlinear.renyi(n, 2) atol=0.1

            @test nl.apen≈2.8339 atol=0.1
            @test nl.sampen≈2.8471 atol=0.1
            @test nl.hurst≈0.3712 atol=0.1
            @test nl.renyi0≈-6.8375 atol=0.1
            @test nl.renyi1≈-6.8399 atol=0.1
            @test nl.renyi2≈-6.8424 atol=0.1
        end
        @testset "rr-interval-healthy-subjects.geometric" begin
            @test g.sd1≈40.6497 atol=0.1
            @test g.sd2≈82.9304 atol=0.1
            @test g.sd2_sd1≈2.0401 atol=0.1
            @test g.sd1_sd2_area≈10590.6196 atol=0.1
            @test g.csi≈2.0401 atol=0.1
            @test g.cvi≈4.7319 atol=0.1
            @test g.ccsi≈2707.0118 atol=0.1
            @test g.ti≈15.3846 atol=0.1
            @test g.tinn≈232 atol=0.1
        end
        @testset "rr-interval-healthy-subjects.geometric_plots" begin
            @test gp.poincare!=nothing
            @test gp.recurrence!=nothing
            @test_throws ErrorException HeartRateVariability.geometric_plots(n, "error")
        end
        @testset "rr-interval-healthy-subjects.windowed_analysis" begin
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.mean))≈934.8574 atol=0.1
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.median))≈937.7632 atol=0.1
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.range))≈266.0 atol=0.1
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.sdnn))≈46.6457 atol=0.1
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.rmssd))≈936.1988 atol=0.1
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.sdsd))≈46.6457 atol=0.1
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.nn50))≈100.0 atol=0.1
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.pnn50))≈99.0099 atol=0.1
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.nn20))≈100.0 atol=0.1
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.pnn20))≈99.0099 atol=0.1
            rr = [2*(n[i]-n[i-1])/(n[i]+n[i-1]) for i in 2:length(n)]
            d = [sqrt((Statistics.mean(rr)-rr[i])^2+(Statistics.mean(rr)-rr[i+1])^2) for i in 1:length(rr)-1]
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.rRR))≈Statistics.median(d)*100 atol=0.1
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.cvsd))≈sqrt(Statistics.mean(diff(n).^2)) / Statistics.mean(n) atol=0.1
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.mean_hr))≈60000/Statistics.mean(n) atol=0.1
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.sd_hr))≈3.3868 atol=0.1
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.max_hr))≈77.4731 atol=0.1
            @test Statistics.mean(Preprocessing.windowed(n; window_size=100, stride=50, f=HeartRateVariability.TimeDomain.min_hr))≈57.0479 atol=0.1
        end
    end

    @testset "Preprocessing" begin
        @testset "Preprocessing.replacing" begin
            @test isequal(Preprocessing.replace_zeros([1, 2, 0, 4, 5, 0, 7]),Float64[1.0, 2.0, NaN, 4.0, 5.0, NaN, 7.0])
            @test isequal(Preprocessing.replace_bio_outliers([400, 500, 200, 2000, 1000, 3000, 1500, 100, 5000]),Float64[400.0, 500.0, NaN, 2000.0, 1000.0, NaN, 1500.0, NaN, NaN])
            N = 100
            a_dist = Float64.(rand(600:1200, N))
            q_high = quantile(a_dist, 0.75)
            q_low = quantile(a_dist, 0.25)
            ridx = [e<q_low || e>q_high for e in a_dist]
            goal = copy(a_dist)
            goal[ridx] .= repeat([NaN], sum(ridx))
           @test isequal(Preprocessing.replace_statistical_outliers(a_dist; low=0.25, high=0.75),goal)
        end
        @testset "Preprocessing.interpolating" begin
            @test isequal(Preprocessing.interpolate(Float64[1, 2, NaN, 4, 5, NaN, 7.0]; method=:constant),Float64[1.0, 2.0, 2.0, 4.0, 5.0, 5.0, 7.0])
            @test isequal(Preprocessing.interpolate(Float64[1, 2, NaN, 4, 5, NaN, 7.0]; method=:linear),Float64[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
            @test isequal(Preprocessing.interpolate(Float64[1, 2, NaN, 4, 5, NaN, 7.0]; method=:quadratic),Float64[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
            @test isequal(Preprocessing.interpolate(Float64[1, 2, NaN, 4, 5, NaN, 7.0]; method=:cubic),Float64[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
        end
        @testset "Preprocessing.windowed" begin
            v = Float64[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            @test isequal(Preprocessing.windowed(v; window_size=3), [[1., 2., 3.], [2., 3., 4.], [3., 4., 5.], [4., 5., 6.], [5., 6., 7.], [6., 7., 8.], [7., 8., 9.], [8., 9., 10.]])
            @test isequal(Preprocessing.windowed(v; window_size=3, stride=3), [[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]])
            @test isequal(Preprocessing.windowed(v; window_size=3, stride=3, f=Statistics.mean), [2., 5., 8.])
        end
    end

    @testset "Types" begin
        @testset "Types.Float64" begin
            N = 1000
            n = Float64.(randn(N) .* 100 .+ 1000)
            td = HeartRateVariability.time_domain(n)
            fd = HeartRateVariability.frequency(n)
            nl = HeartRateVariability.nonlinear(n)
            g = HeartRateVariability.geometric(n)
            gp = HeartRateVariability.geometric_plots(n)

            @testset "Types.Float64.time_domain" begin
                @test td.mean≈1000 atol=10
            end
        end
    end

    @testset "LongMeasurements" begin
        url = "https://physionet.org/files/rr-interval-healthy-subjects/1.0.0/4016.txt"
        hrv = parse.(Float64, filter!(e->e!="", split(String(HTTP.get(url).body), r"[^\d.]")))
        n = hrv
        td = HeartRateVariability.time_domain(n)
        fd = HeartRateVariability.frequency(n)
        # nl = HeartRateVariability.nonlinear(n) # Very slow!
        g = HeartRateVariability.geometric(n)

        @testset "LongMeasurements.time_domain" begin
            @test td.sdann≈49.6634 atol=0.1
        end
    end
end
