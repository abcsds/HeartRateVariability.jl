# Introduction

# Time-Domain Analysis

The time-domain analysis contains the following analysis methods:

###### Mean:
This is the mean value of the RR intervals. It is calculated by summing all NN intervals and then dividing by their number. [Read more](https://en.wikipedia.org/wiki/Mean#Arithmetic_mean_(AM))

###### Median:
This is the median value of the RR intervals. It is calculated by sorting the NN intervals and then selecting the middle value. [Read more](https://en.wikipedia.org/wiki/Median)

###### Range:
This is the difference between the maximum and minimum NN intervals. [Read more](https://en.wikipedia.org/wiki/Range_(statistics))

###### SDNN:
This is the standard deviation of the NN intervals. [Read more](https://en.wikipedia.org/wiki/Heart_rate_variability#Time-domain_methods[36])


###### RMSSD:
This is the root mean square of the differences between successive NN intervals. [Read more](https://en.wikipedia.org/wiki/Heart_rate_variability#Time-domain_methods[36])


###### SDSD:
This is the standard deviation of the differences between successive NN intervals. [Read more](https://en.wikipedia.org/wiki/Heart_rate_variability#Time-domain_methods[36])


###### NN20/NN50:
This is the number of pairs of successive NN intervals that differ by more than 20ms/50ms. [Read more](https://en.wikipedia.org/wiki/Heart_rate_variability#Time-domain_methods[36])


###### pNN20/pNN50:
This is the percentage of pairs of successive NN intervals that differ by more than 20ms/50ms. [Read more](https://en.wikipedia.org/wiki/Heart_rate_variability#Time-domain_methods[36])


###### rRR:
The relative RR intervals are calculated using the equation\
for i=2...n
```math
rr _{i} := \frac{2*(RR_{i}-RR_{i-1})}{RR_{i}+RR_{i-1}}
```
where n is the number of RR intervals.\
The HRV is measured by the median of the euclidean distances of the relative RR intervals to the average of the relative RR intervals. [Read more](https://marcusvollmer.github.io/HRV/files/paper_method.pdf) [^1]

[^1]: Vollmer, M. (2015). A robust, simple and reliable measure of heart rate variability using relative RR intervals. 2015 Computing in Cardiology Conference (CinC), 609–612. https://doi.org/10.1109/CIC.2015.7410984

##### CVSD
The coefficient of variation of successive differences is calculated by dividing the standard deviation of the differences between successive NN intervals by the mean of the NN intervals. That is, the RMSSD divided by the mean NN intervals.

##### MeanHR
The mean heart rate is calculated by dividing 60 by the mean NN intervals.

##### SDHR
The standard deviation of the heart rate is calculated by dividing 60 by the standard deviation of the NN intervals.

##### MaxHR
The maximum heart rate is calculated by dividing 60 by the minimum NN intervals.

##### MinHR
The minimum heart rate is calculated by dividing 60 by the maximum NN intervals.

# Frequency-Domain Analysis

Frequency domain analysis uses a Lomb Scargle Transformation to determine the power spectral density of each frequency domain. The frequency bands are defined as follows:

- **VLF:** very low frequency, from 0.003 to 0.04 Hz

- **LF:** low frequency, from 0.04 to 0.15 Hz

- **HF:** high frequency, from 0.15 to 0.4 Hz

- **LF/HF:** The ratio of LF and HF

- **Total Power:** The sum of VLF, LF and HF

[Read more](https://en.wikipedia.org/wiki/Heart_rate_variability#Frequency-domain_methods[36])

# Nonlinear Analysis

###### Approximate entropy

This is a technique for quantifying the degree of regularity and unpredictability of the RR intervals. [Read more](https://en.wikipedia.org/wiki/Approximate_entropy)

###### Sample entropy
This is a modification of the approximate entropy that is used to assess the complexity of physiological time series signals. [Read more](https://en.wikipedia.org/wiki/Sample_entropy)

###### Hurst exponent
The Hurst exponent is used to measure the long-term memory of time series. [Read more](https://en.wikipedia.org/wiki/Hurst_exponent)

###### Rényi entropy
The renyi entropy is a measure of diversity and forms the basis of the concept of generalized dimensions. [Read more](https://en.wikipedia.org/wiki/R%C3%A9nyi_entropy)

# Geometric Analysis

###### SD1 and SD2
These are the standard deviations of the Poincaré plot along the line of identity and perpendicular to it, respectively. The ratio between them, and area covered by them are also calculated[^2, ^3]. [Read more](https://en.wikipedia.org/wiki/Poincar%C3%A9_plot)

[^2]: Henriques, T. S., Mariani, S., Burykin, A., Rodrigues, F., Silva, T. F., & Goldberger, A. L. (2016). Multiscale Poincaré plots for visualizing the structure of heartbeat time series. BMC Medical Informatics and Decision Making, 16(1), 17. https://doi.org/10.1186/s12911-016-0252-0
[^3]: Tayel, M. B., & AlSaba, E. I. (2015). Poincaré Plot for Heart Rate Variability. 9(9).

###### Cardiac Sympathetic Index
This is a measure of the balance between the sympathetic and parasympathetic nervous systems. It uses the ratio of the SD1 and SD2 features of the Poincaré plot[^4].

[^4]: Jeppesen, J., Beniczky, S., Johansen, P., Sidenius, P., & Fuglsang-Frederiksen, A. (2014). Using Lorenz plot and Cardiac Sympathetic Index of heart rate variability for detecting seizures for patients with epilepsy. 2014 36th Annual International Conference of the IEEE Engineering in Medicine and Biology Society, 4563–4566. https://doi.org/10.1109/EMBC.2014.6944639

###### Poincaré plot
This plot is used to quantify self-similarity in processes. [Read more](https://en.wikipedia.org/wiki/Poincar%C3%A9_plot)

###### Recurrence plot
This plot is used to visualize the periodic nature of a trajectory through a phase space. [Read more](https://en.wikipedia.org/wiki/Recurrence_plot)
