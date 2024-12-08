# Introduction

# Time-Domain Analysis

The time-domain analysis contains the following analysis methods:

###### Mean:
This is the mean value of the RR intervals. It is calculated by summing all NN intervals and then dividing by their number. [Read more](https://en.wikipedia.org/wiki/Mean#Arithmetic_mean_(AM))

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

###### Poincaré plot
This plot is used to quantify self-similarity in processes. [Read more](https://en.wikipedia.org/wiki/Poincar%C3%A9_plot)

###### Recurrence plot
This plot is used to visualize the periodic nature of a trajectory through a phase space. [Read more](https://en.wikipedia.org/wiki/Recurrence_plot)
