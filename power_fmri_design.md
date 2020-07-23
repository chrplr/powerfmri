Tests of the efficiency of various fMRI designs
===============================================

Time-stamp: <2013-01-08 17:06 christophe@pallier.org>

Note: This document is generated from a [Rmarkdown document](power_fmri_design.Rmd) that can be processed with knitr and edited with [rstudio](http://www.rstudio.com)

We present some simulations comparing the power of designs with fixed ISI to designs with varying ISI. We are interested in how precise are the estimates of the parameters associated to each condition.



```
## Loading required package: MASS
```













Basic functions
---------------------

The basic functions are `dirac.comb` (that generates 0/1 activation vector given onsets and durations of events) and `hrf` (that convolves a timeseries by the impulse HRF)


```r
totalduration <- 150
onsets <- seq(1, 120, by = 20)
durations <- 5
par(mfcol = c(2, 1))
plot(dirac.comb(onsets, durations, totalduration))
plot(hrf(onsets, durations, totalduration), type = "l")
```

![plot of chunk basic functions](figure/basic_functions.png) 



For multiple conditions, the function `create_design_matrix` reads a table with columns conditions, onsets and durations.


```r
timings <- data.frame(conditions = c(1, 2, 1, 2, 2, 1, 1, 2), onsets = seq(1, 
    64, by = 8), durations = 1)
X <- create_design_matrix(timings, totalduration = 80)
par(mfcol = c(1, 1))
plot(X[, 1], type = "l", col = "blue")
lines(X[, 2], type = "l", col = "red")
```

![plot of chunk generate design matrix](figure/generate_design_matrix.png) 


To generate paradigms, one can use higher-level function such as `generate_paradigm_fixed_SOA`


```r
a <- generate_paradigm_fixed_SOA(ncond = 5, trialpercond = 5, stimduration = 3, 
    SOA = 6, totalduration = 150)
plot_design_matrix(create_design_matrix(a, totalduration = 150))
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1.png) 


See also `generate_paradigm_fixed_SOA_adding_silences`, `generate_paradigm_varying_SOA`, `generate_paradigm_varying_SOA_null_events`.

Parameters for the simulations
------------------------------


```r

ncond <- 4  # number of conditions
trialpercond <- 20
stimduration <- 1  # in sec
SOA <- 4  # in sec
(totalduration <- (ncond * trialpercond * SOA) + SOA)
```

```
## [1] 324
```


Contrasts of interest


```r
(listcon <- list(linear = (1:ncond) - mean(1:ncond), firstbeta = c(1, rep(0, 
    ncond - 1))))
```

```
## $linear
## [1] -1.5 -0.5  0.5  1.5
## 
## $firstbeta
## [1] 1 0 0 0
```

The number of simulations for each class of designs:


```r
nsim = 10
```


The simulations will use normal noise


```r
normalnoise <- function(npoints) {
    5 * rnorm(npoints)
}
```





Fixed SOA, no silent trial
--------------------------

Let us generate a paradigm with a fixed SOA...


```r
timing <- generate_paradigm_fixed_SOA(ncond, trialpercond, stimduration, SOA, 
    totalduration)
plot_paradigm(timing)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


... and build its associated design matrix


```r
X <- create_design_matrix(timing, totalduration)
plot_design_matrix(X)
```

![plot of chunk generateDesignMat](figure/generateDesignMat.png) 


Now,  for many such designs, we simulate a voxel where the signal increases in a linear fashion with 'condition' (amplitudes=1:5).


```r
o <- simulations(nsim, beta = 1:ncond, listofcontrasts = listcon, normalnoise, 
    generate_paradigm_fixed_SOA, ncond, trialpercond, stimduration, SOA, totalduration)
```



```r
report <- function(outputsim) {
    par(las = 1)
    boxplot(outputsim$estimates, main = "parameter estimates")
    grid()
    boxplot(outputsim$estimatese, main = "standard errors")
    grid()
    for (col in 1:ncol(outputsim$efficiencies)) boxplot(outputsim$efficiencies[, 
        col], main = paste("Efficiency", names(listcon)[col]), horizontal = TRUE)
    grid()
}
```




```r
report(o)
```

![plot of chunk report](figure/report1.png) ![plot of chunk report](figure/report2.png) ![plot of chunk report](figure/report3.png) ![plot of chunk report](figure/report4.png) 


Designs with jitter between trials
----------------------------------

Now, we jitter the SOA between trials and run a similar simulation.


```r
timing <- generate_paradigm_varying_SOA(ncond, trialpercond, stimduration, SOA, 
    jitter = 1, totalduration)
plot_paradigm(timing)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-31.png) 

```r
X <- create_design_matrix(timing, totalduration)
plot_design_matrix(X)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-32.png) 


histogram of SOAs:


```r
hist(diff(timing$onsets))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 




```r
o2 <- simulations(nsim, beta = 1:ncond, listofcontrasts = listcon, normalnoise, 
    generate_paradigm_varying_SOA, ncond, trialpercond, stimduration, 4, SOA, 
    totalduration)
```



```r
report(o2)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-51.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-52.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-53.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-54.png) 




Design with fixed SOA and 'explicit' NULL events
-------------------------------------------------



```r

timing <- generate_paradigm_fixed_SOA_adding_silences(ncond, trialpercond, stimduration, 
    SOA, totalduration)
plot_paradigm(timing)
```

![plot of chunk simulfun_fixed_SOA_and_silences](figure/simulfun_fixed_SOA_and_silences1.png) 

```r
totalduration <- max(timing$onsets) + SOA
X <- create_design_matrix(timing, totalduration)
plot_design_matrix(X)
```

![plot of chunk simulfun_fixed_SOA_and_silences](figure/simulfun_fixed_SOA_and_silences2.png) 

```r

o3 <- simulations(nsim, beta = 1:ncond, listofcontrasts = listcon, normalnoise, 
    generate_paradigm_fixed_SOA_adding_silences, ncond, trialpercond, stimduration, 
    SOA, totalduration)

report(o3)
```

![plot of chunk simulfun_fixed_SOA_and_silences](figure/simulfun_fixed_SOA_and_silences3.png) ![plot of chunk simulfun_fixed_SOA_and_silences](figure/simulfun_fixed_SOA_and_silences4.png) ![plot of chunk simulfun_fixed_SOA_and_silences](figure/simulfun_fixed_SOA_and_silences5.png) ![plot of chunk simulfun_fixed_SOA_and_silences](figure/simulfun_fixed_SOA_and_silences6.png) 




Designs with varying ISI (a la optseq)
--------------------------------------

Here, the same total amount of silence is inserted as NULL events of varying length between the 'real' trials, introducing jitter.


```r
timing <- generate_paradigm_varying_SOA_null_events(ncond, trialpercond, stimduration, 
    SOA, totalduration)
plot_paradigm(timing)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-61.png) 

```r
totalduration <- max(timing$onsets) + SOA
X <- create_design_matrix(timing, totalduration)
plot_design_matrix(X)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-62.png) 


histogram of SOAs:


```r
hist(diff(timing$onsets))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 



```r
o4 <- simulations(nsim, beta = 1:ncond, listofcontrasts = listcon, normalnoise, 
    generate_paradigm_varying_SOA_null_events, ncond, trialpercond, stimduration, 
    SOA, totalduration)

report(o4)
```

![plot of chunk simul_varying_SOA_with_NULL_events](figure/simul_varying_SOA_with_NULL_events1.png) ![plot of chunk simul_varying_SOA_with_NULL_events](figure/simul_varying_SOA_with_NULL_events2.png) ![plot of chunk simul_varying_SOA_with_NULL_events](figure/simul_varying_SOA_with_NULL_events3.png) ![plot of chunk simul_varying_SOA_with_NULL_events](figure/simul_varying_SOA_with_NULL_events4.png) 



Importing designs
-----------------

### Fixed SOA with empty trials


```r
timing <- read.table("fMRI_Order_10.csv", sep = ",", col.names = c("onsets", 
    "conditions", "durations"))
timing <- subset(timing, conditions <= 8)
plot_paradigm(timing)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-81.png) 

```r
(totalduration <- max(timing$onsets) + 10)
```

```
## [1] 770
```

```r
table(timing$conditions)
```

```
## 
## 1 2 3 4 5 6 7 8 
## 8 8 8 8 8 8 8 8
```

```r

X <- create_design_matrix(timing, totalduration)
plot_design_matrix(X)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-82.png) 




```r
betas <- c(1:4, 1:4)
singlesimul(X, betas, list(lin = betas - mean(betas)), normalnoise)
```

```
## $estimates
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 1.3336 2.5316 3.3010 3.8227 0.7448 2.0167 2.9086 3.7563 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 0.2780 0.2780 0.2797 0.2787 0.2779 0.2819 0.2793 0.2797 
## 
## $efficiencies
## [1] 40.36
```



```r
for (csvfile in Sys.glob("fMRI_Order*.csv")) {
    timing <- read.table(csvfile, sep = ",", col.names = c("onsets", "conditions", 
        "durations"))
    print("Totalduration")
    print(totalduration <- max(timing$onsets) + 10)
    print(table(timing$conditions))
    timing <- subset(timing, conditions <= 8)
    X <- create_design_matrix(timing, totalduration)
    betas <- c(1:4, 1:4)
    print(singlesimul(X, betas, list(lin = betas - mean(betas)), normalnoise))
}
```

```
## [1] "Totalduration"
## [1] 770
## 
##  1  2  3  4  5  6  7  8  9 10 11 12 
##  8  8  8  8  8  8  8  8  8  8  8  8 
## $estimates
##   XX1   XX2   XX3   XX4   XX5   XX6   XX7   XX8 
## 1.037 1.988 3.043 3.975 1.331 2.169 3.117 4.147 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 0.2736 0.2736 0.2752 0.2742 0.2735 0.2773 0.2748 0.2752 
## 
## $efficiencies
## [1] 40.36
## 
## [1] "Totalduration"
## [1] 770
## 
##  1  2  3  4  5  6  7  8  9 10 11 12 
##  8  8  8  8  8  8  8  8  8  8  8  8 
## $estimates
##   XX1   XX2   XX3   XX4   XX5   XX6   XX7   XX8 
## 1.062 1.597 2.977 3.953 1.405 1.735 3.117 4.486 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 0.2869 0.2876 0.2854 0.2834 0.2822 0.2845 0.2842 0.2835 
## 
## $efficiencies
## [1] 39.88
## 
## [1] "Totalduration"
## [1] 770
## 
##  1  2  3  4  5  6  7  8  9 10 11 12 
##  8  8  8  8  8  8  8  8  8  8  8  8 
## $estimates
##   XX1   XX2   XX3   XX4   XX5   XX6   XX7   XX8 
## 1.226 1.937 2.759 4.254 1.031 2.206 3.350 3.833 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 0.2847 0.2753 0.2807 0.2756 0.2762 0.2760 0.2782 0.2752 
## 
## $efficiencies
## [1] 40.01
## 
## [1] "Totalduration"
## [1] 770
## 
##  1  2  3  4  5  6  7  8  9 10 11 12 
##  8  8  8  8  8  8  8  8  8  8  8  8 
## $estimates
##   XX1   XX2   XX3   XX4   XX5   XX6   XX7   XX8 
## 1.394 2.370 3.239 3.962 1.447 1.990 3.142 4.466 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 0.2788 0.2775 0.2764 0.2788 0.2783 0.2790 0.2739 0.2760 
## 
## $efficiencies
## [1] 39.72
## 
## [1] "Totalduration"
## [1] 770
## 
##  1  2  3  4  5  6  7  8  9 10 11 12 
##  8  8  8  8  8  8  8  8  8  8  8  8 
## $estimates
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 0.7912 1.9848 2.9499 3.9549 0.8743 1.8224 2.7683 4.1424 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 0.2784 0.2841 0.2823 0.2795 0.2833 0.2814 0.2794 0.2794 
## 
## $efficiencies
## [1] 40.17
## 
## [1] "Totalduration"
## [1] 770
## 
##  1  2  3  4  5  6  7  8  9 10 11 12 
##  8  8  8  8  8  8  8  8  8  8  8  8 
## $estimates
##   XX1   XX2   XX3   XX4   XX5   XX6   XX7   XX8 
## 1.343 1.496 3.144 3.863 1.222 1.870 3.098 4.088 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 0.2779 0.2793 0.2792 0.2753 0.2764 0.2822 0.2774 0.2793 
## 
## $efficiencies
## [1] 39.55
## 
## [1] "Totalduration"
## [1] 770
## 
##  1  2  3  4  5  6  7  8  9 10 11 12 
##  8  8  8  8  8  8  8  8  8  8  8  8 
## $estimates
##   XX1   XX2   XX3   XX4   XX5   XX6   XX7   XX8 
## 1.204 2.137 3.177 3.946 1.179 2.275 3.157 4.017 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 0.2822 0.2794 0.2818 0.2810 0.2815 0.2820 0.2823 0.2799 
## 
## $efficiencies
## [1] 39.77
## 
## [1] "Totalduration"
## [1] 770
## 
##  1  2  3  4  5  6  7  8  9 10 11 12 
##  8  8  8  8  8  8  8  8  8  8  8  8 
## $estimates
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 1.5333 2.2805 3.0908 3.8747 0.9203 1.7688 2.9225 3.6734 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 0.2724 0.2741 0.2715 0.2774 0.2804 0.2737 0.2784 0.2760 
## 
## $efficiencies
## [1] 38.74
## 
## [1] "Totalduration"
## [1] 770
## 
##  1  2  3  4  5  6  7  8  9 10 11 12 
##  8  8  8  8  8  8  8  8  8  8  8  8 
## $estimates
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 0.8325 2.1775 2.7973 4.3606 0.4277 2.3788 2.8196 4.2151 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 0.2791 0.2809 0.2768 0.2800 0.2817 0.2807 0.2793 0.2837 
## 
## $efficiencies
## [1] 39.61
## 
## [1] "Totalduration"
## [1] 770
## 
##  1  2  3  4  5  6  7  8  9 10 11 12 
##  8  8  8  8  8  8  8  8  8  8  8  8 
## $estimates
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 0.7391 1.9882 2.8476 3.8381 1.0623 1.8728 3.0954 3.8900 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8 
## 0.2838 0.2830 0.2820 0.2867 0.2808 0.2829 0.2787 0.2789 
## 
## $efficiencies
## [1] 38.32
```

  
### sequences created by optseq  


```r
timing <- read.table("optseq/simcomp-001.par", col.names = c("onsets", "conditions", 
    "durations", "weight", "condn"))
totalduration <- max(timing$onsets) + 10
print("Totalduration")
```

```
## [1] "Totalduration"
```

```r
print(totalduration <- max(timing$onsets) + 10)
```

```
## [1] 775
```

```r
timing <- subset(timing, conditions != 0)
print(table(timing$conditions))
```

```
## 
##  1  2  3  4  5  6  7  8  9 10 
##  8  8  8  8  8  8  8  8  8  8
```

```r
hist(diff(timing$onsets))
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-111.png) 

```r
plot_paradigm(timing)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-112.png) 

```r

X <- create_design_matrix(timing, totalduration)
plot_design_matrix(X)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-113.png) 



```r
for (csvfile in Sys.glob("optseq/simcomp*.par")) {
    timing <- read.table(csvfile, col.names = c("onsets", "conditions", "durations", 
        "weight", "condn"))
    totalduration <- max(timing$onsets) + 10
    print("Totalduration")
    print(totalduration <- max(timing$onsets) + 10)
    timing <- subset(timing, conditions != 0)
    print(table(timing$conditions))
    hist(diff(timing$onsets))
    X <- create_design_matrix(timing, totalduration)
    betas <- c(1:4, 1:4, 0, 0)
    print(singlesimul(X, betas, list(lin = betas - mean(betas)), normalnoise))
}
```

```
## [1] "Totalduration"
## [1] 775
## 
##  1  2  3  4  5  6  7  8  9 10 
##  8  8  8  8  8  8  8  8  8  8
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-121.png) 

```
## $estimates
##     XX1     XX2     XX3     XX4     XX5     XX6     XX7     XX8     XX9 
##  0.8485  1.9943  3.1519  3.9296  0.9166  2.1482  2.9961  4.4572  0.1053 
##    XX10 
## -0.5974 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8    XX9   XX10 
## 0.2931 0.2946 0.2811 0.2836 0.2902 0.2980 0.2902 0.2891 0.2925 0.2861 
## 
## $efficiencies
## [1] 20
## 
## [1] "Totalduration"
## [1] 772
## 
##  1  2  3  4  5  6  7  8  9 10 
##  8  8  8  8  8  8  8  8  8  8
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-122.png) 

```
## $estimates
##      XX1      XX2      XX3      XX4      XX5      XX6      XX7      XX8 
##  0.84869  1.90286  2.99997  4.05581  0.88969  2.26410  3.09615  3.70067 
##      XX9     XX10 
## -0.15168  0.01529 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8    XX9   XX10 
## 0.2836 0.3033 0.2949 0.2890 0.2816 0.2888 0.2913 0.2948 0.3024 0.2905 
## 
## $efficiencies
## [1] 20.59
## 
## [1] "Totalduration"
## [1] 767.5
## 
##  1  2  3  4  5  6  7  8  9 10 
##  8  8  8  8  8  8  8  8  8  8
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-123.png) 

```
## $estimates
##     XX1     XX2     XX3     XX4     XX5     XX6     XX7     XX8     XX9 
##  0.5226  1.6149  3.0649  4.0112  0.7725  1.0741  2.6850  3.8887 -0.1855 
##    XX10 
## -0.4256 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8    XX9   XX10 
## 0.2932 0.2969 0.3023 0.3030 0.2949 0.3052 0.2881 0.3011 0.2937 0.2955 
## 
## $efficiencies
## [1] 18.53
## 
## [1] "Totalduration"
## [1] 777
## 
##  1  2  3  4  5  6  7  8  9 10 
##  8  8  8  8  8  8  8  8  8  8
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-124.png) 

```
## $estimates
##     XX1     XX2     XX3     XX4     XX5     XX6     XX7     XX8     XX9 
##  0.6530  1.5946  2.8946  3.7627  0.8726  1.6388  2.9115  4.6425  0.3725 
##    XX10 
## -0.1600 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8    XX9   XX10 
## 0.2881 0.2940 0.2869 0.2892 0.2734 0.2920 0.2953 0.2795 0.2668 0.2842 
## 
## $efficiencies
## [1] 20.39
## 
## [1] "Totalduration"
## [1] 777.5
## 
##  1  2  3  4  5  6  7  8  9 10 
##  8  8  8  8  8  8  8  8  8  8
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-125.png) 

```
## $estimates
##     XX1     XX2     XX3     XX4     XX5     XX6     XX7     XX8     XX9 
##  0.7409  1.7667  2.9631  3.3124  0.4462  1.5516  2.6980  3.6188  0.1013 
##    XX10 
## -0.3906 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8    XX9   XX10 
## 0.2934 0.2916 0.2807 0.2895 0.2966 0.2926 0.2982 0.2994 0.2888 0.2936 
## 
## $efficiencies
## [1] 20.49
## 
## [1] "Totalduration"
## [1] 775.5
## 
##  1  2  3  4  5  6  7  8  9 10 
##  8  8  8  8  8  8  8  8  8  8
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-126.png) 

```
## $estimates
##     XX1     XX2     XX3     XX4     XX5     XX6     XX7     XX8     XX9 
##  1.0470  2.1372  3.2034  4.3017  0.9348  2.1841  2.8639  4.2530  0.2046 
##    XX10 
## -0.4097 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8    XX9   XX10 
## 0.2867 0.2829 0.2811 0.2799 0.2861 0.2876 0.2903 0.2941 0.2872 0.2768 
## 
## $efficiencies
## [1] 19.68
## 
## [1] "Totalduration"
## [1] 776
## 
##  1  2  3  4  5  6  7  8  9 10 
##  8  8  8  8  8  8  8  8  8  8
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-127.png) 

```
## $estimates
##      XX1      XX2      XX3      XX4      XX5      XX6      XX7      XX8 
##  1.04559  2.17387  3.43201  3.60060  1.16163  2.35400  3.22337  4.28923 
##      XX9     XX10 
## -0.04203  0.29526 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8    XX9   XX10 
## 0.2866 0.2848 0.2938 0.2938 0.2895 0.2859 0.2884 0.2964 0.2822 0.2906 
## 
## $efficiencies
## [1] 19.9
## 
## [1] "Totalduration"
## [1] 763
## 
##  1  2  3  4  5  6  7  8  9 10 
##  8  8  8  8  8  8  8  8  8  8
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-128.png) 

```
## $estimates
##      XX1      XX2      XX3      XX4      XX5      XX6      XX7      XX8 
##  0.54680  1.98257  3.12512  3.98747  0.80656  1.70044  2.85582  3.93365 
##      XX9     XX10 
## -0.01323  0.27119 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8    XX9   XX10 
## 0.3063 0.3062 0.3077 0.3068 0.3069 0.2999 0.2998 0.3102 0.3000 0.2973 
## 
## $efficiencies
## [1] 19.74
## 
## [1] "Totalduration"
## [1] 776.5
## 
##  1  2  3  4  5  6  7  8  9 10 
##  8  8  8  8  8  8  8  8  8  8
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-129.png) 

```
## $estimates
##     XX1     XX2     XX3     XX4     XX5     XX6     XX7     XX8     XX9 
##  1.1770  2.3127  3.0461  3.7574  1.0375  1.8516  2.8205  3.7224 -0.3937 
##    XX10 
## -0.1947 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8    XX9   XX10 
## 0.2873 0.2850 0.3084 0.2898 0.2841 0.2974 0.2914 0.3000 0.2870 0.2852 
## 
## $efficiencies
## [1] 19.23
## 
## [1] "Totalduration"
## [1] 776.5
## 
##  1  2  3  4  5  6  7  8  9 10 
##  8  8  8  8  8  8  8  8  8  8
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-1210.png) 

```
## $estimates
##      XX1      XX2      XX3      XX4      XX5      XX6      XX7      XX8 
##  1.05743  1.86711  3.00429  4.05282  0.91307  1.85680  2.87120  4.13944 
##      XX9     XX10 
## -0.06797  0.23410 
## 
## $estimatese
##    XX1    XX2    XX3    XX4    XX5    XX6    XX7    XX8    XX9   XX10 
## 0.2825 0.2954 0.2958 0.2966 0.2947 0.2972 0.3001 0.2956 0.2890 0.2949 
## 
## $efficiencies
## [1] 19.25
```

  

