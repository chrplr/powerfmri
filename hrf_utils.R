# A few useful functions
# Time-stamp: <2012-12-28 14:13 christophe@pallier.org>

ihrf <- scan('hrf.dat') # hrf sampled at 10Hz

dirac.comb <- function(onsets, durations, total.duration, sampling.rate=10)
  # given vectors 'onsets' and 'durations' in secs, generate an indicator variable
  # sampled at 10Hz
{
  stopifnot((length(onsets)==length(durations)) | (length(durations)==1))
  if (length(durations)==1) durations<-rep(durations, length(onsets))
  onsets <- ifelse(onsets<0, 0, onsets)
  size <- (total.duration)*sampling.rate
  u <- numeric(size)
  for (i in 1:length(onsets)) 
  {
    i1 <- round(sampling.rate*onsets[i])
    i2 <- round(sampling.rate*(onsets[i]+durations[i]))
    i1 <- min(i1, size)
    i2 <- min(i2, size)
    u[i1:i2] <- 1.0
  }
  u[1:size]
}

# test


hrf <- function(onsets, durations, total.duration, TR=1)
  # convolution of a timpe series by impulse hrf
{
  stopifnot(TR==1) # undersampling not yet implemented
  if (!exists("ihrf")) 
    ihrf <- scan('hrf.dat') # hrf sampled at 10Hz
  x2 <- c(rep(0, length(ihrf)), dirac.comb(onsets, durations, total.duration))
  x3 <- stats::filter(x2, ihrf, sides=1)
  x4 <- x3[-(1:length(ihrf))]
  x4 # here we should do some undersampling if TR>1
}

require(signal)
pink.noise <-  function(n)
  # function adapted from https://ccrma.stanford.edu/~jos/sasp/Example_Synthesis_1_F_Noise.html
{
  B <- c(0.049922035, -0.095993537, 0.050612699, -0.004408786)
  A <- c(1, -2.494956002,   2.017265875,  -0.522189400)
  nT60 <- 1430
  v <- rnorm(n+nT60)
  x <- filter(B,A,v)    # Apply 1/F roll-off to PSD
  x <- x[-(1:nT60+1)]    # Skip transient response
}
