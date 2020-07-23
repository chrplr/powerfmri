#! /usr/bin/env Rscript

# the first argument is input file, in csv format, which must contain 3 named columns: 'onsets','durations','conditions'
args <- commandArgs(trailingOnly = TRUE)
timing <- read.csv(args[1])
totalduration <- max(timing$onsets + 15)
ncond <- length(unique(timing$conditions))

################################
# Misc support functions

source('hrf_utils.R')

generate_hrf_design_mat <- function(timings, totalduration)
{
  with(timings,
{ condnames <- sort(unique(conditions))
  ncond <- length(condnames)
  
  regressor <- vector("list", ncond)
  for (cond in 1:ncond)
  {
    select <- conditions==condnames[cond]
    regressor[[cond]] <- scale(hrf(list(onsets=onsets[select],
                                        durations=durations[select]),
                                   totalduration))
  }
  
  names(regressor) <- condnames
  browser()
  as.matrix(data.frame(regressor,const=1))
}
  )
}

plot_matrix <- function(X)
  # plot the columns of the matrix X 
{ matplot(1:nrow(X), X, type='l', col=1:nrow(X)) }

trac <- function(M) { sum(diag(M)) }
efficiency <- function(mat, contrast)
{
  1/ trac(contrast %*% (chol2inv(crossprod(mat)) %*% contrast))
}


#####################################################
# Main
pdf(paste(args[1],'plots.pdf', sep='_'))

X <- generate_hrf_design_mat(timing, totalduration) 

plot_matrix(X)

#Compute the efficiency (detection power) 1/tr(cT(xTx)-1.c)


linearcontrast <- c((1:ncond)-mean(1:ncond),0)
efficiency(X, linearcontrast)


betas <- c(1:ncond, 0)
y0 <- X %*% betas 
nsim <- 100
estimates <- matrix(nrow=nsim, ncol=ncond+1)
for (sim in 1:nsim)
{
  #y <- y0 + 20*pink.noise(length(y0))
  y <- y0 + rnorm(length(y0), sd=10)
  estimates[sim,] <- coef(lm(y ~ X - 1))
}
plot(y, col='black')
lines(y0, col='red' )
es <- data.frame(estimates=apply(estimates,2,mean),sd=apply(estimates,2,sd))
es
plot(1:nrow(es), es[,1], pch=16, col='black',ylim=c(0,6))
segments(1:nrow(es), es[,1]-es[,2], 1:nrow(es), es[,1]+es[,2])

graphics.off()