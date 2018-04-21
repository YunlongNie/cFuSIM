 # install.packages('fda')
 # install.packages('dplyr')
 # install.packages('DiceKriging')
 # install.packages('NonpModelCheck')
 # install.packages('devtools')
 # library(devtools)
 # install_github('YunlongNie/cFuSIM')

 library(fda)
 library(cFuSIM)
 data(bike_cFuSIM)
 timepts = bike$timepts
 norder=4 ## cubic B-spline
 nbasis=norder+length(timepts)-2; 
 spline_basis=create.bspline.basis(rangeval=c(1,24),nbasis,norder,timepts)
 wull = bike$temp
 xfds=  Data2fd(y=wull%>%t, argvals=bike$timepts)
 y = bike$y
 train_sample = 1:length(y)
 y = y[train_sample]
 xfd = xfds[train_sample]
 res_c = cFuSIM(y, xfd, spline_basis)
 beta_fd = fd(res_c$coefBeta, res_c$basisBeta)
 plot(beta_fd)