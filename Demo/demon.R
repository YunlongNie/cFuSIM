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
 res_c = cFuSIM_index(y, xfd, spline_basis)
 beta_fd = fd(res_c$coefBeta, res_c$basisBeta)
 
 # plot the index function 
 fdagg(beta_fd,ylab="index function", xlab='time')
 score_fit = (res_c$score_fit)
 pred_y = localpoly.reg(score_fit, y, degree.pol = 1, kernel.type = "gaussian",bandwidth = "CV",deriv=0,points=score_fit)
 # plot the fitted integral vs the response 
 plot(x=score_fit, y=y)
 lines(pred_y$predicted[order(score_fit)],x=score_fit[order(score_fit)],col=4)

## compute BIC 
fity = pred_y$predicted
mse  = mean((y - fity)^2)
p = sum(res_c$coefBeta!=0)
bic = length(y)*log(mse)+log(length(y))*(p+1)

