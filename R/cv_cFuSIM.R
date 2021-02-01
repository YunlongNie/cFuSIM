
cv_cFuSIM = function(y,xfds, log10lambda=c(),gamma=c(),bandwidth=c(),spline_basis=xfds$basis, k=3, ...){
  parameters = expand.grid( log10lambda =  log10lambda, gamma=gamma, bandwidth=bandwidth)
  cv_res = c()
  cv_folds = createFolds(length(y), nfolds = k)
  for (j in 1:nrow(parameters)){
    log10lambda = parameters$log10lambda[j]
    gamma = parameters$gamma[j]
    bandwidth = parameters$bandwidth[j]
    cv_errors = c()
    for (i in 1:max(cv_folds)){
      test = which(cv_folds==i)
      train = which(cv_folds!=i)
      ytrain=y[train]
      xfdtrain = xfds[train]
      ytest=y[test]
      xfdtest = xfds[test]
      cv_errors = c(cv_errors,predict_cFuSIM (ytrain, xfdtrain, ytest, xfdtest,spline_basis,
                                              log10lambda=log10lambda,
                                              gamma =gamma,
                                              bandwidth = bandwidth,
                                              ...))
    }
    cv_res = c(cv_res,mean(cv_errors))
  }
  parameters$cv_errors=cv_res
  parameters = parameters%>%arrange(cv_errors)
  return(parameters)
}
