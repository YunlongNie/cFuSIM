predict_cFuSIM = function(ytrain, xfdtrain, ytest=ytrain, xfdtest=xfdtrain, ... ){
  res_cFuSIM = cFuSIM_index(ytrain,xfdtrain,...)
  if(all(res_cFuSIM$coefBeta==0)) {
    print('all coefs are zero!')
    return(mean((ytest-mean(ytrain))^2))
  } else {
    beta_fd = fd(res_cFuSIM$coefBeta, res_cFuSIM$basisBeta)
    score_fit_test = inprod(xfdtest,beta_fd)%>%as.numeric
    score_fit_train = inprod(xfdtrain,beta_fd)%>%as.numeric
    score_fit = score_fit_train
    y = ytrain
    y = y[order(score_fit)]
    score_fit = score_fit[order(score_fit)]
    # plot the index function 
    pred_y = localpoly.reg(score_fit, y, degree.pol = 1, kernel.type = "gaussian",bandwidth = res_cFuSIM$bandwidth,deriv=0)
    gbasis= create.bspline.basis(rangeval=range(c(score_fit_test,score_fit_train)),breaks =c(score_fit_test,score_fit_train)%>%unique%>%sort, norder=3)
    gfd = Data2fd(argvals=score_fit, y=pred_y$predicted, basisobj = gbasis)
    y_pred_test = as.numeric(eval.fd(gfd,score_fit_test))
    return(mean((ytest-y_pred_test )^2))
    
  }
}



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

