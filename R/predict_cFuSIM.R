#' This function uses predict on the test set given the training set and parameters
#' @export predict_cFuSIM
#' @import fda
#' @import dplyr
#' @import ggplot2
#' @import NonpModelCheck
#' @examples
#' \dontrun{
#' }

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
