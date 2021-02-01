#' This function computes the coeffcients of the index function 
#' @param y a vector of response. 
#' @param xfd a fda object contains the functional predictor 
#' @param spline_basis the B-spline basis for the index function 
#' @param lambda a positive number, the sparsity parameter 
#' @export cFuSIM_index
#' @import fda
#' @import dplyr
#' @import slos
#' @import NonpModelCheck
#' @import DiceKriging
#' @examples
#' \dontrun{
#'library(fda)
#'library(cFuSIM)
#'data(bike_cFuSIM)
#'timepts = bike$timepts
#'norder=4 ## cubic B-spline
#'nbasis=norder+length(timepts)-2; 
#'spline_basis=create.bspline.basis(rangeval=c(1,24),nbasis#'norder,timepts)
#'wull = bike$temp
#'xfds=  Data2fd(y=wull%>%t, argvals=bike$timepts)
#'y = bike$y
#'train_sample = 1:length(y)
#'y = y[train_sample]
#'xfd = xfds[train_sample]
#'res_c = cFuSIM_index(y, xfd, spline_basis)
#'beta_fd = fd(res_c$coefBeta, res_c$basisBeta)
#'plot(beta_fd,ylab="index function", xlab='time')
#'fdagg(beta_fd)
#'score_fit = (res_c$score_fit)
#' pred_y = localpoly.reg(score_fit, y, degree.pol = 1, kernel.type = "gaussian",bandwidth = "CV",deriv=0,points#'score_fit)
#'plot(x=score_fit, y=y)
#'lines(pred_y$predicted[order(score_fit)],x=score_fit[order#'score_fit)],col=4)
#' }

# wtemp = function(xbasis0) {
#         xbasis0 %>% knots(, interior = FALSE) %>% unique %>% 
#             data.frame(knot = .) %>% mutate(knotlead = lead(knot)) %>% 
#             dplyr::filter(!is.na(knotlead)) %>% rowwise() %>% 
#             do(temp = eval.penalty(xbasis0, int2Lfd(0), rng = c(.$knot, 
#                 .$knotlead)))
#     }


# pmatf_locploy = function(betac ,ytrain=ytrain,xtrain=xtrain, spline_basis){
# B = inprod(spline_basis,spline_basis,0,0)
# pc_fit = fd(betac,spline_basis)
# alpha = inprod(xtrain,spline_basis)
# score_fit = inprod(xtrain,pc_fit)%>%as.numeric
# loc_fit = localpoly.reg(score_fit, ytrain, degree.pol = 1, kernel.type = "gaussian",bandwidth = "CV",deriv=0)
# fity = loc_fit$predicted
# loc_fit = localpoly.reg(score_fit, ytrain, degree.pol = 1, kernel.type = "gaussian",bandwidth = "CV",deriv=1)
# derfity = loc_fit$predicted
# cc = ytrain - fity-derfity
# Ba = B%*%t(alpha)
# mattemp = matrix(rep(cc, each=spline_basis$nbasis), nrow=spline_basis$nbasis)*Ba
# Pmat = mattemp%*%t(mattemp)
# Pmat 
# }

# W0=  function(beta_j,lambda,bspi,W_m, normthres=10^-2){
# W = inprod(bspi,bspi,1,1)
# R = inprod(bspi,bspi)
# r = inprod(bspi,bspi,2,2)

# range_knots = knots(bspi,interior=FALSE)%>%range

 

#     wi = function(i){
#     zero = NULL       
#     beta_m_norm =  t(as.matrix(beta_j))%*%W_m[[i]]%*%as.matrix(beta_j)%>%as.numeric%>%sqrt
    
#     if(beta_m_norm < normthres)   {

#         zero = c(i:(i+3))
#         res = matrix(0, nrow=bspi$nbasis, ncol=bspi$nbasis)
#     } else{
#     temp1 = SCAD.derivative(sqrt((length(W_m))/((bspi$rangeval[2]-bspi$rangeval[1]))) * beta_m_norm, lambda = lambda)
#     temp2 = beta_m_norm * sqrt(((bspi$rangeval[2]-bspi$rangeval[1]))/(length(W_m)))
#     res = temp1/temp2*W_m[[i]]

#     }
#     return(list(res=res, zero=zero))
#     }
            
# res_list = lapply(1:length(W_m), wi)
# res = lapply(res_list, function(x) x[[1]])
# zero = do.call(c,lapply(res_list, function(x) x[[2]]))%>%unique

# W_j0= Reduce('+', res)*0.5
# return(list( W = W_j0, zero = zero))
# }


cFuSIM_index = function(y, xfd, spline_basis, threshold=1e-5, maxit=150,lambda=1e4, verbose=TRUE){
  lambda = 10^log10lambda
  # set.seed(1001011)
  (betac= rep(1,spline_basis$nbasis))
  B = inprod(spline_basis,spline_basis,0,0)
  betac = betac/sqrt(as.numeric(t(betac)%*%B%*%betac))
  
  GAMMA = inprod(spline_basis,spline_basis,2,2)
  W_m = wtemp(spline_basis )$temp
  thes = 1
  i = 1
  thresb= 3
  while(thes>threshold){
    i=i+1
    if(i>maxit) break
    betac0 = betac
    pmatf = pmatf_locploy(betac,y,xfd,spline_basis,bandwidth="CV")
    yweighted = as.numeric(pmatf$yweighted)
    nrowx = nrow(coef(xfd))
    xfds_weighted = fd(rep(pmatf$weights, each=nrowx)%>%matrix(nrow=nrowx)*coef(xfd), xfd$basis)
    beta_basis = list()
    beta_basis[[1]] = spline_basis
    fit <- try(slos(xfds_weighted,yweighted%>%as.numeric,D=NULL,lambda=lambda,gamma=gamma,intercept=F,beta.basis=beta_basis,cutoff=1e-4,
                          max.iter=1000,tol=1e-12,a=3.7,domain=NULL, tuning = "BIC"), silent=TRUE)
    if(class(fit)=="try-error") 
    {
      betac = rep(0, spline_basis$nbasis) 
      thes = 0
    } else {
      bandwidth = pmatf$bandwidth
      betac = fit$beta[[1]]%>%coef%>%as.numeric
      pc_fit = fit$beta[[1]]
      betac = betac/sqrt(as.numeric(t(betac)%*%B%*%betac))
      pc_fit = fd(betac,spline_basis)
      score_fit = inprod(xfd,pc_fit)%>%as.numeric
      if(as.numeric(inprod(pc_fit))<0) {
        betac= -betac
        pc_fit = fd(betac,spline_basis)
      }
      
      thes  = abs(betac-betac0)%>%max
      thresb =thes
     if(verbose){
      plot(pc_fit) 
      print(thresb)
      }
    }
    
  }
  
  converaged=(thes< threshold)
  if(!all(betac==0)){
    betac = betac/sqrt(as.numeric(t(betac)%*%B%*%betac))
  }
  pc_fit = fd(betac,spline_basis)
  score_fit = inprod(xfd,pc_fit)%>%as.numeric
  list(coefBeta=betac, basisBeta = spline_basis, score_fit = score_fit, Converaged=converaged, threshold=threshold, maxit=maxit, bandwidth = bandwidth)
}

