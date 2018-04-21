#' This function computes the coeffcients of the index function 
#' @param y a vector of response. 
#' @param xfd a fda object contains the functional predictor 
#' @param spline_basis the B-spline basis for the index function 
#' @param lambda a positive number, the sparsity parameter 
#' @export
#' @import fda
#' @import dplyr
#' @examples
#' \dontrun{
#' library(fda)
#' library(cFuSIM)
#' data(bike_cFuSIM)
#' norder=4 ## cubic B-spline
#' nbasis=norder+length(timepts)-2; 
#' spline.basis=create.bspline.basis(rangeval=c(1,24),nbasis,norder,timepts)
#' wull = bike$temp
#' xfds=  Data2fd(y=wull%>%t, argvals=bike$timepts)
#' y = bike$y
#' train_sample = 1:length(y)
#' y = y[train_sample]
#' xfd = xfds[train_sample]
#' res_c = cFuSIM(y, xfd, spline_basis)
#' beta_fd = fd(res_c$coefBeta, res_c$basisBeta)
#' plot(beta_fd)
#' }

cFuSIM = function(y, xfd, spline_basis, threshold=1e-5, maxit=150,lambda=1e4){
(betac= rep(1, spline_basis$nbasis))
thes = 1
i = 1
thresb=3
while(thes>threshold){
i=i+1
if(i>maxit) break
betac0 = betac
Pmat = pmatf_locploy(betac,y,xfd)
Wmat = W0(betac,normthres=normthres,lambda=lambda,bspi=spline_basis)
(nonzero  =setdiff(1:spline_basis$nbasis,Wmat$zero))
Pmat_non = Pmat[nonzero,nonzero]
W_non= Wmat$W[nonzero,nonzero]
betac_non= betac[nonzero] 
B_non = B[nonzero, nonzero]
betact = try(solve(Pmat_non + W_non+gamma*B_non, Pmat_non%*%betac_non), silent=TRUE)
if(class(betact)=="try-error") break
betac = rep(0, spline_basis$nbasis)
(betac[nonzero]= betact)
score_fit = inprod(fd(betac, spline_basis), xfd)
betac = 10*betac/sqrt(sum((score_fit^2)))
pc_fit = fd(betac,spline_basis)
if(as.numeric(inprod(pc_fit))<0) betac= -betac
thes  = abs(betac-betac0)%>%max
thresb =thes
}

converaged=(thes< threshold)

betac = betac/sqrt(as.numeric(t(betac)%*%B%*%betac))
pc_fit = fd(betac,spline_basis)
score_fit = inprod(xfd,pc_fit)%>%as.numeric
list(coefBeta=betac, basisBeta = spline_basis, score_fit = score_fit, Converaged=converaged, threshold=threshold, maxit=maxit)
}

