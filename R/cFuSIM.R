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
#' data(binary_dat);xmat = binary_dat$x;y=binary_dat$y
#' res = sfpcs_binary(xmat,y,npc_select=2,theta=1,xmat_new = xmat)
#' res$fitted
#' mean(y==res$fitted) # fitted accuracy
#' res$fitted # prediction on new data 
#' res$sfpcs # sFPCs 
#' res$beta_fd # coefficient function
#' }
#' 
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

