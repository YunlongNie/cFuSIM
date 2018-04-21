#' @import fda
#' @import dplyr
#' @import NonpModelCheck

pmatf_locploy = function(betac ,ytrain=ytrain,xtrain=xtrain, spline_basis){
B = inprod(spline_basis,spline_basis,0,0)
pc_fit = fd(betac,spline_basis)
alpha = inprod(xtrain,spline_basis)
score_fit = inprod(xtrain,pc_fit)%>%as.numeric
loc_fit = localpoly.reg(score_fit, ytrain, degree.pol = 1, kernel.type = "gaussian",bandwidth = "CV",deriv=0)
fity = loc_fit$predicted
loc_fit = localpoly.reg(score_fit, ytrain, degree.pol = 1, kernel.type = "gaussian",bandwidth = "CV",deriv=1)
derfity = loc_fit$predicted
cc = ytrain - fity-derfity
Ba = B%*%t(alpha)
mattemp = matrix(rep(cc, each=spline_basis$nbasis), nrow=spline_basis$nbasis)*Ba
Pmat = mattemp%*%t(mattemp)
Pmat 
}