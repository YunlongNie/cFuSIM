#' @import fda
#' @import dplyr
#' @import NonpModelCheck

pmatf_locploy = function(betac ,ytrain=ytrain,xtrain=xtrain, spline_basis){
  B = inprod(spline_basis, spline_basis, 0, 0)
  pc_fit = fd(betac, spline_basis)
  Z = inprod(xtrain, spline_basis)
  score_fit = inprod(xtrain, pc_fit) %>% as.numeric
  bandwidth0 = bandwidth
  loc_fit = localpoly.reg(score_fit, ytrain, degree.pol = 1,
                          kernel.type = "gaussian", bandwidth = bandwidth0 , deriv = 0)
  select_bandwith = loc_fit$bandwidth
  fity = loc_fit$predicted
  loc_fit = localpoly.reg(score_fit, ytrain, degree.pol = 1,
                          kernel.type = "gaussian", bandwidth = loc_fit$bandwidth, deriv = 1)
  derfity = loc_fit$predicted
  yweighted = ytrain - fity + derfity*(Z%*%betac)
  list( bandwidth = select_bandwith, weights = -derfity, yweighted = yweighted, Z = Z)
}
