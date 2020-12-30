#' @import fda
#' @import dplyr
#' @import DiceKriging

W0=  function(beta_j,lambda,bspi,W_m, normthres=10^-2){
W = inprod(bspi,bspi,1,1)
R = inprod(bspi,bspi)
r = inprod(bspi,bspi,2,2)

range_knots = fda:::knots.basisfd(bspi,interior=FALSE)%>%range

 

    wi = function(i){
    zero = NULL       
    beta_m_norm =  t(as.matrix(beta_j))%*%W_m[[i]]%*%as.matrix(beta_j)%>%as.numeric%>%sqrt
    
    if(beta_m_norm < normthres)   {

        zero = c(i:(i+3))
        res = matrix(0, nrow=bspi$nbasis, ncol=bspi$nbasis)
    } else{
    temp1 = SCAD.derivative(sqrt((length(W_m))/((bspi$rangeval[2]-bspi$rangeval[1]))) * beta_m_norm, lambda = lambda)
    temp2 = beta_m_norm * sqrt(((bspi$rangeval[2]-bspi$rangeval[1]))/(length(W_m)))
    res = temp1/temp2*W_m[[i]]

    }
    return(list(res=res, zero=zero))
    }
            
res_list = lapply(1:length(W_m), wi)
res = lapply(res_list, function(x) x[[1]])
zero = do.call(c,lapply(res_list, function(x) x[[2]]))%>%unique

W_j0= Reduce('+', res)*0.5
return(list( W = W_j0, zero = zero))
}
