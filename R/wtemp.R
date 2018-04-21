#' @import fda
#' @import dplyr
wtemp = function(xbasis0) {
        xbasis0 %>% knots(, interior = FALSE) %>% unique %>% 
            data.frame(knot = .) %>% mutate(knotlead = lead(knot)) %>% 
            dplyr::filter(!is.na(knotlead)) %>% rowwise() %>% 
            do(temp = eval.penalty(xbasis0, int2Lfd(0), rng = c(.$knot, 
                .$knotlead)))
    }
