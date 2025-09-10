riskdf <- function(fit, s){
  require(magrittr)
  require(tibble)
  coef(fit, s = s) %>% 
    as.matrix %>% 
    as.data.frame %>% 
    tibble::rownames_to_column() %>%
    magrittr::set_colnames(value = c("symbol", "coef")) %>% 
    dplyr::filter(coef != 0)
}