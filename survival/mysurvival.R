mysurvival <- function(time, event, var, data, pvalue = 0.05, ...){
    # loading packages
    require(survival)
    
    fml <- paste("Surv(", time, ", ", event, ") ~", var) %>% as.formula
    # print(fml)
    fit <- coxph(formula = fml, data = data)
    # print(fit)
    
    diff <- survdiff(formula = fml, data = data)
    p <- diff$pvalue
    
    if(p < pvalue){
        return(list(pvalue = p, gene = var))
    }
}
