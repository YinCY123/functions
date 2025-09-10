UniCox <- function(time, event, var, data){
    require(survival)
    require(glmnet)
    y = Surv(data[, time, drop = T], data[, event, drop = T] == 1) 
    x <- data[, var, drop = F] %>% as.matrix
    fit = coxph(y ~ x, data = data)
    fit_summ = summary(fit)
    
    df <- data.frame(
        Variable = var, 
        mean = round(fit_summ$coefficients[, 2, drop = T], 3),
        pvalue = as.character(round(fit_summ$coefficients[, 5, drop = T], 3)), 
        lower = round(fit_summ$conf.int[, 3, drop = T], 3), 
        upper = round(fit_summ$conf.int[, 4, drop = T], 3)
    )
    return(df)
}