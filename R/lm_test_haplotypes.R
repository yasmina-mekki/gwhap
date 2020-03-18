#' lm test haplotypes
#' Perform three tests using R-base lm: 
#' haplotype bloc model test
#' complete model haplotype test
#' variant (single) haplotype model test
#' the common haplotype is removed for two first test
#' Assuming no covariates or Y are residues over covariates for all model
#' No missing values allowed in X
#' Missing values allowed in Y
#'
#' @param X haplotype counts matrix
#' @param Y phenotype or residues matrix
#' @param kind which test to perform. Four values are possible : single, complete, bloc or all
#'
#' @return The resulting test as summary.lm for all haplotypes for 3 test
#' @export
#'
lm_test_haplotypes = function(X, Y, kind='all'){
  
  if(kind != 'single'){
    # get the columns name
    col_orgY = colnames(Y)
    col_orgX = colnames(X)
    
    # get the sum of each columns value and identify the common haplotype
    # Que signifi les nombre 1, 2, 0 pour les haplotypes (dummification) 0, 1 pour présence ou pas mais le 2 ?
    allelcount = unlist(apply(X, 2, sum))
    varL = colnames(X)
    
    # get the index of the max of allel count
    # get haplotypes columns without the common one.
    if (length(varL) > 1) {varL = colnames(X)[setdiff(1:ncol(X), which.max(allelcount))]}
    
    form_Y = paste(col_orgY, collapse = ',')
    form_X = paste(varL, collapse = '+')
    
    # linear regression
    full_Lm = lapply(sprintf("cbind( %s ) ~ %s ", col_orgY, form_X), lm, data = data.frame(cbind(Y, X)))
    sum_Lm  = lapply(full_Lm, summary)
    names(sum_Lm) = col_orgY
  }
  
  if(kind == 'bloc'){return(bloc_test(sum_Lm))}
  if(kind == 'complete'){return(complete_test(sum_Lm))}
  if(kind == 'single'){return(single_test(X, Y))}
  if(kind == 'all'){return(list(block = bloc_test(sum_Lm), complete = complete_test(sum_Lm), variant = single_test(X, Y)))}
  
}


#' haplotype bloc model test
#' association between a given bloc and a phenotype
#' the common haplotype is removed for this test
#' warning : this function should not be used explicitly
#' users must only call the lm_test_haplotypes function and specify bloc as kind parameters value
#' @param sum_Lm summary of a linear model estimated using lm_test_haplotypes function
#'
#' @return Significance of the association
#' @export
#'
bloc_test <- function(sum_Lm){
  return(unlist(lapply(sum_Lm, function(y) {x = y$fstatistic
                                            pf(x[1], x[2], x[3], lower.tail = FALSE)})))
}

#' complete model haplotype test
#' association between each haplotype of a given bloc and a phenotype
#' the common haplotype is removed for this test
#' warning : this function should not be used explicitly
#' users must only call the lm_test_haplotypes function and specify complete as kind parameters value
#' @param sum_Lm summary of a linear model estimated using lm_test_haplotypes function
#'
#' @return Significance of the each association. A two-side p-value of the t-statitic is estimated
#' @export
#'
complete_test <- function(sum_Lm){
  return(lapply(sum_Lm, function(y, varL) {x = coefficients(y)
                                           x[2:nrow(x), ncol(x)]}, varL))
}

#' variant (single) haplotype model test
#' the common haplotype is keeped for this test
#' @param X haplotype counts matrix
#' @param Y phenotype or residues matrix
#'
#' @return Significance of the each association. A two-side p-value of the t-statistic is estimated
#' @export
#'
single_test <- function(X, Y){
  
  # get the columns name
  col_orgY = colnames(Y)
  col_orgX = colnames(X)
  
  # linear regression
  varLM = lapply(col_orgX,
                 function(x) lapply(col_orgY,
                                    function(y) lm(formula = sprintf("cbind( %s ) ~ %s ", y, x), data = data.frame(cbind(Y, X)))
                                    )
                 )
  
  # get the summary of the model trained above
  sum_VLM = lapply(varLM, function(vlm_L) lapply(vlm_L, summary))
  
  # get the p values
  pv_Var = lapply(sum_VLM, function(z){names(z) = col_orgY[1:ncol(Y)]
                                       lapply(z, function(x){y = coefficients(x)
                                                             y[nrow(y), ncol(y)]})
                                       }
                   )
  
  # set the output names
  single_test = setNames(object=(group_by(do.call(rbind,
                                                  lapply(pv_Var, function(x) data.frame(pheno = names(x), pval = unlist(x[as.character(names(x))])))),
                                          pheno) %>% summarise(pval = list(setNames(pval, col_orgX))))$pval,
                         as.character(colnames(Y))
                         )
  
  return(single_test)
}