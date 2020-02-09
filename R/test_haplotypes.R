#' block, complete and variant models using R-base lm.
#' Using lm
#' Assuming no covariates 
#' or Y are residues over covariates for all model
#' Missing values allowed in Y
#' No missing values in hapDF
#' @param Y  phenotype or residues matrix (missing values allowed)
#' @param hapDF  haplotype counts matrix
#' @return \item{results_L}{ The resulting test as summary.lm for all haplotypes for 3 test }
#' @title R-lm based test for haplotypes including block, complete and variant models
#' @export all_lm
all_test_residues = function(hapDF, Y)
{
  {
  
    X = hapDF
    org_rownamesY = rownames(Y)
    org_rownamesX = rownames(X)
    col_orgY = colnames(Y)
    col_orgX = colnames(X)
    allelcount = unlist(apply(X, 2, sum))
    varL = colnames(X)
    if (length(varL) > 1) {
      varL = colnames(X)[setdiff(1:ncol(X), which.max(allelcount))]
    }
    form_Y = paste(col_orgY, collapse = ',')
    form_X = paste(varL, collapse = '+')
    full_Lm = lapply( sprintf("cbind( %s ) ~ %s ", col_orgY, form_X),
                      lm,
                      data = data.frame(cbind(Y, X)))
    sum_Lm = lapply(full_Lm,summary)
    names(sum_Lm) = col_orgY # rownames(full_Lm$coefficients)
    BLOCK_TEST = unlist(lapply(sum_Lm,
                               function(y) {
                                 x = y$fstatistic
                                 pf(x[1], x[2], x[3], lower.tail = FALSE)
                               }))
    COMPLETE_TEST = lapply(sum_Lm,
                           function(y,varL) {
                             x = coefficients(y)
                             x[2:nrow(x), ncol(x)]
                           },varL)
    varLM = lapply(col_orgX, function(x)  
      lapply(col_orgY, function(y)
        lm(
          formula = sprintf("cbind( %s ) ~ %s ", y, x),
          data = data.frame(cbind(Y, X))
        )))
    # names(varLM) = col_orgX
    sum_VLM = lapply(varLM, function(vlm_L) lapply(vlm_L,summary))
    
    
    pv_Var = lapply(sum_VLM,
                    function(z) {
                      names(z) = col_orgY[1:ncol(Y)]
                      lapply(z,
                             function(x) {
                               y = coefficients(x)
                               y[nrow(y), ncol(y)]
                             })
                    })
    VARIANT_TEST = setNames(
      object = (group_by(do.call(
        rbind,
        lapply(pv_Var, function(x)
          data.frame(
            pheno = names(x),
            pval =
              unlist(x[as.character(newColNames)])
          ))
      ),
      pheno) %>% summarise(pval = list(
        setNames(pval, col_orgX)
      )))$pval,
      as.character(newColNames)
    )
  }
  
  results_L = list(block = BLOCK_TEST,
                   complete = COMPLETE_TEST,
                   variant = VARIANT_TEST)
  return(results_L)
  
}
