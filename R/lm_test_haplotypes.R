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
#' @return A list of the results of the three test as summary.lm for all haplotypes
#' @importFrom stats lm
#' @export
#'
lm_test_haplotypes = function(X, Y, kind='all'){

  if(kind != 'single'){

    colnames(Y) = c('phenotype')

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

  final_results = list()

  if(kind == 'bloc' | kind == 'all'){

    # perform the test
    bloc_test_results = bloc_test(sum_Lm)

    # set the results into the right format
    results  = data.frame(strsplit(varL[1], '_')[[1]][2], strsplit(varL[1], '_')[[1]][3], bloc_test_results[[1]], nrow(Y), length(varL))
    names(results) <- c('start', 'end', 'p_value','nb_subjects', 'nb_haplotypes')

    final_results[['bloc']] = results
  }

  if(kind == 'complete' | kind == 'all'){

    # perform the test
    complete_test_results = complete_test(sum_Lm)

    # set the results into the right format
    results  = data.frame()
    for (i in 1:length(complete_test_results$phenotype)){
      results = rbind(results, data.frame(strsplit(names(complete_test_results$phenotype[i]), '_')[[1]][2],
                                          strsplit(names(complete_test_results$phenotype[i]), '_')[[1]][3],
                                          names(complete_test_results$phenotype[i]),
                                          complete_test_results$phenotype[[i]],
                                          nrow(Y),
                                          length(complete_test_results$phenotype)))
    }
    names(results) <- c('start', 'end', 'haplotype','p_value','nb_subjects', 'nb_haplotypes')

    final_results[['complete']] = results
  }

  if(kind == 'single' | kind == 'all'){final_results[['single']] = single_test(X, Y)}

  return(final_results)
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
#' @importFrom stats lm
#' @export
#'
single_test <- function(X, Y){

  # get the columns name
  colnames(Y) = c('phenotype')
  col_orgY = colnames(Y)
  col_orgX = colnames(X)

  # linear regression for multiple phenotypes ...
  varLM = lapply(col_orgX,
                 function(x) lapply(col_orgY,
                                    function(y) lm(formula = sprintf("cbind( %s ) ~ %s ", y, x), data = data.frame(cbind(Y, X)))
                                    )
                 )

  # linear regression for one phenotype
  #varLM = lapply(col_orgX,
  #               function(x) lm(formula = sprintf("cbind( %s ) ~ %s ", col_orgY, x), data = data.frame(cbind(Y, X)))
  #              )

  # get the summary of the model trained above
  sum_VLM = lapply(varLM, function(vlm_L) lapply(vlm_L, summary))
  #sum_VLM = lapply(varLM, function(vlm_L) summary(vlm_L))

  # get the p values
  pv_Var = lapply(sum_VLM, function(z){names(z) = col_orgY[1:ncol(Y)]
                                       lapply(z, function(x){y = coefficients(x)
                                                             y[nrow(y), ncol(y)]})
                                       }
                   )

  # set the output into the right format
  results  = data.frame()
  for (i in 1:length(pv_Var)){
    results = rbind(results, data.frame(strsplit(col_orgX[i], '_')[[1]][2], strsplit(col_orgX[i], '_')[[1]][3], col_orgX[i], pv_Var[[i]]$phenotype, nrow(Y), ncol(X)))
  }
  names(results) <- c('start', 'end', 'haplotype','p_value','nb_subjects', 'nb_haplotypes')

  return(results)
}
