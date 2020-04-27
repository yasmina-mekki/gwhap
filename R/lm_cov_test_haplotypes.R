#' Test haplotypes
#'
#' @description Perform three tests using R-base lm:
#' haplotype bloc model test
#' complete model haplotype test
#' variant (single) haplotype model test
#' the common haplotype is removed for the two first test
#' Assuming no covariates
#' No missing values allowed in X
#' Missing values allowed in Y
#' 
#' @param X0 covariates matrix
#' @param X haplotype counts matrix
#' @param Y phenotype or residues matrix
#' @param kind which test to perform. Four values are possible: single, complete, bloc or all
#' @param verbose silent warning messages. FALSE by default.
#'
#' @return A list of the results of the three test. Each test results is represented by a data frame structure
#' The data frame of the bloc test contains the following information:
#' One bloc per line with the start and the end of the bloc position, p_value, number of subjects and the number of haplotypes of the bloc
#' The data frame of the bloc complete and variant test contains the following information:
#' One haplotype per line with start and the end of the bloc position,
#' the haplotype code, p_value, number of subjects and the number of haplotypes of the bloc
#' @importFrom stats lm
#' @export
#'
lm_cov_test_haplotypes = function(X, Y, X0= null, kind='all', verbose=FALSE){

  # silent warning messages
  if(verbose == TRUE){options(warn=0)} else{options(warn=-1)}
  if(is.null(X0)) return(lm_test_haplotypes(X, Y, kind, verbose))
  else { 

    col_orgY = colnames(Y)
    col_orgX = colnames(X)
    col_cov  =colnames(X0)
    form_X0 = paste(col_cov, collapse = '+')
    cov_Lm = lapply(sprintf("%s ~ %s ", col_orgY, form_X0),
                    lm,
                    data = data.frame(cbind(Y, X0,X)))
    
    if(kind != 'single'){

    # colnames(Y) = c('phenotype')

    # get the columns name
    # get the sum of each columns value and identify the common haplotype
    # Que signifi les nombre 1, 2, 0 pour les haplotypes (dummification) 0, 1 pour présence ou pas mais le 2 ?
    allelcount = unlist(apply(X, 2, sum))
    varL = colnames(X)

    # get the index of the max of allel count
    # get haplotypes columns without the common one.
    if (length(varL) > 1) {varL = colnames(X)[setdiff(1:ncol(X), which.max(allelcount))]}

    #form_Y = paste(col_orgY, collapse = ',')
    form_X = paste(c(form_X0,varL), collapse = '+')
    fupdate= sprintf(". ~ . + %s", paste(varL, collapse = '+'))
    # linear regression
    
    full_Lm = lapply(sprintf("%s  ~ %s ", col_orgY, form_X), 
                     lm, data = data.frame(cbind(Y, X0, X)))
    
    
    
  }

  final_results = list()
  if(kind == 'bloc' | kind == 'all'){
    
    # perform the test
    aov_L=setNames(Map(anova,cov_Lm,full_Lm),col_orgY)
    
    bloc_test_results = lapply(aov_L, function(x)  x[2,6])

    # set the results into the right format
    results  = lapply(bloc_test_results, function(x) data.frame(
                          start=strsplit(col_orgX[1], '_')[[1]][2], 
                          end=strsplit(col_orgX[1], '_')[[1]][3], 
                          p_value=x, 
                          nb_subjects=nrow(Y),
                          nb_haplotypes=ncol(X)))
#    names(results) <- c('start', 'end', 'p_value','nb_subjects', 'nb_haplotypes')

    final_results[['bloc']] = results
  }

  if(kind == 'complete' | kind == 'all'){

    # perform the test
    sum_Lm  = lapply(full_Lm, summary)
    names(sum_Lm) = col_orgY
    complete_test_results = setNames(lapply(sum_Lm,function(x) coefficients(x)[row.names(coefficients(x)) %in% col_orgX,]),col_orgY)

    #if(!is.list(complete_test_results) | length(complete_test_results)<=1){stop("Something went wrong ...")}

    # set the results into the right format
    results  = lapply(complete_test_results,
    function(x) data.frame(start=strsplit(row.names(x), '_')[[1]][2],
               end=strsplit(row.names(x), '_')[[1]][3],
               haplotype=row.names(x),
               p_value=x[,4],
               nb_subjects=nrow(Y),
               nb_haplotypes=ncol(X)))
      
    
    
    #c('phenotype','start', 'end', 'haplotype','p_value','nb_subjects', 'nb_haplotypes')

    final_results[['complete']] = results
  }

  if(kind == 'single' | kind == 'all'){
    
    fupdate_L= sprintf(". ~ . + %s", col_orgX)
    
    aov_L=setNames(Map(anova,cov_Lm,full_Lm),col_orgY)
  
    pv_Var = setNames(lapply(cov_Lm, 
                             function(x) setNames(lapply(
      fupdate_L,function(y) 
        anova(x,
               lm(update(formula(x$terms),y),
           data=cbind(Y, X0,X)),
           test="F")), col_orgX)),
      col_orgY)
    
    
    # set the output into the right format
    results  = lapply(pv_Var,
                      function(x) { pvs=unlist(lapply(x,function(y) y[2,6])) ;
                        data.frame( start=strsplit(names(pvs), '_')[[1]][2],
                                             end=strsplit(names(pvs), '_')[[1]][3],
                                             haplotype=names(pvs),
                                             p_value=pvs,
                                             nb_subjects=nrow(Y),
                                             nb_haplotypes=length(names(pvs)))
                      })
    
    #names(results) <- c('start', 'end', 'haplotype','p_value','nb_subjects', 'nb_haplotypes')
    
      
    final_results[['single']] = results
    
  }

  return(final_results)
  }
}


