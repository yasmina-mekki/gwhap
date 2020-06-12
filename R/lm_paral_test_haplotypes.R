#' Test haplotypes
#'
#' @description Perform three tests using R-base lm:
#' haplotype block model test
#' complete model haplotype test
#' variant (single) haplotype model test
#' the common haplotype is removed for the two first test
#' Assuming no covariates
#' No missing values allowed in X
#' Missing values allowed in Y
#'
#' @param X haplotype counts matrix
#' @param Y phenotype or residues matrix
#' @param kind which test to perform. Four values are possible: single, complete, block or all
#' @param verbose silent warning messages. FALSE by default.
#'
#' @return A list of the results of the three test. Each test results is represented by a data frame structure
#' The data frame of the block test contains the following information:
#' One block per line with the start and the end of the block position, p_value, number of subjects and the number of haplotypes of the block
#' The data frame of the block complete and variant test contains the following information:
#' One haplotype per line with start and the end of the block position,
#' the haplotype code, p_value, number of subjects and the number of haplotypes of the block
#' @importFrom stats lm
#' @export
#'
lm_par_test_haplotypes <- function(X, Y, kind='all', verbose=FALSE){

    # silent warning messages
    if(verbose == TRUE){options(warn=0)} else{options(warn=-1)}

    # this is a vector for all test x vector of pheno
    nsubj = nrow(Y) - unlist(lapply(lapply(Y, is.na), sum)) #nsubj = nrow(Y)

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

        # form_Y not used in formula (see after)
        form_Y = paste(col_orgY, collapse = ',')
        form_X = paste(varL, collapse = '+')

        # linear regression
        # not used since there are NA in the phenotypes
        # could be used to accelerate : only 1 pseudo inverse computations
        # instead of (# of pheno) computations
        # full_Lm = lm(sprintf("cbind( %s ) ~ %s ", form_Y, form_X), data = data.frame(cbind(Y, X)))
        full_Lm = lapply(sprintf("%s  ~ %s ", col_orgY, form_X), lm, data = data.frame(cbind(Y, X)))

        sum_Lm  = lapply(full_Lm, summary)
        full_LM = NULL
        names(sum_Lm) = col_orgY
    }

    final_results = list()

    if(kind == 'block' | kind == 'all'){

        # perform the test
        block_test_results = par_block_test(sum_Lm)

        # set the results into the right format
        block_test_results  = data.frame(
                                strsplit(varL[1], '_')[[1]][2],
                                strsplit(varL[1], '_')[[1]][3], 
                                paste(sep="_",
                                    strsplit(varL[1], '_')[[1]][1], 
                                    strsplit(varL[1], '_')[[1]][2],
                                    strsplit(varL[1], '_')[[1]][3]),
                                block_test_results, nsubj, length(varL))
        names(block_test_results) <- c('start', 'end', 'haplotype','p_value','nb_subjects', 'nb_haplotypes')

        # the rownames becomes the first col (phenotype name)
        block_test_results <- data.frame(phname = row.names(block_test_results),
                                         block_test_results)
        block_test_results['test'] <- 'block'
        rownames(block_test_results) <- c()

        final_results[['block']] = block_test_results
    }

  if(kind == 'complete' | kind == 'all'){

    # perform the test
    complete_test_results = par_complete_test(sum_Lm)
    #if(!is.list(complete_test_results) | length(complete_test_results)<=1){stop("Something went wrong ...")}

    # NAs in the phenotype yield a degenerated haplotype
    if (length(unique(unlist(lapply(complete_test_results, length)))) != 1){
        message(sprintf("Warning: NA in the phenotype yield a degenerated haplotype (%s). NA added in the results",
                  col_orgX[1]))
        # carefully dig in the list data. Some tests on some pheno may yield no results 
        # because of NA pheno data that degenerate the design matrix X
        # brut data.frame casting does nit work.
        #complete_test_results = as.data.frame(complete_test_results)
        fiter=TRUE
        for (i in names(complete_test_results)){
            if (fiter){
                fiter=FALSE
                acc = data.frame(complete_test_results[[i]])
                colnames(acc) = i
                acc = data.frame(haplotype=rownames(acc), acc)
            } else {
                tmp = data.frame(complete_test_results[[i]])
                colnames(tmp) = i
                tmp = data.frame(haplotype=rownames(tmp), tmp)
                acc = merge(acc, tmp, by='haplotype', all=TRUE)
            }
        }
        complete_test_results = acc
    } else {
        complete_test_results = data.frame(complete_test_results)
        complete_test_results = data.frame(haplotype=rownames(complete_test_results), complete_test_results)
    }

    # retain the number of haplotype
    nb_haplotypes=dim(complete_test_results)[1]

    # pivot table from wide to long
    complete_test_results = reshape(complete_test_results, direction="long", 
            times=col_orgY, timevar="phname",
            v.names="p_value",
            varying=col_orgY)
    # rename
    drops = c("id")
    complete_test_results= complete_test_results[,!(names(complete_test_results) %in% drops)]
    
    # add start stop and nb_haplo
    complete_test_results['start'] = as.numeric(sapply(strsplit(as.character(complete_test_results$haplotype), "_"), `[`,2))
    complete_test_results['end'] = as.numeric(sapply(strsplit(as.character(complete_test_results$haplotype), "_"), `[`,3))
    complete_test_results['nb_subjects'] = 0
    complete_test_results['nb_haplotypes'] =length(varL)

    #Set the num of subj by phenotype
    for (n in col_orgY){
        complete_test_results[complete_test_results$phname==n,'nb_subjects'] = nsubj[n]
    }

    # prepare the return
    complete_test_results['test'] <- 'complete'
    rownames(complete_test_results) <- c()
    final_results[['complete']] = complete_test_results
  }

  if(kind == 'single' | kind == 'all'){
    results = par_single_test(X, Y)  # This function may return NULL

    if (!is.null(results)){
        results['test'] <- 'single'
        rownames(results) <- c()
    }

    final_results[['single']] = results
    }

  return(final_results)
}

#' Haplotype block model test
#' @description Association between a given block and a phenotype
#' the common haplotype is removed for this test
#' warning : this function should not be used explicitly
#' users must only call the lm_test_haplotypes function and specify block as kind parameters value
#' @param sum_Lm summary of a linear model estimated using lm_test_haplotypes function
#'
#' @return Significance of the association
#' @export
#'
par_block_test <- function(sum_Lm){
  return(unlist(lapply(sum_Lm, function(y) {x = y$fstatistic
                                            pf(x[1], x[2], x[3], lower.tail = FALSE)[["value"]]})))
}

#' Complete model haplotype test
#' @description Association between each haplotype of a given block and a phenotype
#' the common haplotype is removed for this test
#' warning : this function should not be used explicitly
#' users must only call the lm_test_haplotypes function and specify complete as kind parameters value
#' @param sum_Lm summary of a linear model estimated using lm_test_haplotypes function
#'
#' @return Significance of the each association. A two-side p-value of the t-statitic is estimated
#' @export
#'
par_complete_test <- function(sum_Lm){
  return(lapply(sum_Lm, function(y, varL) {x = coefficients(y)
                                           x[2:nrow(x), ncol(x)]}, varL))
}

#' Variant (single) haplotype model test
#' @description Association between a given block and a phenotype.
#' The common haplotype is keeped for this test
#' @param X haplotype counts matrix
#' @param Y phenotype or residues matrix
#'
#' @return Significance of the each association. A two-side p-value of the t-statistic is estimated
#' @importFrom stats lm
#' @export
#'
par_single_test <- function(X, Y){

  # get the columns name
  #colnames(Y) = c('phenotype')
  col_orgY = colnames(Y)
  col_orgX = colnames(X)

  # nbsubj
  nsubj = nrow(Y) - unlist(lapply(lapply(Y, is.na), sum)) #nsubj = nrow(Y)

  # linear regression for multiple phenotypes ...
  ######################################################################
  # This is clearly the limiting part TO IMPROVE TO FIX
  ######################################################################
  varLM = lapply(col_orgX,
                 function(x) lapply(col_orgY,
                                    function(y) lm(formula = sprintf(" %s  ~ %s ", y, x), data = data.frame(cbind(Y, X)))
                                    )
                 )
  names(varLM)=names(X)


  # linear regression for several phenotypes
  # not usable because of NA in the phenotype
  # varLM = lapply(col_orgX,
  #    function(x) lm(formula = sprintf("cbind( %s ) ~ %s ", col_orgY, x), data = data.frame(cbind(Y, X)))
  #              )

  # get the summary of the model trained above
  sum_VLM = lapply(varLM, function(vlm_L) lapply(vlm_L, summary))

  # get the p values
  pv_Var = lapply(sum_VLM, function(z){names(z) = col_orgY[1:ncol(Y)]
                                       lapply(z, function(x){y = coefficients(x)
                                                             y[nrow(y), ncol(y)]})
                                       }
                   )

  # unexpected return with NULL
  # This NULL does not perturb the do.call(rbind... call'ee.
  if (length(unique(unlist(lapply(sum_VLM, length)))) != 1){
      message(sprintf("Warning: NA in the phenotype yield a degenerated haplotype (%s). This  block is skipped", col_orgX[1]))
      return(NULL)
  }

  results = as.data.frame(t(sapply(pv_Var, function(b) unlist(b))))

  # reshape
  results = reshape(results, direction="long", 
            times=col_orgY, timevar="phname",
            ids=row.names(results),v.names="p_value",
            varying=list(names(results)))
  # rename
  colnames(results)[which(names(results) == "id")] <- "haplotype"
  
  # add suppl.
  results['start'] = as.numeric(sapply(strsplit(results$haplotype, "_"), `[`, 2))
  results['end'] =   as.numeric(sapply(strsplit(results$haplotype, "_"), `[`, 3))
  results['nb_subjects']=unlist(lapply(nsubj, rep, ncol(X)))
  results['nb_haplotypes']=ncol(X)


  return(results)
}

