

#' @title
#' Predict DNA methylation age using CTS clocks
#'
#' @description
#' This is the function for computing DNA methylation age using CTS (cell type
#' specific) clocks. The inputs including DNAm matrix, the CTS clock you want to
#' use, is your DNAm data from bulk tissue sample or sorted cell sample, cell
#' type fraction matrix if you want to use Neu-In/Glia-In/Brain clock, tissue of
#' your DNAm data samples and the number of cores if you want to do parallel computing.
#'
#'
#'
#' @param data.m
#' A DNAm matrix (row: CpGs, column: samples) of the samples you want to get
#' DNAm age predicted by a CTS clock.
#'
#' @param CTSclock
#' Which CTSclock you want to use to predict DNAm age ('Neu-In', 'Neu-Ex',
#' 'Glia-In', 'Glia-Ex', 'Brain', 'Hep', 'Liver').
#'
#' @param dataType
#' Type of the samples in your DNAm data ('bulk' or 'sorted').
#'
#' @param CTF.m
#' Cell type fraction matrix (row: samples in the same order with data.m,
#' column: cell types).
#'
#' @param tissue
#' What tissue are your samples from ('brain', 'otherTissue').
#'
#' @param coreNum
#' The number of cores to use. If set to NULL then it will be the number
#' generated by parallel::detectCores().
#'
#'
#'
#' @return A vecter of the predicted DNAm ages.
#'
#'
#'
#' @import glmnet
#'
#' @import parallel
#'
#'
#'
#' @references Huige Tong, Xiaolong Guo, Qi Luo and Andrew E Teschendorff. 2024. "Cell-type specific epigenetic clocks to quantify biological age at cell-type resolution." Submitted.
#'
#'
#'
#' @export
#'
#'
#'
#' @examples
#' data(MurphyGSE88890)
#' agePred.v = CTSclockAge(beta.m, CTSclock = 'Neu-In', dataType = 'bulk',
#' CTF.m = NULL, tissue = 'brain')
#' plot(phenotype.df$Age, agePred.v)
#'
#' data(ExampleData_Liver)
#' agePred.v = CTSclockAge(Test.m, CTSclock = 'Hep', dataType = 'bulk',
#' CTF.m = NULL, tissue = 'otherTissue')
#' plot(Age, agePred.v)
#'
#'

CTSclockAge = function(data.m, CTSclock = c('Neu-In', 'Neu-Ex', 'Glia-In', 'Glia-Ex', 'Brain', 'Hep', 'Liver'),
                       dataType = c('bulk', 'sorted'), CTF.m = NULL, tissue = c('brain', 'otherTissue'), coreNum = NULL){

  if (is.null(coreNum)){coreNum = detectCores()}

  if (CTSclock %in% c('Neu-In', 'Neu-Ex', 'Glia-In', 'Glia-Ex', 'Brain')){
    data(list = paste0(CTSclock, 'Coef'))

    if (CTSclock %in% c('Neu-Ex', 'Glia-Ex')){
      DNAmAgePred.v = pred(data.m, clockCoef.df)
    }else if(CTSclock %in% c('Neu-In', 'Glia-In', 'Brain')){
      data.m = ProcessData(data.m, dataType = dataType, tissue = tissue, CTF.m = CTF.m, coreNum = coreNum)
      DNAmAgePred.v = pred(data.m, clockCoef.df)
    }


  }else if(CTSclock %in% c('Hep', 'Liver')){
    data(list = paste0(CTSclock, 'Clock'))
    Clock.glm = eval(parse(text = paste0(CTSclock, 'Clock.glm')))
    ClockCpGs.v <- rownames(Clock.glm$beta)
    TestSetCpGs.v <- rownames(data.m)
    # Trim clock. Remember to reload the clock when applying new data sets
    idx <- match(ClockCpGs.v, TestSetCpGs.v)
    Clock.glm$beta <- eval(parse(text = paste0(CTSclock, 'Clock.glm')))$beta[!is.na(idx),, drop = FALSE]
    Clock.glm$dim <- c(sum(!is.na(idx)),1)
    # Trim sample data
    data.m <- data.m[na.omit(idx),]
    DNAmAgePred.v <- as.vector(predict.glmnet(Clock.glm, newx = t(data.m)))
  }

  return(DNAmAgePred.v)
}



#### Auxilliary functions
#' @import HiBED
#'
ProcessData = function(data.m, dataType = c('bulk', 'sorted'), tissue = c('brain', 'otherTissue'), CTF.m = NULL, coreNum = coreNum){

  if(dataType == 'sorted'){
    ## Normalize the data
    dataSD.v = unlist(mclapply(1:nrow(data.m),function(i) sd(data.m[i,]), mc.cores = coreNum))
    dataZ.m = (data.m - rowMeans(data.m))/dataSD.v
    return(dataZ.m)
  }else if(dataType == 'bulk'){
    if (is.null(CTF.m)){
      if(tissue == 'brain'){
        ## Estimate the cell type fractions
        estF.m = HiBED_deconvolution(data.m, h=1)
        estF.m = estF.m/100
        colnames(estF.m) = c('EndoStrom', 'Glia', 'Neu')
        estF.m = estF.m[, c('Neu', 'Glia', 'EndoStrom')]
        CTF.m = as.matrix(estF.m)
      }else if(tissue == 'otherTissue'){
        stop("Cell type fraction matrix (CTF.m) is missing. If you don't have it, you are recommended to use R package EpiSCORE, EpiDISH or some other deconvolution algorithms to estimate the cell type fractions of your samples.")
      }
    }
    ## Adjust for cell type fractions and normalize the data
    print("Processing the data may take some time. Don't worry. :)")
    lm.o = lm(t(data.m) ~ CTF.m)
    res.m = t(lm.o$res)
    resSD.v = unlist(mclapply(1:nrow(res.m),function(i) sd(res.m[i,]), mc.cores = coreNum))
    dataZ.m = (res.m - rowMeans(res.m))/resSD.v
    return(dataZ.m)
  }
}

pred = function(data.m, clockCoef.df){
  modelIntercept = clockCoef.df$coef[1]
  clockCoef.df = clockCoef.df[-1,]
  rownames(clockCoef.df) = clockCoef.df$probe
  clockCoef.df = clockCoef.df[intersect(rownames(data.m), clockCoef.df$probe),, drop = F]
  data.m = data.m[clockCoef.df$probe, , drop = F]
  pred.v = as.vector(clockCoef.df$coef %*% data.m + modelIntercept)
  return(pred.v)
}
