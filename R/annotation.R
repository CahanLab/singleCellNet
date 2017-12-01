# Patrick Cahan (C) 2017
# patrick.cahan@gmail.com


ann_eidToSym<-function# convert a gene symbols to a [single] eid
(syms,
 species='MM'){
  if(species=='MM'){
    library(org.Mm.eg.db);
    egSymbols<-as.list(org.Mm.egSYMBOL);
  }
  else{
    if(species=='HS'){
      library(org.Hs.eg.db);
      egSymbols<-as.list(org.Hs.egSYMBOL);
    }
    else{
      cat("Species ",species," not recognized. Enter MM or HS\n");
      return();
    }
  }
  unique(as.vector(unlist(egSymbols[syms])));
}  

# ksWrapper

#' run KS test given a ranked vector of genes and list of gene sets
#'
#' find gene sets enriched in provided ranked gene list
#'
#' @param queryGenes ranked gene list as a vector
#' @param geneSets named list of gene sets
#'
#' @return data frame of geneSet name, enrichment score, stat, nominal pvalue, holm-corrected pvalue, fdr, ?comma separated character of leading edge genes?
#' 
#' @export
#'
ks.wrapper<-function(
  queryGenes,
  geneSets)
{
 
    # determine the placing of geneSet genes in the query vectot
    q1<-lapply(geneSets, match, queryGenes)   
    ksRes<-lapply(q1, ks.test.2, 1:length(queryGenes))

    nom.p<-unlist(lapply(ksRes, "[[", "p.value"))
    edge<-unlist(lapply(ksRes, "[[", "edge"))
    es<-unlist(lapply(ksRes, "[[", "ES"))

    holms<-p.adjust(nom.p, "holm")
    fdrs<-p.adjust(nom.p, "fdr")

    data.frame(geneSet=names(ksRes), ES=es, p.value=nom.p, holm=holms, fdr=fdrs)
}

# ksWrapperSet

#' run ks.wrapper on list of diffExprs
#'
#' find gene sets enriched in provided ranked gene list
#'
#' @param diffExprList named list of diffExpr results, for example from running gnrAll()
#' @param geneSets named list of gene sets
#' @param sigType defaults to Holm
#' @param sigThresh defaults to 1e-5
#'
#' @return list of  and (2) list of ks.wrapper results
#' 
#' @export
#'
ks.wrapper.set<-function(
    diffExprList,
    geneSets)
{

    ansList<-list()
    qnames<-names(diffExprList)
    for(qname in qnames){
        qdiffRes<-diffExprList[[qname]]
        xgenes<-qdiffRes[order(qdiffRes$cval, decreasing=T),]
        x1<-rownames(xgenes)

        ansList[[qname]]<-ks.wrapper(x1, geneSets)
    }
    ansList
}


# ks.extract

#' extract a matrix of ES scores from rthe result of ks.wrapper.set
#'
#' get a matrix of ES scores for easy display
#'
#' @param enrRes for example from running ks.wrapper.set()
#' @param sigType defaults to Holm
#' @param sigThresh defaults to 1e-5
#'
#' @return atrix of ES for significant geneSets (cols=queryvectrs, rows=genesets),
#' 
#' @export
ks.extract<-function(
    enrRes,
    sigType='holm',
    sigThresh=1e-5)
{

    # first pass, determine the significant gene sets 
    rnames<-names(enrRes)
    sigGS<-vector()
    for(rname in rnames){
        x<-enrRes[[rname]]
        tmpAns<-x[which(x[,sigType]<sigThresh),]
        sigGS<-append(sigGS, as.vector(tmpAns$geneSet))
    }
    sigGS<-unique(sigGS)

    xi<-match(sigGS,x[,"geneSet"])

    ans<-matrix(0, nrow=length(sigGS), ncol=length(rnames))
    rownames(ans)<-sigGS
    colnames(ans)<-rnames
    for(rname in rnames){
        x<-enrRes[[rname]]        
        vals<-x[xi,"ES"]    
        ans[sigGS,rname]<-vals
    }
    ans
}





###
###This function is from https://github.com/franapoli/signed-ks-test
###

ks.test.2 <- function (x, y, ..., alternative = c("two.sided", "less", "greater"), exact = NULL, maxCombSize=10000) 
{
    alternative <- match.arg(alternative)
    DNAME <- deparse(substitute(x))
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 1L) 
        stop("not enough 'x' data")
    PVAL <- NULL
    if (is.numeric(y)) {
        DNAME <- paste(DNAME, "and", deparse(substitute(y)))
        y <- y[!is.na(y)]
        n.x <- as.double(n)
        n.y <- length(y)
        if (n.y < 1L) 
            stop("not enough 'y' data")
        if (is.null(exact)) {
            exact <- (n.x * n.y < maxCombSize)
            if(!exact)
                warning(paste("P-value not computed exactly because",
                              "of combined sample size"))
        }
        METHOD <- "Two-sample Kolmogorov-Smirnov test"
        TIES <- FALSE
        n <- n.x * n.y/(n.x + n.y)
        w <- c(x, y)
        z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
        if (length(unique(w)) < (n.x + n.y)) {
            if (exact) {
                warning("cannot compute exact p-value with ties")
                exact <- FALSE
            }
            else warning("p-value will be approximate in the presence of ties")
            z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
            TIES <- TRUE
        }
        STATISTIC <- switch(alternative, two.sided = max(abs(z)), 
            greater = max(z), less = -min(z))

        edge <- which.max(abs(z))
        ES <- z[edge]
        
        nm_alternative <- switch(alternative, two.sided = "two-sided", 
            less = "the CDF of x lies below that of y", greater = "the CDF of x lies above that of y")
        if (exact && (alternative == "two.sided") && !TIES) 
            PVAL <- 1 - .Call(stats:::C_pSmirnov2x, STATISTIC, n.x, n.y)
    }
    else {
        if (is.character(y)) 
            y <- get(y, mode = "function", envir = parent.frame())
        if (!is.function(y)) 
            stop("'y' must be numeric or a function or a string naming a valid function")
        METHOD <- "One-sample Kolmogorov-Smirnov test"
        TIES <- FALSE
        if (length(unique(x)) < n) {
            warning("ties should not be present for the Kolmogorov-Smirnov test")
            TIES <- TRUE
        }
        if (is.null(exact)) 
            exact <- (n < 100) && !TIES
        x <- y(sort(x), ...) - (0:(n - 1))/n
        STATISTIC <- switch(alternative, two.sided = max(c(x, 
            1/n - x)), greater = max(1/n - x), less = max(x))
        if (exact) {
            PVAL <- 1 - if (alternative == "two.sided")
                result = tryCatch({
                .C(C_pkolmogorov2x, p = as.double(STATISTIC), 
                  as.integer(n), PACKAGE = "stats")$p
                }, warning = function(w) {
                    warning(w)
                }, error = function(e) {
                    .Call(C_pKolmogorov2x, STATISTIC, n)
                }, finally = {
                })

            else {
                pkolmogorov1x <- function(x, n) {
                  if (x <= 0) 
                    return(0)
                  if (x >= 1) 
                    return(1)
                  j <- seq.int(from = 0, to = floor(n * (1 - 
                    x)))
                  1 - x * sum(exp(lchoose(n, j) + (n - j) * log(1 - 
                    x - j/n) + (j - 1) * log(x + j/n)))
                }
                pkolmogorov1x(STATISTIC, n)
            }
        }
        nm_alternative <- switch(alternative, two.sided = "two-sided", 
            less = "the CDF of x lies below the null hypothesis", 
            greater = "the CDF of x lies above the null hypothesis")
    }
    names(STATISTIC) <- switch(alternative, two.sided = "D", 
        greater = "D^+", less = "D^-")
    if (is.null(PVAL)) {
        pkstwo <- function(x, tol = 1e-06) {
            if (is.numeric(x)) 
                x <- as.double(x)
            else stop("argument 'x' must be numeric")
            p <- rep(0, length(x))
            p[is.na(x)] <- NA
            IND <- which(!is.na(x) & (x > 0))
            if (length(IND))
                p[IND] <- tryCatch({
                    tryRes <- .C(stats:::C_pkstwo, length(x[IND]), p = x[IND], 
                             as.double(tol), PACKAGE = "stats")$p
                }, warning = function(w) {
                    warning(w)
                }, error = function(e) {
                    tryRes <- .Call(stats:::C_pKS2, p = x[IND], tol)
                }, finally = {
                })
            p
        }
        PVAL <- ifelse(alternative == "two.sided", 1 - pkstwo(sqrt(n) * 
            STATISTIC), exp(-2 * n * STATISTIC^2))
    }
    PVAL <- min(1, max(0, PVAL))
    RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative, 
        method = METHOD, data.name = DNAME, ES = ES, edge = edge)
    class(RVAL) <- "htest"
    return(RVAL)
}



