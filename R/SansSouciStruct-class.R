#' Create an object of class 'SansSouciStruct'
#' 
#' @param struct A (dyadic) structure
#' @param leaves A list of leaves
#' @param truth An optional numeric vector of $m$ values in ${0,1}$, the status of each null hypothesis (0 means H0 is true, 1 means H1 is true). Typically used in simulations.
#' @return An object of class 'SansSouciStruct'
#' @export
#' @seealso SansSouciDyadic
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @examples
#' s <- 100
#' q <- 7
#' m <- s*2^q
#' 
#' dd <- dyadic.from.window.size(m, s, method = 2)
#' obj <- SansSouciStruct(dd$C, dd$leaf_list)
#' 
SansSouciStruct <- function(struct, leaves, truth = NULL) {
    m <- length(unlist(leaves))
    input = list(struct = struct, 
                 leaves = leaves,
                 m = m)
    input$truth <- truth
    obj <- structure(list(input = input,
                          parameters = NULL,
                          output = NULL), 
                     class = "SansSouciStruct")
    obj
}

#' Create an object of class 'SansSouciStruct' with a dyadic structure
#'
#' @param m A numeric value, the number of hypotheses
#' @param leaf_size An integer value, the number of hypotheses in each leaf.
#'   Defaults to NULL.
#' @param height An integer value, the desired maximal height of the tree.
#'   Defaults to NULL.
#' @return An object of class 'SansSouciStruct'
#' @details If both 'leaf_size' and 'height' are NULL then a binary tree with
#'   maximal height given the total number of hypotheses is built. If both
#'   'leaf_size' and 'height' are non NULL, then 'height' is ignored and a tree
#'   is built with the desired leaf size.
#' @param direction A character value, the direction used for building the tree.
#'   Must be one of  "bottom-up" or "top-down"
#' @param ... Further arguments to be passed to 'SansSouciStruct'
#' @seealso dyadic.from.leave.list dyadic.from.window.size dyadic.from.height
#'   dyadic.from.max.height
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#'   Post hoc false positive control for structured hypotheses. Scandinavian
#'   Journal of Statistics, 47(4), 1114-1148.
#' @export
#' @examples
#' s <- 100
#' q <- 7
#' m <- s*2^q
#'
#' obj <- SansSouciDyadic(m, leaf_size = s, direction = "top-down")
#' obj <- SansSouciDyadic(m, height = 8, direction = "top-down")
#' 
SansSouciDyadic <- function(m, leaf_size = NULL, height = NULL,
                            direction = c("bottom-up", "top-down"),
                            ...) {
    direction <- match.arg(direction)
    method <- switch(direction, 
                     "bottom-up" = 1,
                     "top-down" = 2)
    dd <- NULL
    if (!is.null(leaf_size)) {
        if (!is.null(height)) {
            warning("Both 'leaf_size' and 'height' provided: ignoring 'height'")
        }
        dd <- dyadic.from.window.size(m, s = leaf_size, method)
    } else {
        dd <- dyadic.from.height(m, H = height, method)
    }
    SansSouciStruct(struct = dd$C, leaves = dd$leaf_list, ...)
}


#' Basic methods for class `SansSouciStruct`
#' 
#' @param object An object of class \code{SansSouciStruct}
#' @name SansSouciStruct-methods
#' @examples
#' s <- 100
#' q <- 7
#' m <- s*2^q
#'
#' obj <- SansSouciDyadic(m, leaf_size = s, direction = "top-down")
#' nHyp(obj)
#' nLeaves(obj)
#' mu <- gen.mu.leaves(m = m, K1 = 8, d = 0.9, grouped = TRUE, 
#'   setting = "const", barmu = 3, leaf_list = obj$input$leaves)
#' pvalues <- gen.p.values(m = m, mu = mu, rho = 0)
#' 
#' alpha <- 0.05
#' res <- fit(obj, alpha, pvalues, "DKWM")
#' label(res)
#' str(pValues(res))
#' 

NULL
#> NULL
#' `nHyp`: get the number of hypotheses
#' 
#' @rdname SansSouciStruct-methods
#' @param object An object of class `SansSouciStruct`
#' @export
nHyp.SansSouciStruct <- function(object) {
    object$input$m
}

#' @describeIn all-generics Get the number of leaves
#' @param object An object. See individual methods for specifics
#' @export
nLeaves <- function(object) UseMethod("nLeaves")

#' `nLeaves`: get the number of leaves
#' 
#' @rdname SansSouciStruct-methods
#' @param object An object of class `SansSouciStruct`
#' @export
nLeaves.SansSouciStruct <- function(object) {
    length(object$input$leaves)
}


#' Fit SansSouciStruct object
#' 
#' @param object An object of class `SansSouciStruct`
#' @param alpha Target risk (JER) level
#' @param p.values A vector of length `nHyp(object)`, 
#' @param family A character value describing how the number of true nulls in a set is estimated. Can be either:
#' - "DKWM": estimation by the Dvoretzky-Kiefer-Wolfowitz-Massart inequality (related to the Storey estimator of the proportion of true nulls), valid for independent p-values
#' - "HB": estimation by the Holm-Bonferroni method, always valid
#' - "trivial": dummy estimation as the the size of the set
#' - "Simes": estimation via the Simes inequality, valid for positively-dependent (PRDS) p-values
#' - "Oracle": true number of true null hypotheses Truth" must be available in `object$input$truth`
#' @param flavor A character value which can be 
#' - "tree" (default value): the reference family is the entire tree structure
#' - "partition": the reference family is the partition corresponding to the leaves of the tree
#' @param ... Not used
#' @return A 'fitted' object of class 'SansSouciStruct'. It is a list of three elements
#'  - input: see [SansSouciStruct]
#'  - param: the input parameters, given as a list
#'  - output: a list of two elements
#'    - p.values: the input argument 'p.values'
#'    - ZL: the output of the "zeta function" associated to the input parameter 'family', see e.g. [zeta.DKWM]
#' @details In the particular case where `family=="Simes"` or `family=="Oracle"`, the return value is actually of class `SansSouci` and not `SansSouciStruct`
#' @seealso zeta.DKWM zeta.HB, zeta.tricial
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Dvoretzky, A., Kiefer, J., and Wolfowitz, J. (1956). Asymptotic minimax character of the sample distribution function and of the classical multinomial estimator. The Annals of Mathematical Statistics, pages 642-669.
#' @references Holm, S. A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics 6 (1979), pp. 65-70.
#' @references Massart, P. (1990). The tight constant in the Dvoretzky-Kiefer-Wolfowitz inequality. The Annals of Probability, pages 1269-1283.
#' @export
#' 
#' @examples 
#' s <- 100
#' q <- 7
#' m <- s*2^q
#' obj <- SansSouciDyadic(m, leaf_size = s, direction = "top-down")
#' 
#' mu <- gen.mu.leaves(m = m, K1 = 8, d = 0.9, grouped = TRUE, 
#'   setting = "const", barmu = 3, leaf_list = obj$input$leaves)
#' pvalues <- gen.p.values(m = m, mu = mu, rho = 0)
#' 
#' alpha <- 0.05
#' S1 <- which(mu != 0)
#' 
#' res_DKWM <- fit(obj, alpha, pvalues, "DKWM")
#' predict(res_DKWM, S = S1, what = "FP")
#' 
#' res_Simes <- fit(obj, alpha, pvalues, "Simes")
#' predict(res_Simes, S = S1, what = "FP")
#' 
fit.SansSouciStruct <- function(object, alpha, p.values, 
                                family = c("DKWM", "HB", "trivial", "Simes", "Oracle"), 
                                flavor = c("tree", "partition"), ...) {
    family <- match.arg(family)
    flavor <- match.arg(flavor)
    if (family == "Oracle") {
        truth <- object$input$truth
        if (is.null(truth)) {
            stop("'truth' should be available for 'Oracle'. See ?SansSouci")
        }
    }
    m <- nHyp(object)
    stopifnot(length(p.values) == m)
    params <- list(alpha = alpha, family = family, flavor = flavor)
    if (family %in% c("DKWM", "HB", "trivial")) {
        zeta <- switch(family,
                       "DKWM" = zeta.DKWM,
                       "HB" = zeta.HB,
                       "trivial" = zeta.trivial)
        struct <- object$input$struct
        leaves <- object$input$leaves
        if (flavor == "partition") {
            struct <- struct[length(struct)]
        }
        ZL <- zetas.tree(C = struct, leaf_list = leaves, 
                         method = zeta, pvalues = p.values, alpha = alpha)
        output <- list(p.values = p.values, ZL = ZL)
    } else if (family == "Simes") {
        thr <- t_linear(alpha, seq_len(m), m) ## NB: we are forcing K=m
        output <- list(p.values = p.values, thr = thr)
        class(object) <- "SansSouci" ## so that bound.SansSouci is called...
    } else if (family == "Oracle") {
        thr <- object$input$truth
        output <- list(p.values = p.values, thr = thr)
        class(object) <- "SansSouci" ## so that bound.SansSouci is called...
    }
    object$output <- output
    object$parameters <- params
    object
}

#' `label` Get the label of a post hoc method
#' 
#' @rdname SansSouciStruct-methods
#' @param object An object of class `SansSouciStruct`
#' @export
label.SansSouciStruct <- function(object) {
    param <- object$parameters
    if (is.null(param)) {
        return(NULL)
    }
    lab <- param$family
    flv <- param$flavor
    
    if (!(lab %in% c("Simes", "Oracle"))) {
        lab <- sprintf("%s(%s)", lab, flv)
    } 
    return(lab)
}

#' `pValues`: get p-values
#' 
#' @rdname SansSouciStruct-methods
#' @param object An object of class `SansSouci`
#' @export
pValues.SansSouciStruct <- function(object) {
    object$output$p.values
}

#' @export
predict.SansSouciStruct <- function(object, S = seq_len(nHyp(object)), 
                            what = c("TP", "FDP"), n_out = 1, 
                            order_by_p = TRUE, ...) {
    lab <- label(object)
    if (order_by_p) {
        ordering <- order(pValues(object))
        S <- ordering[S]
    }
    
    if (max(S) > nHyp(object)) {
        stop("'S' is not a subset of hypotheses")
    }
    struct <- object$input$struct
    leaves <- object$input$leaves
    ZL <- object$output$ZL
    
    if (object$param$flavor == "partition") {
        struct <- struct[length(struct)]
    }
    
    idxs <- length(S)
    if (n_out > 1) {
        # idxs <- as.integer(seq(from = 1, to = length(S), length.out = n_out))
        idxs <- 10^seq(from = log10(1), to = log10(length(S)), length.out = n_out)
        idxs <- unique(as.integer(idxs))
    }
    max_FP <- numeric(length(idxs))
    for (ii in seq(along = idxs)) {
        max_FP[ii] <- V.star(S[1:idxs[ii]], C = struct, 
                             ZL = ZL, leaf_list = leaves)
    }
    bounds <- formatBounds(max_FP, idxs = idxs, lab = lab, what = what, all = TRUE)
    if (n_out == 1) {
        bounds <- bounds[, "bound"]
        if (length(bounds) > 1) {
            names(bounds) <- what
        }
    }
    bounds
} 