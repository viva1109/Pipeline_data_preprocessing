  ancom_func<-function(data_table,sig=0.05,mcor=2){
#browser()
    ff<-factor(data_table$y)
    df<-cbind(data.frame(data_table$dataSet),ff)
    ancom.out <- ANCOM2( OTUdat = df, sig = sig,multcorr=mcor)
    return(ancom.out)
  }

ANCOM2<-function (OTUdat, sig = 0.05, multcorr = 3, tau = 0.02, theta = 0.1, repeated = FALSE) 
{
    num_col <- ncol(OTUdat)
    if (repeated == FALSE) {
        colnames(OTUdat)[num_col] <- "Group"
        num_OTU <- ncol(OTUdat) - 1
        sub_drop <- data.frame(nm_drop = "N/A")
        sub_keep <- data.frame(nm_keep = "All subjects")
        colnames(sub_drop) <- "Subjects removed"
        colnames(sub_keep) <- "Subjects retained"
        n_summary <- paste0("No subjects entirely removed (not a repeated-measures design)")
    }
    else {
        colnames(OTUdat)[num_col - 1] <- "Group"
        colnames(OTUdat)[num_col] <- "ID"
        OTUdat$ID <- factor(OTUdat$ID)
        num_OTU <- ncol(OTUdat) - 2
        crossTab <- table(OTUdat$Group, OTUdat$ID) == 0
        id_drop <- apply(crossTab, 2, FUN = function(x) any(x))
        nm_drop <- names(which(id_drop))
        idx_drop <- OTUdat$ID %in% nm_drop
        OTUdat <- OTUdat[idx_drop == FALSE, ]
        if (nrow(OTUdat) == 0) {
            stop("Too many missing values in data, all subjects dropped")
        }
        OTUdat$ID <- droplevels(OTUdat$ID)
        num_dropped <- sum(id_drop)
        num_retain <- length(id_drop) - num_dropped
        sub_drop <- data.frame(nm_drop = paste(nm_drop, collapse = ", "))
        sub_keep <- data.frame(nm_keep = paste(levels(OTUdat$ID), 
            collapse = ", "))
        colnames(sub_drop) <- "Subjects removed"
        colnames(sub_keep) <- "Subjects retained"
        n_summary <- paste0("Analysis used ", num_retain, " subjects (", 
            num_dropped, " were removed due to incomplete data)")
    }
    OTUdat$Group <- factor(OTUdat$Group)
    OTUdat <- data.frame(OTUdat[which(is.na(OTUdat$Group) == 
        FALSE), ], row.names = NULL)
    W.detected <- ancom.detect2(OTUdat, num_OTU, sig, multcorr, 
        ncore = 1)
    W_stat <- W.detected
    if (num_OTU < 10) {
        detected <- colnames(OTUdat)[which(W.detected > num_OTU - 
            1)]
    }
    else {
        if (max(W.detected)/num_OTU >= theta) {
            c.start <- max(W.detected)/num_OTU
            cutoff <- c.start - c(0.05, 0.1, 0.15, 0.2, 0.25)
            prop_cut <- rep(0, length(cutoff))
            for (cut in 1:length(cutoff)) {
                prop_cut[cut] <- length(which(W.detected >= num_OTU * 
                  cutoff[cut]))/length(W.detected)
            }
            del <- rep(0, length(cutoff) - 1)
            for (ii in 1:(length(cutoff) - 1)) {
                del[ii] <- abs(prop_cut[ii] - prop_cut[ii + 1])
            }
            if (del[1] < tau & del[2] < tau & del[3] < tau) {
                nu = cutoff[1]
            }
            else if (del[1] >= tau & del[2] < tau & del[3] < 
                tau) {
                nu = cutoff[2]
            }
            else if (del[2] >= tau & del[3] < tau & del[4] < 
                tau) {
                nu = cutoff[3]
            }
            else {
                nu = cutoff[4]
            }
            up_point <- min(W.detected[which(W.detected >= nu * 
                num_OTU)])
            W.detected[W.detected >= up_point] <- 99999
            W.detected[W.detected < up_point] <- 0
            W.detected[W.detected == 99999] <- 1
            detected <- colnames(OTUdat)[which(W.detected == 
                1)]
        }
        else {
            W.detected <- 0
            detected <- "No significant OTUs detected"
        }
    }
    results <- list(W = W_stat, detected = detected, dframe = OTUdat, 
        repeated = repeated, n_summary = n_summary, sub_drop = sub_drop, 
        sub_keep = sub_keep)
    class(results) <- "ancom"
    return(results)
}

ancom.detect2<-function (otu_data, n_otu, alpha, multcorr, ncore) 
{
    if (ncol(otu_data) == n_otu + 1) {
        Group <- otu_data[, ncol(otu_data)]
        ID <- rep(1, nrow(otu_data))
        repeated <- FALSE
        fformula <- formula("lr ~ Group")
    }
    else if (ncol(otu_data) == n_otu + 2) {
        Group <- otu_data[, ncol(otu_data) - 1]
        ID <- otu_data[, ncol(otu_data)]
        repeated <- TRUE
        fformula <- formula("lr ~ Group | ID")
    }
    else {
        stop("Problem with data. Dataset should contain OTU abundances, groups, \n         and optionally an ID for repeated measures.")
    }
    if (repeated == FALSE) {
        if (length(unique(Group)) == 2) {
            tfun <- wilcox.exact.formula2
        }
        else {
            tfun <- stats::kruskal.test
        }
    }
    else {
        tfun <- stats::friedman.test
    }
    if (FALSE) {
        registerDoParallel(cores = ncore)
        aa <- bb <- NULL
        logratio.mat <- foreach(bb = 1:n_otu, .combine = "rbind", 
            .packages = "foreach") %:% foreach(aa = 1:n_otu, 
            .combine = "c", .packages = "foreach") %dopar% {
            if (aa == bb) {
                p_out <- NA
            }
            else {
                data.pair <- otu_data[, c(aa, bb)]
                lr <- log((1 + as.numeric(data.pair[, 1]))/(1 + 
                  as.numeric(data.pair[, 2])))
                lr_dat <- data.frame(lr = lr, Group = Group, 
                  ID = ID)
                p_out <- tfun(formula = fformula, data = lr_dat)$p.value
            }
            p_out
        }
        rownames(logratio.mat) <- colnames(logratio.mat) <- NULL
    }
    else {
        logratio.mat <- matrix(NA, nrow = n_otu, ncol = n_otu)
        for (ii in 1:(n_otu - 1)) {
            for (jj in (ii + 1):n_otu) {
                data.pair <- otu_data[, c(ii, jj)]
                lr <- log((1 + as.numeric(data.pair[, 1]))/(1 + 
                  as.numeric(data.pair[, 2])))
                lr_dat <- data.frame(lr = lr, Group = Group, 
                  ID = ID)
                logratio.mat[ii, jj] <- tfun(formula = fformula, 
                  data = lr_dat)$p.value
            }
        }
        ind <- lower.tri(logratio.mat)
        logratio.mat[ind] <- t(logratio.mat)[ind]
    }
    logratio.mat[which(is.finite(logratio.mat) == FALSE)] <- 1
    mc.pval <- t(apply(logratio.mat, 1, function(x) {
        s <- p.adjust(x, method = "BH")
        return(s)
    }))
    a <- logratio.mat[upper.tri(logratio.mat, diag = FALSE) == 
        TRUE]
    b <- matrix(0, ncol = n_otu, nrow = n_otu)
    b[upper.tri(b) == T] <- p.adjust(a, method = "BH")
    diag(b) <- NA
    ind.1 <- lower.tri(b)
    b[ind.1] <- t(b)[ind.1]
    if (multcorr == 1) {
        W <- apply(b, 1, function(x) {
            subp <- length(which(x < alpha))
        })
    }
    else if (multcorr == 2) {
        W <- apply(mc.pval, 1, function(x) {
            subp <- length(which(x < alpha))
        })
    }
    else if (multcorr == 3) {
        W <- apply(logratio.mat, 1, function(x) {
            subp <- length(which(x < alpha))
        })
    }
    return(W)
}

wilcox.exact2<-function (x, y = NULL, alternative = c("two.sided", "less", "greater"), 
    mu = 0, paired = FALSE, exact = F, conf.int = FALSE, conf.level = 0.95, 
    ...) {
    alternative <- match.arg(alternative)
    if (!missing(mu) && ((length(mu) > 1) || !is.finite(mu))) 
        stop("mu must be a single number")
    if (conf.int) {
        if (!((length(conf.level) == 1) && is.finite(conf.level) && 
            (conf.level > 0) && (conf.level < 1))) 
            stop("conf.level must be a single number between 0 and 1")
    }
    MIDP <- NULL
    if (!is.numeric(x)) 
        stop("`x' must be numeric")
    if (!is.null(y)) {
        if (!is.numeric(y)) 
            stop("`y' must be numeric")
        DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
        if (paired) {
            if (length(x) != length(y)) 
                stop("x and y must have the same length")
            OK <- complete.cases(x, y)
            x <- x[OK] - y[OK]
            y <- NULL
        }
        else {
            x <- x[is.finite(x)]
            y <- y[is.finite(y)]
        }
    }
    else {
        DNAME <- deparse(substitute(x))
        if (paired) 
            stop("y missing for paired test")
        x <- x[is.finite(x)]
    }
    if (length(x) < 1) 
        stop("not enough (finite) x observations")
    CORRECTION <- 0
    if (is.null(y)) {
        METHOD <- "Exact Wilcoxon signed rank test"
        x <- x - mu
        ZEROES <- any(x == 0)
        if (ZEROES) 
            x <- x[x != 0]
        n <- length(x)
        if (is.null(exact)) 
            exact <- (n < 50)
        r <- rank(abs(x))
        STATISTIC <- sum(r[x > 0])
        names(STATISTIC) <- "V"
        TIES <- (length(r) != length(unique(r)))
        if (exact) {
            PVAL <- switch(alternative, two.sided = pperm(STATISTIC, 
                r, n, alternative = "two.sided", pprob = TRUE), 
                greater = pperm(STATISTIC, r, n, alternative = "greater", 
                  pprob = TRUE), less = pperm(STATISTIC, r, n, 
                  alternative = "less", pprob = TRUE))
            MIDP <- PVAL$PPROB
            PVAL <- PVAL$PVALUE
            if (conf.int && !is.na(x)) {
                x <- x + mu
                alpha <- 1 - conf.level
                diffs <- outer(x, x, "+")
                diffs <- sort(diffs[!lower.tri(diffs)])/2
                if (TIES) {
                  fs <- function(d) {
                    xx <- x - d
                    sum(rank(abs(xx))[xx > 0])
                  }
                  w <- sapply(diffs, fs)
                }
                else {
                  w <- sum(rank(abs(x))):1
                }
                cint <- switch(alternative, two.sided = {
                  qu <- qperm(alpha/2, r, n)
                  ql <- qperm(1 - alpha/2, r, n)
                  if (qu <= min(w)) lci <- max(diffs) else lci <- min(diffs[w <= 
                    qu])
                  if (ql >= max(w)) uci <- min(diffs) else uci <- max(diffs[w > 
                    ql])
                  c(uci, lci)
                }, greater = {
                  ql <- qperm(1 - alpha, r, n)
                  if (ql >= max(w)) uci <- min(diffs) else uci <- max(diffs[w > 
                    ql])
                  c(uci, Inf)
                }, less = {
                  qu <- qperm(alpha, r, n)
                  if (qu <= min(w)) lci <- max(diffs) else lci <- min(diffs[w <= 
                    qu])
                  c(-Inf, lci)
                })
                attr(cint, "conf.level") <- conf.level
                wmean <- sum(r)/2
                ESTIMATE <- mean(c(min(diffs[w <= ceiling(wmean)]), 
                  max(diffs[w > wmean])))
                names(ESTIMATE) <- "(pseudo)median"
            }
        }
        else {
            METHOD <- "Asymptotic Wilcoxon signed rank test"
            wmean <- sum(r)/2
            wvar <- sum(r^2)/4
            PVAL <- pnorm((STATISTIC - wmean)/sqrt(wvar))
            if (alternative == "greater") 
                PVAL <- 1 - PVAL
            if (alternative == "two.sided") 
                PVAL <- 2 * min(PVAL, 1 - PVAL)
            if (conf.int && !is.na(x)) {
                x <- x + mu
                alpha <- 1 - conf.level
                mumin <- min(x)
                mumax <- max(x)
                CORRECTION.CI <- 0
                wdiff <- function(d, zq) {
                  xd <- x - d
                  xd <- xd[xd != 0]
                  nx <- length(xd)
                  dr <- rank(abs(xd))
                  zd <- sum(dr[xd > 0])
                  zd <- (zd - wmean)/sqrt(wvar)
                  zd - zq
                }
                cint <- switch(alternative, two.sided = {
                  l <- uniroot(wdiff, c(mumin, mumax), tol = 1e-04, 
                    zq = qnorm(alpha/2, lower.tail = FALSE))$root
                  u <- uniroot(wdiff, c(mumin, mumax), tol = 1e-04, 
                    zq = qnorm(alpha/2))$root
                  c(l, u)
                }, greater = {
                  l <- uniroot(wdiff, c(mumin, mumax), tol = 1e-04, 
                    zq = qnorm(alpha, lower.tail = FALSE))$root
                  c(l, +Inf)
                }, less = {
                  u <- uniroot(wdiff, c(mumin, mumax), tol = 1e-04, 
                    zq = qnorm(alpha))$root
                  c(-Inf, u)
                })
                attr(cint, "conf.level") <- conf.level
                ESTIMATE <- uniroot(wdiff, c(mumin, mumax), tol = 1e-04, 
                  zq = 0)$root
                names(ESTIMATE) <- "(pseudo)median"
            }
        }
    }
    else {
        if (length(y) < 1) 
            stop("not enough y observations")
        METHOD <- "Exact Wilcoxon rank sum test"
        r <- rank(c(x - mu, y))
        n.x <- length(x)
        n.y <- length(y)
        if (is.null(exact)) 
            exact <- (n.x < 50) && (n.y < 50)
        STATISTIC <- sum(r[seq(along = x)]) - n.x * (n.x + 1)/2
        names(STATISTIC) <- "W"
        TIES <- (length(r) != length(unique(r)))
        if (exact) {
            PVAL <- switch(alternative, two.sided = pperm(STATISTIC + 
                n.x * (n.x + 1)/2, r, n.x, alternative = "two.sided", 
                pprob = TRUE), greater = pperm(STATISTIC + n.x * 
                (n.x + 1)/2, r, n.x, alternative = "greater", 
                pprob = TRUE), less = pperm(STATISTIC + n.x * 
                (n.x + 1)/2, r, n.x, alternative = "less", pprob = TRUE))
            MIDP <- PVAL$PPROB
            PVAL <- PVAL$PVALUE
            if (conf.int) {
                if (mu != 0) 
                  r <- rank(c(x, y))
                alpha <- 1 - conf.level
                diffs <- sort(outer(x, y, "-"))
                if (TIES) {
                  fs <- function(d) sum(rank(c(x - d, y))[seq(along = x)]) - 
                    n.x * (n.x + 1)/2
                  w <- sapply(diffs, fs)
                }
                else {
                  w <- (n.x * n.y):1
                }
                cint <- switch(alternative, two.sided = {
                  qu <- qperm(alpha/2, r, n.x) - n.x * (n.x + 
                    1)/2
                  ql <- qperm(1 - alpha/2, r, n.x) - n.x * (n.x + 
                    1)/2
                  if (qu <= min(w)) lci <- max(diffs) else lci <- min(diffs[w <= 
                    qu])
                  if (ql >= max(w)) uci <- min(diffs) else uci <- max(diffs[w > 
                    ql])
                  c(uci, lci)
                }, greater = {
                  ql <- qperm(1 - alpha, r, n.x) - n.x * (n.x + 
                    1)/2
                  if (ql >= max(w)) uci <- min(diffs) else uci <- max(diffs[w > 
                    ql])
                  c(uci, +Inf)
                }, less = {
                  qu <- qperm(alpha, r, n.x) - n.x * (n.x + 1)/2
                  if (qu <= min(w)) lci <- max(diffs) else lci <- min(diffs[w <= 
                    qu])
                  c(-Inf, lci)
                })
                attr(cint, "conf.level") <- conf.level
                wmean <- n.x/(n.x + n.y) * sum(r) - n.x * (n.x + 
                  1)/2
                ESTIMATE <- mean(c(min(diffs[w <= ceiling(wmean)]), 
                  max(diffs[w > wmean])))
                names(ESTIMATE) <- "difference in location"
            }
        }
        else {
            METHOD <- "Asymptotic Wilcoxon rank sum test"
            N <- n.x + n.y
            wmean <- n.x/N * sum(r)
            wvar <- n.x * n.y/(N * (N - 1)) * sum((r - wmean/n.x)^2)
            PVAL <- pnorm((STATISTIC + n.x * (n.x + 1)/2 - wmean)/sqrt(wvar))
            if (alternative == "greater") 
                PVAL <- 1 - PVAL
            if (alternative == "two.sided") 
                PVAL <- 2 * min(PVAL, 1 - PVAL)
            if (conf.int) {
                alpha <- 1 - conf.level
                mumin <- min(x) - max(y)
                mumax <- max(x) - min(y)
                CORRECTION.CI <- 0
                wdiff <- function(d, zq) {
                  dr <- rank(c(x - d, y))
                  dz <- sum(dr[seq(along = x)])
                  dz <- (dz - wmean)/sqrt(wvar)
                  dz - zq
                }
                cint <- switch(alternative, two.sided = {
                  l <- uniroot(wdiff, c(mumin, mumax), tol = 1e-04, 
                    zq = qnorm(alpha/2, lower.tail = FALSE))$root
                  u <- uniroot(wdiff, c(mumin, mumax), tol = 1e-04, 
                    zq = qnorm(alpha/2))$root
                  c(l, u)
                }, greater = {
                  l <- uniroot(wdiff, c(mumin, mumax), tol = 1e-04, 
                    zq = qnorm(alpha, lower.tail = FALSE))$root
                  c(l, +Inf)
                }, less = {
                  u <- uniroot(wdiff, c(mumin, mumax), tol = 1e-04, 
                    zq = qnorm(alpha))$root
                  c(-Inf, u)
                })
                attr(cint, "conf.level") <- conf.level
                ESTIMATE <- uniroot(wdiff, c(mumin, mumax), tol = 1e-04, 
                  zq = 0)$root
                names(ESTIMATE) <- "difference in location"
            }
        }
    }
    if (!is.null(MIDP)) {
        names(MIDP) <- "point prob"
        RVAL <- list(statistic = STATISTIC, pointprob = MIDP, 
            p.value = PVAL, null.value = c(mu = mu), alternative = alternative, 
            method = METHOD, data.name = DNAME)
    }
    else {
        RVAL <- list(statistic = STATISTIC, p.value = PVAL, null.value = c(mu = mu), 
            alternative = alternative, method = METHOD, data.name = DNAME)
    }
    if (conf.int) {
        RVAL$conf.int <- cint
        RVAL$estimate <- ESTIMATE
    }
    class(RVAL) <- "htest"
    return(RVAL)
}

wilcox.exact.formula2<-function (formula, data, subset, na.action, ...) 
{
    if (missing(formula) || (length(formula) != 3) || (length(attr(terms(formula[-2]), 
        "term.labels")) != 1) || (length(attr(terms(formula[-3]), 
        "term.labels")) != 1)) 
        stop("formula missing or incorrect")
    if (missing(na.action)) 
        na.action <- getOption("na.action")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if (nlevels(g) != 2) 
        stop("grouping factor must have exactly 2 levels")
    DATA <- split(mf[[response]], g)
    names(DATA) <- c("x", "y")
    y <- do.call("wilcox.exact2", c(DATA, list(...)))
    y$data.name <- DNAME
    y
}