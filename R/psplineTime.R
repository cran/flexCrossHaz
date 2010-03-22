psplineTime <-
function (x, df = 4, theta, nterm = 2.5 * df, by.x, pen.diff=2, degree = 3, eps = 0.1, v, crs.point, Bs=FALSE, method, ...){

    if (!missing(theta)) {
        method <- "fixed"
        if (theta <= 0 || theta >= 1)
            stop("Invalid value for theta")
    }
    else if (df == 0 || (!missing(method) && method == "aic")) {
        method <- "aic"
        nterm <- 15
        if (missing(eps))
            eps <- 1e-05
    }
    else {
        method <- "df"
        if (df <= 1)
            stop("Too few degrees of freedom")
        if (df + 1 > nterm)
            stop("`nterm' too small for df=", df)
    }
    xname <- deparse(substitute(by.x))
    keepx <- !is.na(x)
    rx <- range(x[keepx])
    nterm <- round(nterm)
    if (nterm < 3)
        stop("Too few basis functions")
    dx <- (rx[2] - rx[1])/nterm
    knots <- c(rx[1] + dx * ((-degree):(nterm - 1)), rx[2] +
        dx * (0:degree))
    if (all(keepx)) 
        newx <- if(Bs) bs(x,df = nterm, degree= degree, intercept=TRUE) else spline.des(knots, x, degree + 1)$design
    else {
        temp <- if(Bs) bs(x[keepx], df = nterm, degree = degree, intercept=TRUE) else spline.des(knots, x[keepx], degree + 1)$design
        newx <- matrix(NA, length(x), ncol(temp))
        newx[keepx, ] <- temp
    }
        newx <- newx*by.x
    class(newx) <- "coxph.penalty"
    nvar<-ncol(newx)
    dmat<-diff(diag(nvar),diff=pen.diff)
    dmat<-crossprod(dmat)
    if(!missing(v)) dmat<-dmat+diag(I(v<0))*10^{6}
    if(!missing(crs.point)){
        B0<-spline.des(knots, c(min(x),crs.point,max(x)), degree + 1)$design
        dmat<-dmat+tcrossprod(B0[2,])*10^{6}
        }
    xnames <- paste("ps(", xname, ")", 1:nvar, sep = "")
    pfun <- function(coef, theta, n, dmat) {
        if (theta >= 1)
            list(penalty = 100 * (1 - theta), flag = TRUE)
        else {
            if (theta <= 0)
                lambda <- 0
            else lambda <- theta/(1 - theta)
            list(penalty = c(coef %*% dmat %*% coef) * lambda/2,
                first = c(dmat %*% coef) * lambda, second = c(dmat *
                  lambda), flag = FALSE)
        }
    }
    printfun <- function(coef, var, var2, df, history) {
        test1 <- coxph.wtest(var, coef)$test
        xmat <- cbind(1, cbase)
        xsig <- coxph.wtest(var, xmat)$solve
        cmat <- coxph.wtest(t(xmat) %*% xsig, t(xsig))$solve[2,]
        linear <- sum(cmat * coef)
        lvar1 <- c(cmat %*% var %*% cmat)
        lvar2 <- c(cmat %*% var2 %*% cmat)
        test2 <- linear^2/lvar1
        cmat <- rbind(c(linear, sqrt(lvar1), sqrt(lvar2), test2,
            1, 1 - pchisq(test2, 1)), c(NA, NA, NA, test1 - test2,
            df - 1, 1 - pchisq(test1 - test2, max(0.5, df - 1))))
        dimnames(cmat) <- list(c("linear", "nonlin"), NULL)
        nn <- nrow(history$thetas)
        if (length(nn))
            theta <- history$thetas[nn, 1]
        else theta <- history$theta
        list(coef = cmat, history = paste("Theta=", format(theta)))
    }
    cbase <- knots[2:nvar] + (rx[1] - knots[1])
    if (method == "fixed") {
        temp <- list(pfun = pfun, printfun = printfun, pparm = dmat,
            diag = FALSE, cparm = list(theta = theta), varname = xnames,
            cfun = function(parms, iter, old) list(theta = parms$theta,
                done = TRUE))
    }
    else if (method == "df") {
        frailty.controldf<-NULL    #the last
        temp <- list(pfun = pfun, printfun = printfun, diag = FALSE,
            cargs = ("df"), cparm = list(df = df, eps = eps,
                thetas = c(1, 0), dfs = c(1, nterm), guess = 1 -
                  df/nterm, ...), pparm = dmat, varname = xnames,cfun = frailty.controldf)
    }
    else {
    frailty.controlaic<-NULL   #the last
        temp <- list(pfun = pfun, printfun = printfun, pparm = dmat,
            diag = FALSE, cargs = c("neff", "df", "plik"), cparm = list(eps = eps,
                init = c(0.5, 0.95), lower = 0, upper = 1, ...),
            varname = xnames, cfun = frailty.controlaic)
    }
    attributes(newx) <- c(attributes(newx), temp)
    newx
}

