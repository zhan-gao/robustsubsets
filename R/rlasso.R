#' Huber regression with L1 regularization
#' information criteria for tuning parameters
#' @param x
#' @param y
#' @param huber_thres huber threshold parameter; if null, a sequence will be generated;
#' @param var_sel_lasso_lambda Lasso tuning parameter
#' @param ic_select TRUE: use information criteria to choose tuning parameters; FALSE: use the given tuning parameters
#' @param nlambda number of tuning parameters for Lasso
#' @param lambda.min minimum value of tuning parameter for Lasso
#' @param npsi number of tuning parameters for Huber threshold
#' @param psi_min_quantile quantile for minimum Huber threshold
#' @param psi_min_ratio ratio for minimum Huber threshold
#' @return estimated coefficient as a list
#'  \item{b}{estimated coefficient (incldue intercept))
#'  \item{a}{estimated outlier shifter}
#'  \item{d}{inlier indicator}
#'
#' @import hqreg
#'
#' @export

rlasso <- function(x, y, huber_thres = NULL, var_sel_lasso_lambda = NULL, ic_select = TRUE, nlambda = 100, lambda.min = 0.0001, npsi = 100, psi_min_quantile = 0.75, psi_min_ratio = 0.01) {

    if (is.null(huber_thres)) {
        if (ic_select) {
            huber_thres <- get_psi_seq(x, y, nlambda = npsi, lambda_min_quantile = psi_min_quantile, lambda_min_ratio = psi_min_ratio)
        } else {
            huber_thres = 1.345
        }
    }
    if (is.null(var_sel_lasso_lambda)) {
        var_sel_lasso_lambda <- get_lambda_seq(x, y, nlambda = nlambda, lambda_min_ratio = lambda.min)
    }

    if (ic_select) {
        n <- nrow(x)
        p <- ncol(x)
        npsi <- length(huber_thres)
        nlambda <- length(var_sel_lasso_lambda)

        ic_mat <- matrix(NA, nrow = nlambda, ncol = npsi)

        for (i in seq_along(huber_thres)) {
            thres <- huber_thres[i]
            huber_fit <- huber_reg(x, y, huber_thres = thres, var_sel_lasso_lambda = var_sel_lasso_lambda)
            objval <- colSums((y - cbind(1, x) %*% huber_fit$b - huber_fit$a)^2) / 2
            df_vec <- colSums(huber_fit$b != 0) + colSums(huber_fit$a != 0) + 1
            ic_value <- n * log(objval / n) + (df_vec) * log(n)

            ic_mat[, i] <- ic_value
        }

        pos <- which(ic_mat == min(ic_mat), arr.ind = TRUE)
        lambda_hat <- var_sel_lasso_lambda[pos[1]]
        psi_hat <- huber_thres[pos[2]]

        huber_fit_hat <- huber_reg(x, y, huber_thres = psi_hat, var_sel_lasso_lambda = lambda_hat)
        return(c(huber_fit_hat, list(lambda = lambda_hat, huber_thres = psi_hat)))
    } else {
        if (length(huber_thres) != 1) {
            stop("Please provide a single tuning parameter for Huber threshold.")
        } else {
            huber_fit <- huber_reg(x, y, huber_thres = huber_thres, var_sel_lasso_lambda = var_sel_lasso_lambda)
            return(huber_fit)
        }
    }
}
#' Huber regression
#' Direct implementation of Huber regression via CVXR solver
#' Support no intercept case
#'
#' @param x
#' @param y
#' @param intercept
#' @param huber_thres
#'
#' @return estimated coefficient
#'
#' @import CVXR
#'
#' @export

huber_reg_cvxr <- function(x, y, intercept = TRUE, huber_thres = 1.345) {
    # huber_thres = 1.345 suggested by Huber (1981)

    if (intercept) {
        x <- cbind(1, x)
    }

    m <- ncol(x)
    beta <- Variable(m)
    obj <- sum(CVXR::huber(y - x %*% beta, huber_thres))
    problem <- Problem(Minimize(obj))
    prob_data <- get_problem_data(problem, solver = "ECOS")
    ECOS_dims <- ECOS.dims_to_solver_dict(prob_data$data[["dims"]])
    solver_output <- ECOSolveR::ECOS_csolve(
        c = prob_data$data[["c"]],
        G = prob_data$data[["G"]],
        h = prob_data$data[["h"]],
        dims = ECOS_dims,
        A = prob_data$data[["A"]],
        b = prob_data$data[["b"]]
    )
    result <- unpack_results(
        problem, solver_output,
        prob_data$chain, prob_data$inverse_data
    )

    b_hat <- c(result$getValue(beta))
    e_hat <- y - x %*% b_hat
    d_hat <- as.numeric(abs(e_hat) < huber_thres)

    a_hat <- rep(NA, length(y))
    a_hat[abs(e_hat) < huber_thres] <- 0
    a_hat[abs(e_hat) >= huber_thres] <- e_hat[abs(e_hat) >= huber_thres]
    # a_hat[e_hat >= huber_thres] <- e_hat[e_hat >= huber_thres] - huber_thres
    # a_hat[e_hat <= -huber_thres] <- e_hat[e_hat <= -huber_thres] + huber_thres

    return(list(b = b_hat,
                a = a_hat,
                d = d_hat))
}

#' Huber regression based on Yi and Huang (2016)
#'
#' Yi, C. and Huang, J. (2016) Semismooth Newton Coordinate Descent Algorithm
#' for Elastic-Net Penalized Huber Loss Regression and Quantile Regression
#'
#'
#' @param x
#' @param y
#' @param huber_thres
#' @param var_sel_lasso_lambda The tuning parameter for variable selection.
#'      = 0: no regularization on beta; = NULL return results with a sequence of lambda; = a number: use the number as lambda
#'
#' @return estimated coefficient
#'
#' @import hqreg
#'
#' @export
#'
huber_reg <- function(x, y, huber_thres = 1.345, var_sel_lasso_lambda = 0, nlambda = 100, lambda.min = 0.005) {

    n <- nrow(x)
    if (length(var_sel_lasso_lambda) == 1) {
        huber_fit <- hqreg::hqreg(x, y, lambda = c(1, var_sel_lasso_lambda), gamma = huber_thres)
        b_hat <- huber_fit$beta[, 2]
    } else if (length(var_sel_lasso_lambda) == 0) {
        huber_fit <- hqreg::hqreg(x, y, nlambda = nlambda, lambda.min = lambda.min, gamma = huber_thres)
        b_hat <- huber_fit$beta
    } else {
        huber_fit <- hqreg::hqreg(x, y, lambda = var_sel_lasso_lambda, gamma = huber_thres)
        b_hat <- huber_fit$beta
    }
    e_hat <- y - cbind(1, x) %*% b_hat
    d_hat <- abs(e_hat) < huber_thres

    a_hat <- matrix(NA, n, dim(e_hat)[2])
    a_hat[abs(e_hat) < huber_thres] <- 0
    a_hat[abs(e_hat) >= huber_thres] <- e_hat[abs(e_hat) >= huber_thres]


    return(list(b = b_hat,
                a = a_hat,
                d = d_hat))
}
#' The function generate candidate tuning parameters for L1 method (Huber threshold)
#'
#'
#'
get_psi_seq <- function(x, y, intercept = TRUE,
                           nlambda = 100, lambda_min_quantile = 0.75, lambda_min_ratio = 0.01) {
    if (intercept) {
        x <- cbind(1, x)
    }

    b_ols <- solve(t(x) %*% x, t(x) %*% y)
    lambda_max <- max(abs(y - x %*% b_ols))
    if (!is.null(lambda_min_quantile)) lambda_min <- quantile(abs(y - x %*% b_ols), lambda_min_quantile)
    else if (!is.null(lambda_min_ratio)) lambda_min <- lambda_min_ratio * lambda_max
    else stop("No lambda_min available.")
    lambda_seq <- exp(
        seq(log(lambda_max), log(lambda_min), length.out = nlambda)
    )

    return(lambda_seq)

}

sd_n <- function(x){
    sqrt(sum((x - mean(x))^2)/length(x))
}

get_lambda_seq <- function(x, y, scalex = TRUE, nlambda = 100, lambda_min_ratio = 0.0001){

    n <- nrow(x)
    p <- ncol(x)
    y <- as.numeric(y)

    if (scalex) {
        sx <- scale(x, scale = apply(x, 2, sd_n))
        sx <- as.matrix(sx, ncol = p, nrow = n)

        lambda_max <- max(abs(colSums(sx * y))) / n
    } else {
        lambda_max <- max(abs(colSums(x * y))) / n
    }

    lambda_min <- lambda_min_ratio * lambda_max
    lambda_seq <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
    return(lambda_seq)
}

#' LAD
#'
#' @param x
#' @param y
#' @param intercept
#' @param tau quantile
#'
#' @return estimated coefficient
#'
#' @import hqreg
#' @export
#'
lad_reg <- function(x, y, intercept = TRUE, tau = 0.5, lambda = 0) {
    if (intercept) {
        if (length(lambda) == 1) {
            lad_fit <- hqreg::hqreg(
                x,
                y,
                method = "quantile",
                tau = tau,
                lambda = c(1, lambda)
            )
            b_hat <- lad_fit$beta[, 2]
        } else if (length(lambda) == 0) {
            lambda_hat <- hqreg::cv.hqreg(x, y, method = "quantile", tau = tau)$lambda.min
            lad_fit <- hqreg::hqreg(
                x,
                y,
                method = "quantile",
                tau = tau,
                lambda = c(1, lambda_hat)
            )
            b_hat <- lad_fit$beta[, 2]
        } else {
            lad_fit <- hqreg::hqreg(x,
                                    y,
                                    method = "quantile",
                                    tau = tau,
                                    lambda = lambda_hat)
            b_hat <- lad_fit$beta
        }
    } else {
        return(quantreg::rq.fit(x, y, tau)$coefficients)
    }

    return(b_hat)
}


