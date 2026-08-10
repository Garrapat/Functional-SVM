#' Generate a synthetic dataset of noisy curves
#' @param n_curves Number of curves per class
#' @param classes List of functions
#' @param minimum number of points per data entry
#' @param maximum number of points per data entry
#' @return A list of data frames
generate_synthetic_data <- function(n_curves, classes, min_points=25, max_points=50, lower_bound = 0, upper_bound = 1, sd) {
    dataset <- list()
    counter <- 1
    for (j in 1:length(classes)) {
        for (i in 1:n_curves) {
            n_points <- sample(min_points:max_points, 1)
            t <- sort(runif(n_points, lower_bound, upper_bound)) 
            #t <- sort(stratified_runif(n_points, lower_bound, upper_bound)) 
            noise <- rnorm(n_points, mean = 0, sd = sd) # Add noise to the measurements # nolint
            measurement <- classes[[j]](t) + noise
            label = 1
            if (j==1) {
            label = -1 
            }
            dataset[[counter]] <- data.frame(t = t, measurement = measurement, class = j, label = label)
            counter <- counter + 1
        }
    }
  return(dataset)
}


#' Smooth a curve using GCV
#' @param set of points
#' @param the measurements at point t
#' @param basis Basis object from fda
#' @param lambda_values Numeric vector of lambdas in which we want to find best curve
#' @return A list with the best smoothed curve, it's lambda value, log lambda value, and GCV value
smooth_curve_gcv <- function(t, measurement, basis, Lfdobj=2, lambda_values) {
  
  gcv_values <- numeric(length(lambda_values))
  smoothed_curves <- vector("list", length(lambda_values))
    
    for (i in seq_along(lambda_values)) {
      fdparobj <- fdPar(fdobj = basis, Lfdobj = Lfdobj, lambda = lambda_values[i])
      smoothed_curves[[i]] <- smooth.basis(t, measurement, fdparobj)
      gcv_values[i] <- smoothed_curves[[i]]$gcv
    }

    best_i <- which.min(gcv_values)

    return(list(
      best_curve = smoothed_curves[[best_i]],
      best_lambda = lambda_values[best_i],
      best_loglambda = log10(lambda_values[best_i]),
      gcv_values = gcv_values
  ))
}


#' Smooth a dataset using GCV
#' @param dataset of curves as points to be smoothed
#' @param basis object from fda to use as basis
#' @param lambda_values numeric vector of lambdas to test as PSS coefficients
#' @return A list containing the best smoothed curves, their lambda values, and the lambda value base 10 log 
smooth_dataset <- function(dataset, lambda_values = 10^seq(-7, 3, by = 1)) {
    #We now smooth each curve using the best lambda for each curve.
    smoothed_dataset <- vector("list", length(dataset))
    lambda_per_curve <- numeric(length(dataset))
    loglambda_per_curve <- numeric(length(dataset))

    for (i in seq_along(dataset)) {
    smoothed_data <- smooth_curve_gcv(dataset[[i]]$t, dataset[[i]]$measurement, basis, Lfdobj=2, lambda_values = lambda_values)

    smoothed_dataset[[i]] <- list(curve = smoothed_data$best_curve, label= dataset[[i]]$label[[1]])
    lambda_per_curve[i] <- smoothed_data$best_lambda
    loglambda_per_curve[i] <- smoothed_data$best_loglambda
    }
    return(list(smoothed_dataset=smoothed_dataset,lambda_per_curve=lambda_per_curve,loglambda_per_curve=loglambda_per_curve))
}

#Different innerproducts

L2_inner_product <- function(fd1, fd2, n_points = 1000, a=0,b=1) {
  integrand <- function(t) eval.fd(t, fd1) * eval.fd(t, fd2)
  return(integrate(integrand, lower = a, upper = b)$value)
}

L2_inner_product_scratch <- function(fd1, fd2, n_points = 1000, a=0,b=1) {
  t_grid <- seq(a,b, length.out = n_points)
  f1 <- eval.fd(t_grid, fd1)
  f2 <- eval.fd(t_grid, fd2)
  return(sum(diff(t_grid) * (head(f1*f2, -1) + tail(f1*f2, -1)) / 2))
}


L2_inprod <- function(fd1, fd2) {
  return(inprod(fd1, fd2))
}


#' Calculate the matrix of all inner product combinations (Gram matrix)
#' @param Dataset for which to find Gram matrix
#' @param Inner product to use for calculation
#' @return Gram matrix
gram_matrix <- function(dataset, inner_product) {
  n <- length(dataset)
  K <- matrix(0, nrow = n, ncol = n)
  
  for (i in seq_len(n)) {
    for (j in i:n) {
      value <- inner_product(dataset[[i]]$fd, dataset[[j]]$fd)
      K[i, j] <- value
      K[j, i] <- value
    }
  }
  
  return(K)
}



#TODO: SEPARATE B CALCULATION AND CLEAN UP. SPECIFY INPUTS


#' Train the model and return the w and b, parameters
#' @param Dataset which will be used for training
#' @param Penalty coefficient C
#' @param Inner product
#' @param Tolerance for lambda values default is 1e-3 used to take into account numerical error
#' @return w and b, the model parameters for the SVM classifier 
train <- function(dataset, C, inner_product, eps= 1e-3) {
    labels <- sapply(dataset,function(obj) obj$label) #List of labels from dataset
    curves <- lapply(dataset, function(obj) obj$curve) #List of curves from dataset
    basis <- curves[[1]]$fd$basis#Basis used across dataset
    

    K <- gram_matrix(curves, inner_product)
    
    Q <- K * (labels %*% t(labels))


    #Guarantee that Q is positive semidefinete to ensure convexity
    eig_Q <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
    min_eig <- min(eig_Q)
    cat(" min eig ",min_eig)
    if (min_eig < 0) {
        Q <- Q + diag(abs(min_eig) + 1e-8, nrow(Q))
    }


    lambda <- Variable(length(dataset))

    objective <- Maximize(sum(lambda)-(1/2)*quad_form(lambda, Q))
    


    constraints <- list(
        lambda >= 0,
        lambda <= C,
        sum(lambda * labels) == 0
    )

    problem <- Problem(objective,constraints)

    result <- solve(problem)
    
    lambda_value <- as.numeric(result$getValue(lambda))

    #Recover w
    w_fd <- calculate_w(basis, curves, labels, lambda_value)

    #Find free support vectors
    b <- calculate_b(lambda_value, curves, labels, inner_product, w_fd, C)

    
    
    
    scores_from_w <- sapply(seq_along(curves), function(i) {
        inner_product(w_fd, curves[[i]]$fd) + b
    })

    scores_from_K <- as.vector((lambda_value * labels) %*% K) + b

    cat(" error",max(abs(scores_from_w - scores_from_K)))


    classifier <- function(curve) {
        score <- inner_product(w_fd, curve) + b
        ifelse(score < 0, -1, 1)
    }
    return(list(w = w_fd, b = b, classifier = classifier))

}



#' Calculates w parameter
#' @param Basis used to represent the functions
#' @param Curves used for training
#' @param Inner product
#' @param Labels that correspond to the curves
#' @param Lambda value ie solution to the dual problem
#' @return w parameter 
calculate_w <- function(basis, curves, labels, lambda_value) {
    #Recover w coefficients 
    w_coef <- numeric(basis$nbasis)
    for (i in seq_len(basis$nbasis)) {
        ith_coefficients <- sapply(curves, function(obj) obj$fd$coefs[i])
    
        w_coef[i] = sum(lambda_value * labels * ith_coefficients)
    }


    #Turn w into curve
    w_fd <- fd(coef = matrix(w_coef,ncol=1), basisobj= basis)
    return(w_fd)
}



#' Calculates the b parameter
#' @param Lambda value ie solution to the dual problem
#' @param curves
#' @param labels
#' @param Inner product
#' @param w_fd parameter ie function
#' @param Penalty coefficient C
#' @param Tolerance for lambda values default is 1e-3 used to take into account numerical error
#' @return b parameter 
calculate_b <- function(lambda_value, curves, labels, inner_product, w_fd, C, eps= 1e-3) {

    lambda_star_indices <- which(lambda_value>0+eps & lambda_value < C-eps) #Select lambdas we use eps due to numerical error
        
    #Checks if there are free support vectors
    
    if(length(lambda_star_indices) == 0){
        #If they do not exist then we take the average b values
        upper_indices <- which((C - eps < lambda_value & labels == 1) | (lambda_value < 0 + eps & labels == -1))
        b_upper <- sapply(upper_indices, function(i) {1/(labels[i])-inner_product(w_fd,curves[[i]]$fd)})

        lower_indices <- which((C - eps < lambda_value & labels == -1) | (lambda_value < 0 + eps & labels == 1))
        b_lower <- sapply(lower_indices, function(i) {(labels[i])-inner_product(w_fd,curves[[i]]$fd)})

        b = (min(b_upper) + max(b_lower))/2
        print("No free support vectors found, check parameters. Still calculated b but not great")
    }
    #if there are compute normally
    else {   
    b_values <- sapply(lambda_star_indices, function(i){labels[i]-inner_product(w_fd,curves[[i]]$fd)})
    b = mean(b_values)
    }
    return(b)
}

#' Calculates predicted labels for a dataset
#' @param testing data 
#' @param inner_product
#' @return predicted labels
test_model <- function(testing_data, inner_product, model) {
    predicted_labels <- sapply(testing_data, function(obj) {
        score <- inner_product(model$w, obj$curve$fd) + model$b
        ifelse(score<0, -1, 1) 
    })
    return(predicted_labels)
}

#' Calculates accuracy
#' @param testing data 
#' @param prediction
#' @return value between 0, 1 representing accuracy 
calculate_accuracy <- function(testing_data, predicted_labels) {
    testing_labels <- sapply(testing_data, function(obj) obj$label)
    accuracy <- length(which(predicted_labels == testing_labels))/length(testing_data)
    return(accuracy)
}





#kernel attempt
train_ker <- function(dataset, C, ker, eps= 1e-3) {
    labels <- sapply(dataset,function(obj) obj$label) #List of labels from dataset
    curves <- lapply(dataset, function(obj) obj$curve) #List of curves from dataset
    basis <- curves[[1]]$fd$basis#Basis used across dataset
    

    K <- gram_matrix(curves, ker)
    Q <- K * (labels %*% t(labels))

    #Guarantee that Q is positive semidefinete to ensure convexity
    eig_Q <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
    min_eig <- min(eig_Q)
    cat(" min eig ",min_eig)
    if (min_eig < 0) {
        Q <- Q + diag(abs(min_eig) + 1e-8, nrow(Q))
    }


    lambda <- Variable(length(dataset))

    objective <- Maximize(sum(lambda)-(1/2)*quad_form(lambda, Q))
    


    constraints <- list(
        lambda >= 0,
        lambda <= C,
        sum(lambda * labels) == 0
    )

    problem <- Problem(objective,constraints)

    result <- solve(problem)
    
    lambda_value <- as.numeric(result$getValue(lambda))

    #Find free support vectors
    b <- calculate_b(lambda_value, curves, labels, ker, w_fd, C)

    
    
    
    scores_from_w <- sapply(seq_along(curves), function(i) {
        ker(w_fd, curves[[i]]$fd) + b
    })

    scores_from_K <- as.vector((lambda_value * labels) %*% K) + b

    cat(" error",max(abs(scores_from_w - scores_from_K)))


    classifier <- function(curve) {
        score <- ker(w_fd, curve) + b
        ifelse(score < 0, -1, 1)
    }
    return(list(w = w_fd, b = b, classifier = classifier))

}

stratified_runif <- function(n_points, lower_bound = 0, upper_bound = 1) {
    t <- numeric(n_points)
    width <- (upper_bound - lower_bound) / n_points

    for (i in seq_len(n_points)) {
        interval_start <- lower_bound + (i - 1) * width
        interval_end   <- lower_bound + i * width

        t[i] <- stats::runif(
            1,
            min = interval_start,
            max = interval_end
        )
    }

    return(t)
}
