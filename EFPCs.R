source("utils.R")
library(CVXR)
library(fda)
library(fda.usc)
library(e1071)
set.seed(2) 
#Setting enviroment
n_basis <- 11
n_curves <- 100 #number of curves per class 
n_testing <- 20 #number of curves per class for testing

#Number of points per curve
min_number_of_points <- 25
max_number_of_points <- 50


#Domain of function 
lower_bound <- 0
upper_bound <- 1

#Point to point standard deviation in random curves
sd= 0.2

#Define the inner product that will be used
inner_product <- L2_inner_product

#Create basis
basis <- create.fourier.basis(rangeval = c(lower_bound, upper_bound), nbasis = n_basis)

#Slackness coefficient for SVMe
C = 1  



#Define the underlying functions for the two classes we want to classify
# base_class_1 <- function(t) { 1-t^2 }
# base_class_2 <- function(t) { cos((pi/2) * t) }




runExperiment <- function(R, n_curves, n_testing, min_number_of_points, max_number_of_points, lower_bound, upper_bound, sd, classes, inner_product, basis) {

    results <- data.frame(
        repetition = 1:R,
        training_accuracy = NA_real_,
        testing_accuracy = NA_real_,
        b = NA_real_,
        has_free_sv = NA,
        min_eig = NA_real_,
        p = NA_integer_
    )
    # Complex objects: list columns
    results$w <- vector("list", R)
    results$lambda <- vector("list", R)
    results$smoothed_dataset <- vector("list", R)
    
    for (i in 1:R) {
        #Generate synthetic data
        dataset <- generate_synthetic_data(n_curves, classes, min_number_of_points, max_number_of_points, lower_bound, upper_bound, sd)
        smoothed_dataset <- smooth_dataset(dataset)$smoothed_dataset

        #Take the first n_curves of each class for testing 
        testing_data <- c(smoothed_dataset[(n_curves-n_testing+1):n_curves],smoothed_dataset[(2*n_curves-n_testing+1):(2*n_curves)])#class 1 & 2 curves


        #Take the rest for training 
        training_data <- c(smoothed_dataset[1:(n_curves-n_testing)],smoothed_dataset[(n_curves+1):(2*n_curves-n_testing)])

        EFPCs <- get_EFPCs(training_data, variance_explained = 0.85)
        results$p[i] <- EFPCs$p

        training_data_projected <- project_to_EFPCs(training_data, EFPCs, inner_product)
        testing_data_projected <- project_to_EFPCs(testing_data, EFPCs, inner_product)

        #Train the model
        model <- train(training_data_projected, C, inner_product)

        #Test the model
        predicted_labels <- test_model(testing_data_projected, inner_product, model)
        testing_accuracy <- calculate_accuracy(testing_data_projected, predicted_labels)
        cat("ITERATION", i, "\n")

        predicted_labels <- test_model(training_data_projected, inner_product, model)
        training_accuracy <- calculate_accuracy(training_data_projected, predicted_labels)

        results$training_accuracy[i] <- training_accuracy
        results$testing_accuracy[i] <- testing_accuracy
        results$w[[i]] <- model$w
        results$b[[i]] <- model$b
        results$has_free_sv[[i]] <- model$has_free_sv
        results$min_eig[[i]] <- model$min_eig
        results$lambda[[i]] <- model$lambda
        results$smoothed_dataset[[i]] <- smoothed_dataset
        #plot_smoothed_curves_with_truth(
        #  smoothed_dataset,
        #  base_class_1,
        #  base_class_2
        #)
    }
    return(results)

}
base_class_1 <- function(t) {
  sin(2 * pi * t)
}

base_class_2 <- function(t) {
  sin(4 * pi * t)
}

classes <- list(base_class_1, base_class_2)

R <- 100 #NUMBER OF REPETITIONS




result_base <- runExperiment(R, n_curves, n_testing, min_number_of_points, max_number_of_points, lower_bound, upper_bound, sd, classes, inner_product, basis)

a_values <- c(1, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30)
      
delta_phase_values <- c(0,0.05,0.10,0.15,0.20,0.25,0.30)

delta_local_values <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2)

mu_values <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

sigma_values <- c(0.02,0.04,0.06,0.08,0.10,0.15,0.20,0.25,0.30)


parameters <- a_values
experiment_amplitude <- vector("list", length(parameters))
for (k in seq_along(parameters)) {
    i  <- parameters[k]
    base_class_1 <- function(t) {
        sin(2 * pi * t)
    }
    base_class_2 <- function(t) {
        i * sin(2 * pi * t)
    }
    classes <- list(base_class_1, base_class_2)
    result <- runExperiment(R, n_curves, n_testing, min_number_of_points, max_number_of_points, lower_bound, upper_bound, sd, classes, inner_product, basis)
    cat("PARAMETER", i, "\n")
    
    #Save experiment results
    experiment_amplitude[[k]] <- list(
        parameter = i,
        results = result
    )
}

#Phase variation
parameters <- delta_phase_values
experiment_phase <- vector("list", length(parameters))
for (k in seq_along(parameters)) {
    i  <- parameters[k]
    base_class_1 <- function(t) {
        sin(2 * pi * t)
    }
    base_class_2 <- function(t) {
        sin(2 * pi * t+i)
    }
    classes <- list(base_class_1, base_class_2)
    result <- runExperiment(R, n_curves, n_testing, min_number_of_points, max_number_of_points, lower_bound, upper_bound, sd, classes, inner_product, basis)

    #Save experiment results
    experiment_phase[[k]] <- list(
        parameter = i,
        results = result
    )
}

#Delta  variation 
sigma <- 0.05
mu <- 0.5
parameters <- delta_local_values
experiment_delta <- vector("list", length(parameters))
for (k in seq_along(parameters)) {
    i <- parameters[k]
    base_class_1 <- function(t) {
        sin(2 * pi * t)
    }
    base_class_2 <- function(t) {
        sin(2 * pi * t) + (i * exp(-((t - mu)^2) / (2 * sigma^2)))
    }
    classes <- list(base_class_1, base_class_2)
    result <- runExperiment(R, n_curves, n_testing, min_number_of_points, max_number_of_points, lower_bound, upper_bound, sd, classes, inner_product, basis)

    #Save experiment results
    experiment_delta[[k]] <- list(
        parameter = i,
        results = result
    )

}

#Sigma variation
delta <- 0.3
mu <- 0.5
parameters <- sigma_values
experiment_sigma <- vector("list", length(parameters))
for (k in seq_along(parameters)) {
    i <- parameters[k]
    base_class_1 <- function(t) {
        sin(2 * pi * t)
    }
    base_class_2 <- function(t) {
        sin(2 * pi * t) + (delta * exp(-((t - mu)^2) / (2 * i^2)))
    }
    classes <- list(base_class_1, base_class_2)
    result <- runExperiment(R, n_curves, n_testing, min_number_of_points, max_number_of_points, lower_bound, upper_bound, sd, classes, inner_product, basis)

    #Save experiment results
    experiment_sigma[[k]] <- list(
        parameter = i,
        results = result
    )

}

#Mu variation
delta <- 0.3
sigma <- 0.05
parameters <- mu_values
experiment_mu <- vector("list", length(parameters))
for (k in seq_along(parameters)) {
    i <- parameters[k]
    base_class_1 <- function(t) {
        sin(2 * pi * t)
    }
    base_class_2 <- function(t) {
        sin(2 * pi * t) + (delta * exp(-((t - i)^2) / (2 * sigma^2)))
    }
    classes <- list(base_class_1, base_class_2)
    result <- runExperiment(R, n_curves, n_testing, min_number_of_points, max_number_of_points, lower_bound, upper_bound, sd, classes, inner_product, basis)

    #Save experiment results
    experiment_mu[[k]] <- list(
        parameter = i,
        results = result
    )
}


# ============================================================
# SUMMARY OF EXPERIMENT RESULTS
# ============================================================

# Summarise one runExperiment result
summarise_result <- function(results) {

    data.frame(
        training_mean = mean(results$training_accuracy),
        training_sd   = stats::sd(results$training_accuracy),
        testing_mean  = mean(results$testing_accuracy),
        testing_sd    = stats::sd(results$testing_accuracy),
        p = NA_integer_
    )
}


# Summarise an experiment in which a parameter was varied
summarise_experiment <- function(experiment) {
    do.call(
        rbind,
        lapply(experiment, function(x) {

            data.frame(
                parameter = x$parameter,
                training_mean = mean(x$results$training_accuracy),
                training_sd   = sd(x$results$training_accuracy),
                testing_mean  = mean(x$results$testing_accuracy),
                testing_sd    = sd(x$results$testing_accuracy),
                p_mean        = mean(x$results$p),
                p_sd          = sd(x$results$p),
                p_distribution = paste(
                    names(table(x$results$p)),
                    as.integer(table(x$results$p)),
                    sep = ":",
                    collapse = ", "
                )
            )
        })
    )
}

# Basic experiment
summary_base <- summarise_result(result_base)

# Parameter experiments
summary_amplitude <- summarise_experiment(experiment_amplitude)
summary_phase     <- summarise_experiment(experiment_phase)
summary_delta     <- summarise_experiment(experiment_delta)
summary_sigma     <- summarise_experiment(experiment_sigma)
summary_mu        <- summarise_experiment(experiment_mu)


# Print summaries
cat("\n===== BASE EXPERIMENT =====\n")
print(summary_base)

cat("\n===== AMPLITUDE =====\n")
print(summary_amplitude)

cat("\n===== PHASE =====\n")
print(summary_phase)

cat("\n===== LOCAL CHANGE: DELTA =====\n")
print(summary_delta)

cat("\n===== LOCAL CHANGE: SIGMA =====\n")
print(summary_sigma)

cat("\n===== LOCAL CHANGE: MU =====\n")
print(summary_mu)


save_experiment_results <- function(
    result_base,
    experiment_amplitude,
    experiment_phase,
    experiment_delta,
    experiment_sigma,
    experiment_mu,
    file = "experiment_results2_EFPCs.rds"
) {

    results <- list(
        result_base = result_base,
        experiment_amplitude = experiment_amplitude,
        experiment_phase = experiment_phase,
        experiment_delta = experiment_delta,
        experiment_sigma = experiment_sigma,
        experiment_mu = experiment_mu
    )

    saveRDS(results, file = file)

    cat("Results saved in:", file, "\n")
}

save_experiment_results(
    result_base,
    experiment_amplitude,
    experiment_phase,
    experiment_delta,
    experiment_sigma,
    experiment_mu
)
