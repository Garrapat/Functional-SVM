library(CVXR)
library(fda)
library(fda.usc)
library(e1071)
source("utils.R")


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

#Define the underlying functions for the two classes we want to classify
# base_class_1 <- function(t) { 1-t^2 }
# base_class_2 <- function(t) { cos((pi/2) * t) }



base_class_1 <- function(t) {
  sin(2 * pi * t)
}

base_class_2 <- function(t) {
  sin(4 * pi * t)
}

classes <- list(base_class_1, base_class_2)

#Create basis
basis <- create.fourier.basis(rangeval = c(lower_bound, upper_bound), nbasis = n_basis)

#Generate synthetic data
dataset <- generate_synthetic_data(n_curves, classes, min_number_of_points, max_number_of_points, lower_bound, upper_bound, sd)
smoothed_dataset <- smooth_dataset(dataset)$smoothed_dataset

#Take the first n_curves of each class for testing 
testing_data <- c(smoothed_dataset[(n_curves-n_testing+1):n_curves],smoothed_dataset[(2*n_curves-n_testing+1):(2*n_curves)])#class 1 & 2 curves

#Take the rest for training 
training_data <- c(smoothed_dataset[1:(n_curves-n_testing)],smoothed_dataset[(n_curves+1):(2*n_curves-n_testing)])

#Train the model
model <- train(training_data,1, inner_product)

#Test the model
predicted_labels <- test_model(testing_data, inner_product, model)
accuracy <- calculate_accuracy(testing_data, predicted_labels)
print(accuracy)
print("TRAINING")
predicted_labels <- test_model(training_data, inner_product, model)
accuracy <- calculate_accuracy(training_data, predicted_labels)
print(accuracy)


#KERNEL IMPLEMENTATION -> just change innerproduct to \langle phi a, phi b \rangle However need to ensure that 




plot_smoothed_curves_with_truth <- function(
    smoothed_dataset,
    base_class_1,
    base_class_2,
    grid = seq(lower_bound, upper_bound, length.out = 500)
) {
  smoothed_values <- lapply(smoothed_dataset, function(item) {
    eval.fd(grid, item$curve$fd)
  })

  all_values <- c(
    unlist(smoothed_values),
    base_class_1(grid),
    base_class_2(grid)
  )

  plot(
    NA,
    xlim = range(grid),
    ylim = range(all_values, finite = TRUE),
    xlab = "t",
    ylab = "Valor de la Función",
    main = "Smoothed curves and underlying class functions"
  )

  for (i in seq_along(smoothed_dataset)) {
    item <- smoothed_dataset[[i]]

    curve_colour <- if (item$label == -1) {
      rgb(1, 0, 0, alpha = 0.2)
    } else {
      rgb(0, 0, 1, alpha = 0.2)
    }

    lines(
      grid,
      smoothed_values[[i]],
      col = curve_colour,
      lwd = 1
    )
  }

  lines(
    grid,
    base_class_1(grid),
    col = "red",
    lwd = 3,
    lty = 2
  )

  lines(
    grid,
    base_class_2(grid),
    col = "blue",
    lwd = 3,
    lty = 2
  )

  legend(
    "bottomleft",
    legend = c(
      "Curvas de clase 1",
      "Curvas de clase 2",
      expression(paste("Curva subyacente ", C[1])),
      expression(paste("Curva subyacente ", C[2]))
    ),
    col = c(
      rgb(1, 0, 0, alpha = 0.5),
      rgb(0, 0, 1, alpha = 0.5),
      "red",
      "blue"
    ),
    lwd = c(1, 1, 3, 3),
    lty = c(1, 1, 2, 2),
    bty = "n"
  )
}

plot_smoothed_curves_with_truth(
  smoothed_dataset,
  base_class_1,
  base_class_2
)


plot_single_approximation <- function(
    index,
    dataset,
    smoothed_dataset,
    base_class_1,
    base_class_2,
    lower_bound = 0,
    upper_bound = 10,
    n_grid = 500
) {
    raw_curve <- dataset[[index]]
    smooth_curve <- smoothed_dataset[[index]]

    grid <- seq(
        lower_bound,
        upper_bound,
        length.out = n_grid
    )

    # Select the true underlying function from the label
    true_function <- if (smooth_curve$label == -1) {
        base_class_1
    } else {
        base_class_2
    }

    true_values <- true_function(grid)
    fitted_values <- eval.fd(
        grid,
        smooth_curve$curve$fd
    )

    all_y <- c(
        raw_curve$measurement,
        true_values,
        fitted_values
    )

    plot(
        raw_curve$t,
        raw_curve$measurement,
        pch = 16,
        cex = 0.8,
        xlim = c(lower_bound, upper_bound),
        ylim = range(all_y, finite = TRUE),
        xlab = "t",
        ylab = "Valor de la Función",
        main = paste(
            "Curva", "1",
            "- clase",
            ifelse(smooth_curve$label == -1, 1, 2)
        )
    )

    # True underlying function
    lines(
        grid,
        true_values,
        col = "black",
        lwd = 3,
        lty = 2
    )

    # Smoothed approximation
    lines(
        grid,
        fitted_values,
        col = ifelse(
            smooth_curve$label == -1,
            "red",
            "blue"
        ),
        lwd = 3
    )

    legend(
        "bottomleft",
        legend = c(
            "Puntos con ruido",
            "Curva subyacente",
            "Aproximación"
        ),
        pch = c(16, NA, NA),
        lty = c(NA, 2, 1),
        lwd = c(NA, 3, 3),
        col = c(
            "black",
            "black",
            ifelse(
                smooth_curve$label == -1,
                "red",
                "blue"
            )
        ),
        bty = "n"
    )
}

dev.new()

plot_single_approximation(
    index = 1 + n_curves,
    dataset = dataset,
    smoothed_dataset = smoothed_dataset,
    base_class_1 = base_class_1,
    base_class_2 = base_class_2,
    lower_bound = lower_bound,
    upper_bound = upper_bound
)



###
plot_curves_and_weight <- function(
    curves,
    w_fd,
    lower_bound,
    upper_bound,
    n_grid = 500
) {
    stopifnot(inherits(w_fd, "fd"))

    grid <- seq(
        lower_bound,
        upper_bound,
        length.out = n_grid
    )

    curve_values <- lapply(curves, function(obj) {
        if (!inherits(obj$curve$fd, "fd")) {
            stop("A curve does not contain a valid obj$curve$fd object.")
        }

        as.numeric(eval.fd(grid, obj$curve$fd))
    })

    labels <- sapply(curves, function(obj) obj$label)
    w_values <- as.numeric(eval.fd(grid, w_fd))

    dev.new()
    old_par <- par(mfrow = c(2, 1))
    on.exit(par(old_par), add = TRUE)

    plot(
        NA,
        xlim = range(grid),
        ylim = range(unlist(curve_values), finite = TRUE),
        xlab = "t",
        ylab = "x(t)",
        main = "Training curves"
    )

    for (i in seq_along(curves)) {
        curve_colour <- if (labels[i] == -1) {
            rgb(1, 0, 0, alpha = 0.35)
        } else {
            rgb(0, 0, 1, alpha = 0.35)
        }

        lines(
            grid,
            curve_values[[i]],
            col = curve_colour
        )
    }

    plot(
        grid,
        w_values,
        type = "l",
        lwd = 3,
        xlab = "t",
        ylab = "w(t)",
        main = "SVM weight function"
    )

    abline(h = 0, lty = 2)
}

plot_curves_and_weight(
    curves = training_data,
    w_fd = model$w,  # or simply w_fd
    lower_bound = lower_bound,
    upper_bound = upper_bound
)
traceback()
