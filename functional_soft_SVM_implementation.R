library(fda)
library(fda.usc)
library(e1071)



#set.seed(1) 
n_basis <- 30

n_curves <- 10

#Define the underlying functions for the two classes we want to classify
base_class_1 <- function(t) { 1-t^2 }
base_class_2 <- function(t) { cos((pi/2) * t) }



base_class_1 <- function(t) {
  sin(2 * pi * t)
}

base_class_2 <- function(t) {
  sin(4 * pi * t)
}

classes <- list(base_class_1, base_class_2)
#Randomize measurements and time for each curve type
dataset <- list()
smoothed_dataset <- list()
counter <- 1
basis <- create.fourier.basis(rangeval = c(0, 1), nbasis = n_basis)
fdpar <- fdPar(basis, Lfdobj = 2, lambda = 1e-4) #2 for second derivative penalty


for (j in 1:length(classes)) {
    for (i in 1:n_curves) {
        n_points <- sample(25:50, 1)
        t <- sort(runif(n_points, 0, 1)) 
        noise <- rnorm(n_points, mean = 0, sd = 0.5) # Add noise to the measurements # nolint
        measurement <- classes[[j]](t) + noise
        dataset[[counter]] <- data.frame(t = t, measurement = measurement, class = j)

        # Fit the Fourier basis to the noisy measurements
        smoothed_data <- smooth.basis(t, measurement, fdpar)
        smoothed_dataset[[counter]] <- smoothed_data
        counter <- counter + 1
    }
}

#Smooth out the functions using GCV? I need to find best gcv value
loglambda <- seq(-7, 7, by = 1)

for (ilam in 1:length(loglambda)){
    lambda <- 10^loglambda[ilam]
    fdparobj <- fdPar(basis, Lfdobj = 2, lambda = lambda)
    total_gcv <- 0
    total_df <- 0
    
    for (i in seq(1,length(dataset),by=1)){
        smoothed_data_i <- smooth.basis(dataset[[i]]$t, dataset[[i]]$measurement, fdparobj)
        
    }
}





#Split into training and testing sets 80-20

#Train on the noisy data

#Test on the testing set and evaluate performance



curve_id <- n_curves #+ 1 # Select the curve to visualize (e.g., the first curve of class 2)

t_grid <- seq(0, 1, length.out = 1000)

# Clase de la curva elegida: usar solo el primer valor
class_id <- dataset[[curve_id]]$class[1]

# Curva original verdadera
y_true <- classes[[class_id]](t_grid)

# Puntos observados con ruido
t_noisy <- dataset[[curve_id]]$t
y_noisy <- dataset[[curve_id]]$measurement

# Curva suavizada
smooth_fd <- smoothed_dataset[[curve_id]]$fd
y_smooth <- eval.fd(t_grid, smooth_fd)

plot(
  t_grid, y_true,
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "t",
  ylab = "f(t)",
  main = paste("Curve", curve_id, "- Original vs Noisy vs Smoothed"),
  ylim = range(c(y_true, y_noisy, y_smooth))
)

points(
  t_noisy, y_noisy,
  col = "red",
  pch = 16,
  cex = 1
)

lines(
  t_grid, y_smooth,
  col = "darkgreen",
  lwd = 2
)

legend(
  "bottomleft",
  legend = c("Original curve", "Noisy points", "Smoothed curve"),
  col = c("blue", "red", "darkgreen"),
  lty = c(1, NA, 1),
  pch = c(NA, 16, NA),
  lwd = c(2, NA, 2),
  bty = "n"
)