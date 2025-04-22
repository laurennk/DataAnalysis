#load potential energy surface for H2 TDDFT using B3LYP functional
library(tidyverse)
library(ggplot2)
library(ggmagnify)

# Read data from csv files for potential energy surface without cavity and for 
# QED at various lambda values
#
# Returns dataframe
read_in_data <- function(functional) {
  
  print(paste("The functional is: ", functional))
  
  file_name <- paste0("H2_", functional, "_pes.csv")
  
  H2_pes_data <- read.csv(file_name, header=TRUE, sep=",")
  
  lambdas <- c(0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
  
  H2_qed_files <- paste0("H2_", functional, "_qed_", lambdas, ".csv")
  
  H2_qed_data <- lapply(H2_qed_files, read.csv)
  
  H2_data <- c(list(H2_pes_data), H2_qed_data)
  
  return (H2_data %>% reduce(full_join, by="r"))
}

# Create plot of first excited state for no cavity and QED with various values of lambda
#
# Returns plot
generate_plot <- function(data, title, xmin, xmax, ymin, ymax, image_name) {
  
  my_plot <- ggplot(data , aes(r, gs.bwd)) +
    scale_color_gradient2(midpoint = 0.07, low = "green", mid = "blue", high = "blue", name = "Lambda Value") +
    geom_line(data = data, mapping = aes(x = r, y = LP_energy_0.01, color=lambda_0.01)) + 
    geom_line(data = data, mapping = aes(x = r, y = UP_energy_0.01, color=lambda_0.01)) +
    geom_line(data = data, mapping = aes(x = r, y = LP_energy_0.02, color=lambda_0.02)) + 
    geom_line(data = data, mapping = aes(x = r, y = UP_energy_0.02, color=lambda_0.02)) +
    geom_line(data = data, mapping = aes(x = r, y = LP_energy_0.03, color=lambda_0.03)) + 
    geom_line(data = data, mapping = aes(x = r, y = UP_energy_0.03, color=lambda_0.03)) +
    geom_line(data = data, mapping = aes(x = r, y = LP_energy_0.04, color=lambda_0.04)) + 
    geom_line(data = data, mapping = aes(x = r, y = UP_energy_0.04, color=lambda_0.04)) +
    geom_line(data = data, mapping = aes(x = r, y = LP_energy_0.05, color=lambda_0.05)) + 
    geom_line(data = data, mapping = aes(x = r, y = UP_energy_0.05, color=lambda_0.05)) +
    geom_line(data = data, mapping = aes(x = r, y = LP_energy_0.1, color=lambda_0.1)) + 
    geom_line(data = data, mapping = aes(x = r, y = UP_energy_0.1, color=lambda_0.1)) +
    geom_line(data = data, mapping = aes(x = r, y = es.energy)) +
    geom_line(data = data, mapping = aes(x = r, y = energy_0.0)) +
    ggtitle(title) +
    xlab("H-H bond distance (angstrom)") +
    ylab("Energy (hartree)") +
    scale_x_continuous(n.breaks=10) +
    scale_y_continuous(n.breaks=10) +
    coord_cartesian(xlim = c(xmin,xmax), ylim = c(ymin, ymax))
  
  from <- c(xmin = 0.5, xmax = 1.0, ymin = -0.75, ymax = -0.5)
  to <- c(xmin = 0.5, xmax = 1.9, ymin = -0.2, ymax = 0.4)
  
  my_plot <- my_plot + geom_magnify(from = from, to = to)
  
  ggsave(image_name)
  
  return (my_plot)
}

# Fit potential energy well with harmonic and anharmonic potentials and output 
# results to file or console
#
# Returns list of models
fit_ground_state <- function(data, num_points, out_file_name=NULL) {
  
  min_energy <- min(data[,"gs.bwd"])
  min_idx <- as.numeric(which(data$gs.bwd == min_energy))
  
  fitting_data <- data[(min_idx-num_points):(min_idx+num_points),]
  
  harmonic_model <- lm(gs.bwd~poly(r,2,raw=TRUE), data=fitting_data)
  anharmonic_model <- lm(gs.bwd~poly(r,3,raw=TRUE), data=fitting_data)
  
  #
  if (!is.null(out_file_name)) {
    out_file <- file(out_file_name)
    
    write_lines(paste("Harmonic approximation R-squared:   ", summary(harmonic_model)$r.squared,
                      "\nHarmonic force constant (k):        ", summary(harmonic_model)$coef[3],
                      "\nAnharmonic approximation R-squared: ", summary(anharmonic_model)$r.squared,
                      "\nHarmonic force constant (k):        ", summary(anharmonic_model)$coef[3],
                      "\nAnharmonic constant (gamma):        ", summary(anharmonic_model)$coef[4]), out_file)
  }
  
  else {
    print(paste("Harmonic approximation R-squared:   ", summary(harmonic_model)$r.squared))
    print(paste("Harmonic force constant (k):        ", summary(harmonic_model)$coef[3]))
    print(paste("Anharmonic approximation R-squared: ", summary(anharmonic_model)$r.squared))
    print(paste("Harmonic force constant (k):        ", summary(anharmonic_model)$coef[3]))
    print(paste("Anharmonic constant (gamma):        ", summary(anharmonic_model)$coef[4]))
  }
  
  return (list(harmonic_model, anharmonic_model))
}



functionals <- c("B3LYP", "PBE", "wB97X")

for (f in functionals) {
  #Read in data
  H2_data <- read_in_data(f)
  
  #Set bounds for axes
  xmin <- min(H2_data$r)
  xmax <- max(H2_data$r)
  ymin <- -0.9
  ymax <-  0.4
  
  #Create plot of data for each functional
  title = bquote(H[2]~.(paste(f)))
  
  plot_H2 <- generate_plot(H2_data, title, xmin, xmax, ymin, ymax, paste0("H2_", f, "_lambda_plot.pdf"))
  
  #Fit potential energy well with harmonic and anharmonic potentials
  num_points <- 5
  
  models <- fit_ground_state(H2_data, num_points, out_file_name = paste0("H2_", f, ".txt"))
}

#Visually inspect fits for wB97x ground state 
results <- data.frame(r=seq(0.0, 2.0, length = 100))
results$harmonic_model <- predict(models[[1]], results)
results$anharmonic_model <- predict(models[[2]], results)
  
ggplot(H2_data, aes(r, gs.bwd)) +
  geom_point(data = H2_data, mapping = aes(x = r, y = gs.bwd), color="black") +
  geom_line(data = results, mapping = aes(x=r, y=harmonic_model), color="blue") +
  geom_line(data = results, mapping = aes(x=r, y=anharmonic_model), color="red") +
  scale_y_continuous(n.breaks=10) +
  coord_cartesian(xlim = c(0.25, 1.25), ylim = c(-1.18,-1.07)) +
  ggtitle(paste0(f, " Ground State Fit Lines"))

