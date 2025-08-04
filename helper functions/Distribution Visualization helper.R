find_threshold = function(mass, Total_mass = 0.98) {
  sorted_mass = sort(mass, decreasing = TRUE)
  sum = 0
  for (i in c(1:length(mass))) {
    sum = sum + sorted_mass[i]
    if (sum >= Total_mass) {
      return (sorted_mass[i+1])
    }
  }
}

plotter_2D_ALL = function(Total_mass, df_2D, title) {
  threshold = find_threshold(df_2D$prob, Total_mass = Total_mass)
  size_breaks = seq(min(df_2D$prob[df_2D$prob > threshold]), max(df_2D$prob), length.out = 10)
  plot_2D = ggplot(df_2D, aes(x = log(x), y = log(y))) +
    geom_point(aes(size = ifelse(prob < threshold, NA, prob)), color = "blue", alpha = 0.8) +
    scale_size_continuous(range = c(0.2, 5), breaks = size_breaks) +
    theme_minimal() +
    xlim(-5, 2) +
    labs(title = title,
         x = expression(log({tilde(sigma)[iA]}^2)),
         y = expression(log({tilde(sigma)[iB]}^2))) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title = element_text(size = 35),
          axis.text = element_text(size = 28),
          legend.position = 'none')
  return (plot_2D)
}

plotter_2D = function(Total_mass, df_2D, title) {
  threshold = find_threshold(df_2D$prob, Total_mass = Total_mass)
  size_breaks = seq(min(df_2D$prob[df_2D$prob > threshold]), max(df_2D$prob), length.out = 10)
  plot_2D = ggplot(df_2D, aes(x = log(x), y = log(y))) +
    geom_point(aes(size = ifelse(prob < threshold, NA, prob)), color = "blue", alpha = 0.8) +
    scale_size_continuous(range = c(0.2, 5), breaks = size_breaks) +
    theme_minimal() +
    xlim(-2.5, 3) +
    labs(title = title,
         x = expression(log({sigma[iA]}^2)),
         y = expression(log({sigma[iB]}^2))) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title = element_text(size = 35),
          axis.text = element_text(size = 28),
          legend.position = 'none')
  return (plot_2D)
}

plotter_1D = function(df_1D, title) {
  plot_1D = ggplot(df_1D, aes(x = log(x), y = 0)) +
    geom_segment(aes(xend = log(x), yend = prob), size = 1, color = "blue") +
    scale_y_continuous(name = "Density", limits = c(0, max(df_1D$prob))) +
    theme_minimal() +
    labs(title = title,
         x = expression(log(lambda[i]))) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 15))
  return (plot_1D)
}

plotter_1D_F = function(df_1D, title) {
  plot_1D = ggplot(df_1D, aes(x = log(x), y = 0)) +
    geom_segment(aes(xend = log(x), yend = prob), size = 1, color = "blue") +
    scale_y_continuous(name = "Density", limits = c(0, max(df_1D$prob))) +
    theme_minimal() +
    labs(title = title,
         x = expression(log({sigma}^2))) +
    theme(plot.title = element_text(hjust = 0.5))
  return (plot_1D)
}