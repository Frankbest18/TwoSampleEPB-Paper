library('ALL')
library(dplyr)
library(ggplot2)

source('../../EPB main/EPB.R')
source('../../helper functions/Rejection Region helper.R')
source('../../helper functions/Distribution Visualization helper.R')

# Data Preprocessing

data(ALL)
bcell = grep("^B", as.character(ALL$BT))
G1_idx = which(ALL[bcell, ]$mol.biol == 'ALL1/AF4')
G2_idx = which(ALL[bcell, ]$mol.biol == 'E2A/PBX1')
picked = c(G1_idx, G2_idx)
ALL_picked = ALL[, picked]
expr_picked = data.frame(GeneID = rownames(exprs(ALL_picked)), exprs(ALL_picked), row.names = NULL)

write.csv(expr_picked, file = "./Filtered ALL/10ALL1|AF4_vs_5E2A|PBX1.csv", row.names = FALSE)
expr = read.csv("./Filtered ALL/10ALL1|AF4_vs_5E2A|PBX1.csv", header = TRUE, stringsAsFactors = FALSE)
m = nrow(expr)

K1 = 10
K2 = 5
X1 = expr[1:m, 2:11] # K1 = 10 ALL1|AF4
X2 = expr[1:m, 12:16] # K2 = 5 E2A|PBX1

# Run VREPB, DVEPB, Welch, B_F and VE_test

info = information_extractor(X1, X2)

alpha = 0.1
NPMLE_1D_parameter = c(1000, 0, 1.0)
NPMLE_2D_parameter = c(80, 80, 0.01, 1.0)
algorithm_list = c(1,2,3,4,5)

#result = solver(X1, X2, NPMLE_1D_parameter, NPMLE_2D_parameter, algorithm_list)

# Number of discoveries for each method

#length(my_BH(result$VREPB, alpha)) # 320
#length(my_BH(result$DVEPB, alpha)) # 402
#length(my_BH(result$Welch, alpha)) # 168
#length(my_BH(result$B_F, alpha)) # 24
#length(my_BH(result$EV_test, alpha)) # 303

#df_2D_ALL = data.frame(x = result$'2D_grid'[,1], y = result$'2D_grid'[,2], prob = result$'2D_mass')
#write.csv(df_2D_ALL, './Result/2D_ALL.csv')

# Generate the H plot for ALL

df_2D_ALL = read.csv('./Result/2D_ALL.csv')

Total_mass = 0.999999
plot_H_ALL = plotter_2D_ALL(Total_mass, df_2D_ALL, '')
ggsave(filename = './Result/H(var1, var2) of ALL.jpg', plot = plot_H_ALL, width = 12, height = 10, dpi = 500)

# Generate the G plot for ALL

df_1D_ALL = data.frame(x = df_2D_ALL$x/df_2D_ALL$y, prob = df_2D_ALL$prob)

plot_G_ALL = plotter_1D(df_1D_ALL, '')
ggsave(filename = './Result/G(lambda) of ALL.jpg', plot = plot_G_ALL, width = 8, height = 6, dpi = 500)

# Generate the Rejection Region plot

grid = df_2D_ALL[, 2:3]
mass = df_2D_ALL[, 4]
grid_1D = grid[, 1] / grid[, 2]

tau_list = info$S1_list/info$S2_list
log_tau_list = log(tau_list)
log_tau_list = log_tau_list[log_tau_list >= -5]
tau_grid = exp(seq(min(log_tau_list), max(log_tau_list), length.out = 1000))

t_BF_list = c()
for (i in 1:m) {
  X1_i = as.numeric(X1[i,])
  X2_i = as.numeric(X2[i,])
  m1_i = mean(X1_i)
  m2_i = mean(X2_i)
  s1_i = var(X1_i)
  s2_i = var(X2_i)
  t_BF_i = (m1_i - m2_i) / sqrt(s1_i/K1 + s2_i/K2)
  t_BF_list = c(t_BF_list, t_BF_i)
}
t_BF_list = t_BF_list[log(tau_list) >= -5] # Range of sample t_BF to determine range of root solver

## Rejection Region for VREPB

t_VREPB = c()

for (tau in tau_grid) {
  t_VREPB_i = BF_pvalue_VREPB_solver_positive(K1, K2, grid_1D, mass, tau, 0.05)
  t_VREPB = append(t_VREPB, t_VREPB_i)
}

## Rejection Region for B_F
t_BF = c()

for (tau in tau_grid) {
  t_BF_i = BF_pvalue_BF_solver_positive(tau, K1, K2, 0.05)
  t_BF = append(t_BF, t_BF_i)
}

## Rejection Region for Welch

t_Welch = c()

for (tau in tau_grid) {
  t_Welch_i = BF_pvalue_Welch_solver_positive(tau, K1, K2, 0.05)
  t_Welch = append(t_Welch, t_Welch_i)
}

## Rejection Region for VE-test

t_EV = c()

for (tau in tau_grid) {
  t_EV_i = BF_pvalue_EV_solver_positive(tau, K1, K2, 0.05)
  t_EV = append(t_EV, t_EV_i)
}

## Plot rejection region for VREPB, B_F, Welch, EV_test

size = 1

curve_colors = c(
  "VREPB"   = "#0072B2",   # blue
  "B-F"     = "#CC79A7",   # magenta
  "Welch"   = "#E69F00",   # gold
  "EV-test"  = "#D55E00"     # orange
)

plot_2D = function(u1, u2) {
  ggplot(data.frame(u1 = u1, u2 = u2), aes(x = u1, y = u2)) +
    stat_bin2d(bins = 80) +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    labs(x = expression(log(hat(lambda)[i])), y = expression(T[i]^{BF}), fill = "Count") +
    scale_y_continuous(limits = c(-4, 4)) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.major = element_line(colour = "gray80", size=0.5),
      panel.grid.minor = element_line(colour = "gray90", size=0.25),
      legend.position = "top",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      plot.margin = margin(12, 16, 12, 12)
    )
}

p = plot_2D(log_tau_list, t_BF_list) +
  geom_line(data = data.frame(x = log(tau_grid), y = t_VREPB), aes(x = x, y = y, color = "VREPB"), size = size) +
  geom_line(data = data.frame(x = log(tau_grid), y = -t_VREPB), aes(x = x, y = y, color = "VREPB"), size = size) +
  
  geom_line(data = data.frame(x = log(tau_grid), y = t_BF), aes(x = x, y = y, color = "B-F"), size = size) +
  geom_line(data = data.frame(x = log(tau_grid), y = -t_BF), aes(x = x, y = y, color = "B-F"), size = size) +
  
  geom_line(data = data.frame(x = log(tau_grid), y = t_EV), aes(x = x, y = y, color = "EV-test"), size = size) +
  geom_line(data = data.frame(x = log(tau_grid), y = -t_EV), aes(x = x, y = y, color = "EV-test"), size = size) +
  
  geom_line(data = data.frame(x = log(tau_grid), y = t_Welch), aes(x = x, y = y, color = "Welch"), size = size) +
  geom_line(data = data.frame(x = log(tau_grid), y = -t_Welch), aes(x = x, y = y, color = "Welch"), size = size) +
  
  scale_color_manual(
    name = NULL, 
    values = curve_colors,
    breaks = c("VREPB", "B-F", "Welch", "EV-test"),
    labels = c("VREPB", "B-F", "Welch", "EV-test")
  ) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 0.7), 
         color = guide_legend(override.aes = list(size = 2)))

# Save as JPG

ggsave("./Result/Rejection_Region_ALL.jpg", p, width = 7, height = 5.5, dpi = 400)