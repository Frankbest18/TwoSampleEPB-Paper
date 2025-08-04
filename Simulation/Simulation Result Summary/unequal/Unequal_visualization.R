library(RColorBrewer)

Full_FDR_Data = read.csv('./Full Data/Unequal_FDR_Data.csv')
Full_Power_Data = read.csv('./Full Data/Unequal_Power_Data.csv')

method_colors = c(
  EV_test = "#D55E00",   # orange (plus)
  Welch  = "#E69F00",    # gold (triangle)
  B_F    = "#CC79A7",    # magenta (diamond)
  VREPB  = "#0072B2",    # blue (square)
  DVEPB  = "#009E73"     # green (circle)
)

method_pch = c(
  EV_test = 3,     # plus sign
  Welch   = 17,    # triangle
  B_F     = 18,    # diamond
  VREPB   = 15,    # square
  DVEPB   = 16     # circle
)

FDR_plotter = function(n1, length_n2, Full_FDR_Data) {
  n2_range = n1:(n1 + length_n2 - 1)
  dat = Full_FDR_Data[Full_FDR_Data$n1 == n1 & (Full_FDR_Data$n2 %in% n2_range), ]
  
  plot(dat$n2, dat$VREPB, type = 'o', pch = method_pch["VREPB"], col = method_colors["VREPB"], 
       lwd = 3, cex = 1.5, xlab = expression(K[B]), ylab = "FDR",
       ylim = c(0, max(dat[, c("VREPB", "DVEPB", "Welch", "B_F", "EV_test")], 0.13) + 0.05),
       cex.lab = 1.5, cex.axis = 1.5, bty = "n")
  lines(dat$n2, dat$DVEPB, type = 'o', pch = method_pch["DVEPB"], col = method_colors["DVEPB"], lwd = 3, cex = 1.5)
  lines(dat$n2, dat$Welch, type = 'o', pch = method_pch["Welch"], col = method_colors["Welch"], lwd = 3, cex = 1.5)
  lines(dat$n2, dat$B_F, type = 'o', pch = method_pch["B_F"], col = method_colors["B_F"], lwd = 3, cex = 1.5)
  lines(dat$n2, dat$EV_test, type = 'o', pch = method_pch["EV_test"], col = method_colors["EV_test"], lwd = 3, cex = 1.5)
  abline(h = 0.1, lwd = 2, col = "gray40", lty = 2)
}

Power_plotter = function(n1, length_n2, Full_Power_Data) {
  n2_range = n1:(n1 + length_n2 - 1)
  dat = Full_Power_Data[Full_Power_Data$n1 == n1 & (Full_Power_Data$n2 %in% n2_range), ]
  
  plot(dat$n2, dat$VREPB, type = 'o', pch = method_pch["VREPB"], col = method_colors["VREPB"],
       lwd = 3, cex = 1.5, xlab = expression(K[B]), ylab = "Power",
       ylim = c(0, max(dat[, c("VREPB", "DVEPB", "Welch", "B_F", "EV_test")])), 
       cex.lab = 1.5, cex.axis = 1.5, bty = "n")
  lines(dat$n2, dat$DVEPB, type = 'o', pch = method_pch["DVEPB"], col = method_colors["DVEPB"], lwd = 3, cex = 1.5)
  lines(dat$n2, dat$Welch, type = 'o', pch = method_pch["Welch"], col = method_colors["Welch"], lwd = 3, cex = 1.5)
  lines(dat$n2, dat$B_F, type = 'o', pch = method_pch["B_F"], col = method_colors["B_F"], lwd = 3, cex = 1.5)
  lines(dat$n2, dat$EV_test, type = 'o', pch = method_pch["EV_test"], col = method_colors["EV_test"], lwd = 3, cex = 1.5)
}

# Save K_A = 3

png("./Graph/Pub Unequal K_A = 3.png", width = 3400, height = 1400, res = 300, family = "sans")
par(mfrow = c(1,2), mar = c(5.2,5.2,2,2), oma = c(1,1,0,13), mgp = c(3.2,1,0))

FDR_plotter(3, 7, Full_FDR_Data)
Power_plotter(3, 7, Full_Power_Data)

par(xpd = NA) # allow drawing in the outer margin
legend("topright", inset = c(-0.7, 0),
       legend = c("EV-test", "Welch", "B-F", "VREPB", "DVEPB"),
       col = method_colors, pch = method_pch, lwd = 3, pt.cex = 2,
       cex = 1.2, y.intersp = 1.3, bty = "n", xpd = NA)
dev.off()

# Save K_A = 5

png("./Graph/Pub Unequal K_A = 5.png", width = 3400, height = 1400, res = 300, family = "sans")
par(mfrow = c(1,2), mar = c(5.2,5.2,2,2), oma = c(1,1,0,13), mgp = c(3.2,1,0))

FDR_plotter(5, 7, Full_FDR_Data)
Power_plotter(5, 7, Full_Power_Data)

par(xpd = NA) # allow drawing in the outer margin
legend("topright", inset = c(-0.7, 0),
       legend = c("EV-test", "Welch", "B-F", "VREPB", "DVEPB"),
       col = method_colors, pch = method_pch, lwd = 3, pt.cex = 2,
       cex = 1.2, y.intersp = 1.3, bty = "n", xpd = NA)
dev.off()