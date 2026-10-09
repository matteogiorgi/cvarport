# Generates img/contributions.png: for each asset of the current portfolio
# of the synthetic example, its weight next to its share of the
# volatility and of the CVaR95 (Euler contributions, from
# examples/contributions.csv). A 10% position can carry far more than
# 10% of the risk: here the two equity assets carry almost all of it.
#
# Run from the repository root, after examples/analysis.jl.

col_weight <- "grey65"
col_vol <- "#1b9e77"
col_cvar <- "#d95f02"

contrib <- read.csv("examples/contributions.csv")
M <- 100 * t(as.matrix(contrib[, c("weight", "vol", "cvar")]))
colnames(M) <- contrib$asset

# Spectral is the font of the GitHub page (geoteo.net/static/style.css);
# the Cairo device is needed to reference it by family name.
png("img/contributions.png", width = 1800, height = 850, res = 150, type = "cairo")
par(mar = c(3, 4.2, 6.2, 1), family = "Spectral")

mids <- barplot(M,
    beside = TRUE, col = c(col_weight, col_vol, col_cvar), border = NA,
    ylim = c(0, 100), las = 1, ylab = "% of the total", cex.axis = 0.8, cex.names = 0.85, cex.lab = 0.85
)
abline(h = seq(20, 100, by = 20), col = "grey92")
barplot(M, beside = TRUE, col = c(col_weight, col_vol, col_cvar), border = NA, add = TRUE, axes = FALSE, names.arg = rep("", ncol(M)))
text(mids, M, sprintf("%.1f", M), pos = 3, cex = 0.65, offset = 0.3)

mtext("Who carries the risk of the current portfolio", side = 3, line = 3.0, cex = 1.15, font = 2)
mtext("weights against Euler shares of the volatility and of the one-day CVaR95 (synthetic example)", side = 3, line = 1.5, cex = 0.75)

legend("topright",
    legend = c("weight", "share of the volatility", "share of the CVaR95"),
    fill = c(col_weight, col_vol, col_cvar), border = NA, bty = "n", cex = 0.72, inset = c(0.0, 0.02)
)

dev.off()
