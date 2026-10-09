# Generates img/losses.png: the 1260 daily losses of the current portfolio
# of the synthetic example (examples/prices.csv and examples/weights.txt),
# with the VaR95 and CVaR95 that cvarport reports (examples/portfolios.csv)
# and the normal density with the same mean and standard deviation. The
# tail beyond the VaR is what the CVaR averages; the normal curve shows
# how much thinner that tail would be under a Gaussian model.
#
# Run from the repository root, after examples/analysis.jl.

col_body <- "#1b9e77"
col_tail <- "#d95f02"

lighten <- function(col, amount = 0.85) {
    rgb_col <- col2rgb(col) / 255
    rgb_new <- rgb_col + (1 - rgb_col) * amount
    rgb(rgb_new[1, ], rgb_new[2, ], rgb_new[3, ])
}

prices <- read.csv("examples/prices.csv")
P <- as.matrix(prices[, -1])
R <- P[-1, ] / P[-nrow(P), ] - 1

# weights.txt: "TICKER amount" lines, "#" starts a comment, normalized.
lines_w <- trimws(sub("#.*", "", readLines("examples/weights.txt")))
fields <- strsplit(lines_w[lines_w != ""], "[[:space:]]+")
w <- setNames(rep(0, ncol(P)), colnames(P))
for (f in fields) w[f[1]] <- as.numeric(f[2])
w <- w / sum(w)

L <- -drop(R %*% w)
portfolios <- read.csv("examples/portfolios.csv", check.names = FALSE)
cur <- portfolios[portfolios$portfolio == "current", ]

pct <- function(x) 100 * x
breaks <- seq(floor(pct(min(L)) * 4) / 4, ceiling(pct(max(L)) * 4) / 4, by = 0.125)

# Spectral is the font of the GitHub page (geoteo.net/static/style.css);
# the Cairo device is needed to reference it by family name.
png("img/losses.png", width = 1800, height = 850, res = 150, type = "cairo")
par(mar = c(4.2, 4.2, 6.2, 1), family = "Spectral")

h <- hist(pct(L), breaks = breaks, plot = FALSE)
in_tail <- h$mids > pct(cur$var)
plot(h,
    freq = FALSE, col = ifelse(in_tail, lighten(col_tail, 0.35), lighten(col_body, 0.6)),
    border = "white", main = "", xlab = "one-day loss (% of the portfolio, positive = money lost)",
    ylab = "density", las = 1, cex.axis = 0.8, cex.lab = 0.85
)
mtext("Daily losses of the current portfolio", side = 3, line = 3.0, cex = 1.15, font = 2)
mtext("synthetic example: 60% EQTY, 20% BOND, 10% GOLD, 10% SMCP over 1260 trading days", side = 3, line = 1.5, cex = 0.75)

x <- seq(min(breaks), max(breaks), length.out = 400)
lines(x, dnorm(x, mean(pct(L)), sd(pct(L))), col = "grey30", lwd = 1.6, lty = 2)

abline(v = pct(cur$var), col = col_tail, lwd = 2)
abline(v = pct(cur$cvar), col = col_tail, lwd = 2, lty = 4)
top <- par("usr")[4]
text(pct(cur$var), top * 0.92, sprintf("VaR95 = %.2f%%", pct(cur$var)), pos = 2, cex = 0.75, col = col_tail)
text(pct(cur$cvar), top * 0.80, sprintf("CVaR95 = %.2f%%", pct(cur$cvar)), pos = 4, cex = 0.75, col = col_tail)

legend("topleft",
    legend = c("losses below the VaR", "worst 5% of the days (averaged by the CVaR)", "normal with the same mean and sd"),
    fill = c(lighten(col_body, 0.6), lighten(col_tail, 0.35), NA), border = c("white", "white", NA),
    lty = c(NA, NA, 2), lwd = c(NA, NA, 1.6), col = c(NA, NA, "grey30"),
    bty = "n", cex = 0.72, seg.len = 2.2, inset = c(0.0, 0.02)
)

dev.off()
