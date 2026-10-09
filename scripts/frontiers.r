# Generates img/frontiers.png: the mean-variance and mean-CVaR efficient
# frontiers of the synthetic example (examples/frontier.csv), each drawn
# both in the (volatility, mean) plane and in the (CVaR, mean) plane,
# with the reference portfolios of examples/portfolios.csv on top. The
# arrow goes from the current portfolio to the one with the same
# expected return and minimum CVaR: the gain the report quantifies.
#
# Run from the repository root, after examples/analysis.jl.

col_mv <- "#1b9e77"
col_cvar <- "#d95f02"

lighten <- function(col, amount = 0.85) {
    rgb_col <- col2rgb(col) / 255
    rgb_new <- rgb_col + (1 - rgb_col) * amount
    rgb(rgb_new[1, ], rgb_new[2, ], rgb_new[3, ])
}

frontier <- read.csv("examples/frontier.csv")
portfolios <- read.csv("examples/portfolios.csv", check.names = FALSE)
mv <- frontier[frontier$frontier == "mv", ]
cv <- frontier[frontier$frontier == "cvar", ]

# The report stores the Sharpe ratio, not rf: recover rf from the
# tangency portfolio, rf = mean - Sharpe * vol.
tg <- portfolios[portfolios$portfolio == "tangency", ]
rf <- tg$mean - tg$sharpe * tg$vol

pct <- function(x) 100 * x
cur <- portfolios[portfolios$portfolio == "current", ]
best <- portfolios[portfolios$portfolio == "min CVaR, mean >= current", ]
refs <- portfolios[portfolios$portfolio %in% c("min variance", "tangency", "min CVaR"), ]

# Spectral is the font of the GitHub page (geoteo.net/static/style.css);
# the Cairo device is needed to reference it by family name.
png("img/frontiers.png", width = 1800, height = 850, res = 150, type = "cairo")
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 0.6, 1), oma = c(0, 0, 4.4, 0), family = "Spectral")

panel <- function(x_of, xlab, risk_col, show_cml) {
    xs <- c(x_of(mv), x_of(cv), x_of(portfolios))
    plot(NA,
        xlim = c(0, pct(max(xs)) * 1.05), ylim = c(min(0, rf), max(frontier$mean) * 1.08) * 100,
        xlab = xlab, ylab = "annualized mean (%)", bty = "n", las = 1, cex.axis = 0.8, cex.lab = 0.85
    )
    abline(h = 0, col = "grey88")
    if (show_cml) {
        # Capital market line: rf to the tangency portfolio and beyond.
        slope <- (tg$mean - rf) / tg$vol
        abline(a = pct(rf), b = slope, col = "grey60", lty = 3, lwd = 1.2)
        points(0, pct(rf), pch = 21, bg = "white", cex = 1.1)
        text(0, pct(rf), "rf", pos = 4, cex = 0.7, col = "grey40")
    }
    lines(pct(x_of(mv)), pct(mv$mean), col = col_mv, lwd = 2.4)
    lines(pct(x_of(cv)), pct(cv$mean), col = col_cvar, lwd = 2.4, lty = 2)

    arrows(pct(x_of(cur)), pct(cur$mean), pct(x_of(best)) * 1.015, pct(best$mean),
        length = 0.08, col = "grey35", lwd = 1.4
    )
    points(pct(x_of(refs)), pct(refs$mean), pch = 19, cex = 1.2, col = "grey20")
    # min variance and min CVaR almost overlap: label them on opposite sides.
    side <- c("min variance" = 2, "tangency" = 2, "min CVaR" = 4)[refs$portfolio]
    text(pct(x_of(refs)), pct(refs$mean), refs$portfolio, pos = side, cex = 0.68, offset = 0.6)
    points(pct(x_of(cur)), pct(cur$mean), pch = 21, bg = lighten(risk_col, 0.3), col = "grey20", cex = 1.6, lwd = 1.3)
    text(pct(x_of(cur)), pct(cur$mean), "current", pos = 4, cex = 0.72, font = 2, offset = 0.7)
}

panel(function(d) d$vol, "annualized volatility (%)", col_mv, TRUE)
legend("bottomright",
    legend = c("mean-variance frontier", "mean-CVaR frontier", "capital market line"),
    col = c(col_mv, col_cvar, "grey60"), lty = c(1, 2, 3), lwd = c(2.4, 2.4, 1.2),
    bty = "n", cex = 0.72, seg.len = 2.4
)
panel(function(d) d$cvar, "one-day CVaR95 (%)", col_cvar, FALSE)
legend("bottomright",
    legend = c("reference portfolios", "current portfolio", "same mean, minimum CVaR"),
    pch = c(19, 21, NA), pt.bg = c(NA, lighten(col_cvar, 0.3), NA), pt.cex = c(1.2, 1.6, NA),
    lty = c(NA, NA, 1), col = c("grey20", "grey20", "grey35"), bty = "n", cex = 0.72, seg.len = 2.4
)

mtext("Efficient frontiers of the synthetic example", side = 3, line = 2.4, outer = TRUE, cex = 1.15, font = 2)
mtext("long-only, cap 0.7, Ledoit-Wolf covariance, historical CVaR at alpha = 0.95 on 1260 daily scenarios",
    side = 3, line = 0.9, outer = TRUE, cex = 0.75
)

dev.off()
