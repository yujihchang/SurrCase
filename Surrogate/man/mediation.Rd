url <- "https://github.com/yujihchang/SurrCase/raw/main/Surrogate/Data/Wild%20Type%20Live%20Virus.RData"
destfile <- "Wild Type Live Virus.RData"
if (!file.exists(destfile)) {
  download.file(url, destfile, mode = "wb")
}
load("Wild Type Live Virus.RData")

# Ensure data is a data frame with at least two dimensions
if (!is.data.frame(data)) {
  stop("The 'data' object must be a data frame.")
}

# Add random noise to separate ties
ee = runif(30000, 0, 0.5)
data[, 1] = data[, 1] + ee  # Separate ties

# Rename columns for clarity
colnames(data)[1:5] = c("Ts", "D", "xi", "Z", "M")

# Define ppt sequence
ppt = seq(0, max(data$Ts[data$D == 1]), by = 0.01)

# Perform IV analysis
ans_IV = IV(ppt, data, m0 = 0.5, m1 = 2.5, SS = 2.25^0.5, tol = 0.001, PX = c(12, 13, 14), setting = "median")

# Plot results
plot(ppt, ans_IV$WIE, type = "l", main = "IV Analysis Results", xlab = "ppt", ylab = "WIE")
points(ppt, ans_IV$WIE + 1.96 * ans_IV$WIE.v, col = "blue")  # Upper bound
points(ppt, ans_IV$WIE - 1.96 * ans_IV$WIE.v, col = "blue")  # Lower bound
abline(h = 0, col = 2)  # Add a horizontal line at y = 0
