# Set the number of iterations
num_iter <- 1000
# Initialize an empty vector to store the p-values
p_values <- numeric(num_iter)
# Loop over the number of iterations
for (i in 1:num_iter) {
  
  # Generate two vectors of 1000 length with mean 0 and SD 1
  group1 <- rnorm(1000, mean = 0, sd = 1)
  group2 <- rnorm(1000, mean = 0, sd = 1)
  
  # Perform a two-sample t-test and extract the p-value
  result <- t.test(group1, group2)
  p_values[i] <- result$p.value
  
}
# Count the number of p-values below 0.05
num_sig <- sum(p_values < 0.05)
# Print the proportion of significant p-values
cat("Proportion of significant p-values:", num_sig / num_iter, "\n")
