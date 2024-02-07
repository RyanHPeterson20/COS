#both required
library(foreach)
library(parallel)

detectCores()


num_cores <- 4  # You can adjust this number based on your system's capabilities

# Create a cluster with multiple workers
cl <- makeCluster(num_cores)

# Now, the parallel library will utilize the specified number of cores

# Example of parallel processing:
result <- parLapply(cl, 1:num_cores, function(x) {
  # Your parallelized computation here
  # For demonstration, let's just return the square of the input
  return(x^2)
})

# Close the cluster
stopCluster(cl)

# Print the result
print(result)