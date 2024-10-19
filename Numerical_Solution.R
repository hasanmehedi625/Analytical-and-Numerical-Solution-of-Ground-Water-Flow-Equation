# Parameters
T <- 2000  # Transmissivity (m^2/d)
L <- 5000  # Length of the island (m)
R <- 0.002 # Recharge (m/d)
dx <- 50   # Grid spacing (m)

# Analytical solution
x_analytical <- seq(0, L, length.out = 100)
h_analytical <- (R / (2 * T)) * (x_analytical^2 - L^2)

# Numerical solution using finite differences
N <- L / dx # Number of grid points
x_numerical <- seq(0, L, by = dx)

# Initialize head array
h_numerical <- rep(0, length(x_numerical))

# Matrix A and vector b for the system Ah = b
A <- matrix(0, nrow = N + 1, ncol = N + 1)
b <- rep(0, N + 1)

# Fill the matrix A and vector b
for (i in 2:N) {
  A[i, i-1] <- T / dx^2
  A[i, i]   <- -2 * T / dx^2
  A[i, i+1] <- T / dx^2
  b[i] <- -R
}

# Boundary conditions
A[1, 1] <- 1           # dh/dx = 0 at x = 0 (symmetry condition)
A[N+1, N+1] <- 1       # h(L) = 0 at x = L

# Solve the system of equations Ah = b
h_numerical <- solve(A, b)

# Set plot margins
#par(mar = c(2, 4, 2, 2))  # Adjust margins as needed (bottom, left, top, right)
# Plot results
plot(x_analytical, h_analytical, type = "l", col = "blue", lwd = 2,
     xlab = "Distance from center (m)", ylab = "Hydraulic head (m)",
     main = "Head Distribution over Distance")
points(x_numerical, h_numerical, col = "red", pch = 16)
legend("bottomright", legend = c("Analytical Solution", "Numerical Solution (dx=50m)"),
       col = c("blue", "red"), lty = 1, pch = c(NA, 16))
