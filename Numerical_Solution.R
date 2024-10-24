# Parameters
T <- 2000  # Transmissivity (m^2/d)
L <- 5000  # Length of the island (m)
R <- 0.002 # Recharge (m/d)
dx <- 50   # Grid spacing (m)

# Numerical solution using finite differences
N <- L / dx # Number of grid points
x <- seq(0, L, by = dx)

# Initialize head array
h <- rep(0, length(x))

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
h <- solve(A, b)
plot(x, h, type = "p", 
     col = "blue", xlab = "Distance (m)", ylab="Head (m)",
     main = "Hydraulic Head over Disatance")
