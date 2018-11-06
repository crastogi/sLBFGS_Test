# Loss function, gradients and hessian vector products for SVMs with quadratic loss
println("Loading SVM module")

# Computes the likelihood of the function.
# Data matrix columns consist of independent observations
# Data matrix rows consist of independent variables with last row being dependent variable
# Theta is current set of parameters
# Lambda is the regularization parameter
function loss(data::Matrix, theta::Vector, lambda::Float64)
  likelihood = 0.0
  n = size(data, 2)				# Size of dataset. Given by number of columns
  for i = 1:n
    x = data[1:end-1,i]			# Independent variables are everything except the last row
    y = data[end,i]				# Last entry in column is observed value to be predicted
    likelihood += 0.5 * max(0.0, 1 - y * dot(vec(x), theta))^2/n      # Same as Likelihood += 1/(2n)*max(0, 1-y*(X.T))^2, where X.T is the dot product of X and Theta (parameter vector) 
  end
  likelihood += lambda/2 * dot(theta, theta)		# Regularizer. Not re-scaled to data-depth when using stochastic sampling?
  return likelihood
end

# Helper function that computes the gradient for a given point.
# Called grad! because it updates the 'grad' variable
function grad!(datapoint::Vector, theta::Vector, lambda::Float64, grad::Vector)
  # Split individual datapoint 
  x = datapoint[1:end-1]
  y = datapoint[end]
  # This is used to account for the piecewise nature of the max function
  alpha = dot(x, theta)
  if y * alpha < 1.0		# This is the case where max is > 0
  	# grad = grad + L*T - (1-y*X.T)*(y*X) 
    BLAS.axpy!(1.0, lambda * theta -  (1 - y*alpha) * y * x, grad)
  else
  	# grad = grad + L*T
    BLAS.axpy!(1.0, lambda * theta, grad)
  end
 end

# Compute the full gradient over the whole data matrix; calls grad! from above
function full_grad(grad!::Function, data::Matrix, theta::Vector, lambda::Float64)
  n = size(data, 2)						# Size of dataset. Given by the number of columns
  dLikelihood = zeros(length(theta))	# Total gradient (initialize it)
  for i = 1:n
    grad!(vec(data[:,i]), theta, lambda, dLikelihood)
  end
  return dLikelihood
end

# Functions to check gradient and Hessian vector products
# This function computes finite difference gradient
function check_grad(loss::Function, grad!::Function, data::Matrix, theta::Vector, lambda::Float64, eps::Float64)
  n = size(data, 2)		# Number of datapoints
  p = length(theta)		# Number of parameters
  delta = randn(p)		# Compute a unit-normal random vector to find a distance to travel
  # Double-sided FD gradient in direction delta
  deriv_approx = (loss(data, theta+eps*delta, lambda) - loss(data, theta-eps*delta, lambda))/(2*eps)
  # Analytic gradient
  grad = zeros(p)
  for i = 1:n
    grad!(vec(data[:,i]), theta, lambda, grad)
  end
  # Need to adjust analytic gradient to account for 1) the distance delta from above and 2) the number of data points (gradient is not adjusted by n)
  deriv_exact = dot(grad/n, delta)
  # Return the norm of the difference
  return norm(deriv_approx - deriv_exact)
end

# Computes Hessian Vector Product for a single datapoint
function hvp!(datapoint::Vector, theta::Vector, v::Vector, lambda::Float64, result::Vector)
  # Split individual data point
  x = datapoint[1:end-1]
  y = datapoint[end]
  # This is used to account for the piecewise nature of the max function
  alpha = dot(x, theta)
  if y * alpha < 1.0		# The case where max > 0
  	# Here HVP = H(x)*v = (y*X)*(y*X) + L*v
    BLAS.axpy!(1.0, y*y*x * dot(x, v) + lambda * v, result)
  else
    BLAS.axpy!(1.0, lambda * v, result)
  end
end

# Please check with check_gradient if the gradient function is correct before you use this function
# Uses double sided finite difference to check the HVP
function check_hvp(grad!::Function, hvp!::Function, data::Matrix, theta::Vector, lambda::Float64, eps::Float64)
  n = size(data, 2)			# Number of datapoints
  p = length(theta)			# Number of coefficients
  delta = randn(p)			# Compute a unit-normal random vector to find a distance to travel
  hvp_exact = zeros(p)		
  # Computes the exact HVP
  for i = 1:n
    hvp!(vec(data[:,i]), theta, delta, lambda, hvp_exact)
  end
  # Double-sided FD gradient in direction delta
  # H(x)*v ~ [g(x+r*v)-g(x-r*v)]/(2*r)
  hvp_approx = (full_grad(grad!, data, theta+eps*delta, lambda) - full_grad(grad!, data, theta-eps*delta, lambda))/(2*eps)
  return norm(hvp_exact - hvp_approx)
end

# New checking function to see if grad! can be used to approximate hvp! (it can)
function check_hvp_grad(grad!::Function, hvp!::Function, data::Matrix, theta::Vector, lambda::Float64, eps::Float64)
  n = size(data, 2)
  p = length(theta)
  summed_error = 0
  delta = randn(p)
  for i = 1:n
  	# Reset and compute hvp operation
    hvp_exact = zeros(p)
    hvp!(vec(data[:,i]), theta, delta, lambda, hvp_exact)
    # Reset and compute FD hvp
    fd_hvp = zeros(p)
    grad!(vec(data[:,i]), theta-eps*delta, lambda, fd_hvp)
    fd_hvp = -1*fd_hvp
    grad!(vec(data[:,i]), theta+eps*delta, lambda, fd_hvp)
    fd_hvp = fd_hvp/(2*eps)
    summed_error = summed_error + norm(hvp_exact-fd_hvp)
  end
  return summed_error
end

function ad_hoc_check()
  println("Checking gradients and hessian vector products")

  # This creates a 100-point random data matrix with 99 independent variables; set the last row of the data matrix (the class) to a random value 
  data = randn(100, 100)
  data[100,:] = 2.0*randbool(100) - 1

  # Set the theta value, then compute the gradient and hessian norms
  theta = randn(99)
  gnorm = norm(check_grad(loss, grad!, data, theta, 1e-3, 1e-4))
  hnorm = norm(check_hvp(grad!, hvp!, data, theta, 1e-3, 1e-4))
  
  # Ensure that the gradient and hessian norms are below expected error bounds
  println("Gradient Norm: ", gnorm)
  @assert gnorm <= 1e-4
  println("Hessian Vector Product Norm: ", hnorm)
  @assert hnorm <= 1e-4
  
  # Check to see if grad! can be used to approximate hvp!
  inorm = check_hvp_grad(grad!, hvp!, data, theta, 1e-3, 1e-4)  
  println("HVP Gradient Test: ", inorm)
end

# Checks to see if the function definitions are correct
ad_hoc_check()
