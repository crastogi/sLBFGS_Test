# Loss function, gradients and hessian vector products for SVMs with quadratic loss
println("Loading SVM module")

# Data matrix rows consist of rows of independent variables with last entry being dependent variable
function loss(data::Matrix, theta::Vector, lambda::Float64)
  result = 0.0
  n = size(data, 2)
  for i = 1:n
    x = data[1:end-1,i]
    y = data[end,i]
    result += 0.5 * max(0.0, 1 - y * dot(vec(x), theta))^2/n
  end
  result += lambda/2 * dot(theta, theta)
  return result
end

function grad!(datapoint::Vector, theta::Vector, lambda::Float64, grad::Vector)
  x = datapoint[1:end-1]
  y = datapoint[end]
  alpha = dot(x, theta)
  if y * alpha < 1.0
    BLAS.axpy!(1.0, lambda * theta -  (1 - y*alpha) * y * x, grad)
  else
    BLAS.axpy!(1.0, lambda * theta, grad)
  end
end

function hvp!(datapoint::Vector, theta::Vector, v::Vector, lambda::Float64, result::Vector)
  x = datapoint[1:end-1]
  y = datapoint[end]
  alpha = dot(x, theta)
  if y * alpha < 1.0
    BLAS.axpy!(1.0, x * dot(x, v) + lambda * v, result)
  else
    BLAS.axpy!(1.0, lambda * v, result)
  end
end

# Functions to check gradient and Hessian vector products
function check_grad(loss::Function, grad!::Function, data::Matrix, theta::Vector, lambda::Float64, eps::Float64)
  n = size(data, 2)
  p = length(theta)
  delta = randn(p)
  deriv_approx = (loss(data, theta+eps*delta, lambda) - loss(data, theta-eps*delta, lambda))/(2*eps)
  grad = zeros(p)
  for i = 1:n
    grad!(vec(data[:,i]), theta, lambda, grad)
  end
  deriv_exact = dot(grad/n, delta)
  return norm(deriv_approx - deriv_exact)
end

# Compute the full gradient over the whole data matrix
function full_grad(grad!::Function, data::Matrix, theta::Vector, lambda::Float64)
  n = size(data, 2)
  result = zeros(length(theta))
  for i = 1:n
    grad!(vec(data[:,i]), theta, lambda, result)
  end
  return result
end

# plase check with check_gradient if the gradient function is correct before you use this function
function check_hvp(grad!::Function, hvp!::Function, data::Matrix, theta::Vector, lambda::Float64, eps::Float64)
  n = size(data, 2)			# Number of datapoints
  p = length(theta)			# Number of coefficients
  delta = randn(p)			
  hvp_exact = zeros(p)		
  for i = 1:n
    hvp!(vec(data[:,i]), theta, delta, lambda, hvp_exact)
  end
  hvp_approx = (full_grad(grad!, data, theta+eps*delta, lambda) - full_grad(grad!, data, theta-eps*delta, lambda))/(2*eps)
  return norm(hvp_exact - hvp_approx)
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
end

# Checks to see if the function definitions are correct
ad_hoc_check()
