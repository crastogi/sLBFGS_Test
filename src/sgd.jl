# p is number of features
# S is number of epochs
# Plain vanilla sgd (NO mini-batch, with simple annealing scheme)
function sgd(loss::Function, grad!::Function, data::Matrix, theta::Vector, config::Dict, diagnostics::HDF5Group)
  @printf("===== running the stochastic gradient method =====\n")
  n = size(data, 2)							# Size of dataset, given by the number of rows
  theta = copy(theta)						# Make a copy of theta
  dtheta = Array(Float64, length(theta))	# Create a gradient vector
  eta = config["stepsize"]
  S = config["epochs"]
  verbose = true
  T = S*n
  for t in 1:T
    eta_t = eta / sqrt(t)					# Annealing scheme for stepsize e_t = e/sqrt(t)
    fill!(dtheta, 0.0)						# Resets the gradient vector to 0
    index = rand(1:n)						# Randomly selects index of datapoint for next step
    grad!(vec(data[:,index]), theta, config["lambda"], dtheta)		# Computes gradient for selected point
    if verbose && mod(t, n) == 0			# If finished an epoch, print current likelihood to screen and diagnostic object
      l = loss(data, theta, config["lambda"])
      @printf "%4d | %9.4E \n" t/n l
      diagnostics[@sprintf "%04d" t/n] = l
    end
    theta -= eta_t * dtheta					# Update theta with simple update scheme
  end
  return theta								# Return theta (NOTE: this is NOT the same theta as input [thus no 'bang' {sgd!}])
end
