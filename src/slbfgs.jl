@everywhere using DataStructures
using Base.Test
using ArrayViews
using StatsBase

#############################################
### Functions for storing Hessian Updates ###
#############################################
# This defines a data type that stores the history of s, y, and j, and memory size M
type InverseHessian
  # This is the structure that stores s, y, and rhos
  syrhos::Deque{(Vector{Float64}, Vector{Float64}, Float64)}
  # This stores memory depth M
  M::Int
end

# This function is the Constructor of InverseHessian and requires memory size M 
function InverseHessian(M::Int)
  return InverseHessian(Deque{(Vector{Float64}, Vector{Float64}, Float64)}(), M)
end

# This function stores the current values of s and y, and computes and stores the current value of rho
function add!(IH_hist::InverseHessian, s::Vector{Float64}, y::Vector{Float64})
  rho = 1./dot(s, y)						# This is the rho in lbfgs update: rho_j = 1/s_j^T*y_j
  # This adds the current values of s, y, and rho to the InverseHessian data structure
  push!(IH_hist.syrhos, (s, y, rho))
  # In case the memory depth has exceeded
  if length(IH_hist.syrhos) > IH_hist.M
  	# Remove the earliest (in this case the first) item in the collection
    shift!(IH_hist.syrhos)
  end
end

########################################
### Functions for Two-Loop Recursion ###
########################################
# This follows Algorithm 9.1 (Two Loop Recursion) in Nocedal Numerical Optimization.
# It requires the InverseHessian Object and the current gradient, and the pseudo-hessian regularizer delta
function mvp(IH_hist::InverseHessian, v_t::Vector{Float64}, delta::Float64)
  # Ensure that the Hessian update has some history in it
  m = length(IH_hist.syrhos)			# Memory size
  @assert m > 0
  # Clone the input gradient
  q = copy(v_t)
  # Create a vector of alphas with the right memory size
  alphas = zeros(m)
  # The first loop (starts from the latest entry and goes to the first)
  for (i, (s, y, rho)) in enumerate(reverse(collect(IH_hist.syrhos)))
  	# This computes and stores alpha_i = rho_i*s_i*q
    alpha = rho * dot(s, q)
    alphas[m - i + 1] = alpha			# Required for proper indexing
    # Update q with q = q - alpha_i*y_i
    q -= alpha * y
  end
  # Now attempt to compute r. We need this for the last stored iterate
  (s, y, rho) = back(IH_hist.syrhos)
  # gamma_k = s_k*y_k/(y_k*y_k)
  gamma = dot(s, y)/dot(y, y)
  # r = gamma_k*q/(1 + delta*gamma_k); the denominator includes the pseudo-hessian regularization parameter delta
  # NOTE: There is no need to multiply by I here, as that will anyway produce a dot product 
  r = gamma*q/(1 + delta*gamma)
  # Second loop (goes in reverse, starting from the first entry)
  for (i, (s, y, rho)) in enumerate(collect(IH_hist.syrhos))
  	# beta = rho_i*y_i*r
    beta = rho * dot(y, r)
    # r = r + s_i*(alpha_i-beta)
    r += s * (alphas[i] - beta)
  end
  return r
end

# Alternate function to test the two-loop recursion formula (Equation 4)
function test_mvp(IH_hist::InverseHessian, v_t::Vector{Float64}, delta::Float64)
  # Find the latest set of s, y, and rho
  syrhos = collect(IH_hist.syrhos)
  (s, y, rho) = syrhos[end]
  # Find the total number of parameters
  p = length(s)
  # Compute gamma, which is s_r*y_r/(y_r*y_r)
  gamma = dot(s, y)/dot(y, y)
  # Compute IH_hist^(r-M) = gamma*I. Here I is the identity matrix. Below, (1+delta*gamma) is 
  # the correction factor incorporating the regularization of the pseudohessian
  H = gamma / (1 + delta*gamma) * eye(p)
  # Loop over all history stored in the InverseHessian object
  for (s, y, rho) in syrhos
  	# V_j = (I-rho_j*s_j*y_j)
    V = eye(p) - rho * y * s'
    # H_j = V'_j*H_(j-1)*V_j + rho_j*s_j*s_j
    H = V' * H * V + rho * s * s'
  end
  # Compute IH_hist*v_t
  return H*v_t
end

# Checks to see the following functionality: 1) that the InverseHessian object functions 
# as needed; 2) that the 'mvp' function behaves normally using the test_mvp function, 
# therefore verifying that the 2 loop recursion works
function test_inverse_hessian()
  # Create object to store inverse hessian history with memory depth = 3
  self = InverseHessian(3)
  # Populate it for 4 intervals with random vectors. The first entry should be removed
  for i = 1:4
    s = rand(2)
    y = rand(2)
    add!(self, s, y)
  end
  # Create a random gradient vector
  g = rand(2)
  # Compute inverse hessian using two methods: 2 loop recursion and the test function above
  first = mvp(self, g, 0.0)
  second = test_mvp(self, g, 0.0)
  @test norm(first - second) <= 1e-10
end
test_inverse_hessian()

#### FD HVP replacement
function slbfgs_hvp!(datapoint::Vector, theta::Vector, v::Vector, lambda::Float64, result::Vector)
  eps = 1e-3	
  
  fd_hvp = zeros(length(theta))
  grad!(datapoint, theta-eps*v, lambda, fd_hvp)
  fd_hvp = -1*fd_hvp
  grad!(datapoint, theta+eps*v, lambda, fd_hvp)
  fd_hvp = fd_hvp/(2*eps)
  
  BLAS.axpy!(1.0, fd_hvp, result)
end

#########################
### sL-BFGS Main Loop ###
#########################
function slbfgs(loss::Function, grad!::Function, hvp!::Function, data::Matrix, w_0::Vector, config::Dict, diagnostics::HDF5Group)
  ### FIXED CONSTANTS ###
  L = config["hess_period"]					# Number of iterations before the inverse hessian is updated
  M = config["memory_size"]					# This is memory size M
  d = length(w_0)							# Dimensionality
  N = size(data, 2)							# Number of data points
  b_H = config["hess_samples"]				# Batch size for stochastic Hessian samples
  b = convert(Int, b_H/10)					# Batch size for stochastic gradient samples
  m = convert(Int, ceil(N/b))				# Number of iterations to run per epoch
  eta = config["stepsize"]					# Fixed step size parameter eta in x_t+1 = x_t - eta*H_r*v_t
  delta = config["delta"]					# This is an optional Inverse Hessian regularization parameter 
  max_epoch = config["epochs"]				# Maximum number of epochs
    
  ### RUNTIME VARIABLES ###
  r = 0 									# Number of currently computed Hessian correction pairs
  v_t = zeros(d)							# Current iteration's variance reduced gradient estimate
  v_t_prev = zeros(d)						# Previous iteration's variance reduced gradient estimate
  mu_k = zeros(d)							# Full Gradient for Variance-Reduced Gradient
  IH_hist = InverseHessian(M)				# Object that stores the history of the inverse hessian components (s, y, and rho)
  u_r = zeros(d)							# Average of path travelled in an inverse hessian update 
  u_r_prev = zeros(d)  						# Average of path travelled in the previous inverse hessian update
  grad_f_xt = zeros(d)						# Component of variance reduced gradient v_t
  grad_f_wk = zeros(d)						# Component of variance reduced gradient v_t
  w_k = copy(w_0)
  x_t_hist = zeros(m, d)	 

  # Header for verbose output
  println("===== running stochastic lbfgs =====\n")
  println("iter | lbfgs error")
  # Set up HDF5 output
  hdf_loss = g_create(diagnostics, "loss")
  hdf_theta = g_create(diagnostics, "theta")
  
  # Loop over epochs (fixed termination at the end of max_epoch; line 3)
  for k in 1:max_epoch
    # Print epoch number and likelihood at epoch as the result of the iteration 
    l = loss(data, w_k, config["lambda"])
    @printf "%4d | %9.4E \n" k l
    hdf_loss[@sprintf "%04d" k] = l
    hdf_theta[@sprintf "%04d" k] = w_k
   
	# Next, compute full gradient for variance reduction (line 4)
    # \mu_k = 1/N \sum_i^N \nabla f_i(w_k)
    fill!(mu_k, 0.0)						# Fill the variance-reduced gradient with zeros  
    for i = 1:N								# Compute variance-reduced gradient 
      grad!(vec(data[:,i]), w_k, config["lambda"], mu_k)
    end
    mu_k /= N								# Normalize to unit likelihood
	# Need to clone w_k (equivalent to x_0 = w_k, line 5)
    x_t = copy(w_k)
    
    # Perform m iterations before a full gradient computation takes place (line 6) [Ensure t index starts from 1]
    for t in 1:m
      # Compute stochastic gradient estimate: v_t = \nabla f_{i_k}(x_t) - \nabla f_{i_k}(w_k) + \mu_k
	  # Begin by sampling a minibatch for S (line 7)
      copy!(v_t_prev, v_t)
      index = sample(1:N, b, replace = false)
      # Compute stochastic gradient grad_f_xt (line 8)
      fill!(v_t, 0.0)
      fill!(grad_f_xt, 0.0)
      fill!(grad_f_wk, 0.0)
      for i in index
		grad!(vec(data[:,i]), x_t, config["lambda"], grad_f_xt)
		grad!(vec(data[:,i]), w_k, config["lambda"], grad_f_wk)
      end
      grad_f_xt /= b
      grad_f_wk /= b
	  # Compute the reduced variance gradient (line 9)
      v_t += grad_f_xt - grad_f_wk + mu_k
	  # Part of u_r update (line 13)
      u_r += x_t

      # Compute next iteration step (line 10)
      x_t_hist[t,:] = copy(x_t)				# Need to store the history of iteration steps
      if r > 0								# After one Hessian correction, need the two-loop recursion product
        x_t += -eta * mvp(IH_hist, v_t, delta)
      else									# H_0 = I, so update will be simpler
        x_t += -eta*v_t
      end

	  # Check to see if L iterations have passed (triggers hessian update; line 11)
      if mod(t, L) == 0
      	# Increment the number of Hessian correction pairs that have been computed (line 12)
		r += 1
		# Final step in computation of u_r (line 13)
        u_r /= L

        # Use HVP only; start by sampling a minibatch for T (line 14)
        index = sample(1:N, b_H, replace = false)
        # Compute s_r update (line 15)
        s_r = u_r - u_r_prev
        # Compute an estimate for y_r using HVPs [must swap with approximation!] y_r = nabla^2 f_T(u_r)*s_r (line 16)
        y_r = zeros(d)
        for i in index
          slbfgs_hvp!(vec(data[:,i]), u_r, s_r, config["lambda"], y_r)
        end
        y_r /= b_H
        # Store latest values of s_r, y_r updates
        add!(IH_hist, s_r, y_r)

        # Resetting u_r for next evaluation (line 13)
        copy!(u_r_prev, u_r)
        u_r = zeros(d)
      end
    end
    # Set the next w_k value by randomly selecting an iterate from this batch (line 18)
    x_t_hist[m,:] = copy(x_t)
    index = rand(1:m)
    w_k = copy(vec(x_t_hist[index,:]))
    w_k = copy(x_t)
  end
  println(length(w_k))
  println(w_k)
  return w_k
end