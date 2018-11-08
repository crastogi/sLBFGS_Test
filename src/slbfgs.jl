@everywhere using DataStructures
using Base.Test
using ArrayViews

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

#########################
### sL-BFGS Main Loop ###
#########################
function slbfgs(loss::Function, grad!::Function, hvp!::Function, data::Matrix, theta::Vector, config::Dict, diagnostics::HDF5Group)
  ### FIXED CONSTANTS ###
  M = config["memory_size"]							# This is memory size M
  L = config["hess_period"]							# Number of iterations before the inverse hessian is updated
  delta = config["delta"]							# This is an optional Inverse Hessian regularization parameter 
  eta = config["stepsize"]							# Fixed step size parameter eta in x_t+1 = x_t - eta*H_r*v_t
  max_epoch = config["epochs"]						# Maximum number of epochs

  n = size(data, 1)
  m = 2*n
  p = length(theta)
  
  ### RUNTIME VARIABLES ###
  v_t = zeros(p)
  v_t_prev = zeros(p)
  mu_k = zeros(p)
  IH_hist = InverseHessian(M)
  
  thetaold = copy(theta)
  dtheta = zeros(p)
  P = zeros(p)
  
  t = -1 # number of correction pairs currently computed
  stheta = zeros(p)
  sthetaold = zeros(p)

  hdf_loss = g_create(diagnostics, "loss")
  hdf_theta = g_create(diagnostics, "theta")

  # Header for verbose output
  println("===== running stochastic lbfgs =====\n")
  println("iter | lbfgs error")
  # Loop over epochs (fixed termination at the end of max_epoch)
  for s in 1:max_epoch
    # Print epoch number and likelihood at epoch as the result of the iteration 
    l = loss(data, theta, config["lambda"])
    @printf "%4d | %9.4E \n" s l
    hdf_loss[@sprintf "%04d" s] = l
    hdf_theta[@sprintf "%04d" s] = theta
	
	# Need to clone theta (equivalent to x_0 = w_k, line 5)
    ctheta = copy(theta)
    
	# Next, compute full gradient for variance reduction (line 4)
	# WARNING: It is possible that the variance reduction step here does not compute the gradient over the whole dataset
    # \tilde \mu_k = 1/n \sum \nabla f_i(\tilde theta)
    fill!(mu_k, 0.0)			# Fill the variance-reduced gradient with zeros
	# Here n is the number of parameters + 1!!!!    
    for i = 1:n					# Compute variance-reduced gradient 
      grad!(vec(data[:,i]), ctheta, config["lambda"], mu_k)
    end
    mu_k /= n					# Normalize to unit likelihood
    
    # ONLY for first epoch
    if s == 1
      # Do a gradient step to get LBFGS started
      copy!(P, -mu_k/norm(mu_k))
    end
    copy!(theta, ctheta)
    for k in 1:m
      index = rand(1:n)
      copy!(v_t_prev, v_t)
      fill!(v_t, 0.0)
      # z_k = \nabla f_{i_k}(\theta) - \nabla f_{i_k}(\tilde\theta) + \tilde \mu_k
      fill!(dtheta, 0.0)
      grad!(vec(data[:,index]), theta, config["lambda"], dtheta)
      v_t += dtheta
      fill!(dtheta, 0.0)
      grad!(vec(data[:,index]), ctheta, config["lambda"], dtheta)

	  # To compute the reduced variance gradient:
      v_t += mu_k - dtheta

      stheta += theta

      if k > 1 # neccessary?
        if mod(k, L) == 0
          stheta /= L
          t += 1
        end
        if mod(k, L) == 0 && t >= 1
          svec = stheta - sthetaold
          r = zeros(p)
          if config["use_hvp"]
            for i = 1:config["hess_samples"]
              index = rand(1:n)
              hvp!(vec(data[:,index]), stheta, svec, config["lambda"], r)
            end
            r /= config["hess_samples"]
          else
                # delta a regularization parameter
                r = v_t - v_t_prev + delta*svec
            end
            add!(IH_hist, svec, r)
            copy!(sthetaold, stheta)
        end
        if t > 1
          P = -eta * mvp(IH_hist, v_t, delta)
        else
          P = -eta*v_t
        end
      end
      
      # take step
      copy!(thetaold, theta)
      if config["use_gradient_step"]
        theta -= eta*v_t
      else
        theta += P
      end
    end
  end
  return theta
end