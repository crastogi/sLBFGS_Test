@everywhere using DataStructures
using Base.Test
using ArrayViews

#############################################
### Functions for storing Hessian Updates ###
#############################################
# This defines a data type that stores the history of s, y, and j, and memory size M (which is called tau in this code)
type InverseHessian
  # This is the structure that stores s, y, and rhos
  syrhos::Deque{(Vector{Float64}, Vector{Float64}, Float64)}
  # This stores memory depth M
  tau::Int
end

# This function is the Constructor? of the InverseHessian datatype. Requires the memory size to initialize 
function InverseHessian(tau::Int)
  return InverseHessian(Deque{(Vector{Float64}, Vector{Float64}, Float64)}(), tau)
end

# This function stores the current values of s and y, and computes and stores the current value of rho
function add!(self::InverseHessian, s::Vector{Float64}, y::Vector{Float64})
  rho = 1./dot(s, y)						# This is the rho in lbfgs update: rho_j = 1/s_j^T*y_j
  # This adds the current values of s, y, and rho to the InverseHessian data structure
  push!(self.syrhos, (s, y, rho))
  # In case the memory depth has exceeded
  if length(self.syrhos) > self.tau
  	# Remove the earliest (in this case the first) item in the collection
    shift!(self.syrhos)
  end
end

########################################
### Functions for Two-Loop Recursion ###
########################################
# This follows Algorithm 9.1 (Two Loop Recursion) in Nocedal Numerical Optimization.
# It requires the InverseHessian Object and the current gradient, and the pseudo-hessian regularizer delta
function mvp(H_r::InverseHessian, v_t::Vector{Float64}, delta::Float64)
  # Ensure that the Hessian update has some history in it
  m = length(H_r.syrhos)			# Memory size
  @assert m > 0
  # Clone the input gradient
  q = copy(v_t)
  # Create a vector of alphas with the right memory size
  alphas = zeros(m)
  # The first loop (starts from the latest entry and goes to the first)
  for (i, (s, y, rho)) in enumerate(reverse(collect(H_r.syrhos)))
  	# This computes and stores alpha_i = rho_i*s_i*q
    alpha = rho * dot(s, q)
    alphas[m - i + 1] = alpha			# Required for proper indexing
    # Update q with q = q - alpha_i*y_i
    q -= alpha * y
  end
  # Now attempt to compute r. We need this for the last stored iterate
  (s, y, rho) = back(H_r.syrhos)
  # gamma_k = s_k*y_k/(y_k*y_k)
  gamma = dot(s, y)/dot(y, y)
  # r = gamma_k*q/(1 + delta*gamma_k); the denominator includes the pseudo-hessian regularization parameter delta
  # NOTE: There is no need to multiply by I here, as that will anyway produce a dot product 
  r = gamma*q/(1 + delta*gamma)
  # Second loop (goes in reverse, starting from the first entry)
  for (i, (s, y, rho)) in enumerate(collect(H_r.syrhos))
  	# beta = rho_i*y_i*r
    beta = rho * dot(y, r)
    # r = r + s_i*(alpha_i-beta)
    r += s * (alphas[i] - beta)
  end
  return r
end

# Alternate function to test the two-loop recursion formula (Equation 4)
function test_mvp(H_r::InverseHessian, v_t::Vector{Float64}, delta::Float64)
  # Find the latest set of s, y, and rho
  syrhos = collect(H_r.syrhos)
  (s, y, rho) = syrhos[end]
  # Find the total number of parameters
  p = length(s)
  # Compute gamma, which is s_r*y_r/(y_r*y_r)
  gamma = dot(s, y)/dot(y, y)
  # Compute H_r^(r-M) = gamma*I. Here I is the identity matrix. Below, (1+delta*gamma) is 
  # the correction factor incorporating the regularization of the pseudohessian
  H = gamma / (1 + delta*gamma) * eye(p)
  # Loop over all history stored in the InverseHessian object
  for (s, y, rho) in syrhos
  	# V_j = (I-rho_j*s_j*y_j)
    V = eye(p) - rho * y * s'
    # H_j = V'_j*H_(j-1)*V_j + rho_j*s_j*s_j
    H = V' * H * V + rho * s * s'
  end
  # Compute H_r*v_t
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
# S is the number of epochs
function slbfgs(loss::Function, grad!::Function, hvp!::Function, data::Matrix, theta::Vector, config::Dict, diagnostics::HDF5Group)
  tau = config["memory_size"]						# This is memory size M
  delta = config["delta"]							# This is an optional Inverse Hessian regularization parameter 
  eta = config["stepsize"]							# Fixed step size parameter eta in x_t+1 = x_t - eta*H_r*v_t
  S = config["epochs"]								# Number of epochs
  verbose = true
  println("===== running stochastic lbfgs =====\n")
  println("iter | lbfgs error")
  n = size(data, 1)
  m = 2*n
  p = length(theta)
  H = InverseHessian(tau)
  thetaold = copy(theta)
  dtheta = zeros(p)
  mu = zeros(p)
  z = zeros(p)
  zold = zeros(p)
  P = zeros(p)

  t = -1 # number of correction pairs currently computed
  stheta = zeros(p)
  sthetaold = zeros(p)

  hdf_loss = g_create(diagnostics, "loss")
  hdf_theta = g_create(diagnostics, "theta")

  if config["use_gradient_step"]
    println("warning: using gradient steps")
  end

  # Loop over epochs (fixed)
  for s in 1:S
    # Print epoch number and likelihood at epoch as the result of the iteration 
    if verbose
      l = loss(data, theta, config["lambda"])
      @printf "%4d | %9.4E \n" s l
      hdf_loss[@sprintf "%04d" s] = l
      hdf_theta[@sprintf "%04d" s] = theta
    end
	
	# Need to clone theta (equivalent to x_0 = w_k, line 5)
    ctheta = copy(theta)
    
	# Next, compute full gradient for variance reduction (line 4)
    # \tilde \mu = 1/n \sum \nabla f_i(\tilde theta)
    fill!(mu, 0.0)				# Fill the variance-reduced gradient with zeros
    for i = 1:n					# Compute variance-reduced gradient 
      grad!(vec(data[:,i]), ctheta, config["lambda"], mu)
    end
    mu /= n						# Normalize to unit likelihood
    
    # ONLY for first epoch
    if s == 1
      # Do a gradient step to get LBFGS started
      copy!(P, -mu/norm(mu))
    end
    copy!(theta, ctheta)
    for k in 1:m
      index = rand(1:n)
      copy!(zold, z)
      fill!(z, 0.0)
      # z_k = \nabla f_{i_k}(\theta) - \nabla f_{i_k}(\tilde\theta) + \tilde \mu
      fill!(dtheta, 0.0)
      grad!(vec(data[:,index]), theta, config["lambda"], dtheta)
      z += dtheta
      fill!(dtheta, 0.0)
      grad!(vec(data[:,index]), ctheta, config["lambda"], dtheta)
      if config["reduce_grad_var"]
        z -= dtheta
        z += mu
      end

      stheta += theta

      if k > 1 # neccessary?
        if mod(k, config["hess_period"]) == 0
          stheta /= config["hess_period"]
          t += 1
        end
        if mod(k, config["hess_period"]) == 0 && t >= 1
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
                r = z - zold + delta*svec
            end
            add!(H, svec, r)
            copy!(sthetaold, stheta)
        end
        if t > 1
          P = -eta * mvp(H, z, delta)
        else
          P = -eta*z
        end
      end
      # take step
      copy!(thetaold, theta)
      if config["use_gradient_step"]
        theta -= eta*z
      else
        theta += P
      end
    end
  end
  return theta
end
