using Optim
using DataFrames

# Featurize categorical variables
function onehot(column::DataArray)
  categories = unique(column)
  p = length(categories)
  n = size(column, 1)
  result = zeros(n, p)
  for (i, x) in enumerate(column)
    result[i,:] = eye(p)[:,findfirst(categories, x)]
  end
  return result
end

# Featurize categorical variables and unit-normalize numeric ones
function featurize(df)
  X = zeros(size(df, 1), 0)
  for feature in names(df)
    column = df[feature]
    if eltype(column) <: Number
      X = [X (column - mean(column))/std(column)]
    else
      X = [X onehot(column)]
    end
  end
  return X
end

# Turns a vector into a 2 class system
function class(column)
  # Create an output vector the size of the input one
  result = zeros(length(column))
  # Find (and ensure) that there are only 2 categories in the input vector
  categories = unique(column)
  @assert length(categories) == 2
  # Traverses through vector and 
  for (i, x) in enumerate(column)
    if findfirst(categories, x) == 1
      result[i] = -1
    else
      @assert findfirst(categories, x) == 2
      result[i] = 1
    end
  end
  return result
end

##################
### Parse Data ###
##################
df = readtable("adult.data")	# Read in the data (uses function from DataFrames package)
df = df[1:10000,:]				# Truncate to the first 10,000 entries
X = featurize(df[1:end-1])		# Featurize all but the last column of the data from above
y = class(df[end])				# Turn all the last column into class labels
data = cat(2, X, y)				# Create the 'data' matrix by pasting the X and Y values together (cbind in R)
data = data'					# Transpose the data (for what reason I do not know)

##################
### Test Optim ###
##################
# Set the regularization value (lambda) and create a random seed for the parameters (theta)
lambda = 0.001
theta = randn(size(X, 2))		# Parameter vector equal to the number of columns of the data matrix

# Create wrapper functions so that the generic lbfgs optimizer from Optim can use the svm functions (which have already been loaded at this point)
# Likelihood function wrapper
function func(x::Vector)
  return loss(data, x, lambda)
end
# Gradient function wrapper
function gradient!(x::Vector, storage::Vector)
  # Need to copy the return value of full_grad (a vector) into the function-provided vector (storage) for the canned optimize routine below
  copy!(storage, full_grad(grad!, data, x, lambda))
end

# This uses standard L-BFGS to compute the minimum using the optimize package
res = optimize(func, gradient!, theta, method = :l_bfgs, show_trace=true, grtol=1e-5)
# Create some sort of dictionary information to store the optimal (theta star) in some HDF5 file
dataset_info = ["name" => "UCI", "theta_star" => res.minimum]