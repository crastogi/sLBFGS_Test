using Optim
using DataFrames

function featurize(df)
  X = zeros(size(df, 1), 0)
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

function class(column)
  result = zeros(length(column))
  categories = unique(column)
  @assert length(categories) == 2
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

dataset = :adult

# Reads in the data
df = readtable("adult.data")
X = featurize(df[1:end-1])
y = class(df[end])

X = X[1:10000,:]
y = y[1:10000,:]

data = cat(2, X, y)

num_features = size(X, 2)

data = data'

# Theta is random seed?
theta = randn(num_features)
# This is just setting the regularization
lambda = 0.001

function func(x::Vector)
  return loss(data, x, lambda)
end

function gradient!(x::Vector, storage::Vector)
  copy!(storage, full_grad(grad!, data, x, lambda))
end

res = optimize(func, gradient!, theta, method = :l_bfgs, show_trace=true, grtol=1e-5)

dataset_info = ["name" => "UCI", "theta_star" => res.minimum]
