using Random

n = 10  # Replace with the desired number of random numbers

# Generate n normal random numbers
a = randn(n)
b = randn(n)
X = [a b]'
println(size(X))
fit_Z = fit_mle(MvNormal,X)

X = ones(2,n)
fit_Z = fit_mle(MvNormal,X)