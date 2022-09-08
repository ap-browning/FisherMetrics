using Distributions
using FisherMetrics
using Manifolds
using Plots

## Sharp et al. 2022 - Fig. 8b (additive noise)
K = 79.74
σ = 2.301
T = [2.74,6.84,10.95]
N = 10 # Number of data points

# "True" parameters
θtrue = [0.9,0.2]

# Model as a function of θ = [r,c₀] and t
logistic(θ,t) = ((r,c₀) = θ; c₀ * K / (c₀ + (K - c₀)*exp(-r*t)))

# Observation process as a function of ξ = [μ₁,σ₂,μ₂,σ₂,μ₃,σ₃]
obs_process(ξ) = Product(Normal.(ξ[[1,3,5]],ξ[[2,4,6]]))
obs_process(ξ) = Product([
                    Normal(ξ[1:2]...)   # t = t₁
                    Normal(ξ[3:4]...)   # t = t₂
                    Normal(ξ[5:6]...)   # t = t₃
                 ])

# The model acts as a transformation, θ → ξ
model(θ) = [[logistic(θ,t) for t in T]; fill(σ,3)][[1,4,2,5,3,6]]

# Apply the observation process to our model
full_model(θ) = TransformedDistribution(θ,obs_process,model;params=(:r,:c₀))

# Create metric
M = FisherMetricManifold(2,θ -> N * fim(full_model(θ)))

# Sample synthetic data from model (optional)
X = rand(full_model(θtrue),N)


## Produce plots

# Model and (synthetic) data
p1 = plot(t -> logistic(θtrue,t),xlim=(0.0,15.0),label="Model")
scatter!(p1,T,X,c=:red,label="Data")

# Scalar curvature
c₀ = range(0.01,0.5,length=50)
r = range(0.7,1.3,length=51)
Sc = [ricci_curvature(M,[rᵢ,c₀ᵢ],default_basis(M,[rᵢ,c₀ᵢ])) for rᵢ in r, c₀ᵢ in c₀]
p2 = heatmap(c₀,r,Sc,xlabel="c₀",ylabel="r",margin=5Plots.mm)

# Together
fig = plot(p1,p2)