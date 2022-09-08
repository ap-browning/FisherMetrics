using Distributions
using FisherMetrics
using Manifolds
using Plots

## Sharp et al. 2022 - Fig. 8b (additive noise)
K = 79.74
ν = 0.05    # σ / μ
T = [2.74,6.84,10.95]
N = 10 # Number of data points

# "True" parameters
θtrue = [0.9,0.2]

# Model as a function of θ = [r,c₀] and t
logistic(θ,t) = ((r,c₀) = θ; c₀ * K / (c₀ + (K - c₀)*exp(-r*t)))

# Observation process as a function of ξ = [μ₁,ν₂,μ₂,ν₂,μ₃,ν₃]
obs_process_mult(ξ) = Product(MultiplicativeNormal.(ξ[[1,3,5]],ξ[[2,4,6]]))

# The model acts as a transformation, θ → ξ
model_mult(θ) = [[logistic(θ,t) for t in T]; fill(ν,3)][[1,4,2,5,3,6]]

# Apply the observation process to our model
full_model_mult(θ) = TransformedDistribution(θ,obs_process_mult,model_mult;params=(:r,:c₀))

# Create metric
M_mult = FisherMetricManifold(2,θ -> N * fim(full_model_mult(θ)))

# Sample synthetic data from model (NEED TO ADD SAMPLING ALGORITHM FOR TRANSFORMED PRODUCT)
X_mult = hcat([rand(full_model_mult(θtrue).d.v[i],N) for i = 1:3]...)'

## Produce plots

# Model and (synthetic) data
p3 = plot(t -> logistic(θtrue,t),xlim=(0.0,15.0),label="Model")
scatter!(p3,T,X_mult,c=:red,label="Data")

# Scalar curvature (problems with small c₀; i.e., small mean at early time)
Sc_mult = [ricci_curvature(M_mult,[rᵢ,c₀ᵢ],default_basis(M,[rᵢ,c₀ᵢ])) for rᵢ in r, c₀ᵢ in c₀[5:end]]
p4 = heatmap(c₀[5:end],r,Sc_mult,xlabel="c₀",ylabel="r",margin=5Plots.mm)

## All plots
plot(fig,plot(p3,p4),layout=grid(2,1))