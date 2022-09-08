push!(LOAD_PATH,"/Users/Alex/Code/FisherMetrics/src")

using Revise
using Manifolds
using FisherMetrics

using Distributions

# Transform to multiplicative noise:
#   (μ,σ) -> (μ,μ*σ)

f = θ -> [θ[1];θ[1]*θ[2]]

d = TransformedDistribution(rand(2),Normal,f)

d = MultiplicativeNormal(4.0,2.0)

# Additive noise
m1 = FisherMetricManifold(Normal)

μ = range(0.1,2.0,length=20)
σ = range(0.1,2.0,length=20)
r1 = [ricci_curvature(m1,[μᵢ,σᵢ],default_basis(m1,[μᵢ,σᵢ])) for μᵢ in μ, σᵢ in σ]

# Multiplicative noise
m2 = FisherMetricManifold(2,p->fim(MultiplicativeNormal(p...)))

ν = range(0.1,2.0,length=20)    # ν = μ σ ⇒ σ = ν / μ
r2 = [ricci_curvature(m2,[μᵢ,νᵢ],default_basis(m2,[μᵢ,νᵢ])) for μᵢ in μ, νᵢ in ν]

## Plot
