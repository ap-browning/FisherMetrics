using Distributions
using SpecialFunctions
using ForwardDiff

struct TransformedDistribution <: ContinuousUnivariateDistribution
    θ::Vector
    d::Distribution
    f::Function
    TransformedDistribution(θ::Vector,d,f::Function) = new(θ,d(f(θ)...),f)
end
Distributions.params(d::TransformedDistribution) = d.θ
Base.rand(d::TransformedDistribution) = rand(d.d)
Base.rand(d::TransformedDistribution,n::Int) = rand(d.d,n)
MultiplicativeNormal(μ::Number,ν::Number) = TransformedDistribution([μ,ν],Normal,θ -> [θ[1];θ[1]*θ[2]])

function fim(d::Normal)
    μ,σ = params(d)
    [1 / σ^2 0.0; 0.0 2/σ^2]
end
function fim(d::Beta)
    α,β = params(d)
    [trigamma(α) - trigamma(α + β)  -trigamma(α + β);
     -trigamma(α + β)   trigamma(β) - trigamma(α + β)]
end
function fim(d::Gamma)
    α,θ = params(d);
    [polygamma(1,α) θ; θ  α*θ^2]
end
function fim(d::TransformedDistribution)
    J = ForwardDiff.jacobian(d.f,d.θ)
    J' * fim(d.d) * J
end

function fim(d::Product)
    BlockDiagonal([fim(marginal) for marginal in d.v])
end