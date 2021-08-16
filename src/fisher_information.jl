using Distributions
using SpecialFunctions

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


function fim(d::Product)
    BlockDiagonal([fim(marginal) for marginal in d.v])
end