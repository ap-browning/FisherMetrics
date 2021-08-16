# Fisher metric
struct FisherMetric{N} <: AbstractMetric
    I::Function
    dom::Function
end

# Constructors
FisherMetricManifold(N::Int,I::Function,dom::Function=p->true) = MetricManifold(
    Euclidean(N),
    FisherMetric{N}(I,dom)
)
FisherMetricManifold(d,dom::Function=p->true) = MetricManifold(
    Euclidean(length(fieldnames(d))),
    FisherMetric{length(fieldnames(d))}(p -> fim(eval(Expr(:call,d,p...))),dom)
)

# Methods
function local_metric(M::MetricManifold{ℝ,<:AbstractManifold,<:FisherMetric},p,B::InducedBasis)
    metric(M).I(p)
end
function check_point(M::MetricManifold{ℝ,<:AbstractManifold,<:FisherMetric{N}}, p; kwargs...) where {N}
    (size(p)) == (N,) || return DomainError(size(p),"The size of $p is not $((N,)).")
    metric(M).dom(p) || return DomainError(p,"p is not a valid point in the metric domain")
    return nothing
end
function log!(M::MetricManifold{ℝ,<:AbstractManifold,<:FisherMetric},X,p,q)
    # Use shooting
    X .= optimize(X -> norm(q - geodesic(M,p,X)(1.0)), q - p * norm(exp(M,p,X) - p) / norm(q - p)).minimizer
end
function inner(M::MetricManifold{ℝ,<:AbstractManifold,<:FisherMetric},p,X,Y)
    X' * metric(M).I(p) * Y
end

# Default basis (to make things easy)
function default_basis(M,p)
    induced_basis(M,Manifolds.RetractionAtlas(),p,TangentSpace)
end