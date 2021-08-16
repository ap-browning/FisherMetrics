module FisherMetrics

using Manifolds, ManifoldsBase
using Distributions
using Optim
using OrdinaryDiffEq
using BlockDiagonals

import Manifolds: local_metric
import ManifoldsBase: inner, log, log!

# Fisher information matrix for common distributions
include("fisher_information.jl")

# Fisher metric
include("fisher_metrics.jl")

export fim, FisherMetricManifold, default_basis

end