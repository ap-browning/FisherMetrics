# FisherMetrics

<!-- [![Build Status](https://github.com/ap-browning/FisherMetrics.jl/workflows/CI/badge.svg)](https://github.com/ap-browning/FisherMetrics.jl/actions)
[![Coverage](https://codecov.io/gh/ap-browning/FisherMetrics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ap-browning/FisherMetrics.jl) -->

**! Work in progress !**

This is a Julia package that extends [`Manifolds.jl`](https://github.com/JuliaManifolds/Manifolds.jl) to statistical manifolds in [information geometry](https://en.wikipedia.org/wiki/Information_geometry). Essentially, a Euclidean manifold where the local metric is given by the Fisher information matrix.

Example ways to create a manifold
```
    M = FisherMetricManifold(2,p -> fim(Normal(p...)))
    M = FisherMetricManifold(Beta)
```

Once created, one can calculate the distance between two points, `p` and `q`, in addition to calculating the shortest geodesic connecting `p` and `q` (amoung several other functions provided by `Manifolds.jl`).
```
    using Manifolds
    using FisherMetrics
    using Distributions

    M = FisherMetricManifold(Normal)

    p = [1.0,2.0]
    q = [5.0,4.0]

    B = default_basis(M,p)

    dist = distance(M,p,q)
    geo = shortest_geodesic(M,p,q)
```

## Todo:
1. Provide test functions.
2. Provide documentation.
3. Create a more robust boundary value problem solver to compute geodesics (unstable for points that are far away).