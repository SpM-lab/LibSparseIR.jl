function overlap(poly::PiecewiseLegendrePolyVector, f::F) where {F}
    xmin = poly.xmin
    xmax = poly.xmax
    pts = filter(x -> xmin ≤ x ≤ xmax, roots(poly))
    q, _ = quadgk(
        x -> poly(x) * f(x),
        unique!(sort!(vcat(pts, [xmin, xmax])));
        rtol=eps(), order=10, maxevals=10^4)
    q
end

function overlap(poly::PiecewiseLegendrePoly, f::F) where {F}
    xmin = poly.xmin
    xmax = poly.xmax
    pts = filter(x -> xmin ≤ x ≤ xmax, roots(poly))
    q, _ = quadgk(
        x -> poly(x) * f(x),
        unique!(sort!(vcat(pts, [xmin, xmax])));
        rtol=eps(), order=10, maxevals=10^4)
    return q
end

function xmin(poly::PiecewiseLegendrePoly)
    return poly.xmin
end

function xmax(poly::PiecewiseLegendrePoly)
    return poly.xmax
end

function xmin(poly::PiecewiseLegendrePolyVector)
    return poly.xmin
end

function xmax(poly::PiecewiseLegendrePolyVector)
    return poly.xmax
end

function xmin(poly::PiecewiseLegendreFTVector)
    return poly.xmin
end

function xmax(poly::PiecewiseLegendreFTVector)
    return poly.xmax
end
