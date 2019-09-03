#= Product Logarithm Package =#
# Ken Sible | July 5, 2018

""" # Product Logarithm Package
    Product Logarithm (Lambert W Function, Omega Function)
    Function(s): lambertw(z, [tol])
"""
module ProductLog
export lambertw

""" lambertw(z, [tol])

    Solves z = w * exp(w) for w using the Newton-Raphson method
    with the specified tolerance, assuming the principle branch.
"""
function lambertw(z::Number; tol::Real=1e-8)
    z < -1/ℯ && throw(DomainError())
    f(w) = w * exp(w) - z
    df(w) = exp(w) * (w + 1)
    w = z < 1 ? z : log(z) # Initial Estimate
    f(w) == 0 && return w; i = 0
    while abs(f(w)) > tol
        w += -f(w)/df(w)
        (i += 1) > 1000 && break
    end
    return w
end
end