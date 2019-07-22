""" # Numerical Differential Equation Package
    Author: Ken Sible | Last Modified: July 22, 2019
    Exported Function(s)  : odeint, pdeint
    Exported Structure(s) : FiniteDiff
    Internal Function(s)  : RK2, RK4, RKF45
"""
module NumericDE
export FiniteDiff, odeint, pdeint

struct FiniteDiff
    d::Integer
    p::Integer
    s::Vector{<:Integer}
    c::Vector{<:Real}
    h::Real
    function FiniteDiff(d, p, s, h)
        N, F = length(s), factorial(d)
        D = [d == (n - 1) ? 1 : 0 for n = 1:N] .* F
        M = [s[m]^(n - 1) for n = 1:N, m = 1:N]
        new(d, p, s, M\D, h)
    end
end

function (FD::FiniteDiff)(f::Vector, i::Int)
    approx = 0
    for n = 1:length(FD.s)
        approx += FD.c[n]*f[i + FD.s[n]]
    end
    return approx/FD.h^FD.d
end

function (FD::FiniteDiff)(f::Function, x::Real)
    approx = 0
    for n = 1:length(FD.s)
        approx += FD.c[n]*f(x + FD.s[n]*h)
    end
    return approx/FD.h^FD.d
end

function RK2(t::Real, y::Array{<:Real}, f::Function, h::Real, args...)
    p = [0, 2/3]
    q = [0 0; 2/3 0]
    a = [1/4; 3/4]
    
    k1 = f(t, y, args...)
    k2 = f(t + p[2]*h, y + h*(q[2, 1]*k1), args...)
    return y + h*(a[1]*k1 + a[2]*k2)
end

function RK4(t::Real, y::Array{<:Real}, f::Function, h::Real, args...)
    p = [0, 1/2, 1/2, 1]
    q = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0]
    a = [1/6, 1/3, 1/3, 1/6]
    
    k1 = f(t, y, args...)
    k2 = f(t + p[2]*h, y + h*(q[2, 1]*k1), args...)
    k3 = f(t + p[3]*h, y + h*(q[3, 1]*k1 + q[3, 2]*k2), args...)
    k4 = f(t + p[4]*h, y + h*(q[4, 1]*k1 + q[4, 2]*k2 + q[4, 3]*k3), args...)
    return y + h*(a[1]*k1 + a[2]*k2 + a[3]*k3 + a[4]*k4)
end

function RKF45(t::Real, y::Array{<:Real}, f::Function, h::Real, args...)
    p = [0, 1/4, 3/8, 12/13, 1, 1/2]
    q = [0 0 0 0 0 0; 1/4 0 0 0 0 0; 3/32 9/32 0 0 0 0; 1932/2197 -7200/2197 7296/2197 0 0 0;
         439/216 -8 3680/513 -845/4104 0 0; -8/27 2 -3544/2565 1859/4104 -11/40 0]
    a = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55]
    b = [25/216, 0, 1408/2565, 2197/4104, -1/5, 0]
    
    k1 = f(t, y, args...)
    k2 = f(t + p[2]*h, y + h*(q[2, 1]*k1), args...)
    k3 = f(t + p[3]*h, y + h*(q[3, 1]*k1 + q[3, 2]*k2), args...)
    k4 = f(t + p[4]*h, y + h*(q[4, 1]*k1 + q[4, 2]*k2 + q[4, 3]*k3), args...)
    k5 = f(t + p[5]*h, y + h*(q[5, 1]*k1 + q[5, 2]*k2 + q[5, 3]*k3 + q[5, 4]*k4), args...)
    k6 = f(t + p[6]*h, y + h*(q[6, 1]*k1 + q[6, 2]*k2 + q[6, 3]*k3 + q[6, 4]*k4 + q[6, 5]*k5), args...)
    y1 = y + h*(b[1]*k1 + b[2]*k2 + b[3]*k3 + b[4]*k4 + b[5]*k5 + b[6]*k6)
    y2 = y + h*(a[1]*k1 + a[2]*k2 + a[3]*k3 + a[4]*k4 + a[5]*k5 + a[6]*k6)
    return y2, (y2 - y1)
end

function odeint(f::Function, tspan::NTuple{2, Real}, y0::Vector{<:Real}, h::Real, args...;
        tol::Real=1e-8, max_step::Real=Inf, method::Union{Function, Symbol}=:RKF45)
    if method == :RKF45
        i, t, y = 1, [tspan[1]], y0
        while t[i] < tspan[2]
            if (t[i] + h) > tspan[2]
                h = tspan[2] - t[i]
            end
            y_tmp, res = RKF45(t[i], y[:, i], f, h, args...)
            ε = norm(res)
            if ε < tol
                y = hcat(y, y_tmp)
                t = vcat(t, t[i] + h)
                i += 1
            end
            h *= 0.9(tol/ε)^(1/4)
            if h > max_step
                h = max_step
            end
        end
        y = permutedims(y)
    else
        RK = method isa Function ? method : 
            method ∈ (:RK2, :RK4) ? eval(method) : 
            error("Unsupported Runge-Kutta Method")
        t = collect(tspan[1]:h:tspan[2])
        y = zeros(length(t), length(y0))
        y[1, :] = y0 # Initial Condition
        for i = 1:(length(t) - 1)
            y[i + 1, :] = RK(t[i], y[i, :], f, h, args...)
        end
    end
    return t, y
end

function pdeint(f::Function, tspan::NTuple{2, Real}, xspan::NTuple{2, Real}, 
        u0::Array{<:Real}, ht::Real, hx::Real, args...; method::Union{Function, Symbol}=:RK4)
    RK = method isa Function ? method : 
         method ∈ (:RK2, :RK4) ? eval(method) : 
         error("Unsupported Runge-Kutta Method")
    t = collect(tspan[1]:ht:tspan[2])
    x = collect(xspan[1]:hx:xspan[2])
    Nt, Nx = length(t), length(x)
    u = u0 # Initial Condition
    for i = 1:(Nt - 1)
        u = RK(t[i], u, f, ht, args...)
    end
    return x, u
end

end