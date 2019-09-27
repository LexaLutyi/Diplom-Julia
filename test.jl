using DifferentialEquations
using LinearAlgebra

using InteractiveUtils


const agents_amt = 10
const order = 2
const dims_amt = 2
A = randn(dims_amt * agents_amt, dims_amt * agents_amt) - 3I
const tau = 0.0

function f_test(du,u,h,p,t)
    du[:] =p* [zeros(dims_amt * agents_amt,dims_amt * agents_amt) I
               zeros(dims_amt * agents_amt,dims_amt * agents_amt) A] * h(p, t - tau)
end

h(p,t) = randn(order * dims_amt * agents_amt)
u0 = h(nothing, 0)

tspan = (0.0, 4.0)
prob = DDEProblem(f_test,u0,h,tspan,1)
sol = solve(prob)

#]add Plots # You need to install Plots.jl before your first time using it!
using Plots
#plotly() # You can optionally choose a plotting backend

#plot(sol, vars=(0,1))
#plot!(sol, vars=(0,2))
plt = plot(sol)
display(plt)
#println(sol)
"""
@gif for t in range(0.0, length=100, stop=4.0)
    a = sol(t)
    b = reshape(a, (2,agents_amt))
    c,d = [b[1,:], b[2,:]]
    plot(c,d,
        seriestype=:scatter,
        title="My Scatter Plot",
        xlims=(-10,10),
        ylims=(-10,10)
    )
end
"""

"""
@gif for t in range(0.0, length=100, stop=4.0)
    a = sol(t)
    b = reshape(a, (agents_amt, 4))
    c,d,v,w = b[:,1], b[:,2], b[:,3], b[:,4]
    for (i, q) in enumerate(zip(v,w))
        v[i], w[i] = normalize([q[1], q[2]])
    end
    quiver(c,d,
        quiver=(v,w),
        xlims=(-10,10),
        ylims=(-10,10),
        arrow=arrow(5.0)
    )
end
"""

"""
@gif for t in range(0.0, length=100, stop=4.0)
    a = sol(t)
    b = reshape(a, (agents_amt, 4))
    c,d,v,w = b[:,1], b[:,2], b[:,3], b[:,4]
    quiver(c,d,
        quiver=(v,w),
        xlims=(-10,10),
        ylims=(-10,10),
        arrow=arrow(5.0)
    )
end
"""
