using DifferentialEquations
using LinearAlgebra

#using InteractiveUtils

using Random
Random.seed!(2)

include("swarm_topology")
include("agent_dynamic")

Iₙ = Matrix{Number}(I, agents_amt, agents_amt)

R = zeros(agents_amt, agents_amt)
R[1,1] = 1


# τ(t) = floor(10*t) / 10
delay = 0.0
ω = 1000000000000000
# τ(u, t, p) = t - delay
τ(u,p,t) = floor(ω * t) / ω

function swarm_delay_in_control(du, u, h, p, t)
    du[:] = kron(Iₙ, Q) * u +
        kron(L, c*H_obs) * h(p, τ(u, p, t)) +
        kron(L, c*H_agent) * h(p, τ(u, p, t)) +
        kron(R, W) * u
        #(norm(kron(R, W) * u) <= 1 ? kron(R, W) * u : normalize(kron(R, W) * u))
end

function swarm_switch_system(du, u, h, p, t)
    agents_amt, G1 = six_agents_complex()
    agents_amt, G2 = six_agents_with_leader()
    L1 = laplacian_matrix(G1, dir=:in)
    L2 = laplacian_matrix(G2, dir=:in)
    c1 = compute_c(G1)
    c2 = compute_c(G2)
    α(t) = 10
    L, c = mod(floor(α(t) * t), 2) == 0 ? (L1, c1) : (L2, c2)
    du[:] = kron(Iₙ, Q) * u +
        kron(L, c*H_obs) * h(p, τ(u, p, t)) +
        kron(L, c*H_agent) * h(p, τ(u, p, t)) +
        kron(R, W) * u
end

"""
Задаем начальные условия
"""

h(p,t) = randn(2 * agent_order * agents_amt)
"""h(p,t) = Array{Float64}(
    transpose([1 2 0 0 5 6 0 0 9 10 0 0 13 14 0 0 17 18 0 0 21 22 0 0]))
"""
u0 = h(nothing, 0)

"""
Время интегрирования
"""
tstart = 0.
tend = 20.
tspan = (tstart, tend)

prob = DDEProblem(
    swarm_delay_in_control,
    u0,
    h,
    tspan;
    constant_lags = delay==0 ? [] : [delay],
    dependent_lags = [τ]
    )
alg = MethodOfSteps(Tsit5())
sol = solve(prob)

ENV["GKS_ENCODING"]="utf-8"
using Plots
plotly()
 # You can optionally choose a plotting backend()
function draw_sol(sol)
    num_ox = method*agent_order
    for i in 1:div(num_ox, 2)
        ind = num_ox*Array(0:agents_amt-1) .+ i
        agent_ind = Array(transpose(1:6))
        lbl = Array(map(p->"агент $p", agent_ind))
        plt = plot(
            sol,
            vars=ind,
            dpi=1000,
            title="Переменная $i",
            label = lbl
        )
        display(plt)
    end
    x_mid = map(
        p -> kron(ones(agents_amt, 1),
            sum(reshape(p, (num_ox, agents_amt)), dims = 2) / agents_amt
            ),
        sol.u)
    x_error = sol.u - x_mid
    x_sq_error = map(norm, x_error)
    plt = plot(
        sol.t,
        x_sq_error,
        title="ошибка консенсуса"
    )
    display(plt)
    x1 = map(p-> reshape(p, (num_ox, agents_amt)),sol.u)
    x2 = map(p -> p[1:div(num_ox, 2), :] - p[div(num_ox, 2)+1:num_ox, :], x1)
    x_sq_obs_error = map(
        p -> norm(p[1,:]) + norm(p[2, :]),
        x2
    )
    pltq = plot(
        sol.t,
        vcat(x_sq_obs_error...),
        title="ошибка наблюдения"
    )
    display(pltq)

end
draw_sol(sol)
