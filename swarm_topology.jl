using LightGraphs

function six_agents_complex()
    amt = 6
    G = SimpleDiGraph(amt)
    add_edge!(G, 4, 1)
    add_edge!(G, 5, 1)
    add_edge!(G, 6, 1)
    add_edge!(G, 1, 2)
    add_edge!(G, 1, 3)
    add_edge!(G, 1, 4)
    add_edge!(G, 4, 5)
    add_edge!(G, 5, 6)
    return amt, G
end

function two_agents_digraph()
    amt = 2
    G = SimpleDiGraph(amt)
    add_edge!(G, 1, 2)
    return amt, G
end

function six_agents_with_leader()
    amt = 6
    G = SimpleDiGraph(amt)
    add_edge!(G, 1, 2)
    add_edge!(G, 1, 3)
    add_edge!(G, 1, 4)
    add_edge!(G, 2, 3)
    add_edge!(G, 4, 5)
    add_edge!(G, 5, 6)
    add_edge!(G, 3, 6)
    return amt, G
end

# agents_amt, G = six_agents_complex()
agents_amt, G = six_agents_with_leader()


L = laplacian_matrix(G, dir=:in)
function compute_c(G)
    λ = real(laplacian_spectrum(G, dir=:in))
    sort!(λ)
    μ = λ[2:agents_amt]
    c = length(μ) > 0 ? max(1. / (min(μ...)), 1) : 1
end

c = compute_c(G)
c = 10c
