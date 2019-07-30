"""
Динамика одного робота
"""

include("C:\\Users\\liotb\\JuliaProjects\\LMI\\lmi")


" Начальные условия "
A = [0.0  1.0;
     0.0 0.0]

B = [   0.;
        1.]

C = [1. 0.]

function method_2()
    agent_order = 2
    method = 2
    """
    Вычисляем P > 0:
    transpose(A) * P + P * A - 2 * transpose(C)*C < 0
    """
    P = [1. -0.5; -0.5 1.]
    " Настраиваемые параметры системы "
    F = -P^(-1)*transpose(C)
    K = [-1. -2.]
    Q = [A B*K; zeros(agent_order, agent_order) (A + B*K)]

    H_obs = [zeros(agent_order, 2 * agent_order);
        -F*C zeros(agent_order, agent_order)]
    H_agent = [zeros(agent_order, 2 * agent_order);
        zeros(agent_order, agent_order) F*C]
    U = [0 0]
    W = [B*U zeros(agent_order, agent_order);
        zeros(agent_order, agent_order) B*U]
    return agent_order, method, Q, H_obs, H_agent, W, 1
end

function method_1()
    agent_order = 2
    method = 2
    K = compute_K_method_1(A, B)
    println("K = ", K)
    L = compute_L_observer(A, C)
    println("L = ", L)
    Q = [A zeros(agent_order, agent_order); -L*C (A + L*C)]
    H_agent = [zeros(agent_order, agent_order) B*K;
            zeros(agent_order, agent_order) B*K]
    H_obs = zeros(2*agent_order, 2*agent_order)
    U = [-1 0]
    W = [zeros(agent_order, agent_order) B*U;
         zeros(agent_order, agent_order) B*U]
    return agent_order, method, Q, H_obs, H_agent, W, 2
end

function compute_K_method_1(A, B)
    "AP + PAᵀ - 2BBᵀ < 0; P > 0"
    "K = -BᵀP⁻¹"
    Z = lmi_normalize(-transpose(A))
    Q = 2*B*transpose(B)
    R = remove_constant_matrix(Z, Q)
    P = lmi_normalize(Diagonal(ones(Number, 2)))
    V = remove_constant_matrix(P, zeros(Number,size(A)...))
    T = make_standart_lmi_from_system_of_lmi(R, V)
    y = solve_standart_lmi(T, 100)
    x = get_real_x(y)
    X = get_matrix_X(x, 2)
    return -transpose(B) / X
end

function compute_L_observer(A, C)
    "AᵀP + PA - CᵀYᵀ - YC < 0; P > 0"
    "L = P⁻¹Y"

    Z = lmi_normalize(-A)
    x_number = length(Z)
    U = lmi_vector(C)
    append!(Z, U)

    P = lmi_normalize(Matrix(I, 2, 2))
    N = lmi_vector(zeros(size(C)...))
    append!(P, N)

    R = make_standart_lmi_from_system_of_lmi(Z, P)
    T = remove_constant_matrix(R, zeros(2 .* size(A)...))
    y = solve_standart_lmi(T, 100)

    x = get_real_x(cat(y[1:x_number], y[length(T)];dims=1))
    X = get_matrix_X(x, size(A)[1])
    Y = y[x_number+1:length(T)-1]
    L = - X \ Y
end
agent_order, method, Q, H_obs, H_agent, W, obs_type = method_1()
