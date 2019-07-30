using LinearAlgebra

A = [1 0 1; 0 1 1; 1 0 1]

b = [1; 2; 2]

x = qr(A, Val(true)) \ b

if ~isapprox(norm(A * x - b), 0; atol=100*eps(Float64))
    error("Ошибка слишком велика, возможно, система не разрешима")
end
