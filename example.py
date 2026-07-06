from infinite_diff_system import InfiniteOrderDiffSystem
from utils import stringify_poly


input_user_A = "[[1,j*(j-1)],[j*(j-1),0]]"
input_user_P = [[0, 1], [1, 0]]

A = InfiniteOrderDiffSystem(input_user_A, input_user_P)

print("System representation:")
A.print_system()

print("Index of the first non-degenerate matrix: ", A.first_non_degenerate_index())

U = A.undetermined_coefficients_method(5)

print(f"Form of the solution:\n{stringify_poly(U,'u')}")
print(f"Result of substituting the solution into the original system:\n{stringify_poly(A.apply_operator(U),'p')}")