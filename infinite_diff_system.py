import numpy as np
import sympy as sp
import math
from utils import stringify_poly, build_diff_matrix, find_inverse_coefficients
from sympy.solvers.solveset import NonlinearError

j = sp.symbols("j")


class InfiniteOrderDiffSystem:
    def __init__(self, A_array, P_vector):

        matrix_A_data = sp.sympify(A_array)
        A_j = sp.Matrix(matrix_A_data)
        self.A_j = A_j
        self.dimension = A_j.shape[0]
        n_vector = len(P_vector)

        if self.dimension != n_vector:
            raise ValueError(
                f"Невідповідність розмірностей (n): "
                f"Матриця Aj має розмір {self.dimension}x{self.dimension}, "
                f"а вектор P має {n_vector} компонент(и). "
            )

        A_func = sp.lambdify(j, A_j, "sympy")

        self.A_func = A_func
        self.P = np.array(P_vector)
        self.degree = len(P_vector[0])
        pass

    def print_system(self):
        sp.pprint(f"A_j = {self.A_j}\n")
        print(stringify_poly(self.P, "p"))

    def is_bijective_operator(self):
        if sp.det(self.A_func(0)) != 0:
            return True

        return False

    def first_non_degenerate_index(self):

        det_expr = self.A_j.det()

        roots = sp.solve(det_expr, j)

        integer_roots = [r for r in roots if r.is_integer and r >= 0]

        d = 0
        while d in integer_roots:
            d += 1
        if (self.A_func(d)).det() == 0 and d == 0:
            return None

        return d

    def Broggi_method(self):

        if not self.is_bijective_operator():
            return None

        s = self.degree
        n = self.dimension
        A_j = self.A_func
        P = self.P

        A = np.zeros((s, n, n))
        B = np.zeros((s, n, n))
        res = np.zeros((n, s))

        D = build_diff_matrix(s)

        for i in range(s):
            A[i, :, :] = np.array(A_j(i)).astype(np.float64)
            B[i, :, :] = find_inverse_coefficients(A, B, i, n)

        res += B[0, :, :] @ P
        for i in range(1, s):
            res += B[i, :, :] @ (P @ np.linalg.matrix_power(D, i))

        return res

    def degenerate_Broggi_method(self):

        A_j = self.A_func
        s = self.degree
        n = self.dimension
        P = self.P
        k = self.first_non_degenerate_index()

        for j in range(k):
            if not np.array_equal(sp.matrix2numpy(A_j(j)), np.zeros_like(A_j(j))):
                return None

        A = np.zeros((s, n, n))
        B = np.zeros((s, n, n))
        sol_system = np.zeros((n, s))

        for i in range(s):
            A[i, :, :] = np.array(A_j(i + k)).astype(np.float64)
            B[i, :, :] = find_inverse_coefficients(A, B, i, n)

        for i in range(s):
            for j in range(s - i):
                sol_system[:, i] += (
                    B[j, :, :]
                    @ P[:, i + j]
                    * math.factorial(i + j)
                    / math.factorial(i + k)
                )

        res = np.zeros((n, s + k))
        res[:, k:] = sol_system

        c_vars = sp.Matrix(
            [
                [sp.Symbol(f"c{j}{i + 1}") / math.factorial(j) for i in range(n)]
                for j in range(k)
            ]
        )

        output = sp.Matrix(res)
        output[:, :k] = sp.Matrix(c_vars).T
        return sp.nsimplify(output)

    def undetermined_coefficients_method(self, m):
        s = self.degree - 1
        n = self.dimension
        A_func = self.A_func
        deg_max = m

        for t in range(s, deg_max + 1):
            P_current = self.P.T.tolist()

            try:
                u_vars = [
                    [sp.Symbol(f"c{j}{i + 1}") for i in range(n)] for j in range(t + 1)
                ]
                all_symbols = [sym for sublist in u_vars for sym in sublist]
                U = [sp.Matrix(u_vars[i]) for i in range(t + 1)]

                A = [A_func(j) for j in range(t + 1)]

                while len(P_current) <= t:
                    P_current.append([0] * n)

                all_equations = []

                for k in range(t + 1):
                    lhs_degree_k = -sp.Matrix(P_current[k])

                    for j in range(t - k + 1):
                        coeff = math.perm(k + j, j)
                        lhs_degree_k += coeff * A[j] * U[k + j]

                    all_equations.extend(lhs_degree_k[:])

                solution_set = sp.linsolve(all_equations, all_symbols)

                if solution_set and solution_set != sp.EmptySet:
                    solution = list(solution_set)[0]
                    return sp.Matrix(
                        [[solution[j * n + i] for j in range(t + 1)] for i in range(n)]
                    )
                else:
                    continue

            except NonlinearError as e:
                continue
            except Exception as e:
                continue

        return None

    def apply_operator(self, U):
        n, t = U.shape
        A_j_func = self.A_func

        res = sp.zeros(n, t)

        D = sp.Matrix(build_diff_matrix(t))

        current_u_diff = sp.Matrix(U)

        for j in range(t):
            Aj = sp.Matrix(A_j_func(j))

            res += Aj * current_u_diff

            current_u_diff = current_u_diff * D

        return res

    def is_solution(self, U):

        ##TODO: Need end this method. 

        return 

