import numpy as np
import sympy as sp


def build_diff_matrix(n):
    return np.diag(np.arange(1, n), k=-1)


def find_inverse_coefficients(A, B, j, n):
    res = np.zeros((n, n))
    if j == 0:
        return np.linalg.inv(A[0, :, :])
    i = 1
    while i <= j:
        res += A[i, :, :] @ B[j - i, :, :]
        i += 1
    return -B[0, :, :] @ res


def stringify_poly(coef_Matrix, letter):
    col, row = coef_Matrix.shape
    type_input = type(coef_Matrix)
    x = sp.symbols("x")
    results = []

    if type_input == np.ndarray:
        mask = np.abs(coef_Matrix) < 1e-10
        zeros = np.zeros_like(coef_Matrix)
        coef_Matrix[mask] = zeros[mask]

        for i in range(col):
            poly_expr = sum(
                [sp.nsimplify(coef_Matrix[i, j]) * x**j for j in range(row)]
            )

            poly_expr = poly_expr
            results.append(f"{letter}_{i + 1} = {poly_expr}")

        return "\n".join(results)

    else:
        clean_matrix = coef_Matrix.applyfunc(
            lambda x: x if (sp.sympify(x).free_symbols or abs(x) > 10e-10) else 0
        )

        for i in range(col):
            poly_expr = 0
            for j in range(row):
                coef = clean_matrix[i, j]
                if coef == 0:
                    continue

                # if coef.is_Number:
                #    coef = round(float(coef), 6)
                #    if coef.is_integer():
                #        coef = int(coef)

                poly_expr += coef * x**j

            poly_expr = sp.simplify(poly_expr)
            results.append(f"{letter}_{i + 1} = {poly_expr}")

        return "\n".join(results)
