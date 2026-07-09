# ♾️ infinite-order-lde-solver 📐

**A Python library for finding polynomial solutions to systems of infinite-order linear differential equations.**

This repository contains a symbolic and numerical solver built with Python, `NumPy` and `SymPy`. It is designed to find exact polynomial solutions $u(x)$ for systems of infinite-order linear differential equations of the form:

$$\sum_{j=0}^{\infty}A_{j}u^{(j)}(x)=P(x)$$

where $P(x)$ is a given polynomial vector, and $A_j$ are square matrix coefficients. 

---

## 🔬 Academic Background

This software was developed as the practical component of my **Bachelor's Qualification Thesis** in Applied Mathematics. 
The repository includes the full theoretical justification, mathematical proofs, and presentations (currently available in Ukrainian).
* 📄 `Кваліфікаційна робота.pdf` — The full Bachelor's thesis detailing the theory of polynomial solutions for infinite-order LDE systems.
* 📊 `Презентація Кваліфікаційної роботи.pdf` — Defense presentation slides.

---

## 🛠 Core API & Features

The solver is built around the main `InfiniteOrderDiffSystem` class, which defines the system through its matrix coefficients $A_j$ and the right-hand side vector $P(x)$.

### Solution Methods
The class provides three distinct mathematical approaches to find the polynomial solution:

* `undetermined_coefficients_method(t)`
    A universal algebraic approach that searches for a polynomial solution of a degree not exceeding $t$ by constructing and solving a system of linear algebraic equations.
    
* `Broggi_method()`
    An analytical method that finds the unique polynomial solution for the non-degenerate case, where the leading matrix is invertible ($\det A_0 \neq 0$).

* `degenerate_Broggi_method()`
    A generalized method that finds the general solution for the highly degenerate case, specifically when the initial matrices are zero ($A_0 = A_1 = \dots = A_{k-1} = 0$) but a subsequent matrix is non-degenerate ($\det A_k \neq 0$).

### Utilities
* `stringify_poly(coef_Matrix, letter)`
    A formatting utility that converts raw SymPy/NumPy coefficient matrices into human-readable, beautifully formatted polynomial equations using any specified variable letter.

---

## 🚀 Quick Start Example

Here is a basic example of defining a system and finding its solution using the undetermined coefficients method (see `example.py` for more details):

```python
from infinite_diff_system import InfiniteOrderDiffSystem
from utils import stringify_poly

# Define matrix A_j (as a string evaluated by SymPy) and polynomial vector P
input_user_A = "[[1, j*(j-1)], [j*(j-1), 0]]"
input_user_P = [[0, 1], [1, 0]]

# Initialize the system
A = InfiniteOrderDiffSystem(input_user_A, input_user_P)

# Find solution of degree <= 5
U = A.undetermined_coefficients_method(5)

# Print human-readable polynomial
print(stringify_poly(U, 'u'))
```

---

## 💻 Tech Stack
* Python 3

* SymPy (for symbolic mathematics, symbolic matrix operations, and solving algebraic systems)

* NumPy (for numerical arrays and matrix manipulations)

---

## 📜 License
Distributed under the MIT License. See LICENSE for more information.