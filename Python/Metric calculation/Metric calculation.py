import sympy as sp
from einsteinpy.symbolic import (
    MetricTensor,
    ChristoffelSymbols,
    RiemannCurvatureTensor,
    RicciTensor,
    RicciScalar,
    EinsteinTensor
)

# ---------------------------------------------------------
# General geometry routine
# Given a metric, this function computes:
# - Christoffel symbols
# - Riemann curvature tensor
# - Ricci tensor
# - Ricci scalar
# - Einstein tensor
# ---------------------------------------------------------
def run_geometry(metric, coords, name=""):
    Gamma = ChristoffelSymbols.from_metric(metric)
    Riemann = RiemannCurvatureTensor.from_christoffels(Gamma)
    Ricci = RicciTensor.from_riemann(Riemann)
    R_scalar = RicciScalar.from_riccitensor(Ricci, metric)
    Einstein = EinsteinTensor.from_metric(metric)
    return Gamma, Riemann, Ricci, R_scalar, Einstein


# ---------------------------------------------------------
# 1. Coordinate system and Schwarzschild metric definition
# ---------------------------------------------------------
# Coordinates: (t, r, theta, phi)
t, r, theta, phi = sp.symbols('t r theta phi')

# Mass parameter (geometrized units: G = c = 1)
M = sp.Symbol('M')
rs = 2 * M  # Schwarzschild radius

# Schwarzschild metric tensor in standard coordinates
g = [
    [-(1 - rs/r),          0,              0,                      0],
    [0,                    1/(1 - rs/r),   0,                      0],
    [0,                    0,              r**2,                   0],
    [0,                    0,              0,        r**2 * sp.sin(theta)**2]
]

# Create the metric tensor object
metric = MetricTensor(g, (t, r, theta, phi))


# ---------------------------------------------------------
# 2. Compute geometric quantities from the metric
# ---------------------------------------------------------
Gamma, Riemann, Ricci, R_scalar, Einstein = run_geometry(
    metric,
    (t, r, theta, phi)
)


# ---------------------------------------------------------
# 3. Helper function for LaTeX output of tensors
# ---------------------------------------------------------
def write_tensor_math(tensor, name, symmetry="symmetric"):
    """
    Generates LaTeX expressions for tensor components.

    - If all components vanish, prints: name_{μν} = 0
    - Otherwise, prints only independent nonzero components
      (taking tensor symmetry into account)
    """
    dim = 4
    has_nonzero = False

    # Check whether the tensor has any nonzero component
    for i in range(dim):
        for j in range(dim):
            if symmetry == "symmetric" and j < i:
                continue
            if sp.simplify(tensor[i, j]) != 0:
                has_nonzero = True
                break
        if has_nonzero:
            break

    if not has_nonzero:
        return [f"${name}_{{\\mu\\nu}} = 0$\\\\\n"]

    # Collect and print only nonzero independent components
    lines = []
    for i in range(dim):
        for j in range(dim):
            if symmetry == "symmetric" and j < i:
                continue
            val = sp.simplify(tensor[i, j])
            if val != 0:
                lines.append(
                    f"${name}_{{{i}{j}}} = {sp.latex(val)}$\\\\\n"
                )
    return lines


# ---------------------------------------------------------
# 4. Write results to a LaTeX file
# ---------------------------------------------------------
output_file = "schwarzschild_tensors.tex"

with open(output_file, "w", encoding="utf-8") as f:
    # LaTeX document header
    f.write(r"""
\documentclass{article}
\usepackage{amsmath}
\usepackage{setspace}
\usepackage{geometry}
\geometry{a4paper, margin=1in}
\setstretch{2}             
\begin{document}
""")

    f.write("\\section*{Schwarzschild Tensor Calculations}\n")

    # Metric tensor
    f.write("\\subsection*{Metric Tensor}\n")
    f.write("$g_{\\mu\\nu} = " + sp.latex(metric.tensor()) + "$\n\n")

    # Christoffel symbols
    f.write("\\subsection*{Christoffel Symbols}\n")
    dim = 4
    for a in range(dim):
        for b in range(dim):
            for c in range(b, dim):  # Γ^a_{bc} = Γ^a_{cb}
                val = sp.simplify(Gamma[a, b, c])
                if val != 0:
                    f.write(
                        f"$\\Gamma^{a}_{{{b}{c}}} = {sp.latex(val)}$\\\\\n"
                    )

    # Riemann curvature tensor
    f.write("\\subsection*{Riemann Tensor}\n")
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(c + 1, dim):  # antisymmetry in last indices
                    val = sp.simplify(Riemann[a, b, c, d])
                    if val != 0:
                        f.write(
                            f"$R^{a}_{{{b}{c}{d}}} = {sp.latex(val)}$\\\\\n"
                        )

    # Ricci tensor
    f.write("\\subsection*{Ricci Tensor}\n")
    for line in write_tensor_math(Ricci, "R", symmetry="symmetric"):
        f.write(line)

    # Ricci scalar
    f.write("\\subsection*{Ricci Scalar}\n")
    f.write(f"$R = {sp.latex(R_scalar.expr)}$\\\\\n\n")

    # Einstein tensor
    f.write("\\subsection*{Einstein Tensor}\n")
    for line in write_tensor_math(Einstein, "G", symmetry="symmetric"):
        f.write(line)

    # End of document
    f.write(r"\end{document}")

print("✓ DONE! File written to schwarzschild_tensors.tex")
