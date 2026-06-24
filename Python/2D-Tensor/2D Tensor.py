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
# 1. Coordinate system and metric definition
# ---------------------------------------------------------
x, y = sp.symbols('x y')
alpha = sp.Symbol('alpha', real=True, positive=True)

# 2D hyperbolic metric: ds^2 = α^2 / y^2 (dx^2 + dy^2)
g = [
    [alpha**2 / y**2, 0],
    [0, alpha**2 / y**2]
]

metric = MetricTensor(g, (x, y))

# ---------------------------------------------------------
# 2. General geometry routine
# Computes Christoffel symbols, Riemann, Ricci, Ricci scalar, Einstein tensor
# ---------------------------------------------------------
Gamma = ChristoffelSymbols.from_metric(metric)
Riemann = RiemannCurvatureTensor.from_christoffels(Gamma)
Ricci = RicciTensor.from_riemann(Riemann)
R_scalar = RicciScalar.from_riccitensor(Ricci, metric)
Einstein = EinsteinTensor.from_metric(metric)

# ---------------------------------------------------------
# 3. Helper function for LaTeX output of tensors
# ---------------------------------------------------------
def write_tensor_math(tensor, name, symmetry="symmetric"):
    T = sp.Matrix(tensor.tensor())  # تبدیل به Sympy Matrix
    dim = T.shape[0]
    lines = []
    for i in range(dim):
        for j in range(dim):
            if symmetry == "symmetric" and j < i:
                continue
            val = sp.simplify(T[i, j])
            if val != 0:
                lines.append(f"${name}_{{{i}{j}}} = {sp.latex(val)}$\\\\\n")
    if not lines:
        lines.append(f"${name}_{{\\mu\\nu}} = 0$\\\\\n")
    return lines

# ---------------------------------------------------------
# 4. Write results to a LaTeX file
# ---------------------------------------------------------
output_file = "hyperbolic2d_tensors.tex"

with open(output_file, "w", encoding="utf-8") as f:
    f.write(r"""
\documentclass{article}
\usepackage{amsmath}
\usepackage{setspace}
\usepackage{geometry}
\geometry{a4paper, margin=1in}
\setstretch{2}
\begin{document}
""")

    f.write("\\section*{2D Hyperbolic Tensor Calculations}\n")

    # Metric tensor
    f.write("\\subsection*{Metric Tensor}\n")
    f.write("$g_{\\mu\\nu} = " + sp.latex(metric.tensor()) + "$\n\n")

    # Christoffel symbols
    f.write("\\subsection*{Christoffel Symbols}\n")
    dim = 2
    Gamma_tensor = Gamma.tensor()
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                val = sp.simplify(Gamma_tensor[a, b, c])
                if val != 0:
                    f.write(f"$\\Gamma^{a}_{{{b}{c}}} = {sp.latex(val)}$\\\\\n")

    # Riemann curvature tensor
    f.write("\\subsection*{Riemann Tensor}\n")
    Riemann_tensor = Riemann.tensor()
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    val = sp.simplify(Riemann_tensor[a, b, c, d])
                    if val != 0:
                        f.write(f"$R^{a}_{{{b}{c}{d}}} = {sp.latex(val)}$\\\\\n")

    # Ricci tensor
    f.write("\\subsection*{Ricci Tensor}\n")
    for line in write_tensor_math(Ricci, "R", symmetry="symmetric"):
        f.write(line)

    # Ricci scalar
    f.write("\\subsection*{Ricci Scalar}\n")
    f.write(f"$R = {sp.latex(sp.simplify(R_scalar.expr))}$\\\\\n\n")

    # Einstein tensor
    f.write("\\subsection*{Einstein Tensor}\n")
    for line in write_tensor_math(Einstein, "G", symmetry="symmetric"):
        f.write(line)

    # End document
    f.write(r"\end{document}")

print("✓ DONE! File written to hyperbolic2d_tensors.tex")
