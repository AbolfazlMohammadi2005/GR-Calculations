import sympy as sp
from einsteinpy.symbolic import (
    MetricTensor,
    ChristoffelSymbols,
    RiemannCurvatureTensor,
    RicciTensor,
    RicciScalar,
    EinsteinTensor
)

# =========================================================
# Static Spherically Symmetric Space-Time Geometry
#
# ds² = -e^(2Φ(r))dt² + e^(2Λ(r))dr²
#       + r²dθ² + r²sin²θ dφ²
#
# Computes:
#   - Christoffel Symbols
#   - Riemann Tensor
#   - Ricci Tensor
#   - Ricci Scalar
#   - Einstein Tensor
# =========================================================


# ---------------------------------------------------------
# Geometry Routine
# ---------------------------------------------------------
def run_geometry(metric):
    Gamma = ChristoffelSymbols.from_metric(metric)
    Riemann = RiemannCurvatureTensor.from_christoffels(Gamma)
    Ricci = RicciTensor.from_riemann(Riemann)
    R_scalar = RicciScalar.from_riccitensor(Ricci, metric)
    Einstein = EinsteinTensor.from_metric(metric)

    return Gamma, Riemann, Ricci, R_scalar, Einstein


# ---------------------------------------------------------
# Coordinates
# ---------------------------------------------------------
t, r, theta, phi = sp.symbols(
    "t r theta phi",
    real=True
)


# ---------------------------------------------------------
# Metric Functions
# ---------------------------------------------------------
Phi = sp.Function("Phi")(r)
Lambda = sp.Function("Lambda")(r)


# ---------------------------------------------------------
# Static Spherically Symmetric Metric
# ---------------------------------------------------------
g = [
    [-sp.exp(2 * Phi), 0, 0, 0],

    [0, sp.exp(2 * Lambda), 0, 0],

    [0, 0, r**2, 0],

    [0, 0, 0, r**2 * sp.sin(theta)**2]
]

metric = MetricTensor(
    g,
    (t, r, theta, phi)
)


# ---------------------------------------------------------
# Compute Geometry
# ---------------------------------------------------------
Gamma, Riemann, Ricci, R_scalar, Einstein = run_geometry(
    metric
)


# ---------------------------------------------------------
# Tensor LaTeX Writer
# ---------------------------------------------------------
def write_tensor_math(
    tensor,
    name,
    symmetry="symmetric"
):
    dim = 4
    lines = []

    has_nonzero = False

    for i in range(dim):
        for j in range(dim):

            if symmetry == "symmetric" and j < i:
                continue

            val = sp.simplify(
                sp.nsimplify(tensor[i, j])
            )

            if val != 0:
                has_nonzero = True
                break

        if has_nonzero:
            break

    if not has_nonzero:
        return [f"${name}_{{\\mu\\nu}} = 0$\\\\\n"]

    for i in range(dim):
        for j in range(dim):

            if symmetry == "symmetric" and j < i:
                continue

            val = sp.simplify(
                sp.nsimplify(tensor[i, j])
            )

            if val != 0:
                lines.append(
                    f"${name}_{{{i}{j}}}"
                    f" = {sp.latex(val)}$\\\\\n"
                )

    return lines


# ---------------------------------------------------------
# Output File
# ---------------------------------------------------------
output_file = "static_spherical_geometry.tex"


with open(output_file, "w", encoding="utf-8") as f:

    f.write(r"""
\documentclass{article}
\usepackage{amsmath}
\usepackage{geometry}
\usepackage{setspace}

\geometry{a4paper, margin=1in}
\setstretch{1.5}

\begin{document}
""")

    # -----------------------------------------------------
    # Title
    # -----------------------------------------------------
    f.write(
        "\\section*"
        "{Static Spherically Symmetric Space-Time}\n"
    )

    f.write(
        "Metric:\n"
        "\\[\n"
        "ds^2="
        "-e^{2\\Phi(r)}dt^2"
        "+e^{2\\Lambda(r)}dr^2"
        "+r^2d\\theta^2"
        "+r^2\\sin^2\\theta d\\phi^2\n"
        "\\]\n"
    )

    # -----------------------------------------------------
    # Metric
    # -----------------------------------------------------
    f.write("\\subsection*{Metric Tensor}\n")

    f.write(
        "$g_{\\mu\\nu}="
        + sp.latex(metric.tensor())
        + "$\n\n"
    )

    # -----------------------------------------------------
    # Christoffel Symbols
    # -----------------------------------------------------
    f.write(
        "\\subsection*{Christoffel Symbols}\n"
    )

    dim = 4

    for a in range(dim):
        for b in range(dim):
            for c in range(b, dim):

                val = sp.simplify(
                    sp.nsimplify(
                        Gamma[a, b, c]
                    )
                )

                if val != 0:
                    f.write(
                        f"$\\Gamma^{a}_{{{b}{c}}}"
                        f"={sp.latex(val)}$\\\\\n"
                    )

    # -----------------------------------------------------
    # Riemann Tensor
    # -----------------------------------------------------
    f.write(
        "\\subsection*{Riemann Curvature Tensor}\n"
    )

    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(c + 1, dim):

                    val = sp.simplify(
                        sp.nsimplify(
                            Riemann[a, b, c, d]
                        )
                    )

                    if val != 0:
                        f.write(
                            f"$R^{a}_{{{b}{c}{d}}}"
                            f"={sp.latex(val)}$\\\\\n"
                        )

    # -----------------------------------------------------
    # Ricci Tensor
    # -----------------------------------------------------
    f.write("\\subsection*{Ricci Tensor}\n")

    for line in write_tensor_math(
        Ricci,
        "R",
        symmetry="symmetric"
    ):
        f.write(line)

    # -----------------------------------------------------
    # Ricci Scalar
    # -----------------------------------------------------
    f.write("\\subsection*{Ricci Scalar}\n")

    f.write(
        "$R="
        + sp.latex(
            sp.simplify(
                sp.nsimplify(
                    R_scalar.expr
                )
            )
        )
        + "$\\\\\n\n"
    )

    # -----------------------------------------------------
    # Einstein Tensor
    # -----------------------------------------------------
    f.write("\\subsection*{Einstein Tensor}\n")

    for line in write_tensor_math(
        Einstein,
        "G",
        symmetry="symmetric"
    ):
        f.write(line)

    f.write("\n\\end{document}")


print(
    f"✓ DONE! File written to {output_file}"
)