"""
Calculate and plot q-shell averaged intensity from one merged HKL file.


Usage: 

cd project/path/phytochrome2026/dataPhyto_reduction_light/results_metric

python /path/to/shell_intensity.py \
    ../outputs/upto1ps_test01/dataPhyto_upto1ps_C12_LPF_DRL_SCL_AVG.hkl

* Please change the DEFAULT_INPUT_FILE/path parameters 
* and Primitive orthorhombic unit cell parameters below if needed
* Note that: All relative paths are interpreted 
* relative to the directory containing shell_intensity.py.




Directory convention
--------------------
The script is expected to live in:

    dataPhyto_reduction_light/results_metric/shell_intensity.py

A test-run folder may be passed relative to the script location, for example:

    python shell_intensity.py ../outputs/upto1ps_test01

The script then reads:

    ../outputs/upto1ps_test01/
        dataPhyto_upto1ps_C12_LPF_DRL_SCL_AVG.hkl

and writes, by default:

    results_metric/upto1ps_test01/
        dataPhyto_upto1ps_C12_LPF_DRL_SCL_AVG_q_shell.png
        dataPhyto_upto1ps_C12_LPF_DRL_SCL_AVG_q_shell.txt

You may also pass the .hkl file itself:

    python shell_intensity.py \
        ../outputs/upto1ps_test01/dataPhyto_upto1ps_C12_LPF_DRL_SCL_AVG.hkl

Important: relative input and output paths are resolved relative to the folder
containing this script, not relative to the terminal's current directory.

Expected input columns
----------------------
Only the first four whitespace-separated columns are used:

    h  k  l  intensity  [other columns ...]

Any fifth or later column is ignored.

Crystal convention
------------------
Primitive orthorhombic unit cell:

    a = 54.22 A
    b = 115.78 A
    c = 117.08 A
    alpha = beta = gamma = 90 degrees

The script follows q = 1/d, not q = 2*pi/d:

    q = sqrt((h/a)^2 + (k/b)^2 + (l/c)^2)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

# Use a non-interactive backend so the script works on HPC compute nodes.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


###############################################################################
# Paths and default parameters
###############################################################################

# Absolute path to the folder containing this Python script.
# All relative input/output paths below are anchored here.
SCRIPT_DIR = Path(__file__).resolve().parent

DEFAULT_TEST_RUN = Path("../outputs/upto1ps_test01")
DEFAULT_HKL_NAME = "dataPhyto_upto1ps_C12_LPF_DRL_SCL_AVG.hkl"

# Output filenames are derived automatically from the actual input file name.

# Primitive orthorhombic unit cell, in Angstrom.
DEFAULT_A = 54.22
DEFAULT_B = 115.78
DEFAULT_C = 117.08
DEFAULT_N_SHELLS = 50


###############################################################################
# Path helpers
###############################################################################


def resolve_from_script(path: Path) -> Path:
    """Resolve a path relative to the folder containing this script."""
    path = path.expanduser()
    if not path.is_absolute():
        path = SCRIPT_DIR / path
    return path.resolve()



def resolve_input_file(input_path: Path, hkl_name: str) -> Path:
    """Resolve either a test-run directory or a direct HKL file path."""
    resolved = resolve_from_script(input_path)

    if resolved.is_dir():
        resolved = resolved / hkl_name

    if not resolved.is_file():
        raise FileNotFoundError(
            "Input HKL file not found.\n"
            f"Resolved path: {resolved}\n"
            "Pass either a test-run folder or the full .hkl file path."
        )

    return resolved


###############################################################################
# Input and preprocessing functions
###############################################################################


def read_hkl_intensity(file_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read h, k, l, and intensity from the first four columns.

    Blank lines, comment lines, text headers, malformed lines, and columns after
    intensity are ignored. Miller indices must be integer-valued.
    """
    rows: list[tuple[int, int, int, float]] = []

    with file_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()

            if not stripped or stripped.startswith(("#", ";")):
                continue

            fields = stripped.split()
            if len(fields) < 4:
                continue

            try:
                h_float, k_float, l_float = map(float, fields[:3])
                intensity = float(fields[3])
            except ValueError:
                # Skip text headers, for example: h k l I sigma(I)
                continue

            h = int(round(h_float))
            k = int(round(k_float))
            l = int(round(l_float))

            if not np.allclose(
                [h_float, k_float, l_float],
                [h, k, l],
                rtol=0.0,
                atol=1.0e-8,
            ):
                raise ValueError(
                    f"Non-integer Miller index in {file_path} "
                    f"at line {line_number}."
                )

            rows.append((h, k, l, intensity))

    if not rows:
        raise ValueError(
            f"No valid h, k, l, intensity rows were found in {file_path}."
        )

    data = np.asarray(rows, dtype=float)
    hkl = data[:, :3].astype(np.int64)
    intensity = data[:, 3]
    return hkl, intensity


def average_duplicate_reflections(
    hkl: np.ndarray,
    intensity: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Average finite intensities for repeated (h, k, l) rows."""
    unique_hkl, inverse = np.unique(hkl, axis=0, return_inverse=True)

    finite = np.isfinite(intensity)
    sums = np.bincount(
        inverse[finite],
        weights=intensity[finite],
        minlength=unique_hkl.shape[0],
    )
    counts = np.bincount(
        inverse[finite],
        minlength=unique_hkl.shape[0],
    )

    averaged_intensity = np.full(unique_hkl.shape[0], np.nan, dtype=float)
    has_data = counts > 0
    averaged_intensity[has_data] = sums[has_data] / counts[has_data]

    return unique_hkl, averaged_intensity


###############################################################################
# Orthorhombic q and q-shell calculations
###############################################################################


def calculate_q_orthorhombic(
    hkl: np.ndarray,
    a: float,
    b: float,
    c: float,
) -> np.ndarray:
    """Calculate q = 1/d for an orthorhombic unit cell."""
    h = hkl[:, 0].astype(float)
    k = hkl[:, 1].astype(float)
    l = hkl[:, 2].astype(float)

    q_squared = (h / a) ** 2 + (k / b) ** 2 + (l / c) ** 2
    return np.sqrt(q_squared)


def q_shell_average(
    q: np.ndarray,
    intensity: np.ndarray,
    q_edges: np.ndarray,
    omit_zero: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Return mean intensity and reflection count in every q-shell."""
    n_shells = q_edges.size - 1

    # Internal edges follow [left, right); q_max is included in the final shell.
    shell_index = np.searchsorted(q_edges, q, side="right") - 1
    shell_index[q == q_edges[-1]] = n_shells - 1

    valid = (
        np.isfinite(q)
        & np.isfinite(intensity)
        & (shell_index >= 0)
        & (shell_index < n_shells)
    )
    if omit_zero:
        valid &= intensity != 0.0

    valid_index = shell_index[valid]
    valid_intensity = intensity[valid]

    shell_count = np.bincount(valid_index, minlength=n_shells).astype(np.int64)
    shell_sum = np.bincount(
        valid_index,
        weights=valid_intensity,
        minlength=n_shells,
    )

    shell_mean = np.full(n_shells, np.nan, dtype=float)
    nonempty = shell_count > 0
    shell_mean[nonempty] = shell_sum[nonempty] / shell_count[nonempty]

    return shell_mean, shell_count


def standardize(values: np.ndarray) -> np.ndarray:
    """Return z-scores using MATLAB-compatible sample standard deviation."""
    if values.size < 2:
        raise ValueError("At least two valid q-shells are required to standardize.")

    standard_deviation = np.std(values, ddof=1)
    if np.isclose(standard_deviation, 0.0):
        raise ValueError("Cannot standardize a constant q-shell intensity curve.")

    return (values - np.mean(values)) / standard_deviation


###############################################################################
# Command-line interface and main analysis
###############################################################################


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate and plot q-shell averaged intensity from one HKL file "
            "using an orthorhombic unit cell. Relative paths are resolved "
            "relative to the folder containing this script."
        )
    )
    parser.add_argument(
        "input_path",
        nargs="?",
        type=Path,
        default=DEFAULT_TEST_RUN,
        help=(
            "Test-run folder or direct .hkl file path. Relative paths are "
            "resolved from the script folder. Default: "
            f"{DEFAULT_TEST_RUN}"
        ),
    )
    parser.add_argument(
        "--file-name",
        default=DEFAULT_HKL_NAME,
        help=(
            "HKL filename used when input_path is a directory; default: "
            f"{DEFAULT_HKL_NAME}"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output directory. Relative paths are resolved from the script "
            "folder. Default: <script_folder>/<test_run_name>"
        ),
    )
    parser.add_argument(
        "--output-figure",
        type=Path,
        default=None,
        help=(
            "Explicit PNG path. Relative paths are resolved from the script "
            "folder. Overrides --output-dir for the figure."
        ),
    )
    parser.add_argument(
        "--output-data",
        type=Path,
        default=None,
        help=(
            "Explicit text-output path. Relative paths are resolved from the "
            "script folder. Overrides --output-dir for the data."
        ),
    )
    parser.add_argument(
        "--n-shells",
        type=int,
        default=DEFAULT_N_SHELLS,
        help=f"Number of q-shells; default: {DEFAULT_N_SHELLS}",
    )
    parser.add_argument(
        "--a",
        type=float,
        default=DEFAULT_A,
        help=f"Orthorhombic lattice a in Angstrom; default: {DEFAULT_A}",
    )
    parser.add_argument(
        "--b",
        type=float,
        default=DEFAULT_B,
        help=f"Orthorhombic lattice b in Angstrom; default: {DEFAULT_B}",
    )
    parser.add_argument(
        "--c",
        type=float,
        default=DEFAULT_C,
        help=f"Orthorhombic lattice c in Angstrom; default: {DEFAULT_C}",
    )
    parser.add_argument(
        "--standardize",
        action="store_true",
        help="Plot z-scored shell means instead of raw shell means",
    )
    parser.add_argument(
        "--include-zero",
        action="store_true",
        help="Include exact zero intensities in shell averages",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.n_shells < 1:
        raise ValueError("--n-shells must be at least 1.")
    if args.a <= 0 or args.b <= 0 or args.c <= 0:
        raise ValueError("Lattice constants a, b, and c must be positive.")
    input_file = resolve_input_file(args.input_path, args.file_name)
    test_run_name = input_file.parent.name

    # Derive output names from the actual input filename by default.
    # Example:
    #   dataPhyto_upto1ps_C12_LPF_DRL_SCL_AVG.hkl
    # becomes:
    #   dataPhyto_upto1ps_C12_LPF_DRL_SCL_AVG_q_shell.png/.txt
    output_stem = f"{input_file.stem}_q_shell"

    # Default output structure:
    #   <script folder>/<test run name>/<input_stem>_q_shell.png
    #   <script folder>/<test run name>/<input_stem>_q_shell.txt
    if args.output_dir is None:
        output_dir = SCRIPT_DIR / test_run_name
    else:
        output_dir = resolve_from_script(args.output_dir)

    if args.output_figure is None:
        output_figure = output_dir / f"{output_stem}.png"
    else:
        output_figure = resolve_from_script(args.output_figure)

    if args.output_data is None:
        output_data = output_dir / f"{output_stem}.txt"
    else:
        output_data = resolve_from_script(args.output_data)

    # Read h, k, l, and intensity. A fifth or later column is ignored.
    hkl, intensity = read_hkl_intensity(input_file)
    n_rows = hkl.shape[0]

    # Combine repeated Miller indices before q-shell averaging.
    hkl, intensity = average_duplicate_reflections(hkl, intensity)
    n_unique = hkl.shape[0]

    # Orthorhombic q = 1/d; centering P requires no formula change.
    q = calculate_q_orthorhombic(hkl, args.a, args.b, args.c)

    finite_q = np.isfinite(q)
    q = q[finite_q]
    intensity = intensity[finite_q]

    if q.size == 0:
        raise ValueError("No finite q values were calculated.")
    if np.isclose(q.min(), q.max()):
        raise ValueError("All reflections have the same q value.")

    # N shells require N + 1 edges. Plot each shell at its midpoint.
    q_edges = np.linspace(q.min(), q.max(), args.n_shells + 1)
    q_centers = 0.5 * (q_edges[:-1] + q_edges[1:])

    shell_mean, shell_count = q_shell_average(
        q,
        intensity,
        q_edges,
        omit_zero=not args.include_zero,
    )

    # Keep only shells that contain at least one valid intensity.
    valid_shells = np.isfinite(shell_mean)
    q_centers = q_centers[valid_shells]
    shell_mean = shell_mean[valid_shells]
    shell_count = shell_count[valid_shells]

    if q_centers.size == 0:
        raise ValueError("No q-shell contains a valid intensity value.")

    output_figure.parent.mkdir(parents=True, exist_ok=True)
    output_data.parent.mkdir(parents=True, exist_ok=True)

    # Save unstandardized shell means for later analysis.
    output_table = np.column_stack((q_centers, shell_mean, shell_count))
    np.savetxt(
        output_data,
        output_table,
        fmt=["%.10e", "%.10e", "%d"],
        header="q_center_A^-1  mean_intensity  n_reflections",
    )

    if args.standardize:
        y_values = standardize(shell_mean)
        y_label = r"standardized $\langle I(q) \rangle$"
        title_suffix = "standardized shell mean"
    else:
        y_values = shell_mean
        y_label = r"$\langle I(q) \rangle$"
        title_suffix = "raw shell mean"

    # Draw and save one q-shell averaged intensity curve.
    figure, axis = plt.subplots(figsize=(8, 6))
    axis.plot(q_centers, y_values, ":o", linewidth=1, markersize=6)
    axis.set_xlabel(r"$q$ [$\mathrm{\AA}^{-1}$]", fontsize=14)
    axis.set_ylabel(y_label, fontsize=14)
    axis.set_title(f"{test_run_name}: {input_file.stem}\n({title_suffix})", fontsize=12)
    axis.grid(True, alpha=0.3)
    figure.tight_layout()
    figure.savefig(output_figure, dpi=300, bbox_inches="tight")
    plt.close(figure)

    print(f"Script directory: {SCRIPT_DIR}")
    print(f"Test run: {test_run_name}")
    print(f"Input file: {input_file}")
    print(f"Rows read: {n_rows}")
    print(f"Unique reflections: {n_unique}")
    print("Lattice type: primitive orthorhombic (P)")
    print(f"Unit cell: a={args.a} A, b={args.b} A, c={args.c} A")
    print("q convention: q = 1/d")
    print(f"Valid q-shells: {q_centers.size} / {args.n_shells}")
    print(f"Zero intensities omitted: {not args.include_zero}")
    print(f"Figure saved to: {output_figure}")
    print(f"Shell data saved to: {output_data}")


if __name__ == "__main__":
    main()
