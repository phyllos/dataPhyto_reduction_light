#!/usr/bin/env python3
"""
Compare two merged HKL datasets using their unique common reflections.

This is a Python rewrite of the main workflow in the original MATLAB files:

    [C, id1, id2] = intersect(HKL_1, HKL_2, 'rows');

The script:

1. Reads h, k, l, intensity from the first four columns of two HKL files.
2. Collapses duplicate (h, k, l) rows within each file by averaging their
   finite intensities.
3. Finds the sorted, unique intersection of the two HKL arrays.
4. Aligns the two intensity arrays using the returned intersection indices.
5. Calculates the correlation coefficient between aligned intensities.
6. Calculates q = 1/d for a primitive orthorhombic unit cell.
7. Calculates q-shell mean intensities for the common reflections.
8. Saves the common reflections, q-shell table, and comparison figure.

Relative input paths are resolved from the terminal's current working
folder. The default output folder is created beside this script.

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

Example
-------
Run from the project root:

    python results_metric/shell_intensity_compare.py \
        inputs/DATA_upto1ps/upto1ps_AVG_of_All.hkl \
        outputs/upto1ps_test01/dataPhyto_upto1ps_C12_LPF_DRL_SCL_AVG.hkl \
        --label1 original \
        --label2 python
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
# Default parameters
###############################################################################

SCRIPT_DIR = Path(__file__).resolve().parent

# Primitive orthorhombic unit cell, in Angstrom.
DEFAULT_A = 54.22
DEFAULT_B = 115.78
DEFAULT_C = 117.08
DEFAULT_N_SHELLS = 50


###############################################################################
# Input and preprocessing
###############################################################################


def resolve_input_file(path: Path) -> Path:
    """Resolve an input path relative to the current terminal directory."""
    resolved = path.expanduser().resolve()

    if not resolved.is_file():
        raise FileNotFoundError(
            "Input HKL file not found.\n"
            f"Resolved path: {resolved}"
        )

    return resolved


def resolve_output_dir(path: Path | None, default_name: str) -> Path:
    """Resolve an output folder.

    An explicitly supplied relative path is resolved from the current terminal
    directory. When omitted, output is written beside this script.
    """
    if path is None:
        return (SCRIPT_DIR / default_name).resolve()

    return path.expanduser().resolve()


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
                # Skip text headers such as: h k l I sigma(I)
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
    """Return unique HKLs and mean finite intensity for each reflection.

    ``np.unique(..., axis=0)`` sorts the unique HKL rows lexicographically.
    Repeated reflections are combined before the two datasets are intersected.
    """
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

    mean_intensity = np.full(unique_hkl.shape[0], np.nan, dtype=float)
    has_data = counts > 0
    mean_intensity[has_data] = sums[has_data] / counts[has_data]

    return unique_hkl, mean_intensity


###############################################################################
# MATLAB-style unique row intersection
###############################################################################


def hkl_as_structured(hkl: np.ndarray) -> np.ndarray:
    """Represent an (N, 3) integer HKL array as a one-dimensional record array."""
    contiguous = np.ascontiguousarray(hkl, dtype=np.int64)
    dtype = np.dtype(
        [
            ("h", np.int64),
            ("k", np.int64),
            ("l", np.int64),
        ]
    )
    return contiguous.view(dtype).reshape(-1)


def intersect_unique_hkl(
    hkl1: np.ndarray,
    hkl2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return sorted unique common HKLs and indices into both input arrays.

    This corresponds to MATLAB:

        [C, id1, id2] = intersect(HKL_1, HKL_2, 'rows');

    Both input arrays must already contain unique HKL rows.
    """
    structured1 = hkl_as_structured(hkl1)
    structured2 = hkl_as_structured(hkl2)

    _, id1, id2 = np.intersect1d(
        structured1,
        structured2,
        assume_unique=True,
        return_indices=True,
    )

    common_hkl = hkl1[id1]

    # Cross-check equivalent to MATLAB's isequal(...).
    if not np.array_equal(common_hkl, hkl2[id2]):
        raise RuntimeError("Internal HKL intersection alignment check failed.")

    return common_hkl, id1, id2


###############################################################################
# Correlation coefficient
###############################################################################


def corr_coeff(x0: np.ndarray, x1: np.ndarray) -> float:
    """Return the correlation coefficient between two aligned vectors.

    The formula follows corrCoeff.m. Non-finite pairs are removed together so
    the reflection alignment is preserved.
    """
    x0 = np.asarray(x0, dtype=float).reshape(-1, order="F")
    x1 = np.asarray(x1, dtype=float).reshape(-1, order="F")

    if x0.size != x1.size:
        raise ValueError("The two intensity arrays must have equal length.")

    valid = np.isfinite(x0) & np.isfinite(x1)
    x0 = x0[valid]
    x1 = x1[valid]

    if x0.size < 2:
        raise ValueError(
            "At least two finite common intensity pairs are required "
            "for correlation."
        )

    x0 = x0 - np.mean(x0)
    x1 = x1 - np.mean(x1)

    denominator = np.sqrt(np.dot(x0, x0) * np.dot(x1, x1))
    if np.isclose(denominator, 0.0):
        raise ValueError(
            "Correlation is undefined because at least one aligned "
            "intensity array is constant."
        )

    return float(np.dot(x0, x1) / denominator)


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
    """Return mean intensity and reflection count in every q-shell.

    Shells use [left, right), except that q_max is included in the final shell.
    Exact zeros are omitted by default to reproduce MATLAB ``nonzeros()``.
    """
    n_shells = q_edges.size - 1

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

    shell_count = np.bincount(
        valid_index,
        minlength=n_shells,
    ).astype(np.int64)

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
        raise ValueError("At least two valid q-shells are required.")

    standard_deviation = np.std(values, ddof=1)
    if np.isclose(standard_deviation, 0.0):
        raise ValueError("Cannot standardize a constant q-shell curve.")

    return (values - np.mean(values)) / standard_deviation


###############################################################################
# Command-line interface
###############################################################################


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Intersect two merged HKL datasets, align their intensities, "
            "calculate correlation, and compare q-shell mean intensities."
        )
    )

    parser.add_argument(
        "file1",
        type=Path,
        help="First HKL file; relative paths use the terminal working directory",
    )
    parser.add_argument(
        "file2",
        type=Path,
        help="Second HKL file; relative paths use the terminal working directory",
    )
    parser.add_argument(
        "--label1",
        default=None,
        help="Plot label for the first dataset; default: first filename stem",
    )
    parser.add_argument(
        "--label2",
        default=None,
        help="Plot label for the second dataset; default: second filename stem",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output folder. Relative paths use the terminal working directory. "
            "Default: <script_folder>/<file1_parent>_vs_<file2_parent>"
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
        "--plot-mode",
        choices=("standardized", "raw"),
        default="standardized",
        help=(
            "Plot standardized shell means, as in the MATLAB script, "
            "or raw shell means; default: standardized"
        ),
    )
    parser.add_argument(
        "--include-zero",
        action="store_true",
        help="Include exact zero intensities in q-shell averages",
    )

    return parser.parse_args()


###############################################################################
# Main analysis
###############################################################################


def main() -> None:
    args = parse_args()

    if args.n_shells < 1:
        raise ValueError("--n-shells must be at least 1.")
    if args.a <= 0 or args.b <= 0 or args.c <= 0:
        raise ValueError("Lattice constants a, b, and c must be positive.")

    file1 = resolve_input_file(args.file1)
    file2 = resolve_input_file(args.file2)

    label1 = args.label1 or file1.stem
    label2 = args.label2 or file2.stem

    default_output_name = f"{file1.parent.name}_vs_{file2.parent.name}"
    output_dir = resolve_output_dir(args.output_dir, default_output_name)
    output_dir.mkdir(parents=True, exist_ok=True)

    compare_stem = f"{file1.stem}_vs_{file2.stem}"

    common_output = output_dir / f"{compare_stem}_common_reflections.txt"
    shell_output = output_dir / f"{compare_stem}_q_shell.txt"
    figure_output = output_dir / f"{compare_stem}_q_shell.png"

    ###########################################################################
    # Read and make each source unique
    ###########################################################################

    hkl1_raw, intensity1_raw = read_hkl_intensity(file1)
    hkl2_raw, intensity2_raw = read_hkl_intensity(file2)

    hkl1, intensity1 = average_duplicate_reflections(
        hkl1_raw,
        intensity1_raw,
    )
    hkl2, intensity2 = average_duplicate_reflections(
        hkl2_raw,
        intensity2_raw,
    )

    ###########################################################################
    # MATLAB equivalent:
    #
    # [C, id1, id2] = intersect(HKL_1, HKL_2, 'rows');
    ###########################################################################

    common_hkl, id1, id2 = intersect_unique_hkl(hkl1, hkl2)

    if common_hkl.shape[0] == 0:
        raise ValueError("The two files have no common (h, k, l) reflections.")

    common_intensity1 = intensity1[id1]
    common_intensity2 = intensity2[id2]

    # Save the unique common reflections and aligned intensities.
    common_table = np.column_stack(
        (
            common_hkl,
            common_intensity1,
            common_intensity2,
        )
    )
    np.savetxt(
        common_output,
        common_table,
        fmt=["%d", "%d", "%d", "%.10e", "%.10e"],
        header=(
            f"h  k  l  intensity_{label1}  intensity_{label2}\n"
            "Unique common reflections, sorted lexicographically by h, k, l"
        ),
    )

    correlation = corr_coeff(
        common_intensity1,
        common_intensity2,
    )

    ###########################################################################
    # Calculate q-shell means using only the common reflections
    ###########################################################################

    q = calculate_q_orthorhombic(
        common_hkl,
        args.a,
        args.b,
        args.c,
    )

    if q.size == 0:
        raise ValueError("No q values were calculated.")
    if np.isclose(q.min(), q.max()):
        raise ValueError("All common reflections have the same q value.")

    q_edges = np.linspace(
        q.min(),
        q.max(),
        args.n_shells + 1,
    )
    q_centers = 0.5 * (
        q_edges[:-1]
        + q_edges[1:]
    )

    shell_mean1, shell_count1 = q_shell_average(
        q,
        common_intensity1,
        q_edges,
        omit_zero=not args.include_zero,
    )
    shell_mean2, shell_count2 = q_shell_average(
        q,
        common_intensity2,
        q_edges,
        omit_zero=not args.include_zero,
    )

    # Match the MATLAB logic: retain only shells valid in both datasets.
    valid_shells = (
        np.isfinite(shell_mean1)
        & np.isfinite(shell_mean2)
    )

    q_centers = q_centers[valid_shells]
    shell_mean1 = shell_mean1[valid_shells]
    shell_mean2 = shell_mean2[valid_shells]
    shell_count1 = shell_count1[valid_shells]
    shell_count2 = shell_count2[valid_shells]

    if q_centers.size == 0:
        raise ValueError(
            "No q-shell contains a valid intensity in both datasets."
        )

    # Save raw shell averages and counts, regardless of plot mode.
    shell_table = np.column_stack(
        (
            q_centers,
            shell_mean1,
            shell_count1,
            shell_mean2,
            shell_count2,
        )
    )
    np.savetxt(
        shell_output,
        shell_table,
        fmt=["%.10e", "%.10e", "%d", "%.10e", "%d"],
        header=(
            f"q_center_A^-1  mean_{label1}  count_{label1}  "
            f"mean_{label2}  count_{label2}"
        ),
    )

    ###########################################################################
    # Standardize and plot, matching the original MATLAB comparison
    ###########################################################################

    if args.plot_mode == "standardized":
        plot_values1 = standardize(shell_mean1)
        plot_values2 = standardize(shell_mean2)
        y_label = r"standardized $\langle I(q) \rangle$"
        title_suffix = "both rescaled by mean and sample standard deviation"
    else:
        plot_values1 = shell_mean1
        plot_values2 = shell_mean2
        y_label = r"$\langle I(q) \rangle$"
        title_suffix = "raw shell means"

    figure, axis = plt.subplots(figsize=(8, 6))

    axis.plot(
        q_centers,
        plot_values2,
        ":s",
        linewidth=1,
        markersize=6,
        label=label2,
    )
    axis.plot(
        q_centers,
        plot_values1,
        ":o",
        linewidth=1,
        markersize=6,
        label=label1,
    )

    axis.set_xlabel(r"$q$ [$\mathrm{\AA}^{-1}$]", fontsize=14)
    axis.set_ylabel(y_label, fontsize=14)
    axis.set_title(
        f"{label1} vs {label2}\n"
        f"{title_suffix}; r = {correlation:.6f}",
        fontsize=12,
    )
    axis.legend()
    axis.grid(True, alpha=0.3)

    figure.tight_layout()
    figure.savefig(
        figure_output,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(figure)

    ###########################################################################
    # Summary
    ###########################################################################

    print(f"First input: {file1}")
    print(f"Second input: {file2}")
    print(
        f"First dataset: {hkl1_raw.shape[0]} rows, "
        f"{hkl1.shape[0]} unique reflections"
    )
    print(
        f"Second dataset: {hkl2_raw.shape[0]} rows, "
        f"{hkl2.shape[0]} unique reflections"
    )
    print(f"Unique common reflections: {common_hkl.shape[0]}")
    print("Intersection alignment check: True")
    print(f"Correlation coefficient: {correlation:.10f}")
    print("Lattice type: primitive orthorhombic (P)")
    print(f"Unit cell: a={args.a} A, b={args.b} A, c={args.c} A")
    print("q convention: q = 1/d")
    print(f"Valid common q-shells: {q_centers.size} / {args.n_shells}")
    print(f"Zero intensities omitted from q-shell means: {not args.include_zero}")
    print(f"Common reflections saved to: {common_output}")
    print(f"Q-shell table saved to: {shell_output}")
    print(f"Figure saved to: {figure_output}")


if __name__ == "__main__":
    main()
