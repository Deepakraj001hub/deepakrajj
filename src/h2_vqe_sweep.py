"""Run and record VQE at the published O'Malley H2 bond distances."""

import csv
import json
from pathlib import Path
from typing import Any

from h2_table_i_data import TABLE_I
from h2_vqe import run_h2_vqe


CSV_COLUMNS = (
    "bond_distance_angstrom",
    "vqe_energy_hartree",
    "exact_energy_hartree",
    "absolute_error_hartree",
    "optimizer",
    "optimizer_configuration",
    "optimizer_evaluations",
    "optimizer_max_evaluations",
    "optimizer_success",
    "optimizer_status",
    "optimizer_message",
    "maximum_evaluations_reached",
    "run_error",
)


def _result_row(distance: float) -> dict[str, Any]:
    try:
        result = run_h2_vqe(distance)
    except Exception as exc:
        return {
            "bond_distance_angstrom": distance,
            "vqe_energy_hartree": "",
            "exact_energy_hartree": "",
            "absolute_error_hartree": "",
            "optimizer": "COBYLA",
            "optimizer_configuration": "maxiter=2000, tol=1e-5, rhobeg=1.0",
            "optimizer_evaluations": "",
            "optimizer_max_evaluations": 2000,
            "optimizer_success": False,
            "optimizer_status": "",
            "optimizer_message": "",
            "maximum_evaluations_reached": False,
            "run_error": f"{type(exc).__name__}: {exc}",
        }

    return {
        "bond_distance_angstrom": result.bond_distance_angstrom,
        "vqe_energy_hartree": result.vqe_energy,
        "exact_energy_hartree": result.exact_energy,
        "absolute_error_hartree": result.absolute_error,
        "optimizer": result.optimizer,
        "optimizer_configuration": result.optimizer_configuration,
        "optimizer_evaluations": result.optimizer_evaluations,
        "optimizer_max_evaluations": result.optimizer_max_evaluations,
        "optimizer_success": result.optimizer_success,
        "optimizer_status": result.optimizer_status,
        "optimizer_message": result.optimizer_message,
        "maximum_evaluations_reached": result.maximum_evaluations_reached,
        "run_error": "",
    }


def _summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    complete_rows = [
        row for row in rows
        if row["vqe_energy_hartree"] != "" and row["exact_energy_hartree"] != ""
    ]
    errors = [float(row["absolute_error_hartree"]) for row in complete_rows]
    max_error = max(complete_rows, key=lambda row: float(row["absolute_error_hartree"])) if complete_rows else None
    min_exact = min(complete_rows, key=lambda row: float(row["exact_energy_hartree"])) if complete_rows else None
    min_vqe = min(complete_rows, key=lambda row: float(row["vqe_energy_hartree"])) if complete_rows else None

    return {
        "source_point_count": len(TABLE_I),
        "completed_point_count": len(complete_rows),
        "run_error_count": len(rows) - len(complete_rows),
        "optimizer_success_count": sum(bool(row["optimizer_success"]) for row in rows),
        "optimizer_failure_count": sum(not bool(row["optimizer_success"]) for row in rows),
        "maximum_evaluation_termination_count": sum(
            bool(row["maximum_evaluations_reached"]) for row in rows
        ),
        "max_absolute_error_hartree": (
            float(max_error["absolute_error_hartree"]) if max_error else None
        ),
        "max_absolute_error_distance_angstrom": (
            float(max_error["bond_distance_angstrom"]) if max_error else None
        ),
        "mean_absolute_error_hartree": sum(errors) / len(errors) if errors else None,
        "minimum_exact_energy_distance_angstrom": (
            float(min_exact["bond_distance_angstrom"]) if min_exact else None
        ),
        "minimum_exact_energy_hartree": (
            float(min_exact["exact_energy_hartree"]) if min_exact else None
        ),
        "minimum_vqe_energy_distance_angstrom": (
            float(min_vqe["bond_distance_angstrom"]) if min_vqe else None
        ),
        "minimum_vqe_energy_hartree": (
            float(min_vqe["vqe_energy_hartree"]) if min_vqe else None
        ),
        "optimizer": "COBYLA",
        "optimizer_configuration": "maxiter=2000, tol=1e-5, rhobeg=1.0",
        "source_distances_angstrom": [row[0] for row in TABLE_I],
        "run_errors": [
            {
                "bond_distance_angstrom": row["bond_distance_angstrom"],
                "error": row["run_error"],
            }
            for row in rows if row["run_error"]
        ],
    }


def run_sweep(output_directory: Path | None = None) -> tuple[Path, Path]:
    """Run exactly the 54 published distances and write incremental results."""
    root = Path(__file__).resolve().parents[1]
    output_directory = output_directory or root / "results"
    output_directory.mkdir(parents=True, exist_ok=True)
    csv_path = output_directory / "h2_vqe_sweep.csv"
    summary_path = output_directory / "h2_vqe_summary.json"

    rows: list[dict[str, Any]] = []
    with csv_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        csv_file.flush()

        for index, (distance, *_) in enumerate(TABLE_I, start=1):
            row = _result_row(distance)
            rows.append(row)
            writer.writerow(row)
            csv_file.flush()
            print(
                f"[{index}/{len(TABLE_I)}] R={distance:.2f} Å "
                f"success={row['optimizer_success']} "
                f"evaluations={row['optimizer_evaluations']} "
                f"error={row['absolute_error_hartree']} "
                f"run_error={row['run_error']}",
                flush=True,
            )

    summary = _summarize(rows)
    summary_path.write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return csv_path, summary_path


if __name__ == "__main__":
    csv_output, summary_output = run_sweep()
    print(f"CSV: {csv_output}")
    print(f"Summary: {summary_output}")
