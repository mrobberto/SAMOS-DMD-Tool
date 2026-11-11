from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd

from .config import SAMOSGeometry, TargetToolConfig
from .ranking import prepare_ranked_table


@dataclass(frozen=True)
class OptimizationSummary:
    n_input: int
    n_in_field: int
    n_assigned: int
    assigned_by_mask: dict[int, int]
    assigned_by_priority: dict[int, int]


class MaskOptimizer:
    """Randomized greedy SAMOS slit-mask optimizer.

    The algorithm is intentionally close to the original prototype: for each
    mask it tries multiple randomized orderings of the top unassigned priority
    class, then greedily accepts targets whose spectral footprints do not
    overlap in an occupancy grid.
    """

    def __init__(self, geometry: SAMOSGeometry, config: TargetToolConfig):
        self.geometry = geometry
        self.config = config

    def optimize(self, df: pd.DataFrame) -> tuple[pd.DataFrame, OptimizationSummary]:
        work = df.copy()
        if "in_field" in work.columns:
            work = work[work["in_field"]].copy()
        work = prepare_ranked_table(work, priority_column="priority")

        if self.config.max_targets is not None:
            work = work.iloc[: self.config.max_targets].copy()

        rng = np.random.default_rng(self.config.random_seed)

        for mask_index in range(self.config.n_masks):
            best_table = None
            best_score = (-1, -1)  # priority-1 assigned, total assigned in this mask

            for _trial in range(self.config.n_trials):
                trial = self._shuffle_current_top_priority(work, rng)
                accepted = self._fill_one_mask(trial, mask_index)

                in_this = accepted["assigned_mask"].eq(mask_index)
                n_priority1 = int((in_this & accepted["priority"].eq(1)).sum())
                n_total = int(in_this.sum())
                score = (n_priority1, n_total)

                if score >= best_score:
                    best_score = score
                    best_table = accepted

            if best_table is None:
                break
            work = best_table.sort_values(["priority", "RA"], ascending=[True, False]).reset_index(drop=True)

        assigned = work[work["assigned_mask"] >= 0]
        summary = OptimizationSummary(
            n_input=int(len(df)),
            n_in_field=int(len(work)),
            n_assigned=int(len(assigned)),
            assigned_by_mask=assigned.groupby("assigned_mask").size().astype(int).to_dict(),
            assigned_by_priority=assigned.groupby("priority").size().astype(int).to_dict(),
        )
        return work, summary

    def _shuffle_current_top_priority(self, df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
        out = df.copy().reset_index(drop=True)
        unassigned = out["assigned_mask"].eq(-1)
        if not unassigned.any():
            return out
        top_priority = int(out.loc[unassigned, "priority"].min())
        # Pandas 2.x is strict about assigning heterogeneous row blocks back
        # into columns with narrow dtypes inherited from Astropy/Gaia tables
        # (for example int16 columns).  Instead of replacing row values in
        # place, build a new row order and reindex.  This preserves column
        # dtypes and avoids lossy setitem/upcast errors.
        idx = out.index[unassigned & out["priority"].eq(top_priority)].to_numpy()
        shuffled = idx.copy()
        rng.shuffle(shuffled)

        order = np.arange(len(out))
        order[idx] = shuffled
        return out.iloc[order].reset_index(drop=True)

    def _fill_one_mask(self, df: pd.DataFrame, mask_index: int) -> pd.DataFrame:
        g = self.geometry
        out = df.copy().reset_index(drop=True)
        occupancy = np.zeros((g.detector_ny, g.detector_nx), dtype=np.uint8)

        for j, row in out.iterrows():
            if int(row["assigned_mask"]) != -1:
                continue

            x_pix = float(row["x_pix"])
            y_pix = float(row["y_pix"])

            if abs(x_pix) > g.half_field_x_pix or abs(y_pix) > g.half_field_y_pix:
                out.at[j, "rejection_reason"] = "outside_field"
                continue

            # Exclusion rectangle: broad footprint used only for collision check.
            exclusion_width = g.exclusion_width_pix
            x0 = int(round(x_pix - exclusion_width / 2 + g.x_center_pix))
            x1 = int(round(x_pix + exclusion_width / 2 + g.x_center_pix + 1))
            y0 = int(round(y_pix - g.spectrum_length_pix / 2 + g.y_center_pix))
            y1 = int(round(y_pix + g.spectrum_length_pix / 2 + g.y_center_pix))

            x0 = max(x0, 0)
            x1 = min(x1, g.detector_nx)
            y0 = max(y0, 0)
            y1 = min(y1, g.detector_ny)
            if x1 <= x0 or y1 <= y0:
                out.at[j, "rejection_reason"] = "empty_footprint"
                continue

            if occupancy[y0:y1, x0:x1].any():
                out.at[j, "rejection_reason"] = "spectral_overlap"
                continue

            # Accepted slit footprint: narrower than exclusion rectangle, as in
            # the original script where two guard pixels are removed on each side.
            slit_x0 = min(max(x0 + 2, 0), g.detector_nx)
            slit_x1 = min(max(x1 - 2, 0), g.detector_nx)
            if slit_x1 <= slit_x0:
                slit_x0, slit_x1 = x0, x1

            occupancy[y0:y1, slit_x0:slit_x1] = 1
            out.at[j, "assigned_mask"] = mask_index
            out.at[j, "rejection_reason"] = ""

        return out
