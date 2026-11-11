from __future__ import annotations

import numpy as np
import pandas as pd


def assign_priority_bins(
    df: pd.DataFrame,
    sort_column: str,
    n_bins: int = 5,
    lower_is_better: bool = True,
    output_column: str = "priority",
) -> pd.DataFrame:
    """Assign integer priority bins from a sortable catalog column.

    Priority 1 is the best/highest priority. For magnitudes, use
    lower_is_better=True. Rows with missing sort_column values are assigned to
    the lowest priority bin.
    """
    if sort_column not in df.columns:
        raise ValueError(f"Cannot assign priorities: missing column {sort_column!r}")

    out = df.copy()
    values = pd.to_numeric(out[sort_column], errors="coerce")
    valid = values.notna()
    out[output_column] = n_bins

    order = values[valid].sort_values(ascending=lower_is_better).index.to_numpy()
    if len(order) == 0:
        return out

    # Split into roughly equal ranked groups. np.array_split naturally handles
    # remainders and small catalogs.
    for i, idx in enumerate(np.array_split(order, n_bins), start=1):
        out.loc[idx, output_column] = i

    return out


def prepare_ranked_table(
    df: pd.DataFrame,
    priority_column: str = "priority",
    ra_column: str = "RA",
) -> pd.DataFrame:
    """Normalize ranking columns and initialize assignment bookkeeping."""
    out = df.copy()
    if priority_column not in out.columns:
        raise ValueError(f"Missing priority column {priority_column!r}")
    out[priority_column] = pd.to_numeric(out[priority_column], errors="coerce").fillna(5).astype(int)
    if "target_id" not in out.columns:
        out["target_id"] = np.arange(len(out), dtype=int)
    out["assigned_mask"] = -1
    out["rejection_reason"] = ""
    out = out.sort_values([priority_column, ra_column], ascending=[True, False]).reset_index(drop=True)
    return out
