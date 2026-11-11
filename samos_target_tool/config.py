from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SAMOSGeometry:
    """Instrument/field constants used by the slit-mask allocator.

    The defaults reproduce the assumptions in the original SAMOS_GS_V5.7.py
    prototype: a 3 arcmin square DMD field sampled at 0.18 arcsec per mirror,
    mapped onto an occupancy grid with X across the field and Y along the
    dispersion direction.
    """

    field_width_arcsec: float = 180.0
    field_height_arcsec: float = 180.0
    image_scale_arcsec_per_pix: float = 0.18

    # Occupancy grid used for collision testing.
    detector_nx: int = 1100
    detector_ny: int = 4096
    x_center_pix: int = 550
    y_center_pix: int = 2048

    # Slit and spacing in DMD/image pixels.
    slit_width_pix: int = 3
    slit_length_pix: int = 9
    min_spacing_pix: int = 2

    # Approximate spectrum length along the detector/dispersion direction.
    spectrum_length_pix: int = 4000

    @property
    def half_field_x_pix(self) -> float:
        return 0.5 * self.field_width_arcsec / self.image_scale_arcsec_per_pix

    @property
    def half_field_y_pix(self) -> float:
        return 0.5 * self.field_height_arcsec / self.image_scale_arcsec_per_pix

    @property
    def exclusion_width_pix(self) -> int:
        """Center-to-center exclusion width for spectra in the cross-dispersion axis."""
        return int(self.slit_length_pix + self.min_spacing_pix + 2)


@dataclass(frozen=True)
class TargetToolConfig:
    n_masks: int = 5
    n_trials: int = 20
    random_seed: int = 42
    max_targets: int | None = 1000
    priority_column: str | None = None
    sort_column: str | None = "Gmag"
    lower_sort_is_better: bool = True
