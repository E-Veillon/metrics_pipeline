"""Compute Earth Mover's Distance (or Wasserstein-1 Distance) metric."""

from scipy.stats import wasserstein_distance

from pymatgen.core import Structure

from .metric_base import Metric


class EMD(Metric):
    """
    Compute Earth Mover's Distance (or Wasserstein-1 Distance) metric.

    Definition
    ----------
    Measure similarity between 1D distributions of one of the properties of
    reference and computed structures.

    Reference
    ---------
    Xie, T., Fu, X., Ganea, O., Barzilay, R., & Jaakkola, T. (2022).
    Crystal Diffusion Variational Autoencoder for Periodic Material Generation.
    Bulletin of the American Physical Society, 67.
    """
    def __init__(
        self,
        structures: list[Structure],
        ref_structs: list[Structure],
        computed_property: str
    ) -> None:
        """
        Compute Earth Mover's Distance (or Wasserstein-1 Distance) metric.

        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Earth Mover's Distance metric on.

        ref_structs: list[Structure]
            Known structures to use as reference distribution.

        computed_property: str
            Structures property to compare. The property must be stored under this name
            as a numerical value (int or float) in structures `properties` attribute.
        """
        super().__init__(structures)

        assert all(struct.properties.get(computed_property) is not None for struct in self.structures)
        assert all(struct.properties.get(computed_property) is not None for struct in ref_structs)

        self.ref_structs = ref_structs
        self.computed_property = computed_property

        self._compute()

    def _compute(self) -> None:
        computed_values = [struct.properties[self.computed_property] for struct in self.structures]
        ref_values = [struct.properties[self.computed_property] for struct in self.ref_structs]
        self._distance = wasserstein_distance(ref_values, computed_values)

    @property
    def computed_distance(self) -> float:
        """Get computed EMD value."""
        return self._distance

    def write_result(self, filename: str, decimals: int = 6) -> None:
        text = "===== Earth Mover's Distance Results ====="
        text += f"Total computed structures:  {len(self)}"
        text += f"Total reference structures: {len(self.ref_structs)}"
        text += f"EMD ({self.computed_property}): {self.computed_distance:.{decimals}f}"

        with open(filename, "wt", encoding="utf-8") as fp:
            fp.write("\n".join(text))