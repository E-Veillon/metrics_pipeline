"""Compute Earth Mover's Distance (or Wasserstein-1 Distance) metric."""

import warnings

from scipy.stats import wasserstein_distance

from pymatgen.core import Structure

from .metric_base import Metric, StructureDistribution


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
        transform: StructureDistribution | None = None,
        computed_property: str | None = None,
        **kwargs
    ) -> None:
        """
        Compute Earth Mover's Distance (or Wasserstein-1 Distance) metric.

        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Earth Mover's Distance metric on.

        ref_structs: list[Structure]
            Known structures to use as reference distribution.

        transform: StructureDistribution, optional
            Any callable taking a list of Structure objects and eventual keyword arguments
            and returning a numpy array representation of the structures 1D distribution
            of the computed property. If not given, `computed_property` must be given.

        computed_property: str, optional
            Structures property to compare. The property must be stored under this name
            as a numerical value (int or float) in structures `properties` attribute.
            If not given, `transform` must be given.

        kwargs: Any
            Additional keyword arguments to pass to `transform` (if used).

        Notes
        -----
        - If both `transform` and `computed_property` are given, `computed_property`
        is used in priority to save computations, and if it is not sufficient
        (i.e. not all structures have a valid value for the property), `transform`
        is used instead and stored property is ignored.
        """
        super().__init__(structures)

        def has_property(struct: Structure) -> bool:
            return struct.properties.get(computed_property) is not None

        if transform is None and computed_property is None:
            raise ValueError("Either 'transform' or 'computed_property' must be given.")

        if computed_property is not None:
            try:
                assert all(has_property(struct) for struct in self.structures), KeyError(
                    "Some computed structures do not have the property "
                    f"{computed_property!r} defined."
                )
                assert all(has_property(struct) for struct in ref_structs), KeyError(
                    "Some reference structures do not have the property "
                    f"{computed_property!r} defined."
                )
            except KeyError as exc:
                if transform is None:
                    raise exc
                else:
                    warn_msg = str(exc) + " Provided 'transform' callable will be used instead."
                    warnings.warn(warn_msg)
                    computed_property = None

        self.ref_structs = ref_structs
        self.transform = transform
        self.computed_property = computed_property
        self.kwargs = kwargs

        self._compute()

    def _compute(self) -> None:
        if self.computed_property is None:
            assert self.transform is not None, RuntimeError("Type checker assertion.")
            computed_values = self.transform(self.structures, **self.kwargs)
            ref_values = self.transform(self.ref_structs, **self.kwargs)

        else:
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