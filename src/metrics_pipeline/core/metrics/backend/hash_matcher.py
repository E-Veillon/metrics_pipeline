"""Group Pymatgen structure-like objects from their composition."""

import typing as tp
import itertools as itt
import functools as ft

from pymatgen.core import Structure, IStructure, Composition, SiteCollection
from pymatgen.analysis.phase_diagram import PDEntry, Entry


HasComposition = (
    Composition | SiteCollection | Structure | IStructure | Entry | PDEntry
)

def hash_compositions(comp: HasComposition, by: tp.Literal["system", "formula"] = "system") -> int:
    """
    Generate a hash from the composition of a compatible pymatgen object, 
    i.e. a Composition object or an object having a `composition` attribute returning a
    Composition object.

    Parameters
    ----------
    comp: Entry|SiteCollection|Composition
        A pymatgen object containing a composition formula.

    by: "system", "formula"
        How to hash compositions:

        - "system" (default) returns same hash value for compositions containing the
        same elements, i.e. in the same chemical system. For example, a composition of formula
        "Fe2O3" will have the same hash as another one of formula "Fe3O4", but not one such as
        "LiFeO3".

        - "formula" returns same hash value for compositions having the exact same formula.

    Returns
    -------
    int
        The hash value.
    """
    if isinstance(comp, (Entry, SiteCollection)):
        comp = comp.composition

    match by:
        case "system":
            return hash(comp)
        case "formula":
            return hash(comp.reduced_formula)
        case str():
            raise ValueError(f"'by' only supports 'system' and 'formula', got {by!r}.")
        case _:
            raise TypeError(f"'by' expected a type 'str', got {type(by).__name__!r}.")


def group_compositions(
    comps: list[HasComposition],
    by: tp.Literal["system", "formula"] = "system"
) -> list[list[HasComposition]]:
    """
    Group Composition objects or objects having a `composition` attribute using the
    `hash_compositions` hash function (see the `by` argument for hashing algorithm).
    Note that objects from different classes but having their respective associated
    composition equal will be grouped together anyway.

    Parameters
    ----------
    comps: list[Composition|Entry|SiteCollection]
        The sequence of objects to group using their composition.

    by: "system", "formula"
        How to group compositions:

        - "system" (default) groups together compositions containing the
        same elements, i.e. in the same chemical system. For example, a composition of formula
        "Fe2O3" will be grouped with another one of formula "Fe3O4", but not with "LiFeO3".

        - "formula" groups together compositions having the exact same formula.

    Returns
    -------
    list[list[Composition|Entry|SiteCollection]]
        A list containing lists of objects with equivalent compositions.
    """
    assert all(isinstance(obj, Composition) or hasattr(obj, "composition") for obj in comps)

    hasher = ft.partial(hash_compositions, by=by)
    sorted_comps = sorted(comps, key=hasher)
    return [
        list(grouped)
        for _, grouped in itt.groupby(sorted_comps, hasher)
    ]
