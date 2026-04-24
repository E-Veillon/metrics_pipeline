"""Utility functions to process phase diagram entries."""

import typing as tp

from pymatgen.core import Element, Species, DummySpecies, Composition
from pymatgen.analysis.phase_diagram import PDEntry

from src.utils import flatten


def get_elements_from_entries(entries: tp.Sequence[PDEntry]) -> list[Element | Species | DummySpecies]:
    """
    Get a list of all unique elements present in a sequence of PDEntry objects.

    Parameters
    ----------
    entries: Sequence[PDEntry]
        The entries to get unique elements from.

    Returns
    -------
    list[Element | Species | DummySpecies]
        The list of found Element objects.
    """
    if not entries:
        return []

    if len(entries) == 1:
        return entries[0].elements

    return list(set(flatten([entry.elements for entry in entries])))


def get_lacking_elts_entries(
    entries: tp.Sequence[PDEntry],
    ref_elts: tp.Sequence[Element] | set[Element]
) -> list[PDEntry]:
    """
    Check elemental entries with respect to given reference elements,
    then initialize lacking elemental entries with an energy of 0.0 eV.

    Parameter
    ---------
    entries: Sequence[PDEntry]
        Entries to check elemental entries in.

    ref_elts: Sequence[Element] | set[Element]
        Reference Element objects defining the wanted chemical space.

    Returns
    -------
    list[PDEntry]
        A list of auto-defined elemental entries.
    """
    if not entries:
        user_elt_entries = []
    else:
        user_elt_entries = list(
            filter(
                lambda entry: entry.is_element and entry.elements[0] in ref_elts,
                entries
            )
        )
    user_defined_elts = get_elements_from_entries(user_elt_entries)
    lacking_elts = list(
        filter(
            lambda elt: elt not in user_defined_elts,
            ref_elts
        )
    )
    if not lacking_elts: # No lacking entry
        return []

    auto_defined_elts_entries = [
        PDEntry(
            composition=Composition(str(elt)),
            energy=0.0,
            name=elt.symbol,
            attribute="element_ref",
        ) for elt in lacking_elts
    ]

    return auto_defined_elts_entries