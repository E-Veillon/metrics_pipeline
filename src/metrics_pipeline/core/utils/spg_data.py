"""
Crystallographic spacegroup symmetry classification data.
"""

import typing_extensions as tpe
from collections import OrderedDict

ALL_CRYSTAL_FAMILIES = ("triclinic", "monoclinic", "orthorhombic", "tetragonal", "hexagonal", "cubic")
ALL_CRYSTAL_SYSTEMS = ("triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic")
ALL_POINT_GROUPS = (
    # Triclinic
    "1", "-1",
    # Monoclinic
    "2", "m", "2/m",
    # Orthorhombic
    "222", "mm2", "mmm",
    # Tetragonal
    "4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm",
    # Trigonal
    "3", "-3", "32", "3m", "-3m",
    # Hexagonal
    "6", "-6", "6/m", "622", "6mm", "-62m", "6/mmm",
    # Cubic
    "23", "m-3", "432", "-43m", "m-3m"
)
ALL_SPACEGROUPS = (
    # Triclinic
    "P1", "P-1",
    # Monoclinic PG "2"
    "P2", "P2_1", "C2",
    # Monoclinic PG "m"
    "Pm", "Pc", "Cm", "Cc",
    # Monoclinic PG "2/m"
    "P2/m", "P2_1/m", "C2/m", "P2/c", "P2_1/c", "C2/c",
    # Orthorhombic PG "222"
    "P222", "P222_1", "P2_12_12", "P2_12_12_1", "C222_1", "C222", "F222", "I222", "I2_12_12_1",
    # Orthorhombic PG "mm2"
    "Pmm2", "Pmc2_1", "Pcc2", "Pma2", "Pca2_1", "Pnc2", "Pmn2_1", "Pba2", "Pna2_1", "Pnn2",
    "Cmm2", "Cmc2_1", "Ccc2", "Amm2", "Aem2", "Ama2", "Aea2", "Fmm2", "Fdd2", "Imm2", "Iba2",
    "Ima2",
    # Orthorhombic PG "mmm"
    "Pmmm", "Pnnn", "Pccm", "Pban", "Pmma", "Pnna", "Pmna", "Pcca", "Pbam", "Pccn", "Pbcm", "Pnnm",
    "Pmmn", "Pbcn", "Pbca", "Pnma", "Cmcm", "Cmce", "Cmmm", "Cccm", "Cmme", "Ccce", "Fmmm", "Fddd",
    "Immm", "Ibam", "Ibca", "Imma",
    # Tetragonal PG "4"
    "P4", "P4_1", "P4_2", "P4_3", "I4", "I4_1",
    # Tetragonal PG "-4"
    "P-4", "I-4",
    # Tetragonal PG "4/m"
    "P4/m", "P4_2/m", "P4/n", "P4_2/n", "I4/m", "I4_1/a",
    # Tetragonal PG "422"
    "P422", "P42_12", "P4_122", "P4_12_12", "P4_222", "P4_22_12", "P4_322", "P4_32_12", "I422",
    "I4_122",
    # Tetragonal PG "4mm"
    "P4mm", "P4bm", "P4_2cm", "P4_2nm", "P4cc", "P4nc", "P4_2mc", "P4_2bc", "I4mm", "I4cm",
    "I4_1md", "I4_1cd",
    # Tetragonal PG "-42m"
    "P-42m", "P-42c", "P-42_1m", "P-42_1c", "P-4m2", "P-4c2", "P-4b2", "P-4n2", "I-4m2", "I-4c2",
    "I-42m", "I-42d",
    # Tetragonal PG "4/mmm"
    "P4/mmm", "P4/mcc", "P4/nbm", "P4/nnc", "P4/mbm", "P4/mnc", "P4/nmm",
    "P4/ncc", "P4_2/mmc", "P4_2/mcm", "P4_2/nbc", "P4_2/nnm", "P4_2/mbc", "P4_2/mnm", "P4_2/nmc",
    "P4_2/ncm", "I4/mmm", "I4/mcm", "I4_1/amd", "I4_1/acd",
    # Trigonal PG "3"
    "P3", "P3_1", "P3_2", "R3",
    # Trigonal PG "-3"
    "P-3", "R-3",
    # Trigonal PG "32"
    "P312", "P321", "P3_112", "P3_121", "P3_212", "P3_221", "R32",
    # Trigonal PG "3m"
    "P3m1", "P31m", "P3c1", "P31c", "R3m", "R3c",
    # Trigonal PG "-3m"
    "P-31m", "P-31c", "P-3m1", "P-3c1", "R-3m", "R-3c",
    # Hexagonal PG "6"
    "P6", "P6_1", "P6_5", "P6_2", "P6_4", "P6_3",
    # Hexagonal PG "-6"
    "P-6",
    # Hexagonal PG "6/m"
    "P6/m", "P6_3/m",
    # Hexagonal PG "622"
    "P622", "P6_122", "P6_522", "P6_222", "P6_422", "P6_322",
    # Hexagonal PG "6mm"
    "P6mm", "P6cc", "P6_3cm", "P6_3mc",
    # Hexagonal PG "-6m2"
    "P-6m2", "P-6c2", "P-62m", "P-62c",
    # Hexagonal PG "6/mmm"
    "P6/mmm", "P6/mcc", "P6_3/mcm", "P6_3/mmc",
    # Cubic PG "23"
    "P23", "F23", "I23", "P2_13", "I2_13",
    # Cubic PG "m-3"
    "Pm-3", "Pn-3", "Fm-3", "Fd-3", "Im-3", "Pa-3", "Ia-3",
    # Cubic PG "432"
    "P432", "P4_232", "F432", "F4_132", "I432", "P4_332", "P4_132", "I4_132",
    # Cubic PG "-43m"
    "P-43m", "F-43m", "I-43m", "P-43n", "F-43c", "I-43d",
    # Cubic PG "m-3m"
    "Pm-3m", "Pn-3n", "Pm-3n", "Pn-3m", "Fm-3m", "Fm-3c", "Fd-3m", "Fd-3c", "Im-3m", "Ia-3d"
)
ALL_SYMMETRY_CLASSES = {
    "CRYSTAL_FAMILIES": ALL_CRYSTAL_FAMILIES,
    "CRYSTAL_SYSTEMS": ALL_CRYSTAL_SYSTEMS,
    "POINT_GROUPS": ALL_POINT_GROUPS,
    "SPACEGROUPS": ALL_SPACEGROUPS,
}
SPG_NUM_TO_PG: OrderedDict[int, str] = OrderedDict({
    # Triclinic
    1: "1", 2: "-1",
    # Monoclinic
    **dict.fromkeys(list(range(3, 6)), "2"),
    **dict.fromkeys(list(range(6, 10)), "m"),
    **dict.fromkeys(list(range(10, 16)), "2/m"),
    # Orthorhombic
    **dict.fromkeys(list(range(16, 25)), "222"),
    **dict.fromkeys(list(range(25, 47)), "mm2"),
    **dict.fromkeys(list(range(47, 75)), "mmm"),
    # Tetragonal
    **dict.fromkeys(list(range(75, 81)), "4"),
    **dict.fromkeys(list(range(81, 83)), "-4"),
    **dict.fromkeys(list(range(83, 89)), "4/m"),
    **dict.fromkeys(list(range(89, 99)), "422"),
    **dict.fromkeys(list(range(99, 111)), "4mm"),
    **dict.fromkeys(list(range(111, 123)), "-42m"),
    **dict.fromkeys(list(range(123, 143)), "4/mmm"),
    # Trigonal
    **dict.fromkeys(list(range(143, 147)), "3"),
    **dict.fromkeys(list(range(147, 149)), "-3"),
    **dict.fromkeys(list(range(149, 156)), "32"),
    **dict.fromkeys(list(range(156, 162)), "3m"),
    **dict.fromkeys(list(range(162, 168)), "-3m"),
    # Hexagonal
    **dict.fromkeys(list(range(168, 174)), "6"),
    **dict.fromkeys(list(range(174, 175)), "-6"),
    **dict.fromkeys(list(range(175, 177)), "6/m"),
    **dict.fromkeys(list(range(177, 183)), "622"),
    **dict.fromkeys(list(range(183, 187)), "6mm"),
    **dict.fromkeys(list(range(187, 191)), "-6m2"),
    **dict.fromkeys(list(range(191, 195)), "6/mmm"),
    # Cubic
    **dict.fromkeys(list(range(195, 200)), "23"),
    **dict.fromkeys(list(range(200, 207)), "m-3"),
    **dict.fromkeys(list(range(207, 215)), "432"),
    **dict.fromkeys(list(range(215, 221)), "-43m"),
    **dict.fromkeys(list(range(221, 231)), "m-3m"),
})
PG_TO_SYSTEM: OrderedDict[str, str] = OrderedDict({
    **dict.fromkeys(("1", "-1"), "triclinic"),
    **dict.fromkeys(("2", "m", "2/m"), "monoclinic"),
    **dict.fromkeys(("222", "mm2", "mmm"), "orthorhombic"),
    **dict.fromkeys(("4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm"), "tetragonal"),
    **dict.fromkeys(("3", "-3", "32", "3m", "-3m"), "trigonal"),
    **dict.fromkeys(("6", "-6", "6/m", "622", "6mm", "-62m", "6/mmm"), "hexagonal"),
    **dict.fromkeys(("23", "m-3", "432", "-43m", "m-3m"), "cubic"),
})

class Spacegroup:
    """
    A crystallographic spacegroup with its classification data.

    Attributes
    ----------
    int_number: int
        The spacegroup international number, between 0 and 230,
        with 0 representing an undetermined spacegroup.
    symbol: str
        The spacegroup Hermann-Mauguin symbol.
    point_group: str
        The point group symbol corresponding to this spacegroup.
    crystal_system: str
        The crystal system corresponding to this spacegroup.
    crystal_family: str
        The crystal family corresponding to this spacegroup
        (same as crystal system except trigonal system is of hexagonal family).
    """
    def __init__(self, int_number: int, *, is_uncomputable: bool = False) -> None:
        """
        A crystallographic spacegroup with its classification data.

        Parameters
        ----------
        int_number: int
            The spacegroup international number, between 0 and 230,
            with 0 representing an undetermined spacegroup.
        is_uncomputable: bool
            Whether the spacegroup is uncomputable. Defaults to False.
        """
        if not isinstance(int_number, int):
            raise TypeError(f"Spacegroup number must be an integer, got {type(int_number).__name__!r}.")
        if not 0 <= int_number <= 230:
            raise ValueError(f"Spacegroup number must be between 0 and 230, got {int_number}.")

        self._int_number = int_number
        self._is_uncomputable = is_uncomputable

    @property
    def int_number(self) -> int:
        """Get the international number of the spacegroup."""
        return self._int_number

    @int_number.setter
    def int_number(self, value: int) -> None:
        """Set the international number of the spacegroup."""
        if not isinstance(value, int):
            raise TypeError(f"Spacegroup number must be an integer, got {type(value).__name__!r}.")
        if self.is_uncomputable:
            raise ValueError(
                "This spacegroup is set as uncomputable and cannot be modified unless"
                "'is_uncomputable' is set to False."
            )
        if not 0 <= value <= 230:
            raise ValueError(f"Spacegroup number must be between 0 and 230, got {value}.")
        self._int_number = value

    @property
    def symbol(self) -> str:
        """Get the Hermann-Mauguin symbol of the spacegroup."""
        if self.int_number == 0:
            return "Undetermined"
        return ALL_SPACEGROUPS[self.int_number - 1]

    @property
    def point_group(self) -> str:
        """Get the point group corresponding to this spacegroup."""
        if self.int_number == 0:
            return "Undetermined"
        return SPG_NUM_TO_PG[self.int_number]
    
    @property
    def crystal_system(self) -> str:
        """Get the crystal system corresponding to this spacegroup."""
        if self.int_number == 0:
            return "Undetermined"
        return PG_TO_SYSTEM[self.point_group]
    
    @property
    def crystal_family(self) -> str:
        """Get the crystal family corresponding to this spacegroup."""
        if self.crystal_system == "trigonal":
            return "hexagonal"
        return self.crystal_system
    
    @property
    def is_uncomputable(self) -> bool:
        """
        Whether the spacegroup is uncomputable
        (i.e., SpacegroupAnalyzer could not determine a spacegroup).
        """
        return self._is_uncomputable
    
    @is_uncomputable.setter
    def is_uncomputable(self, value: bool) -> None:
        """
        Set the spacegroup as uncomputable (i.e., SpacegroupAnalyzer could not determine
        a spacegroup).
        Setting it to `True` automatically freezes the spacegroup number to 0 (undetermined).
        """
        if not isinstance(value, bool):
            raise TypeError(f"is_uncomputable must be a boolean, got {type(value).__name__!r}.")
        self._is_uncomputable = value
        if value:
            self._int_number = 0

    def as_dict(self, verbose: bool = False) -> dict:
        """
        Return a dictionary representation of the spacegroup.
        
        Parameters
        ----------
        verbose: bool
            Whether to add classification details to the dict.
            The default (False) only saves necessary data to rebuild the object.
        """
        dct = {
            "int_number": self.int_number,
            "is_uncomputable": self.is_uncomputable
        }
        if verbose:
            dct.update({
                "symbol": self.symbol,
                "point_group": self.point_group,
                "crystal_system": self.crystal_system,
                "crystal_family": self.crystal_family
            })
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> tpe.Self:
        """Create a Spacegroup object from a dictionary representation."""
        return cls(int_number=dct["int_number"], is_uncomputable=dct["is_uncomputable"])