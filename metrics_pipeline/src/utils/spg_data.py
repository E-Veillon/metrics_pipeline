#!/usr/bin/python
"""
Spacegroup symmetry classification data structures.
"""

SPG_NUM_TO_PG = {
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
}

PG_TO_SYSTEM = {
    **dict.fromkeys(("1", "-1"), "triclinic"),
    **dict.fromkeys(("2", "m", "2/m"), "monoclinic"),
    **dict.fromkeys(("222", "mm2", "mmm"), "orthorhombic"),
    **dict.fromkeys(("4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm"), "tetragonal"),
    **dict.fromkeys(("3", "-3", "32", "3m", "-3m"), "trigonal"),
    **dict.fromkeys(("6", "-6", "6/m", "622", "6mm", "-62m", "6/mmm"), "hexagonal"),
    **dict.fromkeys(("23", "m-3", "432", "-43m", "m-3m"), "cubic"),
}
