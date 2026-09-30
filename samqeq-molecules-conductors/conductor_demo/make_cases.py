#!/usr/bin/env python3
"""Referee item 5: write the conductor demonstration as two regression cases in the release clone
(internal/tests/cases/conductor_water_{samqeq,qeq}). Single solves on the fixed Au(111) + 81 TIP4P-FQ film of
experiment I: samQEq (per-molecule constraints, self-energy of width 0.5 A on gold) at eta_Au = 0.05 eV, the
conductor limit, and conventional QEq (one global constraint, no self-energy) at 5.172 eV, where its operator is
still definite. Usage: python3 make_cases.py <clone>/internal/tests/cases"""
import sys
from pathlib import Path
import run_exp1

HERE = Path(__file__).resolve().parent
HEAD = {
 "samqeq": """# Au(111) slab (240 atoms) carrying 81 TIP4P-FQ water molecules, one charge solve. samQEq: one neutrality
# constraint per water molecule and one for the slab, and the finite-width self-energy (width 0.5 A) on the
# gold, at a gold hardness of 0.05 eV (the conductor limit). The operator stays positive definite; the water
# molecules stay neutral and the slab carries the image charge of the film.
""",
 "qeq": """# Au(111) slab (240 atoms) carrying 81 TIP4P-FQ water molecules, one charge solve. Conventional QEq: one
# global neutrality constraint and no self-energy, at the Thomas-Fermi gold hardness of 5.172 eV, where the
# operator is still definite. Charge moves between the water and the metal and between molecules; below a
# gold hardness of about 3.5 eV this operator is indefinite.
""",
}
CASES = {"samqeq": ("samqeq", "0.05"), "qeq": ("qeq", "5.172")}

out = Path(sys.argv[1])
for name, (arm, eta) in CASES.items():
    cons, gs = run_exp1.ARMS[arm]
    d = out / f"conductor_water_{name}"; d.mkdir(parents=True, exist_ok=True)
    deck = run_exp1.DECK.format(arm=arm, eta=eta, nometal="",
        cons="set             group all mol 1   # one global neutrality constraint" if cons == "global" else
             "# per-molecule neutrality: each water is its own fragment, the slab is one fragment",
        gself="fix_modify      chg gself on width 0.50 types 1" if gs else "", lmin="")
    deck = deck.replace("../../../data.exp1", "data.lmp").replace("{{", "{").replace("}}", "}")
    deck = deck.replace("write_dump      all custom charges.dump id mol type q x y z modify sort id format line "
                        "\"%d %d %d %.12g %.8f %.8f %.8f\"",
                        "write_dump      all custom charges.dump id type q modify sort id format line \"%d %d %.12g\"")
    lines = [l for l in deck.splitlines()[1:] if l.strip() or True]
    body = "\n".join(lines)
    while "\n\n\n" in body: body = body.replace("\n\n\n", "\n\n")
    (d / "in.lammps").write_text(HEAD[name] + body.lstrip("\n") + "\n")
    (d / "param").write_text(run_exp1.PARAM.format(eta=eta))
    data = (HERE / "data.exp1").read_text().splitlines()
    data[0] = "Au(111) 5-layer slab + 81 TIP4P-FQ water molecules"
    (d / "data.lmp").write_text("\n".join(data) + "\n")
    print("wrote", d)
