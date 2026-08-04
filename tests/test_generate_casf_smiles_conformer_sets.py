from __future__ import annotations

import pandas as pd
from rdkit import Chem

from casf_benchmark.generation import conformer_sets as generator




def mol_with_conformer() -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
    mol.RemoveAllConformers()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for atom_idx in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(atom_idx, (float(atom_idx), 0.0, 0.0))
    mol.AddConformer(conf, assignId=True)
    return mol


def test_posebusters_geometry_passes_isolates_runtime_failure(monkeypatch):
    calls = {"single": 0}

    def fake_run(records, reference_mol, max_workers, chunk_size, energy_num_threads):
        if len(records) > 1:
            raise RuntimeError("batch failed")
        calls["single"] += 1
        if calls["single"] == 2:
            raise RuntimeError("Cannot normalize a zero length vector")
        return generator.PoseBustersResult([True], [{"identity": True}])

    monkeypatch.setattr(generator, "_run_posebusters_records", fake_run)
    mols = [mol_with_conformer(), mol_with_conformer(), mol_with_conformer()]

    result = generator.posebusters_geometry_passes(
        mols,
        mols[0],
        max_workers=3,
        energy_num_threads=1,
    )

    assert result.passes == [True, False, True]
    assert result.check_rows[1]["posebusters_runtime_error"] is False
    assert result.check_rows[1]["posebusters_exception_RuntimeError"] is False
