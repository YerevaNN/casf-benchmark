"""Reproduce a Kabsch-based false negative using the audit or repository matcher.

Run with Python >=3.10 and NumPy installed. This only tests synthetic feature
points; it neither changes nor estimates the error rate of the molecular runs.
"""

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

source = Path(__file__).resolve().parent / 'evidence/src/casf_benchmark/analysis/pharmacophore.py'
if not source.is_file():
    source = Path(__file__).resolve().parents[3] / 'src/casf_benchmark/analysis/pharmacophore.py'
spec = importlib.util.spec_from_file_location('audit_pharmacophore', source)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)

query = np.array([[-5., -5., 0.], [5., -5., 0.], [-5., 5., 0.], [5., 5., 0.]])
candidate = query.copy()
candidate[0, 0] += .9
types = ['HydrogenDonor', 'HydrogenAcceptor', 'Hydrophobic', 'PositiveIon']
rotation, translation, _ = module.proper_kabsch(candidate, query)
residuals = np.linalg.norm(candidate @ rotation + translation - query, axis=1)
feasible_translation = np.array([-.45, 0., 0.])
feasible_residuals = np.linalg.norm(candidate + feasible_translation - query, axis=1)
matched = module.exact_feature_correspondence_match(query, types, candidate, types, .5)
assert np.all(feasible_residuals < .5)
assert residuals.max() > .5002
assert not matched
print(json.dumps({
    'status': 'reproduced_false_negative',
    'scope': 'Synthetic typed feature points, not molecular-panel error-rate evidence',
    'radius_angstrom': .5,
    'numerical_tolerance_angstrom': .0002,
    'types': types,
    'query_coordinates': query.tolist(),
    'candidate_coordinates': candidate.tolist(),
    'feasible_translation': feasible_translation.tolist(),
    'feasible_residuals': feasible_residuals.tolist(),
    'kabsch_rotation': rotation.tolist(),
    'kabsch_translation': translation.tolist(),
    'kabsch_residuals': residuals.tolist(),
    'implemented_matcher_hit': matched,
}, indent=2))
