"""Source-only precursor balance audit; no optimization or phenotype values."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import cobra

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.map_quinone_pathway import met_record, reaction_record

BASE = 'results/ppnp_repair_2026_09_06/runs/main/both_forward/model.xml.gz'
BASE_HASH = 'b053283adf1257675b771439e99bb0dcaf2ee044c89e8df34096888c26d7f812'
PRECURSOR = ['4hbz_c', '4hbz_p', '4hbz_e', '4hbzcoa_c']
TAIL = ['3ophb_c', '2oph_c', '2ohph_c', '2omph_c', '2ombzl_c',
        '2ommbl_c', '2omhmbl_c', 'q8_c', 'q8h2_c']
COUMARATE = ['T4hcinnm_c', 'T4hcinnm_p', 'T4hcinnm_e', 'coucoa_c', '4hbald_c']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def combined_rows(model, metabolites):
    weights = {m: Fraction(1) for m in metabolites}
    result = {}
    for r in model.reactions:
        total = sum(Fraction.from_float(float(c)) * weights.get(m.id, 0) for m, c in r.metabolites.items())
        if total:
            result[r.id] = {'coefficient': float(total), 'exact': str(total), 'bounds': list(r.bounds)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('Preserve prior evidence; use a fresh output file')
    paths = [BASE, 'results/ppnp_repair_2026_09_06/runs/main/both_forward/gene_map.json',
             'data/reference/fitness_browser_media_bigg.tsv',
             'scripts/map_quinone_precursor.py', 'scripts/map_quinone_pathway.py']
    hashes = {p: sha(ROOT / p) for p in paths}
    if hashes[BASE] != BASE_HASH:
        raise ValueError('The declared PpnP parent changed')
    model = cobra.io.read_sbml_model(str(ROOT / BASE))
    mapping = json.loads((ROOT / paths[1]).read_text())['model_to_browser']
    ids = sorted(set(PRECURSOR + TAIL + COUMARATE))
    records = {r.id: reaction_record(r) for mid in ids for r in model.metabolites.get_by_id(mid).reactions}
    epsilon = -model.reactions.Growth.get_coefficient('q8h2_c')
    rows = {name: {'metabolites': mids, 'nonzero_rows': combined_rows(model, mids)} for name, mids in
            [('precursor', PRECURSOR), ('tail', TAIL), ('coumarate', COUMARATE),
             ('total', PRECURSOR + TAIL + COUMARATE)]}
    expected = {'CHRPL': 1.0, 'EX_4hbz_e': -1.0, 'EX_T4hcinnm_e': -1.0,
                '4HBHYOX': -1.0, 'sink_2ohph_c': -1.0, 'Growth': -epsilon}
    actual = {k: v['coefficient'] for k, v in rows['total']['nonzero_rows'].items()}
    if actual != expected:
        raise ValueError(f'Unexpected boundary in declared precursor rows: {actual}')
    for rid in ('4HBHYOX', 'sink_2ohph_c', 'CHRPL'):
        if model.reactions.get_by_id(rid).lower_bound != 0:
            raise ValueError('The necessary supply inequality requires nonnegative disposal fluxes')
    aliases = {}
    for locus in ('PP_5317', 'PP_5318', 'PP_1376', 'PP_1121', 'PP_4198'):
        hits = [g for g, value in mapping.items() if value == locus]
        if len(hits) != 1:
            raise ValueError(f'Ambiguous or absent diagnostic gene: {locus}')
        aliases[locus] = {'model_gene': hits[0], 'reactions': sorted(r.id for r in model.genes.get_by_id(hits[0]).reactions)}
    result = {'created_utc': datetime.now(timezone.utc).isoformat(), 'base_model': {'path': BASE, 'sha256': BASE_HASH, 'id': model.id},
        'evaluation_role': 'development diagnostic design; not an intervention or validation',
        'input_sha256': hashes, 'gene_aliases': aliases, 'quinone_biomass_coefficient_mmol_per_gDW': epsilon,
        'metabolites': {mid: met_record(model.metabolites.get_by_id(mid)) for mid in ids},
        'reactions': dict(sorted(records.items())), 'row_combinations': rows,
        'necessary_balance': 'CHRPL - EX_4hbz_e - EX_T4hcinnm_e = epsilon*Growth + 4HBHYOX + sink_2ohph_c',
        'implication': 'With CHRPL disabled and no coumarate import, net 4HBZ uptake must be at least epsilon*Growth. Equality is an LP feasibility question, not guaranteed by this row sum.',
        'transport_note': 'Both pre-existing 4HBZtex and UHBZ1t_pp are reversible. Secretion directions are permitted in the saved model, not experimentally established secretion rates.',
        'chemical_note': 'This is a declared equal-weight combination of stored rows, not atom-mapped proof of a chemical conserved moiety.',
        'units_note': 'Exchange bounds are fluxes per dry biomass per time. They are not extracellular concentrations; no concentration is inferred without a biomass/time balance.',
        'optimization_performed': False, 'numeric_fitness_values_accessed': False, 'model_intervention_performed': False}
    if hashes != {p: sha(ROOT / p) for p in paths}:
        raise ValueError('Source inputs changed during extraction')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'exact_precursor_balance_verified': True, 'epsilon': epsilon, 'recorded_reactions': len(records)}))


if __name__ == '__main__':
    main()
