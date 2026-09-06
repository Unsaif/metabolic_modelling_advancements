"""Source-only exhaustive iJN1463 EX identifier/metabolite correspondence audit."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import cobra

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
MODEL = 'models/bigg/iJN1463.xml'
EXPECTED = 'd573833328ffae0dfa752a1fa3262ed939ed5862288beab287fca30d0fefb4a1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=BASE / 'ijn1463_exchange_alias_audit.json')
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('Preserve the existing source audit; use a fresh output file')
    def forbidden(*args, **kwargs):
        raise AssertionError('Source audit must not optimize')
    cobra.Model.optimize = forbidden
    cobra.Model.slim_optimize = forbidden
    inputs = [MODEL, str(Path(__file__).resolve().relative_to(ROOT)),
              str((BASE / 'medium_structure_inventory.json').relative_to(ROOT)),
              'data/reference/fitness_browser_media_bigg.tsv']
    before = {path: sha(ROOT / path) for path in inputs}
    if before[MODEL] != EXPECTED:
        raise ValueError('Declared source model changed')
    model = cobra.io.read_sbml_model(str(ROOT / MODEL))
    exchanges = sorted((r for r in model.reactions if r.id.startswith('EX_')), key=lambda r: r.id)
    mismatches, all_records = [], []
    for reaction in exchanges:
        record = {'reaction': reaction.id, 'name': reaction.name, 'bounds': list(reaction.bounds),
                  'gpr': reaction.gene_reaction_rule, 'annotation': reaction.annotation,
                  'metabolites': [{'id': m.id, 'coefficient': c, 'name': m.name, 'compartment': m.compartment,
                      'formula': m.formula, 'charge': m.charge, 'annotation': m.annotation}
                      for m, c in reaction.metabolites.items()]}
        all_records.append({'reaction': reaction.id,
                            'metabolites': {m.id: c for m, c in reaction.metabolites.items()}})
        if len(reaction.metabolites) != 1 or next(iter(reaction.metabolites)).id != reaction.id[3:]:
            mismatches.append(record)
    inventory = json.loads((BASE / 'medium_structure_inventory.json').read_text())
    parent = next(item for item in inventory['models'] if item['label'] == 'Putida_iJN1463')
    carbon = {exchange for condition in parent['conditions'] for exchange in condition['exchanges']}
    with (ROOT / 'data/reference/fitness_browser_media_bigg.tsv').open() as handle:
        media = [r for r in csv.DictReader(handle, delimiter='\t') if r['media'] in parent['media']]
    pairs = []
    for record in mismatches:
        if len(record['metabolites']) != 1:
            raise ValueError('Unexpected multi-metabolite exchange; requires separate review')
        metabolite = record['metabolites'][0]['id']
        candidates = {record['reaction'], 'EX_' + metabolite}
        components = {record['reaction'][3:-2], metabolite[:-2]}
        pairs.append({'reaction': record['reaction'], 'metabolite': metabolite,
                      'alternative_canonical_reaction_exists': 'EX_' + metabolite in model.reactions,
                      'reaction_stem_metabolite_exists': record['reaction'][3:] in model.metabolites,
                      'other_exchanges_on_same_metabolite': [r.id for r in exchanges if r.id != record['reaction'] and
                                                           model.metabolites.get_by_id(metabolite) in r.metabolites],
                      'requested_carbon_aliases': sorted(candidates & carbon),
                      'base_medium_component_matches': {r['media']: sorted(components & set(r['bigg_ids'].split(';')))
                                                        for r in media}})
    result = {'scope': 'Exhaustive stored EX ID-to-metabolite comparison; no optimization, fitness parsing or model amendment.',
              'model': MODEL, 'model_sha256': EXPECTED, 'cobra_version': cobra.__version__,
              'n_EX_reactions': len(exchanges), 'n_id_metabolite_mismatches': len(mismatches),
              'all_EX_single_negative_unit': all(len(r.metabolites) == 1 and next(iter(r.metabolites.values())) == -1 for r in exchanges),
              'mismatches': mismatches, 'identity_and_medium_context': pairs,
              'all_checked_EX_id_metabolite_pairs': all_records, 'input_sha256': before,
              'source_conclusion': 'EX_AEP_e is the only mismatch: its stored substrate is 2ameph_e (2-aminoethylphosphonate). A reaction-specific identity adapter can preserve the existing reaction ID and chemistry; this is not a nutrient addition or evidence of transport activity.',
              'scope_refinement': 'The initial inventory checked one-metabolite negative-unit shapes but did not test reaction-stem identity. The main preparation guard abstained on this unsupported representation before any curated-model solve. Earlier source files and the failed primary attempt are preserved.'}
    if before != {path: sha(ROOT / path) for path in inputs}:
        raise ValueError('Source changed during audit')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')
    print(json.dumps({'n_EX': len(exchanges), 'n_mismatches': len(mismatches), 'mismatch_ids': [r['reaction'] for r in mismatches]}))


if __name__ == '__main__':
    main()
