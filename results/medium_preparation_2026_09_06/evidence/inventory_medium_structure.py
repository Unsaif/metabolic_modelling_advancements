"""Source-only inventories and parent loaders for declared development models.

Public API: build_parent(label) -> (model, provenance); metadata_conditions(org)
-> list of mapped Condition objects. Neither helper optimizes or reads numerical
fitness values. CLI additionally inventories inherited medium completion in memory.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import logging
from pathlib import Path
import sys
from types import SimpleNamespace

import cobra
import pandas as pd
from cobra.medium import find_external_compartment

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
sys.path.insert(0, str(ROOT))
from gembench.fitness_browser import base_medium, carbon_source_conditions
from gembench.gene_mapping import build_gene_map
from gembench.patches import apply_gpr_patches, apply_model_patches, apply_universe_patches, load_patch_files
from gembench.protocols.carbon_fitness_generic import GenericParams, complete_medium_transport

LABELS = ('Putida', 'Btheta', 'MR1', 'Smeli', 'Putida_iJN1463')
PUTIDA_PARENT = 'results/ppnp_repair_2026_09_06/runs/main/both_forward/model.xml.gz'
PUTIDA_HASH = 'b053283adf1257675b771439e99bb0dcaf2ee044c89e8df34096888c26d7f812'
PATCHES = ['data/reference/universe_patches_v0.1.json', 'data/reference/model_patches_v0.4.json',
           'data/reference/gpr_patches_v0.2.json', 'data/reference/gpr_patches_v0.3.json',
           'data/reference/gpr_patches_v0.4.json']
CODE = ['gembench/patches.py', 'gembench/gene_mapping.py', 'gembench/fitness_browser.py',
        'gembench/media.py', 'gembench/protocols/carbon_fitness_generic.py',
        str(Path(__file__).resolve().relative_to(ROOT))]


def require(test, message):
    if not test:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hashes(paths):
    return {str(path): sha(ROOT / path) for path in sorted(set(paths))}


def metadata_paths(org):
    require(org in {'Putida', 'Btheta', 'MR1', 'Smeli'}, 'Only declared development metadata')
    prefix = f'data/fitness_browser/{org}/'
    return [prefix + 'genes.tsv', prefix + 'experiments.tsv', prefix + 'fit_logratios.tsv',
            'data/reference/fitness_browser_carbon_sources_bigg.tsv',
            'data/reference/fitness_browser_media_bigg.tsv']


def metadata_conditions(org):
    """Read experimental metadata and the fitness table HEADER ONLY."""
    paths = metadata_paths(org)
    experiments = pd.read_table(ROOT / paths[1], dtype=str, keep_default_na=False)
    with (ROOT / paths[2]).open() as handle:
        header = next(csv.reader(handle, delimiter='\t'))
    meta = {'orgId', 'locusId', 'sysName', 'geneName', 'desc'}
    columns = [field.split(' ')[0] for field in header if field not in meta]
    org_metadata = SimpleNamespace(experiments=experiments, fitness=SimpleNamespace(columns=columns))
    return [condition for condition in carbon_source_conditions(org_metadata) if condition.bigg_ids]


def build_parent(label):
    """Reproduce a source model before completion; never call the benchmark."""
    require(label in LABELS, f'Undeclared model label: {label}')
    inputs, applied = list(CODE), []
    if label == 'Putida':
        path = PUTIDA_PARENT
        require(sha(ROOT / path) == PUTIDA_HASH, 'Saved Putida parent changed')
        mode = 'exact_saved_PpnP_both_forward_model'
        inputs += [path, 'results/ppnp_repair_2026_09_06/runs/main/manifest.json']
        model = cobra.io.read_sbml_model(str(ROOT / path))
    elif label == 'Putida_iJN1463':
        path = 'models/bigg/iJN1463.xml'
        mode = 'unmodified_local_curated_structural_reference_already_exposed_Putida'
        inputs.append(path)
        model = cobra.io.read_sbml_model(str(ROOT / path))
    else:
        path = f'models/gapfilled/{label}.xml.gz'
        invocation = f'results/development_sprint_2026_09_06/runs/{label}_cycle7/invocation.json'
        command = json.loads((ROOT / invocation).read_text())['command']
        for flag, expected in [('--orgs', label), ('--variant', 'gapfilled'),
                               ('--patch', PATCHES[0]), ('--model-patch', PATCHES[1]),
                               ('--gpr-patch', ','.join(PATCHES[2:]))]:
            require(command[command.index(flag) + 1] == expected, 'Cycle7 recipe mismatch')
        require('--complete-medium-transport' in command and '--no-drop-rich' in command,
                'Unexpected cycle7 preparation flags')
        inputs += [path, invocation, f'data/fitness_browser/{label}/genes.tsv',
                   f'data/genpept/{label}_genpept_map.tsv', *PATCHES,
                   'results/development_sprint_2026_09_06/prespecified_manifest.json']
        with gzip.open(ROOT / path, 'rt') as handle:
            model = cobra.io.read_sbml_model(handle)
        genes = pd.read_table(ROOT / f'data/fitness_browser/{label}/genes.tsv', dtype=str, keep_default_na=False)
        universe = json.loads((ROOT / PATCHES[0]).read_text())
        applied += apply_universe_patches(model, label, universe)
        gene_map = build_gene_map(label, [g.id for g in model.genes], set(genes.sysName))
        changes = apply_model_patches(model, label, json.loads((ROOT / PATCHES[1]).read_text()), gene_map, verbose=False)
        applied += changes
        if changes:
            gene_map = build_gene_map(label, [g.id for g in model.genes], set(genes.sysName))
        applied += apply_gpr_patches(model, label, gene_map,
                                     load_patch_files([str(ROOT / p) for p in PATCHES[2:]]), verbose=False)
        mode = 'reconstructed_cycle7_before_medium_completion_no_serialized_cycle7_SBML_available'
    provenance = {'label': label, 'mode': mode, 'source_path': path, 'source_sha256': sha(ROOT / path),
                  'applied_historical_changes': applied, 'input_sha256': hashes(inputs),
                  'excluded_operations': ['optimization', 'energy probe', 'knockout simulation', 'fitness value parsing',
                                          'metric calculation', 'medium completion', 'new patch selection']}
    return model, provenance


def reaction_record(reaction):
    return {'id': reaction.id, 'stoichiometry': {m.id: c for m, c in sorted(reaction.metabolites.items(), key=lambda x: x[0].id)},
            'bounds': list(reaction.bounds), 'compartments': sorted(reaction.compartments),
            'gpr': reaction.gene_reaction_rule, 'sbo': reaction.annotation.get('sbo')}


def inventory(model, include_records=True):
    try:
        external = find_external_compartment(model)
        external_error = None
    except (RuntimeError, ValueError) as error:
        external, external_error = None, str(error)
    inferred = {}
    for role in ('exchanges', 'demands', 'sinks'):
        try:
            inferred[role] = sorted(r.id for r in getattr(model, role))
        except (RuntimeError, ValueError):
            inferred[role] = None
    boundary = [r for r in model.reactions if r.boundary]
    external_by_id = [m for m in model.metabolites if m.id.endswith('_e')]
    suffixes = {}
    for suffix in ('_c', '_p', '_e'):
        counts = Counter(m.compartment for m in model.metabolites if m.id.endswith(suffix))
        suffixes[suffix] = dict(sorted(counts.items()))
    shapes = Counter()
    for r in model.reactions:
        values = list(r.metabolites.values())
        if values and (all(c < 0 for c in values) or all(c > 0 for c in values)):
            prefix = 'EX' if r.id.startswith('EX_') else 'DM' if r.id.startswith('DM_') else 'sink' if r.id.lower().startswith(('sink', 'sk_', 'sn_')) else 'other'
            shapes[(prefix, len(values), 'negative' if values[0] < 0 else 'positive')] += 1
    ex_ids = {r.id for r in model.reactions if r.id.startswith('EX_')}
    unrecognized = ex_ids - set(inferred['exchanges'] or [])
    result = {'model_id': model.id, 'n_reactions': len(model.reactions), 'n_metabolites': len(model.metabolites),
        'compartments': model.compartments, 'metabolite_count_by_compartment': dict(sorted(Counter(m.compartment for m in model.metabolites).items())),
        'suffix_assignment_counts': suffixes, 'inferred_external_compartment': external, 'inference_error': external_error,
        'boundary_counts_by_compartment': dict(sorted(Counter(next(iter(r.compartments)) for r in boundary).items())),
        'inferred_role_counts': {key: len(value) if value is not None else None for key, value in inferred.items()},
        'single_sided_shapes': [{'id_prefix_class': p, 'n_metabolites': n, 'coefficient_sign': sign, 'count': count}
                                 for (p, n, sign), count in sorted(shapes.items())],
        'external_suffix_metabolites_outside_inferred_compartment': [
            {'id': m.id, 'compartment': m.compartment, 'formula': m.formula, 'charge': m.charge}
            for m in sorted(external_by_id, key=lambda m: m.id) if m.compartment != external],
        'EX_ids_not_in_cobra_exchanges': [reaction_record(model.reactions.get_by_id(r)) for r in sorted(unrecognized)],
        'nonstandard_exchange_shapes': [reaction_record(r) for r in model.reactions if r.id.startswith('EX_') and
                                          (len(r.metabolites) != 1 or next(iter(r.metabolites.values())) != -1)],
        'internal_or_missing_compartment_EX': [reaction_record(r) for r in model.reactions if r.id.startswith('EX_') and
                                              (external is None or external not in r.compartments)]}
    if include_records:
        result['all_single_metabolite_boundary_records'] = [reaction_record(r) for r in sorted(boundary, key=lambda r: r.id)]
    return result


def completion_requirements(model, media):
    requested = set().union(*(base_medium(name).uptakes for name in media))
    exclude = set(GenericParams().medium_completion_exclude)
    records = []
    for exchange in sorted(requested):
        extracellular, stem = exchange[3:], exchange[3:-2]
        cytoplasm, carrier = stem + '_c', 'MEDt_' + stem
        e = model.metabolites.get_by_id(extracellular) if extracellular in model.metabolites else None
        c = model.metabolites.get_by_id(cytoplasm) if cytoplasm in model.metabolites else None
        conflicts = []
        if c is not None and e is not None:
            for field in ('formula', 'charge'):
                old, new = getattr(e, field), getattr(c, field)
                if old is not None and new is not None and old != new:
                    conflicts.append({'field': field, 'extracellular': old, 'cytoplasmic': new})
        reason = 'already_has_exchange' if exchange in model.reactions else 'excluded_by_declared_policy' if stem in exclude else 'no_cytoplasmic_metabolite' if c is None else 'would_add_exchange_and_inward_carrier'
        records.append({'exchange': exchange, 'decision': reason, 'cytoplasmic_metabolite': cytoplasm,
                        'cytoplasmic_compartment': c.compartment if c is not None else None,
                        'extracellular_metabolite_exists': e is not None, 'extracellular_compartment': e.compartment if e is not None else None,
                        'carrier_id_exists': carrier in model.reactions, 'known_chemical_metadata_conflicts': conflicts})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=BASE / 'medium_structure_inventory.json')
    args = parser.parse_args()
    require(not args.out.exists(), 'Refusing to overwrite inventory')
    logging.getLogger('cobra').setLevel(logging.ERROR)
    def no_optimization(*args, **kwargs):
        raise AssertionError('This inventory must not optimize')
    # CLI guard only: importable loaders remain usable by a separately frozen runner.
    cobra.Model.optimize = no_optimization
    cobra.Model.slim_optimize = no_optimization
    entries, inputs = [], set(CODE)
    for label in LABELS:
        org = 'Putida' if label == 'Putida_iJN1463' else label
        model, provenance = build_parent(label)
        inputs.update(provenance['input_sha256'])
        inputs.update(metadata_paths(org))
        conditions = metadata_conditions(org)
        media = sorted({condition.media for condition in conditions})
        before = inventory(model)
        requirements = completion_requirements(model, media)
        prepared = model.copy()
        for exchange in prepared.exchanges:
            exchange.bounds = (0, 1000)
        added = complete_medium_transport(prepared, media, GenericParams().medium_completion_exclude)
        additions = [reaction_record(r) for r in prepared.reactions if r.id not in model.reactions]
        entries.append({'label': label, 'provenance': provenance, 'before_completion': before,
                        'mapped_condition_count': len(conditions), 'media': media,
                        'conditions': [{'key': c.key, 'name': c.name, 'media': c.media, 'exchanges': c.exchanges} for c in conditions],
                        'completion_requirements': requirements, 'actual_legacy_added_exchange_ids': added,
                        'actual_legacy_added_reactions': additions, 'after_completion': inventory(prepared, False)})
    initial_hashes = hashes(inputs)
    result = {'scope': 'Source structure and in-memory historical preparation only; no optimization or numeric fitness parsing.',
              'cobra_version': cobra.__version__, 'models': entries, 'input_sha256': initial_hashes,
              'Keio': {'status': 'unavailable', 'reason': 'No iML1515 SBML in clone or the three expected original model paths; no new download or surrogate organism used.'},
              'fitness_access': 'Only first-line experiment column names are parsed from fit_logratios.tsv. Exact file bytes are hashed for provenance, without parsing their numerical values.',
              'panel_proposal': {'role': 'predeclared development differential preparation test; no biological validation claim',
                  'real_models': list(LABELS), 'conditions': 'Every mapped condition listed for each parent, with complete-medium union identical to the historical protocol.',
                  'before_any_real_solve': 'Freeze parent source/reconstruction recipes, metadata, old and new helper code, explicit compartment mapping and completion policy.',
                  'structural_gate': 'Compare stoichiometry, GPRs, formula, charge, all bounds and objective exactly for every condition; permit only declared compartment-label and display-name normalization.',
                  'numerical_followup': 'Only after that gate, predeclare paired WT checks with GLPK and HiGHS; no mutant fitness rescore or changes selected from outcomes.',
                  'adversarial_fixtures': ['canonical c/p/e and C_c/C_p/C_e with equivalent stoichiometry',
                    'already existing extracellular metabolite under a second recognized label',
                    'ambiguous external-compartment tie and explicit override',
                    'single-positive source versus negative exchange and malformed multi-metabolite EX',
                    'internal demand/sink boundary preservation',
                    'pre-existing carrier/exchange ID with conflicting stoichiometry or direction',
                    'known formula/charge conflict between an existing external and internal metabolite',
                    'sequential medium reuse, fresh-copy equivalence, idempotence, atomic failure and input-model immutability']},
              'interpretation': 'Different labels do not imply different biology. A suffix or inferred external compartment alone is insufficient to prove semantic identity; records identify what an explicit helper must validate.'}
    require(initial_hashes == hashes(inputs), 'Inputs changed during inventory')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')
    print(json.dumps({'output': str(args.out), 'parents': len(entries),
                      'conditions_per_parent': {r['label']: r['mapped_condition_count'] for r in entries},
                      'legacy_added_exchanges': {r['label']: r['actual_legacy_added_exchange_ids'] for r in entries}}))


if __name__ == '__main__':
    main()
