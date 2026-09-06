"""Predeclared development test of a homology-supported PpnP hypothesis.

Start from the saved quinone-demand model; never amend historical patches.
Source reversibility is overridden in the primary arms because the primary
experimental supplement contradicts itself about reverse catalysis.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path
import sys
import time

import cobra
import highspy
import numpy as np
import swiglpk

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gembench.cards import BenchmarkCard, ModelProvenance, carbon_fitness_leakage, now, sha256_of, software_versions
from gembench.cofactor import pool_balance_certificate, probe_metabolite_production
from gembench.comparison import aligned_pairs, compare_runs, load_run
from gembench.fitness_browser import carbon_source_conditions, load_organism
from gembench.protocols import carbon_fitness_generic as P
from gembench.study import freeze_study
from scripts.map_quinone_pathway import met_record, reaction_record
from scripts.run_carbon_fitness_generic import gene_map_for, per_condition, score
from scripts.run_quinone_repair import (input_path, metadata_organism, physical_checks,
    read_json, resolve_aliases, verify_inputs, write_json)
from scripts.verify_quinone_repair_solvers import reconstruct

PLAN = ROOT / 'data/studies/ppnp_repair_v1.json'
ARM_IDS = ['parent', 'adenosine_forward', 'inosine_forward', 'both_forward',
           'parent_rhcys_closed', 'both_forward_rhcys_closed', 'both_reversible']
BENCHMARK_IDS = ['parent', 'both_forward', 'both_forward_rhcys_closed', 'both_reversible']
REACTION_IDS = {'adenosine': 'PUNP1', 'inosine': 'PUNP5'}
RIBOSYL_POOL = {mid: 1.0 for mid in ('ahcys_c', 'rhcys_c', 'adn_c', 'ins_c', 'rib__D_c')}


def freeze_inputs(root, recipe_path, recipe):
    registry = read_json(root / 'data/studies/exposure_registry_v1.json')
    old = read_json(input_path(root, recipe['parent_manifest']))
    paths = {f['path'] for f in old['files']}
    candidate = read_json(input_path(root, recipe['candidate_path']))
    paths.update(candidate['provenance']['input_sha256'])
    paths.update(str(p.relative_to(root)) for p in (root / 'gembench').rglob('*.py'))
    paths.update(['scripts/run_ppnp_repair.py', 'scripts/verify_quinone_repair_solvers.py',
                  'scripts/diagnose_quinone_producibility.py', 'tests/test_ppnp_repair.py'])
    paths.update(str(p.relative_to(root)) for p in (root / recipe['evidence_directory']).rglob('*') if p.is_file())
    paths.update([str(recipe_path.resolve().relative_to(root)), recipe['base_model'],
                  recipe['parent_manifest'], recipe['sequence_evidence'], recipe['study_definition'],
                  'requirements-audit.txt', 'data/studies/exposure_registry_v1.json'])
    return freeze_study(root, {
        'schema_version': 1, 'study_id': recipe['study_id'], 'evaluation_role': 'development',
        'development_organisms': registry['development_organisms'],
        'quarantined_organisms': registry['quarantined_organisms'], 'paths': sorted(paths),
        'recipe': recipe, 'runtime_versions': {**software_versions(), 'glpk': swiglpk.glp_version(),
                                             'highs': highspy.Highs().version()},
        'note': 'Frozen before real model interventions/solves. Hypothesis arose after previous development outcomes; this is not preregistration or independent validation.'})


def validate_candidate(base, source, candidate):
    """Validate all chemistry/metadata before any mutation of the parent."""
    records = candidate['reaction_records']
    if [r['id'] for r in records] != ['PUNP1', 'PUNP5']:
        raise ValueError('Require exactly the two canonical source reactions')
    if [reaction_record(source.reactions.get_by_id(r['id'])) for r in records] != records:
        raise ValueError('Source reaction records changed')
    involved = {m.id: m for r in records for m in source.reactions.get_by_id(r['id']).metabolites}
    if {mid: met_record(m) for mid, m in involved.items()} != candidate['all_involved_metabolites']:
        raise ValueError('Source metabolite metadata changed')
    expected = candidate['expected_target_metabolites']
    if {mid: met_record(base.metabolites.get_by_id(mid)) for mid in expected} != expected:
        raise ValueError('Parent metabolite metadata changed')
    proposals = candidate['proposed_reactions']
    if [p['id'] for p in proposals] != ['PUNP1', 'PUNP5']:
        raise ValueError('Unexpected proposed reaction inventory')
    if candidate['target_metabolite_map'] != {mid: mid for mid in involved} or set(expected) != set(involved):
        raise ValueError('This experiment requires identical existing cytosolic metabolite identities')
    for proposal in proposals:
        r = source.reactions.get_by_id(proposal['source_id'])
        if r.boundary or r.check_mass_balance() or any(not m.formula or m.charge is None for m in r.metabolites):
            raise ValueError('Require complete, balanced source biochemistry')
        if proposal['source_id'] != proposal['id'] or proposal['id'] != REACTION_IDS[proposal['substrate']]:
            raise ValueError('Substrate/source identity mismatch')
        if proposal['metabolites'] != {m.id: c for m, c in r.metabolites.items()}:
            raise ValueError('Proposed stoichiometry differs from source')
        if r.id in base.reactions:
            raise ValueError('Candidate reaction already exists')
        for m in r.metabolites:
            target = base.metabolites.get_by_id(m.id)
            allowed = (m.charge, target.charge) == ((-2, 0) if m.id in {'pi_c', 'r1p_c'} else (m.charge, m.charge))
            if not allowed or m.formula != target.formula or m.compartment != target.compartment:
                raise ValueError(f'Undeclared metabolite conflict: {m.id}')
        for existing in base.reactions:
            row = {m.id: c for m, c in existing.metabolites.items()}
            if row == proposal['metabolites'] or row == {k: -v for k, v in proposal['metabolites'].items()}:
                raise ValueError(f'Equivalent chemistry already exists: {existing.id}')


def apply_arm(base, source, candidate, mapping, configuration, browser_names):
    validate_candidate(base, source, candidate)
    if configuration['direction'] not in {'forward', 'reversible'}:
        raise ValueError('Unknown direction policy')
    if len(set(configuration['substrates'])) != len(configuration['substrates']) or set(configuration['substrates']) - set(REACTION_IDS):
        raise ValueError('Unexpected/duplicate substrate inventory')
    aliases = resolve_aliases(base, mapping, ['PP_4248'])
    model = base.copy()
    before_biomass = reaction_record(model.reactions.Growth)
    added = []
    for substrate in configuration['substrates']:
        rid = REACTION_IDS[substrate]
        original = source.reactions.get_by_id(rid)
        r = cobra.Reaction(rid, name=original.name, subsystem=original.subsystem,
                           lower_bound=0 if configuration['direction'] == 'forward' else -1000,
                           upper_bound=1000)
        r.add_metabolites({model.metabolites.get_by_id(m.id): c for m, c in original.metabolites.items()})
        r.gene_reaction_rule = aliases['PP_4248']
        r.annotation = dict(original.annotation)
        model.add_reactions([r])
        if r.check_mass_balance():
            raise ValueError('Added reaction fails target metadata balance')
        added.append(reaction_record(r))
    closed = []
    if configuration['close_rhcys']:
        closed.append({'reaction': 'RHCYS', 'before': list(model.reactions.RHCYS.bounds), 'after': [0, 0]})
        model.reactions.RHCYS.bounds = (0, 0)
    if reaction_record(model.reactions.Growth) != before_biomass:
        raise ValueError('Biomass demand must remain exactly unchanged')
    model.id = f'{base.id}__ppnp_{configuration["id"]}'
    gm = gene_map_for('Putida', model, browser_names, 'gapfilled')
    if added and gm.model_to_browser.get(aliases['PP_4248']) != 'PP_4248':
        raise ValueError('Candidate gene failed unique organism mapping')
    return model, gm, {'configuration': configuration, 'source_reactions': [r['id'] for r in added],
                      'added_reaction_records': added, 'closed_reactions': closed, 'gene_aliases': aliases,
                      'evidence_grade': 'KT2440 function by homology; E. coli substrate assays',
                      'target_charge_note': 'Preserved historical zero placeholders; source Pi and r1p charges are -2.',
                      'biomass_unchanged': True}


def witness(model, recipe):
    solution = model.optimize()
    if solution.status != 'optimal' or solution.objective_value is None or not np.isfinite(solution.objective_value):
        raise RuntimeError(f'Failed/nonfinite mechanistic solve: {solution.status}')
    flux = solution.fluxes
    if not np.isfinite(flux.to_numpy()).all():
        raise RuntimeError('Nonfinite witness')
    residual = max(abs(sum(c * flux[r.id] for r in m.reactions for met, c in r.metabolites.items() if met is m)) for m in model.metabolites)
    bound_error = max(max(r.lower_bound - flux[r.id], flux[r.id] - r.upper_bound, 0) for r in model.reactions)
    if max(residual, bound_error) > recipe['physical_settings']['residual_tolerance']:
        raise RuntimeError('Witness violates stoichiometry/bounds')
    observed = {rid: float(flux[rid]) for rid in recipe['witness_reactions'] if rid in model.reactions}
    return {'status': 'optimal', 'growth': float(solution.objective_value),
            'mass_balance_max_abs': float(residual), 'bound_violation_max': float(bound_error),
            'selected_fluxes': observed, 'witness_note': 'One optimal solution; selected fluxes are not unique or measured.'}


def mechanism_checks(model, mapping, physical, conditions, recipe):
    model = model.copy()
    reconstruct(model, physical, recipe, conditions)
    aliases = resolve_aliases(model, mapping, recipe['targeted_gene_loci'])
    cases = []
    for specification in recipe['mechanistic_cases']:
        record = {'id': specification['id'], 'loci': specification['loci'],
                  'close_added_reactions': specification.get('close_added_reactions', False)}
        absent = [l for l in specification['loci'] if aliases[l] not in model.genes or not model.genes.get_by_id(aliases[l]).reactions]
        if absent:
            record.update(status='unrepresented_gene_or_function', growth=None, unrepresented=absent)
        else:
            with model:
                for locus in specification['loci']:
                    model.genes.get_by_id(aliases[locus]).knock_out()
                if record['close_added_reactions']:
                    for rid in REACTION_IDS.values():
                        if rid in model.reactions:
                            model.reactions.get_by_id(rid).bounds = (0, 0)
                record.update(witness(model, recipe))
                if specification.get('probe_quinone'):
                    record['quinone_production'] = probe_metabolite_production(model, list(recipe['pool_weights']),
                        capacity=recipe['physical_settings']['production_capacity'],
                        threshold=recipe['physical_settings']['production_threshold'])
        cases.append(record)
    return {'cases': cases, 'ribosyl_pool_certificate': pool_balance_certificate(model, RIBOSYL_POOL),
            'role': 'Diagnostic model predictions only. Missing genes are not interpreted as dispensable.',
            'ppm_note': 'PPM deletion need not block rescue: reverse PNP + RNMK + NMNN can convert r1p to r5p.'}


def save_benchmark(directory, model, mapping, intervention, fb, conditions, recipe, manifest):
    params = P.GenericParams(**recipe['protocol'])
    result = P.run(model, fb, conditions, mapping, params)
    grows = np.isfinite(result.wt_growth) & (result.wt_growth >= params.growth_threshold)
    results = {'condition_level': {'n_conditions_mapped': len(result.conditions), 'n_conditions_wt_grows': int(grows.sum())},
               'gene_level_conditions_where_wt_grows': score(result, grows), 'counts': result.counts}
    card = BenchmarkCard(benchmark='carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Putida', created=now(),
        dataset_provenance=fb.provenance,
        model=ModelProvenance(model_id=model.id, file=str(directory / 'model.xml.gz'),
            source='Exact saved provisional quinone-demand parent plus declared PpnP/direction/RHCYS intervention',
            sha256=sha256_of(directory / 'model.xml.gz'), n_reactions=len(model.reactions), n_metabolites=len(model.metabolites), n_genes=len(model.genes)),
        protocol={'params': asdict(params), 'evaluation_role': 'development', 'intervention': intervention,
                  'study_fingerprint': manifest['content_fingerprint']},
        leakage=carbon_fitness_leakage('Putida', 'gapfilled', patched=True, medium_completion=True),
        results=results, warnings=recipe['caveats'])
    write_json(directory / 'card.json', asdict(card))
    with (directory / 'card.md').open('x') as stream:
        stream.write('\n'.join(line.rstrip() for line in card.to_markdown().splitlines()) + '\n')
    with (directory / 'matrices.npz').open('xb') as stream:
        np.savez_compressed(stream, sim_growth=result.sim_growth, wt_growth=result.wt_growth, fitness=result.fitness,
            model_genes=np.array(result.model_genes), browser_genes=np.array(result.browser_genes), conditions=np.array([c.key for c in result.conditions]))
    result.condition_table().to_csv(directory / 'conditions.tsv', sep='\t', index=False, mode='x')
    per_condition(result).to_csv(directory / 'per_condition_metrics.tsv', sep='\t', index=False, mode='x')
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, default=PLAN)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--freeze-only', action='store_true')
    parser.add_argument('--physical-only', action='store_true')
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('Preserve every attempt; select a fresh output directory')
    args.out.mkdir(parents=True)
    started = time.time()
    try:
        recipe = read_json(args.plan)
        if recipe['evaluation_role'] != 'development' or recipe['organism'] != 'Putida':
            raise ValueError('Only the declared Putida development study is supported')
        configs = recipe['configurations']
        if [c['id'] for c in configs] != ARM_IDS or [c['id'] for c in configs if c['run_benchmark']] != BENCHMARK_IDS:
            raise ValueError('Unexpected study arms/benchmark inventory')
        verify_inputs(ROOT, read_json(input_path(ROOT, recipe['parent_manifest'])))
        candidate = read_json(input_path(ROOT, recipe['candidate_path']))
        for path, digest in candidate['provenance']['input_sha256'].items():
            if sha256_of(input_path(ROOT, path)) != digest:
                raise ValueError(f'Candidate evidence input changed: {path}')
        manifest = freeze_inputs(ROOT, args.plan, recipe)
        write_json(args.out / 'manifest.json', manifest)
        if args.freeze_only:
            return
        verify_inputs(ROOT, manifest)
        base_path = input_path(ROOT, recipe['base_model'])
        source_path = input_path(ROOT, candidate['source_model']['path'])
        if sha256_of(base_path) != recipe['base_model_sha256'] or sha256_of(source_path) != candidate['source_model']['sha256']:
            raise ValueError('Source/parent model hash changed')
        base = cobra.io.read_sbml_model(str(base_path))
        source = cobra.io.read_sbml_model(str(source_path))
        if source.id != candidate['source_model']['model_id']:
            raise ValueError('Source model identity changed')
        metadata = metadata_organism(ROOT)
        browser_names, fitness_rows = set(metadata.genes.sysName), set(metadata.fitness.index)
        conditions = carbon_source_conditions(metadata)
        mapping = gene_map_for('Putida', base, browser_names, 'gapfilled')
        arms = {}
        for config in configs:
            verify_inputs(ROOT, manifest)
            model, gm, intervention = apply_arm(base, source, candidate, mapping, config, browser_names)
            arms[config['id']] = (model, gm, intervention)
            directory = args.out / config['id']
            directory.mkdir()
            write_json(directory / 'intervention.json', intervention)
            write_json(directory / 'gene_map.json', asdict(gm))
            cobra.io.write_sbml_model(model, str(directory / 'model.xml.gz'))
            physical = physical_checks(model, gm, conditions, recipe, fitness_rows)
            write_json(directory / 'physical.json', physical)
            write_json(directory / 'mechanism.json', mechanism_checks(model, gm, physical, conditions, recipe))
            print('Physical checks completed:', config['id'], physical['wild_type'], flush=True)
        if args.physical_only:
            write_json(args.out / 'input_verification_after_run.json', verify_inputs(ROOT, manifest))
            return
        fb = load_organism('Putida')
        records = []
        for config in configs:
            if config['run_benchmark']:
                verify_inputs(ROOT, manifest)
                model, gm, intervention = arms[config['id']]
                record = save_benchmark(args.out / config['id'], model, gm, intervention, fb, conditions, recipe, manifest)
                records.append({'arm': config['id'], **record})
                print('Benchmark completed:', config['id'], flush=True)
        parent = load_run(args.out / 'parent')
        comparisons, changes = [], []
        for record in records[1:]:
            other = load_run(args.out / record['arm'])
            comparisons.append({'arm': record['arm'], **compare_runs(parent, other,
                n_boot=recipe['comparison']['resamples'], seed=recipe['comparison']['seed'])})
            sa, sb, fit, genes, common, _ = aligned_pairs(parent, other)
            threshold = recipe['protocol']['growth_threshold']
            different = np.isfinite(fit) & ((sa >= threshold) != (sb >= threshold))
            for i, j in zip(*np.where(different)):
                changes.append({'arm': record['arm'], 'gene': genes[i], 'condition': common[j],
                                'before': sa[i, j], 'after': sb[i, j], 'fitness': fit[i, j]})
        write_json(args.out / 'changed_predictions.json', changes)
        write_json(args.out / 'input_verification_after_run.json', verify_inputs(ROOT, manifest))
        write_json(args.out / 'summary.json', {'study': recipe, 'study_fingerprint': manifest['content_fingerprint'],
            'runs': records, 'comparisons': comparisons, 'n_changed_gene_condition_predictions': len(changes),
            'seconds': time.time() - started})
    except Exception as error:
        write_json(args.out / 'failure.json', {'created': now(), 'type': type(error).__name__, 'message': str(error),
            'seconds': time.time() - started, 'note': 'Partial outputs preserved. Refreeze reviewed changes in a fresh output directory.'})
        raise


if __name__ == '__main__':
    main()
