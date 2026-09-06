"""Frozen development experiment for an exact-source provisional quinone repair.

No old patches are modified. Every invocation requires a new output directory;
failed and partial runs are preserved. Physical-only mode does not load numeric
fitness values. Neither arm selection nor reaction transfer uses fitness scores.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import gzip
import json
import math
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import cobra
from cobra.manipulation.modify import rename_genes
import highspy
import numpy as np
import pandas as pd
import swiglpk

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gembench.cards import BenchmarkCard, ModelProvenance, carbon_fitness_leakage, now, sha256_of, software_versions  # noqa: E402
from gembench.checks import energy_from_nothing  # noqa: E402
from gembench.cofactor import pool_balance_certificate, probe_metabolite_production  # noqa: E402
from gembench.comparison import aligned_pairs, compare_runs, load_run  # noqa: E402
from gembench.fitness_browser import base_medium, carbon_source_conditions, load_organism  # noqa: E402
from gembench.media import apply_medium  # noqa: E402
from gembench.protocols import carbon_fitness_generic as P  # noqa: E402
from gembench.study import freeze_study, verify_study  # noqa: E402
from scripts.map_quinone_pathway import CANDIDATES, met_record, reaction_record  # noqa: E402
from scripts.run_carbon_fitness_generic import gene_map_for, per_condition, score  # noqa: E402
from scripts.run_quinone_biomass_sensitivity import prepare  # noqa: E402

PLAN = ROOT / 'data/studies/quinone_repair_v1.json'
ARM_IDS = ['baseline', 'demand_only', 'curated_path_only', 'curated_path_template',
           'curated_path_curated', 'sequence_path_template']


def read_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f'Duplicate JSON key: {key}')
            result[key] = value
        return result
    def invalid(value):
        raise ValueError(f'Nonfinite JSON value: {value}')
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs, parse_constant=invalid)


def json_safe(value):
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return None if isinstance(value, float) and not math.isfinite(value) else value


def write_json(path, value):
    with Path(path).open('x') as stream:
        json.dump(json_safe(value), stream, indent=2, allow_nan=False)
        stream.write('\n')


def input_path(root, relative):
    path = Path(relative)
    if path.is_absolute() or '..' in path.parts or not path.parts:
        raise ValueError(f'Unsafe input path: {relative}')
    resolved = (root / path).resolve(strict=True)
    if not resolved.is_relative_to(root.resolve()) or not resolved.is_file():
        raise ValueError(f'Input must be a file within the repository: {relative}')
    return resolved


def freeze_inputs(root, recipe_path, recipe, candidate):
    registry_path = root / 'data/studies/exposure_registry_v1.json'
    registry = read_json(registry_path)
    paths = {str(recipe_path.resolve().relative_to(root)), 'scripts/run_quinone_repair.py',
             'scripts/run_quinone_biomass_sensitivity.py', 'scripts/run_carbon_fitness_generic.py',
             'scripts/map_quinone_pathway.py', 'data/studies/exposure_registry_v1.json',
             recipe['candidate_path'], recipe['gene_evidence_path'], 'requirements-audit.txt',
             recipe['transfer_metadata_policy']['basis_path'],
             'models/gapfilled/Putida.xml.gz', 'models/gapfilled/Putida_gapfill.json',
             candidate['source_model']['path'], 'data/genpept/Putida_genpept_map.tsv',
             'data/reference/universe_patches_v0.1.json', 'data/reference/model_patches_v0.4.json'}
    paths.update(candidate['provenance']['input_sha256'])
    paths.update(f'data/reference/gpr_patches_v0.{i}.json' for i in (2, 3, 4))
    paths.update(str(p.relative_to(root)) for p in (root / 'gembench').rglob('*.py'))
    paths.update(str(p.relative_to(root)) for p in (root / 'data/fitness_browser/Putida').glob('*') if p.is_file())
    paths.update(str(p.relative_to(root)) for p in (root / 'data/reference').glob('*.tsv'))
    frozen_plan = {'schema_version': 1, 'study_id': recipe['study_id'], 'evaluation_role': 'development',
                   'development_organisms': registry['development_organisms'],
                   'quarantined_organisms': registry['quarantined_organisms'], 'paths': sorted(paths),
                   'recipe': recipe, 'runtime_versions': {**software_versions(), 'glpk': swiglpk.glp_version(),
                                                        'highs': highspy.Highs().version()},
                   'note': 'Frozen before model preparation or interventions; hashing does not establish independence.'}
    return freeze_study(root, frozen_plan)


def verify_inputs(root, manifest):
    report = verify_study(root, manifest)
    if not report['valid']:
        raise ValueError(f'Frozen inputs changed: {report}')
    return report


def resolve_aliases(model, mapping, loci):
    """Use each existing identity exactly once; never create a duplicate locus."""
    result = {}
    for locus in sorted(set(loci)):
        hits = sorted(gid for gid, value in mapping.model_to_browser.items() if value == locus)
        if len(hits) > 1:
            raise ValueError(f'Ambiguous model gene identity for {locus}: {hits}')
        target = hits[0] if hits else locus
        if target in model.genes and mapping.model_to_browser.get(target) != locus:
            raise ValueError(f'Gene alias collision for {locus}: {target}')
        result[locus] = target
    if len(set(result.values())) != len(result):
        raise ValueError('Multiple source genes resolve to the same target identifier')
    return result


def transfer_candidate(model, source, candidate, mapping, metadata_policy=None):
    """Validate the complete transfer before copying any source reactions."""
    ids = candidate['reaction_ids']
    if ids != CANDIDATES or [r['id'] for r in candidate['reaction_records']] != ids:
        raise ValueError('The candidate must name the five declared source reactions exactly once')
    if any(rid in model.reactions for rid in ids):
        raise ValueError('Candidate reaction already exists; refusing to overwrite')
    reactions = [source.reactions.get_by_id(rid) for rid in ids]
    if any(r.boundary for r in reactions):
        raise ValueError('A biochemical transfer cannot contain an artificial source or sink')
    if [reaction_record(r) for r in reactions] != candidate['reaction_records']:
        raise ValueError('Candidate reaction records differ from the exact source model')
    if any(r.check_mass_balance() or any(not m.formula or m.charge is None for m in r.metabolites) for r in reactions):
        raise ValueError('Source reactions require complete formulas/charges and elemental/charge balance')
    involved = {m.id: m for reaction in reactions for m in reaction.metabolites}
    if {mid: met_record(m) for mid, m in involved.items()} != candidate['all_involved_metabolites']:
        raise ValueError('Candidate metabolite records differ from the exact source model')
    new_mets = {mid: met_record(m) for mid, m in involved.items() if mid not in model.metabolites}
    if new_mets != candidate['new_metabolites']:
        raise ValueError('Expected new metabolites differ from the prepared model')
    policy = metadata_policy or {}
    compartments = policy.get('compartment_aliases', {})
    expected_charges = policy.get('expected_shared_charge_differences', {})
    observed_charges, shared_metadata = {}, []
    for mid, met in involved.items():
        if mid in model.metabolites:
            target = model.metabolites.get_by_id(mid)
            if target.formula != met.formula or target.compartment != compartments.get(met.compartment, met.compartment):
                raise ValueError(f'Metabolite formula/compartment conflict: {mid}')
            if target.charge != met.charge:
                observed_charges[mid] = {'source': met.charge, 'target': target.charge}
            shared_metadata.append({'metabolite': mid, 'source': met_record(met), 'preserved_target': met_record(target)})
    if observed_charges != expected_charges:
        raise ValueError(f'Shared charge differences disagree with the explicit recipe: {observed_charges}')
    aliases = resolve_aliases(model, mapping, {g.id for r in reactions for g in r.genes})
    if aliases != candidate['gene_aliases']:
        raise ValueError('Declared source gene aliases disagree with the prepared model mapping')
    identity = candidate['gene_identity_records']
    if (len(identity) != len(aliases) or {r['source_locus'] for r in identity} != set(aliases) or
            any(r['identity_mapping_ambiguous'] or r['target_model_gene_id'] != aliases[r['source_locus']] for r in identity)):
        raise ValueError('Candidate has incomplete or ambiguous gene identity evidence')
    holder = cobra.Model('exact_source_transfer')
    holder.add_reactions([r.copy() for r in reactions])
    rename_genes(holder, aliases)
    for mid in new_mets:
        met = holder.metabolites.get_by_id(mid)
        met.compartment = compartments.get(met.compartment, met.compartment)
    model.add_reactions([r.copy() for r in holder.reactions])
    return {'source_reaction_ids': ids, 'gene_aliases': aliases, 'new_metabolites': sorted(new_mets),
            'metadata_policy': policy, 'shared_metabolite_metadata': shared_metadata,
            'preserved_charge_differences': observed_charges,
            'new_metabolite_records': {mid: met_record(model.metabolites.get_by_id(mid)) for mid in new_mets},
            'applied_reactions': [reaction_record(model.reactions.get_by_id(rid)) for rid in ids],
            'source_balance': {r.id: r.check_mass_balance() for r in reactions},
            'balance_caveat': 'Source balance is checked separately. Preserved target charge placeholders do not establish target chemical correctness.',
            'caveats': candidate['uncertainties']}


def apply_configuration(base, source, candidate, mapping, config, coefficients, browser_names, metadata_policy=None):
    model = base.copy()
    record = {'configuration': config, 'transfer': None, 'gpr_changes': []}
    if config['add_pathway']:
        record['transfer'] = transfer_candidate(model, source, candidate, mapping, metadata_policy)
    biomass = model.reactions.get_by_id('Growth')
    if any(biomass.metabolites.get(model.metabolites.get_by_id(mid), 0) != 0 for mid in ('q8h2_c', 'mql8_c')):
        raise ValueError('Expected a baseline with no biomass quinone demand')
    coefficient = coefficients[config['biomass_coefficient_source']]
    if coefficient:
        if not math.isfinite(coefficient) or coefficient >= 0:
            raise ValueError('Biomass demand coefficient must be finite and negative')
        biomass.add_metabolites({model.metabolites.get_by_id('q8h2_c'): coefficient})
    record['applied_biomass_coefficient'] = coefficient
    gm = gene_map_for('Putida', model, browser_names, 'gapfilled')
    for rid, source_rule in config['gpr_overrides'].items():
        reaction = model.reactions.get_by_id(rid)
        holder = cobra.Model('gene_rule_override')
        probe = cobra.Reaction('rule')
        probe.gene_reaction_rule = source_rule
        holder.add_reactions([probe])
        loci = {g.id for g in probe.genes}
        if not loci <= browser_names:
            raise ValueError(f'Override genes absent from Putida annotation: {loci - browser_names}')
        aliases = resolve_aliases(model, gm, loci)
        rename_genes(holder, aliases)
        before = reaction.gene_reaction_rule
        reaction.gene_reaction_rule = probe.gene_reaction_rule
        record['gpr_changes'].append({'reaction': rid, 'before': before, 'after': reaction.gene_reaction_rule,
                                      'source_locus_rule': source_rule, 'gene_aliases': aliases})
        gm = gene_map_for('Putida', model, browser_names, 'gapfilled')
    model.id += '__' + config['id']
    return model, gm, record


def lp_signature(model):
    objective = cobra.util.solver.linear_reaction_coefficients(model)
    return {'direction': model.objective_direction,
            'reactions': {r.id: (r.bounds, {m.id: c for m, c in r.metabolites.items()}, objective.get(r, 0)) for r in model.reactions},
            'constraints': {c.name: (c.lb, c.ub, str(c.expression)) for c in model.constraints}}


def verify_sequence_equivalence(reference, variant, reference_map, variant_map, fitness_rows):
    if lp_signature(reference) != lp_signature(variant):
        raise ValueError('Sequence arm differs in stoichiometry, bounds, objective, or solver constraints')
    changes = []
    for before in reference.reactions:
        after = variant.reactions.get_by_id(before.id)
        if before.gene_reaction_rule != after.gene_reaction_rule:
            loci = set()
            for reaction, mapping in ((before, reference_map), (after, variant_map)):
                if any(g.id not in mapping.model_to_browser for g in reaction.genes):
                    raise ValueError('An altered GPR has an unmapped gene; benchmark omission cannot be justified')
                loci.update(mapping.model_to_browser[g.id] for g in reaction.genes)
            if loci & fitness_rows:
                raise ValueError(f'Sequence GPR changes involve benchmark-covered genes: {sorted(loci & fitness_rows)}')
            changes.append({'reaction': before.id, 'loci': sorted(loci)})
    return {'stoichiometry_bounds_objective_constraints_identical': True,
            'changed_gprs_use_exclusively_unscored_loci': True, 'changed_reactions': changes,
            'benchmark_omitted_for_sequence_arm': True,
            'interpretation': 'Identical LP for wild type and every covered single-gene deletion; targeted unscored predictions remain distinct hypotheses.'}


def metadata_organism(root):
    folder = root / 'data/fitness_browser/Putida'
    genes = pd.read_table(folder / 'genes.tsv', dtype=str, keep_default_na=False)
    experiments = pd.read_table(folder / 'experiments.tsv', dtype=str, keep_default_na=False)
    meta_cols = ['orgId', 'locusId', 'sysName', 'geneName', 'desc']
    metadata = pd.read_table(folder / 'fit_logratios.tsv', usecols=meta_cols, dtype=str, keep_default_na=False)
    rows = metadata.sysName.where(metadata.sysName != '', metadata.locusId)
    if not rows.is_unique:
        raise ValueError('Ambiguous duplicate fitness gene identifiers')
    columns = pd.read_table(folder / 'fit_logratios.tsv', nrows=0).columns
    fitness = SimpleNamespace(index=rows, columns=[c.split(' ')[0] for c in columns if c not in meta_cols])
    return SimpleNamespace(org_id='Putida', genes=genes, experiments=experiments, fitness=fitness)


def growth_solution(model):
    solution = model.optimize()
    if solution.status != 'optimal' or solution.objective_value is None or not math.isfinite(solution.objective_value):
        raise RuntimeError(f'Growth solve failed: {solution.status}')
    if not np.isfinite(solution.fluxes.to_numpy()).all():
        raise RuntimeError('Growth solve returned nonfinite fluxes')
    return {'status': solution.status, 'growth': float(solution.objective_value)}


def targeted_knockouts(model, mapping, loci, fitness_rows, *, wild_type=None, growth_threshold=1e-3):
    knockouts = []
    wild_type = growth_solution(model) if wild_type is None else wild_type
    wt_grows = wild_type['growth'] >= growth_threshold
    aliases = resolve_aliases(model, mapping, loci)
    for locus in loci:
        gid = aliases[locus]
        record = {'locus': locus, 'model_gene_id': gid if gid in model.genes else None,
                  'exported_fitness_row_present': locus in fitness_rows, 'role': 'model prediction only; no experimental validation',
                  'wild_type_grows_at_threshold': wt_grows, 'growth_threshold': growth_threshold,
                  'wild_type_growth': wild_type['growth'],
                  'interpretation': 'Conditional model prediction only' if wt_grows else 'Uninterpretable for gene essentiality because wild type does not grow'}
        if gid not in model.genes:
            record.update(status='unrepresented_gene', growth=None, associated_reactions=[])
        else:
            gene = model.genes.get_by_id(gid)
            record['associated_reactions'] = sorted(r.id for r in gene.reactions)
            if not gene.reactions:
                record.update(status='unrepresented_function', growth=None)
            else:
                with model:
                    gene.knock_out()
                    record.update(growth_solution(model))
        knockouts.append(record)
    return knockouts


def physical_checks(model, mapping, conditions, recipe, fitness_rows):
    model = model.copy()
    params = P.GenericParams(**recipe['protocol'])
    settings = recipe['physical_settings']
    model.solver = params.solver
    model.solver.configuration.tolerances.feasibility = settings['feasibility_tolerance']
    for reaction in model.exchanges:
        reaction.bounds = (0, 1000)
    for gid in params.knockout_genes:
        model.genes.get_by_id(gid).knock_out()
    added = P.complete_medium_transport(model, sorted({c.media for c in conditions if c.bigg_ids}), params.medium_completion_exclude) if params.complete_medium_transport else []
    egc = energy_from_nothing(model)
    if any(value > settings['egc_threshold'] for value in egc.values()):
        raise ValueError(f'Energy from nothing after intervention: {egc}')
    selected = recipe['physical_condition']
    condition = next(c for c in conditions if c.name == selected['name'] and c.media == selected['media'])
    missing = apply_medium(model, base_medium(condition.media), close_all=True)
    for exchange in condition.exchanges:
        model.reactions.get_by_id(exchange).lower_bound = params.carbon_uptake
    certificate = pool_balance_certificate(model, recipe['pool_weights'])
    wt = growth_solution(model)
    with model:
        model.reactions.Growth.lower_bound = 0
        production = probe_metabolite_production(model, list(recipe['pool_weights']), capacity=settings['production_capacity'], threshold=settings['production_threshold'])
    knockouts = targeted_knockouts(model, mapping, recipe['targeted_gene_loci'], fitness_rows,
                                  wild_type=wt, growth_threshold=params.growth_threshold)
    return {'condition': condition.key, 'settings': settings, 'solver': model.solver.interface.__name__,
            'effective_feasibility_tolerance': model.solver.configuration.tolerances.feasibility,
            'objective': {'expression': str(model.objective.expression), 'direction': model.objective_direction},
            'growth_bounds': list(model.reactions.Growth.bounds),
            'atp_maintenance_bounds': list(model.reactions.ATPM.bounds) if 'ATPM' in model.reactions else None,
            'maintenance_note': ('No positive ATP-maintenance lower bound is imposed.'
                                 if 'ATPM' in model.reactions and model.reactions.ATPM.lower_bound == 0
                                 else 'ATPM bounds are reported as stored; no maintenance requirement was added.'),
            'mandatory_flux_bounds': {r.id: list(r.bounds) for r in model.reactions if r.lower_bound > 0 or r.upper_bound < 0},
            'final_exchange_bounds': {r.id: list(r.bounds) for r in sorted(model.exchanges, key=lambda r: r.id)},
            'declared_carbon_source_ids': condition.bigg_ids, 'declared_carbon_exchanges': condition.exchanges,
            'declared_carbon_uptake': params.carbon_uptake,
            'medium_completion_added': added, 'missing_medium_components': missing,
            'energy_from_nothing': egc, 'tested_energy_currencies': sorted(egc),
            'pool_certificate': certificate, 'wild_type': wt, 'production_at_growth_lower_bound_zero': production,
            'targeted_gene_predictions': knockouts}


def save_benchmark(directory, model, mapping, intervention, fb, conditions, recipe, manifest):
    params = P.GenericParams(**recipe['protocol'])
    result = P.run(model, fb, conditions, mapping, params)
    grows = np.isfinite(result.wt_growth) & (result.wt_growth >= params.growth_threshold)
    results = {'condition_level': {'n_conditions_mapped': len(result.conditions), 'n_conditions_wt_grows': int(grows.sum())},
               'gene_level_conditions_where_wt_grows': score(result, grows), 'counts': result.counts}
    card = BenchmarkCard(benchmark='carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Putida', created=now(),
        dataset_provenance=fb.provenance,
        model=ModelProvenance(model_id=model.id, file=str(directory / 'model.xml.gz'),
                             source='Frozen v0.4 prepared model plus explicitly declared provisional intervention',
                             sha256=sha256_of(directory / 'model.xml.gz'), n_reactions=len(model.reactions),
                             n_metabolites=len(model.metabolites), n_genes=len(model.genes)),
        protocol={'params': asdict(params), 'evaluation_role': 'development', 'intervention': intervention,
                  'study_fingerprint': manifest['content_fingerprint']},
        leakage=carbon_fitness_leakage('Putida', 'gapfilled', patched=True, medium_completion=True), results=results,
        warnings=recipe['caveats'])
    write_json(directory / 'card.json', asdict(card))
    with (directory / 'card.md').open('x') as stream:
        stream.write(card.to_markdown())
    with (directory / 'matrices.npz').open('xb') as stream:
        np.savez_compressed(stream, sim_growth=result.sim_growth, wt_growth=result.wt_growth, fitness=result.fitness,
                            model_genes=np.array(result.model_genes), browser_genes=np.array(result.browser_genes),
                            conditions=np.array([c.key for c in result.conditions]))
    result.condition_table().to_csv(directory / 'conditions.tsv', sep='\t', index=False, mode='x')
    per_condition(result).to_csv(directory / 'per_condition_metrics.tsv', sep='\t', index=False, mode='x')
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, default=PLAN)
    parser.add_argument('--out', type=Path, default=ROOT / 'results/quinone_repair_2026_09_06/runs/main')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--freeze-only', action='store_true')
    mode.add_argument('--physical-only', action='store_true')
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(f'Output exists; preserve this attempt and choose a new --out directory: {args.out}')
    args.out.mkdir(parents=True)
    started = time.time()
    try:
        recipe = read_json(args.plan)
        if recipe['evaluation_role'] != 'development' or recipe['organism'] != 'Putida' or recipe['base_model_patch_version'] != '0.4':
            raise ValueError('This runner supports the declared Putida development v0.4 baseline only')
        if [c['id'] for c in recipe['configurations']] != ARM_IDS:
            raise ValueError('Unexpected arm inventory or order')
        candidate = read_json(input_path(ROOT, recipe['candidate_path']))
        manifest = freeze_inputs(ROOT, args.plan, recipe, candidate)
        write_json(args.out / 'manifest.json', manifest)
        if args.freeze_only:
            return
        verify_inputs(ROOT, manifest)
        for path, digest in candidate['provenance']['input_sha256'].items():
            if sha256_of(input_path(ROOT, path)) != digest:
                raise ValueError(f'Candidate evidence input changed: {path}')
        source_path = input_path(ROOT, candidate['source_model']['path'])
        if sha256_of(source_path) != candidate['source_model']['sha256']:
            raise ValueError('Source model hash differs from candidate')
        source = cobra.io.read_sbml_model(str(source_path))
        if source.id != candidate['source_model']['model_id']:
            raise ValueError('Source model identity differs from candidate')
        with gzip.open(ROOT / 'models/gapfilled/Putida.xml.gz', 'rt') as stream:
            template = cobra.io.read_sbml_model(stream)
        coefficients = {'none': 0.0, 'template': template.reactions.Growth.get_coefficient('mql8_c'),
                        'curated': source.reactions.BIOMASS_KT2440_WT3.get_coefficient('q8h2_c')}
        metadata = metadata_organism(ROOT)
        conditions = carbon_source_conditions(metadata)
        browser_names, fitness_rows = set(metadata.genes.sysName), set(metadata.fitness.index)
        base, mapping, _, _ = prepare(metadata)
        arms = {c['id']: apply_configuration(base, source, candidate, mapping, c, coefficients, browser_names, recipe['transfer_metadata_policy'])
                for c in recipe['configurations']}
        a, am, _ = arms['curated_path_template']
        b, bm, _ = arms['sequence_path_template']
        equivalence = verify_sequence_equivalence(a, b, am, bm, fitness_rows)
        if recipe['configurations'][-1]['run_benchmark'] or not all(c['run_benchmark'] for c in recipe['configurations'][:-1]):
            raise ValueError('Expected five benchmark arms with only the verified sequence arm omitted')
        write_json(args.out / 'sequence_benchmark_equivalence.json', equivalence)
        physical = []
        for config in recipe['configurations']:
            verify_inputs(ROOT, manifest)
            model, gm, intervention = arms[config['id']]
            directory = args.out / config['id']
            directory.mkdir()
            write_json(directory / 'intervention.json', intervention)
            write_json(directory / 'gene_map.json', asdict(gm))
            cobra.io.write_sbml_model(model, str(directory / 'model.xml.gz'))
            report = physical_checks(model, gm, conditions, recipe, fitness_rows)
            write_json(directory / 'physical.json', report)
            physical.append({'arm': config['id'], 'wild_type': report['wild_type'], 'production': report['production_at_growth_lower_bound_zero']})
            print('Physical checks completed:', config['id'], report['wild_type'], flush=True)
        write_json(args.out / 'physical_summary.json', {'arms': physical, 'coefficients': coefficients, 'role': 'physical model predictions, not biological validation'})
        if args.physical_only:
            write_json(args.out / 'input_verification_after_run.json', verify_inputs(ROOT, manifest))
            return
        fb = load_organism('Putida')
        records = []
        for config in recipe['configurations']:
            if not config['run_benchmark']:
                continue
            verify_inputs(ROOT, manifest)
            model, gm, intervention = arms[config['id']]
            record = save_benchmark(args.out / config['id'], model, gm, intervention, fb, conditions, recipe, manifest)
            records.append({'arm': config['id'], **record})
        baseline = load_run(args.out / 'baseline')
        comparisons, changes = [], []
        for record in records[1:]:
            other = load_run(args.out / record['arm'])
            comparisons.append({'arm': record['arm'], **compare_runs(baseline, other, n_boot=recipe['comparison']['resamples'], seed=recipe['comparison']['seed'])})
            sa, sb, fit, genes, common, _ = aligned_pairs(baseline, other)
            threshold = recipe['protocol']['growth_threshold']
            different = np.isfinite(fit) & ((sa >= threshold) != (sb >= threshold))
            for i, j in zip(*np.where(different)):
                changes.append({'arm': record['arm'], 'gene': genes[i], 'condition': common[j],
                                'before': sa[i, j], 'after': sb[i, j], 'fitness': fit[i, j]})
        write_json(args.out / 'changed_predictions.json', changes)
        write_json(args.out / 'input_verification_after_run.json', verify_inputs(ROOT, manifest))
        write_json(args.out / 'summary.json', {'study': recipe, 'study_fingerprint': manifest['content_fingerprint'],
                   'runs': records, 'comparisons': comparisons, 'physical_arms': physical,
                   'n_changed_gene_condition_predictions': len(changes), 'seconds': time.time() - started})
    except Exception as error:
        write_json(args.out / 'failure.json', {'created': now(), 'type': type(error).__name__, 'message': str(error),
                   'seconds': time.time() - started, 'note': 'Partial outputs preserved; use a new directory for a reviewed retry.'})
        raise


if __name__ == '__main__':
    main()
