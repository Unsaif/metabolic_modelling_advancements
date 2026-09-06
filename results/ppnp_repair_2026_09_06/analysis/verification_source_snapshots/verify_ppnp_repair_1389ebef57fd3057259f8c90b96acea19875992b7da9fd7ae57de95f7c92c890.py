"""Verify saved PpnP models and benchmark artifacts without running simulations.

Only the earlier independent verifier supplies metric/reading helpers. Neither
the experiment runner nor production scoring/comparison functions are imported.
Raw Putida fitness is re-averaged; native and paired gene-bootstrap statistics
are recomputed. A fresh JSON report preserves failure and input/source hashes.
Physical numerical re-solves belong to the separate solver verification.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
from pathlib import Path
import platform
import re
import sys
import time

import cobra
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import verify_quinone_repair as V  # noqa: E402

ARM_IDS = ['parent', 'adenosine_forward', 'inosine_forward', 'both_forward',
           'parent_rhcys_closed', 'both_forward_rhcys_closed', 'both_reversible']
BENCHMARK_IDS = ['parent', 'both_forward', 'both_forward_rhcys_closed', 'both_reversible']
REACTIONS = {'adenosine': 'PUNP1', 'inosine': 'PUNP5'}
STOICHIOMETRY = {
    'PUNP1': {'adn_c': -1.0, 'pi_c': -1.0, 'ade_c': 1.0, 'r1p_c': 1.0},
    'PUNP5': {'ins_c': -1.0, 'pi_c': -1.0, 'hxan_c': 1.0, 'r1p_c': 1.0},
}
ARM_FILES = ('model.xml.gz', 'gene_map.json', 'intervention.json', 'physical.json', 'mechanism.json')
BENCHMARK_FILES = ('card.json', 'card.md', 'matrices.npz', 'conditions.tsv', 'per_condition_metrics.tsv')


def exact(actual, expected, label):
    """Model coefficients/bounds are exact; metric comparisons use V.equal."""
    V.require(actual == expected, f'{label}: exact declared model data differ')


def local_file(root, relative):
    path = Path(relative)
    V.require(not path.is_absolute() and '..' not in path.parts, f'Unsafe repository path: {relative}')
    resolved = (root / path).resolve()
    V.require(resolved.is_relative_to(root) and resolved.is_file(), f'Missing/unsafe repository input: {relative}')
    return resolved


def reaction_equal(actual, expected, label):
    a, b = V.reaction_record(actual), V.reaction_record(expected)
    for record in (a, b):
        record.pop('equation')
        record.pop('mass_charge_imbalance')  # deterministic stoichiometry checked directly
        record['gene_reaction_rule'] = V.canonical_rule(record['gene_reaction_rule'])
    exact(a, b, label)
    exact(actual.subsystem, expected.subsystem, label + ' subsystem')


def mapping_from_annotations(model, table, browser_names):
    """Resolve identities from raw accession/locus metadata, not saved maps."""
    by_version, by_accession, old = {}, {}, {}
    for row in table.itertuples(index=False):
        if row.version in by_version:
            exact(by_version[row.version], row.locus_tags, 'Conflicting versioned gene identity')
        by_version[row.version] = row.locus_tags
        accession = row.version.split('.')[0]
        if accession in by_accession:
            exact(by_accession[accession], row.locus_tags, 'Conflicting unversioned gene identity')
        by_accession[accession] = row.locus_tags
        old[row.version] = row.old_locus_tags
    result = {}
    for gene in model.genes:
        gid = gene.id
        if gid in browser_names:
            result[gid] = gid
            continue
        accession = gid.removeprefix('G_')
        if '.' not in accession:
            accession = re.sub(r'_(\d+)$', r'.\1', accession)
        tags = by_version.get(accession, by_accession.get(accession.split('.')[0]))
        if tags is None:
            continue
        hits = set()
        for tag in (tags + ';' + old.get(accession, '')).split(';'):
            if tag:
                hits.update({tag, tag.replace('_', ''), tag.replace('_', '', 1)} & browser_names)
        V.require(len(hits) <= 1, f'Ambiguous Browser identity for {gid}: {sorted(hits)}')
        if hits:
            result[gid] = hits.pop()
    return result


def verify_candidate(root, recipe, candidate):
    base_path = local_file(root, recipe['base_model'])
    V.equal(V.sha(base_path), recipe['base_model_sha256'], 'Exact parent model hash')
    V.equal(candidate['baseline_model']['sha256'], recipe['base_model_sha256'], 'Candidate parent hash')
    V.equal(candidate['baseline_model']['path'], recipe['base_model'], 'Candidate parent path')
    source_path = local_file(root, candidate['source_model']['path'])
    V.equal(V.sha(source_path), candidate['source_model']['sha256'], 'Candidate source hash')
    for relative, digest in candidate['provenance']['input_sha256'].items():
        V.equal(V.sha(local_file(root, relative)), digest, 'Candidate provenance ' + relative)
    base = cobra.io.read_sbml_model(str(base_path))
    source = cobra.io.read_sbml_model(str(source_path))
    exact(base.id, candidate['baseline_model']['model_id'], 'Candidate parent identity')
    exact(source.id, candidate['source_model']['model_id'], 'Candidate source identity')
    exact([r['id'] for r in candidate['reaction_records']], list(STOICHIOMETRY), 'Source reaction inventory')
    exact([r['id'] for r in candidate['proposed_reactions']], list(STOICHIOMETRY), 'Proposed reaction inventory')
    for record, proposal in zip(candidate['reaction_records'], candidate['proposed_reactions']):
        rid = record['id']
        actual = V.reaction_record(source.reactions.get_by_id(rid))
        V.compare_reaction(actual, record, 'Source ' + rid)
        exact(actual['metabolites'], STOICHIOMETRY[rid], 'Canonical phosphorolysis ' + rid)
        exact(proposal['metabolites'], STOICHIOMETRY[rid], 'Proposed phosphorolysis ' + rid)
        exact(proposal['source_id'], rid, 'Proposed source ID')
        exact(REACTIONS[proposal['substrate']], rid, 'Proposed substrate ID')
        V.require(not actual['mass_charge_imbalance'] and not actual['missing_formula_or_charge'], 'Incomplete/unbalanced source chemistry')
        V.require(rid not in base.reactions, 'Candidate already present in parent')
    mids = set().union(*(set(row) for row in STOICHIOMETRY.values()))
    exact({mid: V.met_record(source.metabolites.get_by_id(mid)) for mid in mids}, candidate['all_involved_metabolites'], 'Source metabolite metadata')
    exact({mid: V.met_record(base.metabolites.get_by_id(mid)) for mid in mids}, candidate['expected_target_metabolites'], 'Parent metabolite metadata')
    exact(candidate['target_metabolite_map'], {mid: mid for mid in mids}, 'Existing metabolite identities')
    differences = {}
    for mid in mids:
        s, b = source.metabolites.get_by_id(mid), base.metabolites.get_by_id(mid)
        exact((s.formula, s.compartment), (b.formula, b.compartment), 'Source/parent chemistry ' + mid)
        if s.charge != b.charge:
            differences[mid] = {'source': s.charge, 'target': b.charge}
    exact(differences, {'pi_c': {'source': -2, 'target': 0}, 'r1p_c': {'source': -2, 'target': 0}}, 'Declared historical charge placeholders')
    exact(candidate['metadata_policy']['expected_charge_differences'], differences, 'Candidate charge policy')
    return base, source


def verify_declared_arm(base, model, source, config, alias):
    """Check the loaded model against a parent copy with only declared edits."""
    V.require(config['direction'] in ('forward', 'reversible'), 'Unsupported direction policy')
    V.require(len(set(config['substrates'])) == len(config['substrates']) and set(config['substrates']) <= set(REACTIONS), 'Invalid substrate inventory')
    added = [REACTIONS[s] for s in config['substrates']]
    exact({r.id for r in model.reactions}, {r.id for r in base.reactions} | set(added), 'Reaction inventory ' + config['id'])
    exact({m.id for m in model.metabolites}, {m.id for m in base.metabolites}, 'Unchanged metabolite inventory')
    exact({g.id for g in model.genes}, {g.id for g in base.genes} | ({alias} if added else set()), 'Gene inventory')
    exact(model.id, base.id + '__ppnp_' + config['id'], 'Model arm identity')
    exact((model.name, model.annotation), (base.name, base.annotation), 'Preserved model metadata')
    exact(model.objective_direction, base.objective_direction, 'Objective direction')
    exact({r.id: c for r, c in cobra.util.solver.linear_reaction_coefficients(model).items()},
          {r.id: c for r, c in cobra.util.solver.linear_reaction_coefficients(base).items()}, 'Objective coefficients')
    for metabolite in base.metabolites:
        exact(V.met_record(model.metabolites.get_by_id(metabolite.id)), V.met_record(metabolite), 'Preserved metabolite ' + metabolite.id)
    for gene in base.genes:
        actual = model.genes.get_by_id(gene.id)
        exact((actual.name, actual.annotation), (gene.name, gene.annotation), 'Preserved gene metadata ' + gene.id)
    for reaction in base.reactions:
        expected = reaction.copy()
        if reaction.id == 'RHCYS' and config['close_rhcys']:
            expected.bounds = (0, 0)
        reaction_equal(model.reactions.get_by_id(reaction.id), expected, 'Preserved reaction ' + reaction.id)
    for rid in added:
        original = source.reactions.get_by_id(rid)
        expected = original.copy()
        expected.bounds = (0 if config['direction'] == 'forward' else -1000, 1000)
        expected.gene_reaction_rule = alias
        expected.annotation = dict(original.annotation)
        reaction_equal(model.reactions.get_by_id(rid), expected, 'Added reaction ' + rid)
    return added


def verify_models(root, directory, recipe, candidate, genes):
    base, source = verify_candidate(root, recipe, candidate)
    table = pd.read_table(root / 'data/genpept/Putida_genpept_map.tsv', header=None,
        names=['version', 'locus_tags', 'old_locus_tags', 'gene_names', 'coded_by', 'definition'], dtype=str, keep_default_na=False)
    browser_names = set(genes.sysName)
    base_map = mapping_from_annotations(base, table, browser_names)
    hits = [gid for gid, locus in base_map.items() if locus == 'PP_4248']
    V.require(len(hits) <= 1, 'Ambiguous parent PpnP gene identity')
    alias = hits[0] if hits else 'PP_4248'
    exact(candidate['gene_aliases'], {'PP_4248': alias}, 'Candidate gene alias')
    exact(candidate['proposed_gene_rule'], 'PP_4248', 'Candidate single-locus hypothesis')
    reports = {}
    for config in recipe['configurations']:
        arm = config['id']
        model = cobra.io.read_sbml_model(str(directory / arm / 'model.xml.gz'))
        added = verify_declared_arm(base, model, source, config, alias)
        intervention = V.read_json(directory / arm / 'intervention.json')
        exact(intervention['configuration'], config, arm + ' saved configuration')
        exact(intervention['source_reactions'], added, arm + ' declared additions')
        exact(intervention['gene_aliases'], {'PP_4248': alias}, arm + ' declared aliases')
        V.require(intervention['biomass_unchanged'] is True, arm + ' missing biomass invariant')
        records = intervention['added_reaction_records']
        exact([r['id'] for r in records], added, arm + ' added record inventory')
        for record in records:
            V.compare_reaction(V.reaction_record(model.reactions.get_by_id(record['id'])), record, arm + ' added reaction record')
        closure = [{'reaction': 'RHCYS', 'before': list(base.reactions.RHCYS.bounds), 'after': [0, 0]}] if config['close_rhcys'] else []
        exact(intervention['closed_reactions'], closure, arm + ' RHCYS closure record')
        mapping = V.read_json(directory / arm / 'gene_map.json')
        exact(mapping['org_id'], 'Putida', arm + ' mapping organism')
        expected_map = mapping_from_annotations(model, table, browser_names)
        exact(mapping['model_to_browser'], expected_map, arm + ' independently resolved gene map')
        exact(set(mapping['unmapped_model_genes']), {g.id for g in model.genes} - set(expected_map), arm + ' unmapped genes')
        V.require(len(mapping['unmapped_model_genes']) == len(set(mapping['unmapped_model_genes'])), arm + ' duplicate unmapped gene')
        physical = V.read_json(directory / arm / 'physical.json')
        exact(physical['settings'], recipe['physical_settings'], arm + ' physical settings')
        exact(physical['condition'], ' | '.join(recipe['physical_condition'][k] for k in ('name', 'media')), arm + ' physical condition')
        mechanism = V.read_json(directory / arm / 'mechanism.json')
        exact([r['id'] for r in mechanism['cases']], [r['id'] for r in recipe['mechanistic_cases']], arm + ' mechanism case inventory')
        reports[arm] = {'model': model, 'mapping': mapping, 'intervention': intervention, 'physical': physical}
    return reports


def verify_completion(report, recipe, conditions, media):
    model = report['model']
    exclusion = recipe['protocol'].get('medium_completion_exclude', ['pnto__R', 'fol', 'hco3'])
    components = {c for name in {r['media'] for r in conditions.values() if r['bigg_ids']}
                  for c in V.medium_components(media[name])}
    added = sorted('EX_' + c + '_e' for c in components if 'EX_' + c + '_e' not in model.reactions
                   and c + '_c' in model.metabolites and c not in exclusion)
    if not recipe['protocol']['complete_medium_transport']:
        added = []
    exact(report['physical']['medium_completion_added'], added, model.id + ' completion inventory')


def coverage_change(parent, other):
    new = sorted(set(other['browser_genes']) - set(parent['browser_genes']))
    lost = sorted(set(parent['browser_genes']) - set(other['browser_genes']))
    rows = []
    for locus in new:
        index = list(other['browser_genes']).index(locus)
        finite = np.isfinite(other['fitness'][index]) & np.isfinite(other['sim_growth'][index])
        grows = other['wt_growth'] >= other['st']
        rows.append({'locus': locus, 'n_finite_pairs': int(finite.sum()),
                     'n_finite_pairs_where_wt_grows': int((finite & grows).sum()),
                     'excluded_from_common_gene_comparison': True})
    return {'new_browser_genes': new, 'lost_browser_genes': lost, 'new_gene_coverage': rows}


def verify(root, directory):
    V.require(not (directory / 'failure.json').exists(), 'Main run has a saved failure marker')
    manifest = V.verify_manifest(root, directory)
    recipe = manifest['plan']['recipe']
    exact([c['id'] for c in recipe['configurations']], ARM_IDS, 'Seven-arm declaration')
    exact([c['id'] for c in recipe['configurations'] if c['run_benchmark']], BENCHMARK_IDS, 'Four scored arms')
    exact((recipe['organism'], recipe['evaluation_role']), ('Putida', 'development'), 'Study scope')
    actual_arms = [p.name for p in directory.iterdir() if p.is_dir() and any((p / f).exists() for f in ARM_FILES + BENCHMARK_FILES)]
    exact(sorted(actual_arms), sorted(ARM_IDS), 'Saved arm directory inventory')
    fitness, experiments, genes = V.read_fitness(root)
    conditions, media = V.condition_inventory(root, fitness, experiments)
    candidate = V.read_json(local_file(root, recipe['candidate_path']))
    reports = verify_models(root, directory, recipe, candidate, genes)
    for report in reports.values():
        verify_completion(report, recipe, conditions, media)
    summary = V.read_json(directory / 'summary.json')
    exact(summary['study'], recipe, 'Summary recipe')
    exact(summary['study_fingerprint'], manifest['content_fingerprint'], 'Summary fingerprint')
    exact([r['arm'] for r in summary['runs']], BENCHMARK_IDS, 'Summary benchmark inventory')
    runs, native = {}, {}
    for config in recipe['configurations']:
        arm = config['id']
        if not config['run_benchmark']:
            V.require(not any((directory / arm / f).exists() for f in BENCHMARK_FILES), arm + ' contains an undeclared benchmark')
            continue
        run, recomputed = V.verify_benchmark(directory / arm, reports[arm], recipe, manifest, fitness, experiments, conditions)
        exact(run['card']['benchmark'], 'carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Putida', arm + ' benchmark identity')
        record = next(r for r in summary['runs'] if r['arm'] == arm)
        V.equal({k: v for k, v in record.items() if k != 'arm'}, run['card']['results'], arm + ' summary/card results')
        runs[arm], native[arm] = run, recomputed
        print('Verified model, raw replicate means, coverage and native statistics:', arm, flush=True)
    comparisons, changes, coverage = {}, [], {}
    exact([r['arm'] for r in summary['comparisons']], BENCHMARK_IDS[1:], 'Paired comparison inventory')
    for record in summary['comparisons']:
        arm = record['arm']
        computed = V.paired_comparison(runs['parent'], runs[arm], resamples=recipe['comparison']['resamples'], seed=recipe['comparison']['seed'])
        for key, value in computed.items():
            V.equal(record[key], value, arm + ' paired statistics ' + key)
        exact(record['evaluation_role'], 'retrospective_development', arm + ' comparison role')
        comparisons[arm] = computed
        changes.extend(V.changed_predictions(runs['parent'], runs[arm], arm))
        coverage[arm] = coverage_change(runs['parent'], runs[arm])
        expected_new = ['PP_4248'] if 'PP_4248' in fitness.index and 'PP_4248' not in runs['parent']['browser_genes'] else []
        exact(coverage[arm]['new_browser_genes'], expected_new, arm + ' expected new gene coverage')
        exact(coverage[arm]['lost_browser_genes'], [], arm + ' lost gene coverage')
    V.equal(V.read_json(directory / 'changed_predictions.json'), changes, 'Every changed prediction')
    exact(summary['n_changed_gene_condition_predictions'], len(changes), 'Changed prediction count')
    return {'study_fingerprint': manifest['content_fingerprint'], 'verified_model_arms': ARM_IDS,
            'verified_benchmark_arms': BENCHMARK_IDS, 'native_statistics': native,
            'paired_comparisons': comparisons, 'coverage_changes': coverage,
            'fitness_coverage': {'annotated_genes': len(genes), 'exported_fitness_rows': len(fitness)},
            'changed_prediction_count': len(changes),
            'model_inventories': {a: {'reactions': len(r['model'].reactions), 'metabolites': len(r['model'].metabolites), 'genes': len(r['model'].genes)} for a, r in reports.items()},
            'scope': 'Exact declared parent/model changes; raw-data, metric and artifact reproducibility on exposed development data. No optimizer called; physical optima and physiological validity are not established.'}


def input_records(root, directory):
    paths = {directory / f for f in ('manifest.json', 'summary.json', 'changed_predictions.json', 'input_verification_after_run.json')}
    paths.update(directory / arm / f for arm in ARM_IDS for f in ARM_FILES)
    paths.update(directory / arm / f for arm in BENCHMARK_IDS for f in BENCHMARK_FILES)
    manifest = V.read_json(directory / 'manifest.json')
    paths.update(local_file(root, r['path']) for r in manifest['files'])
    candidate = V.read_json(local_file(root, manifest['plan']['recipe']['candidate_path']))
    paths.update(local_file(root, p) for p in candidate['provenance']['input_sha256'])
    paths.update((Path(__file__).resolve(), Path(V.__file__).resolve()))
    return [{'path': str(p), 'size_bytes': p.stat().st_size, 'sha256': V.sha(p)} for p in sorted(paths)]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True, help='Fresh JSON report; an existing report is never overwritten')
    args = parser.parse_args(argv)
    directory, output = args.run_dir.resolve(), args.out.resolve()
    V.require(not output.exists(), 'Report already exists; preserve it and use a fresh --out')
    started = time.time()
    report = {'created_utc': datetime.now(timezone.utc).isoformat(), 'passed': False,
              'run_directory': str(directory), 'source_sha256': {str(p): V.sha(p) for p in (Path(__file__).resolve(), Path(V.__file__).resolve())},
              'runtime': {'python': platform.python_version(), 'platform': platform.platform(),
                          **{n: importlib.metadata.version(n) for n in ('numpy', 'pandas', 'cobra', 'python-libsbml', 'optlang')}},
              'model_path_policy': 'Card checkout prefixes are historical; arm/filename suffix and actual exported model hash/identity remain mandatory.'}
    try:
        records = input_records(ROOT, directory)
        report['inputs'] = records
        report.update(verify(ROOT, directory))
        exact(input_records(ROOT, directory), records, 'Inputs changed during verification')
        report['passed'] = True
    except Exception as error:
        report['failure'] = {'type': type(error).__name__, 'message': str(error)}
    report['seconds'] = time.time() - started
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x') as stream:
        json.dump(V.safe(report), stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'passed': report['passed'], 'report': str(output), 'failure': report.get('failure')}))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
