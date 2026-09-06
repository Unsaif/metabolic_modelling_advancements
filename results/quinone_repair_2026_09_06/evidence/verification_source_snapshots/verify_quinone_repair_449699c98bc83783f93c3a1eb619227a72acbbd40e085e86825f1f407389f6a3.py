"""Independently verify a completed six-arm quinone-repair development run.

Reads saved artifacts and frozen inputs; never runs an optimizer or imports the
runner, production scoring, comparison, or cofactor helpers. Confusion, MCC,
ranking metrics, gene bootstrap, alignment, and pool arithmetic are implemented
here. Physical fluxes are checked for internal consistency, not independently
re-solved. A new report is written exclusively; any disagreement exits nonzero.
"""
from __future__ import annotations

import argparse
import ast
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import re
import sys
import time

import cobra
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ARM_IDS = ['baseline', 'demand_only', 'curated_path_only', 'curated_path_template',
           'curated_path_curated', 'sequence_path_template']
REACTION_IDS = ['OHPHM', 'OMPHHX', 'OMBZLM', 'OMMBLHX', 'DMQMT']
METRICS = ['aucpr_bernstein', 'aucpr_standard', 'auroc_standard', 'mcc', 'balanced_accuracy', 'accuracy']


class VerificationError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise VerificationError(message)


def read_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, f'Duplicate JSON key in {path}: {key}')
            result[key] = value
        return result
    def invalid(value):
        raise VerificationError(f'Nonfinite JSON literal in {path}: {value}')
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs, parse_constant=invalid)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def safe(value):
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, dict):
        return {str(k): safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe(v) for v in value]
    return None if isinstance(value, float) and not math.isfinite(value) else value


def equal(actual, expected, label, *, atol=1e-10):
    actual, expected = safe(actual), safe(expected)
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and actual.keys() == expected.keys(), f'{label}: dictionary fields differ')
        for key in expected:
            equal(actual[key], expected[key], f'{label}.{key}', atol=atol)
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(actual) == len(expected), f'{label}: list lengths differ')
        for index, (a, e) in enumerate(zip(actual, expected)):
            equal(a, e, f'{label}[{index}]', atol=atol)
    elif isinstance(expected, (float, int)) and not isinstance(expected, bool):
        require(isinstance(actual, (float, int)) and not isinstance(actual, bool) and
                math.isclose(actual, expected, rel_tol=1e-12, abs_tol=atol), f'{label}: {actual!r} != {expected!r}')
    else:
        require(actual == expected, f'{label}: {actual!r} != {expected!r}')


def confusion(sim, fit, st=0.001, ft=-2.0):
    sim, fit = np.asarray(sim, dtype=float), np.asarray(fit, dtype=float)
    require(sim.shape == fit.shape, 'Confusion inputs have different shapes')
    valid = np.isfinite(sim) & np.isfinite(fit)
    predicted, observed = sim[valid] >= st, fit[valid] >= ft
    return {'tp': int(np.count_nonzero(predicted & observed)),
            'tn': int(np.count_nonzero(~predicted & ~observed)),
            'fp': int(np.count_nonzero(predicted & ~observed)),
            'fn': int(np.count_nonzero(~predicted & observed))}


def mcc_from_counts(counts):
    tp, tn, fp, fn = (counts[k] for k in ('tp', 'tn', 'fp', 'fn'))
    if tp + tn + fp + fn == 0:
        return float('nan')
    denominator = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    # The benchmark uses the conventional zero for a nonempty degenerate
    # confusion table. No observations is a distinct, undefined quantity.
    return (tp * tn - fp * fn) / math.sqrt(denominator) if denominator else 0.0


def rank_metrics(labels, scores):
    """Tie-grouped PR trapezoid area, average precision, and ROC area."""
    if not len(labels) or np.all(labels == labels[0]):
        return (float('nan'),) * 3
    order = np.argsort(-scores, kind='stable')
    ordered_scores, positives = scores[order], labels[order].astype(int)
    endpoints = np.r_[np.flatnonzero(ordered_scores[:-1] != ordered_scores[1:]), len(scores) - 1]
    tp = np.cumsum(positives)[endpoints].astype(float)
    selected = endpoints + 1
    precision = np.r_[1.0, tp / selected]
    recall = np.r_[0.0, tp / positives.sum()]
    fpr = np.r_[0.0, (selected - tp) / (len(labels) - positives.sum())]
    return (float(np.sum(np.diff(recall) * (precision[1:] + precision[:-1]) / 2)),
            float(np.sum(np.diff(recall) * precision[1:])),
            float(np.sum(np.diff(fpr) * (recall[1:] + recall[:-1]) / 2)))


def metric_points(sim, fit, st, ft):
    counts = confusion(sim, fit, st, ft)
    n = sum(counts.values())
    pg, ng = counts['tp'] + counts['fn'], counts['tn'] + counts['fp']
    finite = np.isfinite(sim) & np.isfinite(fit)
    s, f = sim[finite], fit[finite]
    bernstein = rank_metrics(s < st, -f)[0]
    _, ap, roc = rank_metrics(f < ft, -s)
    return {'aucpr_bernstein': bernstein, 'aucpr_standard': ap, 'auroc_standard': roc,
            'mcc': mcc_from_counts(counts),
            'balanced_accuracy': (counts['tp'] / pg + counts['tn'] / ng) / 2 if pg and ng else float('nan'),
            'accuracy': (counts['tp'] + counts['tn']) / n if n else float('nan')}


def interval(values):
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    return np.percentile(finite, [2.5, 97.5]).tolist() if len(finite) else [float('nan')] * 2


def native_score(run, *, resamples=500, seed=0):
    st, ft = run['st'], run['ft']
    grows = np.isfinite(run['wt_growth']) & (run['wt_growth'] >= st)
    sim, fit = run['sim_growth'][:, grows], run['fitness'][:, grows]
    rows = np.isfinite(fit).any(axis=1)
    sim, fit = sim[rows], fit[rows]
    result = {'n_genes': len(sim), 'n_conditions': int(grows.sum()),
              'n_gene_condition_pairs': int((np.isfinite(sim) & np.isfinite(fit)).sum()),
              'n_missing_fitness': int((~np.isfinite(fit)).sum()),
              'n_nonfinite_simulation': int((~np.isfinite(sim)).sum())}
    if not sim.size:
        return result
    points = metric_points(sim, fit, st, ft)
    draws = {name: [] for name in METRICS}
    rng = np.random.default_rng(seed)
    for _ in range(resamples):
        indices = rng.integers(0, len(sim), len(sim))
        sampled = metric_points(sim[indices], fit[indices], st, ft)
        for name, value in sampled.items():
            draws[name].append(value)
    result.update({name: {'point': points[name], 'ci95': interval(draws[name])} for name in METRICS})
    result['confusion'] = confusion(sim, fit, st, ft)
    return result


def align(a, b):
    require(a['st'] == b['st'] and a['ft'] == b['ft'], 'Comparison thresholds differ')
    genes = sorted(set(a['browser_genes']) & set(b['browser_genes']))
    common = sorted(set(a['conditions']) & set(b['conditions']))
    gi = [{v: i for i, v in enumerate(r['browser_genes'])} for r in (a, b)]
    ci = [{v: i for i, v in enumerate(r['conditions'])} for r in (a, b)]
    grows = [condition for condition in common if all(np.isfinite(r['wt_growth'][ix[condition]]) and
             r['wt_growth'][ix[condition]] >= a['st'] for r, ix in zip((a, b), ci))]
    arrays = []
    for r, g, c in zip((a, b), gi, ci):
        indices = np.ix_([g[x] for x in genes], [c[x] for x in grows])
        arrays.append((r['sim_growth'][indices], r['fitness'][indices]))
    (sa, fa), (sb, fb) = arrays
    require(np.array_equal(fa, fb, equal_nan=True), 'Aligned experimental observations disagree')
    valid = np.isfinite(sa) & np.isfinite(sb) & np.isfinite(fa)
    return sa, sb, np.where(valid, fa, np.nan), genes, grows, common


def paired_comparison(a, b, *, resamples=500, seed=0):
    sa, sb, fit, genes, grows, common = align(a, b)
    counts = {label: confusion(s, fit, a['st'], a['ft']) for label, s in [('A', sa), ('B', sb)]}
    scores = {label: mcc_from_counts(c) for label, c in counts.items()}
    rng, differences = np.random.default_rng(seed), []
    for _ in range(resamples if genes else 0):
        selected = rng.integers(0, len(genes), len(genes))
        ma = mcc_from_counts(confusion(sa[selected], fit[selected], a['st'], a['ft']))
        mb = mcc_from_counts(confusion(sb[selected], fit[selected], a['st'], a['ft']))
        differences.append(mb - ma)
    return {'n_common_genes': len(genes), 'n_common_conditions': len(common),
            'n_conditions_both_grow': len(grows), 'conditions_both_grow': grows,
            'n_shared_finite_pairs': int(np.isfinite(fit).sum()),
            'wt_grows': {label: int((np.isfinite(r['wt_growth']) & (r['wt_growth'] >= a['st'])).sum())
                         for label, r in [('A', a), ('B', b)]},
            'mcc': scores, 'confusion': counts,
            'paired_mcc_difference_B_minus_A': {'point': scores['B'] - scores['A'], 'ci95': interval(differences),
                'finite_bootstrap_samples': int(np.isfinite(differences).sum()), 'requested_bootstrap_samples': resamples}}


def changed_predictions(a, b, arm):
    sa, sb, fit, genes, grows, _ = align(a, b)
    changed = np.isfinite(sa) & np.isfinite(sb) & np.isfinite(fit) & ((sa >= a['st']) != (sb >= a['st']))
    return [{'arm': arm, 'gene': genes[i], 'condition': grows[j], 'before': float(sa[i, j]),
             'after': float(sb[i, j]), 'fitness': float(fit[i, j])} for i, j in zip(*np.where(changed))]


def met_record(met):
    return {'id': met.id, 'name': met.name, 'formula': met.formula, 'charge': met.charge,
            'compartment': met.compartment, 'annotation': met.annotation}


def reaction_record(reaction):
    # The parser and SBML reader are shared third-party infrastructure; the
    # experiment's extraction and transfer helpers are deliberately not used.
    balance = {}
    for met, coefficient in reaction.metabolites.items():
        for element, amount in (met.elements or {}).items():
            balance[element] = balance.get(element, 0.0) + coefficient * amount
        if met.charge is not None:
            balance['charge'] = balance.get('charge', 0.0) + coefficient * met.charge
    return {'id': reaction.id, 'name': reaction.name, 'equation': reaction.reaction,
            'metabolites': {m.id: c for m, c in reaction.metabolites.items()},
            'lower_bound': reaction.lower_bound, 'upper_bound': reaction.upper_bound,
            'gene_reaction_rule': reaction.gene_reaction_rule, 'gene_ids': sorted(g.id for g in reaction.genes),
            'annotation': reaction.annotation, 'mass_charge_imbalance': {k: v for k, v in balance.items() if abs(v) > 1e-12},
            'missing_formula_or_charge': sorted(m.id for m in reaction.metabolites if not m.formula or m.charge is None)}


def renamed_rule(rule, aliases):
    return re.sub(r'\b[A-Za-z_][A-Za-z_0-9]*\b', lambda match: aliases.get(match.group(), match.group()), rule)


def canonical_rule(rule):
    if not rule:
        return ()
    def walk(node):
        if isinstance(node, ast.Name):
            return ('gene', node.id)
        if isinstance(node, ast.BoolOp) and isinstance(node.op, (ast.And, ast.Or)):
            return (type(node.op).__name__, tuple(sorted(walk(x) for x in node.values)))
        raise VerificationError('Unsupported gene-rule syntax')
    return walk(cobra.core.gene.GPR.from_string(rule).body)


def compare_reaction(actual, expected, label):
    actual, expected = dict(actual), dict(expected)
    # Equation text and ordering of a Boolean OR are serialization details;
    # stoichiometry and the complete Boolean expression remain mandatory.
    actual.pop('equation', None)
    expected.pop('equation', None)
    equal(canonical_rule(actual.pop('gene_reaction_rule')), canonical_rule(expected.pop('gene_reaction_rule')), label + '.gpr')
    equal(actual, expected, label)


def pool_certificate(model, weights):
    rows = []
    for reaction in sorted(model.reactions, key=lambda r: r.id):
        terms = {m.id: float(c) for m, c in reaction.metabolites.items() if m.id in weights and c}
        if not terms:
            continue
        total = sum((Fraction.from_float(float(weights[mid])) * Fraction.from_float(c) for mid, c in terms.items()), Fraction(0))
        rows.append({'reaction': reaction.id, 'pool_terms': terms, 'net_pool_coefficient': float(total),
                     'exact_numerator': str(total.numerator), 'exact_denominator': str(total.denominator),
                     'bounds': list(reaction.bounds), 'boundary': len(reaction.metabolites) == 1,
                     'can_supply_by_direction': bool(total > 0 and reaction.upper_bound > 0 or total < 0 and reaction.lower_bound < 0),
                     'can_consume_by_direction': bool(total < 0 and reaction.upper_bound > 0 or total > 0 and reaction.lower_bound < 0)})
    return {'pool_weights': weights, 'exactly_conserved_in_stoichiometry': all(r['exact_numerator'] == '0' for r in rows),
            'n_incident_reactions': len(rows), 'reactions': rows,
            'direction_compatible_supply_reactions': [r['reaction'] for r in rows if r['can_supply_by_direction']]}


def verify_manifest(root, directory):
    manifest = read_json(directory / 'manifest.json')
    files = manifest['files']
    require(len({r['path'] for r in files}) == len(files), 'Duplicate frozen input paths')
    equal(sorted(r['path'] for r in files), sorted(manifest['plan']['paths']), 'Manifest input inventory')
    payload = {key: manifest[key] for key in ('schema_version', 'manifest_type', 'plan')}
    payload['files'] = sorted(files, key=lambda record: record['path'])
    fingerprint = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False).encode()).hexdigest()
    equal(manifest['content_fingerprint'], fingerprint, 'Manifest fingerprint')
    for record in files:
        path = (root / record['path']).resolve()
        require(path.is_relative_to(root) and path.is_file(), f'Unsafe/missing frozen input: {record["path"]}')
        equal(path.stat().st_size, record['size_bytes'], f'{record["path"]} byte count')
        equal(sha(path), record['sha256'], f'{record["path"]} frozen hash')
    after = read_json(directory / 'input_verification_after_run.json')
    require(after['valid'] and after['checked_files'] == len(files), 'Runner did not verify every frozen input after completion')
    equal(after['content_fingerprint'], fingerprint, 'Post-run fingerprint')
    return manifest


def read_fitness(root):
    directory = root / 'data/fitness_browser/Putida'
    raw = pd.read_table(directory / 'fit_logratios.tsv', keep_default_na=False)
    metadata_columns = {'orgId', 'locusId', 'sysName', 'geneName', 'desc'}
    identifiers = raw.sysName.where(raw.sysName != '', raw.locusId).astype(str)
    require(identifiers.is_unique, 'Duplicate exported fitness gene identifiers')
    values = raw.drop(columns=list(metadata_columns)).apply(pd.to_numeric, errors='raise')
    values.columns = [column.split(' ')[0] for column in values.columns]
    values.index = identifiers
    require(values.columns.is_unique, 'Duplicate exported fitness experiment identifiers')
    experiments = pd.read_table(directory / 'experiments.tsv', dtype=str, keep_default_na=False)
    genes = pd.read_table(directory / 'genes.tsv', dtype=str, keep_default_na=False)
    return values, experiments, genes


def condition_inventory(root, fitness, experiments):
    carbon = pd.read_table(root / 'data/reference/fitness_browser_carbon_sources_bigg.tsv', dtype=str, keep_default_na=False)
    media = pd.read_table(root / 'data/reference/fitness_browser_media_bigg.tsv', dtype=str, keep_default_na=False)
    require(carbon.condition.is_unique and media.media.is_unique, 'Ambiguous condition/medium mapping')
    carbon = carbon.set_index('condition').to_dict('index')
    media = media.set_index('media').to_dict('index')
    conditions = {}
    for row in experiments.itertuples(index=False):
        if (row.expGroup != 'carbon source' or row.condition_2 not in ('', 'Dimethyl Sulfoxide') or
                row.media not in media or row.expName not in fitness.columns):
            continue
        key = row.condition_1 + ' | ' + row.media
        if key not in conditions:
            mapping = carbon.get(row.condition_1, {'bigg_ids': '', 'confidence': 'none'})
            conditions[key] = {'name': row.condition_1, 'media': row.media,
                               'bigg_ids': mapping['bigg_ids'], 'confidence': mapping['confidence'], 'experiments': []}
        conditions[key]['experiments'].append(row.expName)
    return conditions, media


def medium_components(record):
    components = [value for value in record['bigg_ids'].split(';') if value]
    if record['aerobic'].strip().lower() in ('yes', 'true', '1') and 'o2' not in components:
        components.append('o2')
    return components


def verify_models(root, directory, recipe, candidate, fitness_rows):
    source_path = root / candidate['source_model']['path']
    equal(sha(source_path), candidate['source_model']['sha256'], 'Candidate source-model hash')
    source = cobra.io.read_sbml_model(str(source_path))
    equal(source.id, candidate['source_model']['model_id'], 'Candidate source-model identity')
    equal(candidate['reaction_ids'], REACTION_IDS, 'Candidate reaction inventory')
    equal([record['id'] for record in candidate['reaction_records']], REACTION_IDS, 'Candidate reaction record inventory')
    for record in candidate['reaction_records']:
        compare_reaction(reaction_record(source.reactions.get_by_id(record['id'])), record, 'Candidate source ' + record['id'])
        require(not record['mass_charge_imbalance'] and not record['missing_formula_or_charge'], 'Candidate source chemistry is incomplete or unbalanced')
    source_involved = {met.id: met_record(met) for rid in REACTION_IDS for met in source.reactions.get_by_id(rid).metabolites}
    equal(source_involved,
          candidate['all_involved_metabolites'], 'Candidate source metabolite metadata')
    for path, digest in candidate['provenance']['input_sha256'].items():
        equal(sha(root / path), digest, 'Candidate evidence hash ' + path)
    models = {arm: cobra.io.read_sbml_model(str(directory / arm / 'model.xml.gz')) for arm in ARM_IDS}
    base = models['baseline']
    base_map = read_json(directory / 'baseline/gene_map.json')['model_to_browser']
    template = cobra.io.read_sbml_model(str(root / 'models/gapfilled/Putida.xml.gz'))
    coefficients = {'none': 0.0, 'template': template.reactions.Growth.get_coefficient('mql8_c'),
                    'curated': source.reactions.BIOMASS_KT2440_WT3.get_coefficient('q8h2_c')}
    for met in ('q8h2_c', 'mql8_c'):
        equal(base.reactions.Growth.metabolites.get(base.metabolites.get_by_id(met), 0.0), 0.0, 'Baseline quinone demand ' + met)
    new_mets = candidate['new_metabolites']
    equal(sorted(new_mets), sorted(set(candidate['all_involved_metabolites']) - {m.id for m in base.metabolites}), 'Candidate new-metabolite inventory')
    observed_charge_differences = {}
    for mid in set(candidate['all_involved_metabolites']) - set(new_mets):
        original = source.metabolites.get_by_id(mid)
        target = base.metabolites.get_by_id(mid)
        equal(target.formula, original.formula, 'Shared source/target formula ' + mid)
        equal(target.compartment, recipe['transfer_metadata_policy']['compartment_aliases'].get(original.compartment, original.compartment), 'Declared shared compartment alias ' + mid)
        if target.charge != original.charge:
            observed_charge_differences[mid] = {'source': original.charge, 'target': target.charge}
    equal(observed_charge_differences, recipe['transfer_metadata_policy']['expected_shared_charge_differences'], 'Actual shared charge metadata differences')
    alias = candidate['gene_aliases']
    for locus, target in alias.items():
        hits = sorted(gid for gid, name in base_map.items() if name == locus)
        require(len(hits) <= 1, f'Ambiguous existing source locus: {locus}')
        equal(target, hits[0] if hits else locus, 'Candidate alias ' + locus)
    identities = candidate['gene_identity_records']
    equal(sorted(r['source_locus'] for r in identities), sorted(alias), 'Candidate identity inventory')
    for record in identities:
        require(not record['identity_mapping_ambiguous'], 'Candidate declares an ambiguous identity')
        equal(record['target_model_gene_id'], alias[record['source_locus']], 'Candidate identity alias')
    reports = {}
    for config in recipe['configurations']:
        arm, add = config['id'], config['add_pathway']
        model = models[arm]
        intervention = read_json(directory / arm / 'intervention.json')
        equal(intervention['configuration'], config, arm + ' configuration')
        equal(intervention['applied_biomass_coefficient'], coefficients[config['biomass_coefficient_source']], arm + ' demand coefficient')
        expected_reactions = {r.id for r in base.reactions} | (set(REACTION_IDS) if add else set())
        equal(sorted(r.id for r in model.reactions), sorted(expected_reactions), arm + ' reaction inventory')
        equal(sorted(m.id for m in model.metabolites), sorted({m.id for m in base.metabolites} | (set(new_mets) if add else set())), arm + ' metabolite inventory')
        equal(model.objective_direction, base.objective_direction, arm + ' objective direction')
        equal({r.id: float(c) for r, c in cobra.util.solver.linear_reaction_coefficients(model).items()},
              {r.id: float(c) for r, c in cobra.util.solver.linear_reaction_coefficients(base).items()}, arm + ' objective')
        for met in base.metabolites:
            equal(met_record(model.metabolites.get_by_id(met.id)), met_record(met), arm + ' preserved metadata ' + met.id)
        changes = {record['reaction']: record for record in intervention['gpr_changes']}
        equal(sorted(changes), sorted(config['gpr_overrides']), arm + ' GPR change inventory')
        for reaction in base.reactions:
            actual = reaction_record(model.reactions.get_by_id(reaction.id))
            expected = reaction_record(reaction)
            if reaction.id == 'Growth':
                coefficient = coefficients[config['biomass_coefficient_source']]
                if coefficient:
                    expected['metabolites']['q8h2_c'] = coefficient
                # Biomass is intentionally a demand; compare its chemistry and
                # metadata directly, without asserting a balance it lacks.
                actual.pop('mass_charge_imbalance'); expected.pop('mass_charge_imbalance')
            if reaction.id in changes:
                expected['gene_reaction_rule'] = changes[reaction.id]['after']
                expected['gene_ids'] = sorted(g.id for g in model.reactions.get_by_id(reaction.id).genes)
            compare_reaction(actual, expected, arm + ' existing reaction ' + reaction.id)
        transfer = intervention['transfer']
        require((transfer is not None) == add, arm + ' transfer presence differs from plan')
        if add:
            equal(transfer['source_reaction_ids'], REACTION_IDS, arm + ' transferred reactions')
            equal(transfer['gene_aliases'], alias, arm + ' transferred aliases')
            equal(transfer['metadata_policy'], recipe['transfer_metadata_policy'], arm + ' metadata policy')
            equal(transfer['new_metabolites'], sorted(new_mets), arm + ' new metabolites')
            equal(transfer['source_balance'], {rid: {} for rid in REACTION_IDS}, arm + ' source balance')
            equal(transfer['preserved_charge_differences'], recipe['transfer_metadata_policy']['expected_shared_charge_differences'], arm + ' preserved charge differences')
            for mid, record in new_mets.items():
                expected = dict(record)
                expected['compartment'] = recipe['transfer_metadata_policy']['compartment_aliases'].get(record['compartment'], record['compartment'])
                equal(met_record(model.metabolites.get_by_id(mid)), expected, arm + ' normalized new metabolite ' + mid)
                equal(transfer['new_metabolite_records'][mid], expected, arm + ' saved new metabolite ' + mid)
            shared = {item['metabolite']: item for item in transfer['shared_metabolite_metadata']}
            equal(sorted(shared), sorted(set(candidate['all_involved_metabolites']) - set(new_mets)), arm + ' shared metadata inventory')
            for mid, record in shared.items():
                equal(record['source'], candidate['all_involved_metabolites'][mid], arm + ' source metadata ' + mid)
                equal(record['preserved_target'], met_record(base.metabolites.get_by_id(mid)), arm + ' target metadata ' + mid)
            applied = {record['id']: record for record in transfer['applied_reactions']}
            equal(sorted(applied), sorted(REACTION_IDS), arm + ' applied reaction record inventory')
            for source_record in candidate['reaction_records']:
                rid = source_record['id']
                expected = dict(source_record)
                expected['gene_reaction_rule'] = renamed_rule(source_record['gene_reaction_rule'], alias)
                expected['gene_ids'] = sorted(alias[gid] for gid in source_record['gene_ids'])
                actual = reaction_record(model.reactions.get_by_id(rid))
                expected['mass_charge_imbalance'] = applied[rid]['mass_charge_imbalance']
                compare_reaction(applied[rid], expected, arm + ' transfer record ' + rid)
                if rid in changes:
                    expected['gene_reaction_rule'] = changes[rid]['after']
                    expected['gene_ids'] = sorted(g.id for g in model.reactions.get_by_id(rid).genes)
                compare_reaction(actual, expected, arm + ' transferred reaction ' + rid)
        mapping = read_json(directory / arm / 'gene_map.json')
        require(mapping['org_id'] == 'Putida', arm + ' gene map organism differs')
        for target, locus in base_map.items():
            equal(mapping['model_to_browser'].get(target), locus, arm + ' preserved baseline gene identity ' + target)
        for rid, change in changes.items():
            before = (next(r for r in transfer['applied_reactions'] if r['id'] == rid)['gene_reaction_rule']
                      if rid in REACTION_IDS else base.reactions.get_by_id(rid).gene_reaction_rule)
            equal(canonical_rule(change['before']), canonical_rule(before), arm + ' GPR before ' + rid)
            equal(change['source_locus_rule'], config['gpr_overrides'][rid], arm + ' source-locus override ' + rid)
            expected = renamed_rule(config['gpr_overrides'][rid], change['gene_aliases'])
            equal(canonical_rule(change['after']), canonical_rule(expected), arm + ' aliased override ' + rid)
            for locus, target in change['gene_aliases'].items():
                original_hits = sorted(gid for gid, value in base_map.items() if value == locus)
                expected_target = original_hits[0] if original_hits else alias.get(locus, locus)
                equal(target, expected_target, arm + ' override reuses original identity ' + locus)
                equal(mapping['model_to_browser'].get(target), locus, arm + ' override alias identity ' + locus)
                equal(sorted(gid for gid, name in mapping['model_to_browser'].items() if name == locus), [target], arm + ' unique override identity ' + locus)
        model_gene_ids = {g.id for g in model.genes}
        require(set(mapping['model_to_browser']) <= model_gene_ids, arm + ' map names genes absent from model')
        equal(sorted(mapping['unmapped_model_genes']), sorted(model_gene_ids - set(mapping['model_to_browser'])), arm + ' unmapped gene inventory')
        expected_genes = {g.id for g in base.genes}
        if add:
            expected_genes.update(alias.values())
        for change in changes.values():
            expected_genes.update(change['gene_aliases'].values())
        equal(sorted(model_gene_ids), sorted(expected_genes), arm + ' model gene inventory (including retained orphans)')
        reports[arm] = {'model': model, 'mapping': mapping, 'intervention': intervention,
                        'n_reactions': len(model.reactions), 'n_metabolites': len(model.metabolites), 'n_genes': len(model.genes)}
    reference, variant = models['curated_path_template'], models['sequence_path_template']
    changed = []
    for reaction in reference.reactions:
        other = variant.reactions.get_by_id(reaction.id)
        if canonical_rule(reaction.gene_reaction_rule) != canonical_rule(other.gene_reaction_rule):
            loci = set()
            for r, arm in ((reaction, 'curated_path_template'), (other, 'sequence_path_template')):
                loci.update(reports[arm]['mapping']['model_to_browser'][g.id] for g in r.genes)
            require(not loci & fitness_rows, 'Sequence skip changes a benchmark-covered locus')
            changed.append({'reaction': reaction.id, 'loci': sorted(loci)})
        equal({m.id: c for m, c in reaction.metabolites.items()}, {m.id: c for m, c in other.metabolites.items()}, 'Sequence identical stoichiometry ' + reaction.id)
        equal(reaction.bounds, other.bounds, 'Sequence identical bounds ' + reaction.id)
    equivalence = read_json(directory / 'sequence_benchmark_equivalence.json')
    equal(equivalence['changed_reactions'], changed, 'Sequence changed reaction inventory')
    for key in ('stoichiometry_bounds_objective_constraints_identical', 'changed_gprs_use_exclusively_unscored_loci', 'benchmark_omitted_for_sequence_arm'):
        require(equivalence[key] is True, 'Sequence skip check not affirmative: ' + key)
    return reports, coefficients


def verify_physical(directory, recipe, reports, fitness_rows, coefficients, conditions=None, media=None):
    physical_arms = []
    for arm, report in reports.items():
        model, mapping = report['model'], report['mapping']['model_to_browser']
        physical = read_json(directory / arm / 'physical.json')
        equal(physical['condition'], ' | '.join(recipe['physical_condition'][k] for k in ('name', 'media')), arm + ' physical condition')
        equal(physical['settings'], recipe['physical_settings'], arm + ' physical settings')
        certificate = pool_certificate(model, recipe['pool_weights'])
        for key, value in certificate.items():
            equal(physical['pool_certificate'][key], value, arm + ' pool certificate ' + key)
        wt = physical['wild_type']
        require(wt['status'] == 'optimal' and isinstance(wt['growth'], (int, float)) and math.isfinite(wt['growth']) and wt['growth'] >= -1e-6, arm + ' invalid physical WT status/value')
        if not certificate['direction_compatible_supply_reactions'] and any(row['reaction'] == 'Growth' and row['net_pool_coefficient'] < 0 for row in certificate['reactions']):
            require(wt['growth'] < recipe['protocol']['growth_threshold'], arm + ' growing WT contradicts pool demand with no supply')
        egc = physical['energy_from_nothing']
        equal(physical['tested_energy_currencies'], sorted(egc), arm + ' EGC currency inventory')
        required_currency_mets = {'atp': {'atp_c', 'h2o_c', 'adp_c', 'pi_c', 'h_c'},
                                 'nadh': {'nadh_c', 'nad_c', 'h_c'}, 'nadph': {'nadph_c', 'nadp_c', 'h_c'},
                                 'q8h2': {'q8h2_c', 'q8_c', 'h_c'}, 'h_p': {'h_p', 'h_c'}}
        available = {m.id for m in model.metabolites}
        equal(sorted(egc), sorted(key for key, mets in required_currency_mets.items() if mets <= available), arm + ' all available EGC currencies tested')
        require(egc and all(isinstance(v, (float, int)) and math.isfinite(v) and -1e-6 <= v <= recipe['physical_settings']['egc_threshold'] for v in egc.values()), arm + ' invalid/positive EGC result')
        if conditions is not None:
            excluded = recipe['protocol'].get('medium_completion_exclude', ['pnto__R', 'fol', 'hco3'])
            names = {record['media'] for record in conditions.values() if record['bigg_ids']}
            components = {component for name in names for component in medium_components(media[name])}
            added = sorted(f'EX_{component}_e' for component in components if f'EX_{component}_e' not in model.reactions and
                           f'{component}_c' in model.metabolites and component not in excluded)
            if not recipe['protocol']['complete_medium_transport']:
                added = []
            equal(physical['medium_completion_added'], added, arm + ' declared medium completion inventory')
            selected = medium_components(media[recipe['physical_condition']['media']])
            missing = [f'EX_{component}_e' for component in selected if f'EX_{component}_e' not in model.reactions and f'EX_{component}_e' not in added]
            equal(physical['missing_medium_components'], missing, arm + ' missing physical medium inventory')
        production = physical['production_at_growth_lower_bound_zero']
        equal([row['metabolite'] for row in production], list(recipe['pool_weights']), arm + ' production probe inventory')
        for row in production:
            rate = row['maximum_demand_flux']
            require(row['status'] == 'optimal' and isinstance(rate, (int, float)) and math.isfinite(rate), arm + ' invalid production solve')
            equal(row['threshold'], recipe['physical_settings']['production_threshold'], arm + ' probe threshold')
            equal(row['probe_capacity'], recipe['physical_settings']['production_capacity'], arm + ' probe capacity')
            require(-1e-6 <= rate <= row['probe_capacity'] + 1e-6, arm + ' production exceeds probe bounds')
            equal(row['producible_at_threshold'], rate >= row['threshold'], arm + ' production threshold classification')
            equal(row['capacity_reached'], rate >= row['probe_capacity'] or math.isclose(rate, row['probe_capacity'], rel_tol=1e-9, abs_tol=0), arm + ' production capacity classification')
            if not certificate['direction_compatible_supply_reactions']:
                require(rate < row['threshold'], arm + ' positive production contradicts exact pool certificate')
        knockouts = physical['targeted_gene_predictions']
        equal([row['locus'] for row in knockouts], recipe['targeted_gene_loci'], arm + ' targeted gene inventory')
        for row in knockouts:
            hits = sorted(gid for gid, locus in mapping.items() if locus == row['locus'])
            require(len(hits) <= 1, arm + ' ambiguous targeted gene identity')
            gid = hits[0] if hits else row['locus']
            equal(row['model_gene_id'], gid if gid in model.genes else None, arm + ' targeted model identity')
            associated = sorted(r.id for r in model.genes.get_by_id(gid).reactions) if gid in model.genes else []
            equal(row['associated_reactions'], associated, arm + ' targeted reaction inventory ' + row['locus'])
            equal(row['exported_fitness_row_present'], row['locus'] in fitness_rows, arm + ' targeted fitness coverage')
            equal(row['wild_type_growth'], wt['growth'], arm + ' targeted WT growth')
            equal(row['growth_threshold'], recipe['protocol']['growth_threshold'], arm + ' targeted growth threshold')
            grows = wt['growth'] >= recipe['protocol']['growth_threshold']
            equal(row['wild_type_grows_at_threshold'], grows, arm + ' targeted WT classification')
            require(grows or 'Uninterpretable for gene essentiality' in row['interpretation'], arm + ' no-growing WT interpreted as essentiality evidence')
            if not associated:
                equal(row['status'], 'unrepresented_function' if gid in model.genes else 'unrepresented_gene', arm + ' absent/orphan status')
                require(row['growth'] is None, arm + ' absent/orphan gene has a fabricated growth value')
            else:
                require(row['status'] == 'optimal' and isinstance(row['growth'], (int, float)) and math.isfinite(row['growth']), arm + ' invalid targeted solve')
                require(-1e-6 <= row['growth'] <= wt['growth'] + 1e-5, arm + ' knockout growth exceeds WT feasible optimum')
        report['physical'] = physical
        physical_arms.append({'arm': arm, 'wild_type': wt, 'production': production})
    stored = read_json(directory / 'physical_summary.json')
    equal(stored['arms'], physical_arms, 'Physical summary arms')
    equal(stored['coefficients'], coefficients, 'Physical summary biomass coefficients')
    return physical_arms


def load_benchmark(directory):
    card = read_json(directory / 'card.json')
    with np.load(directory / 'matrices.npz', allow_pickle=False) as archive:
        expected = {'sim_growth', 'fitness', 'wt_growth', 'browser_genes', 'model_genes', 'conditions'}
        equal(sorted(archive.files), sorted(expected), directory.name + ' matrix inventory')
        run = {key: archive[key].copy() for key in expected}
    for key in ('browser_genes', 'model_genes', 'conditions'):
        require(run[key].ndim == 1 and len(set(run[key])) == len(run[key]), directory.name + ' duplicate/invalid ' + key)
    shape = (len(run['browser_genes']), len(run['conditions']))
    require(run['sim_growth'].shape == run['fitness'].shape == shape and run['wt_growth'].shape == (shape[1],) and len(run['model_genes']) == shape[0], directory.name + ' mismatched matrix axes')
    run.update(st=card['protocol']['params']['growth_threshold'], ft=card['protocol']['params']['fitness_threshold'], card=card)
    return run


def verify_benchmark(directory, report, recipe, manifest, fitness, experiments, conditions=None):
    arm, model = directory.name, report['model']
    run = load_benchmark(directory)
    card = run['card']
    for key, value in recipe['protocol'].items():
        equal(card['protocol']['params'][key], value, arm + ' protocol ' + key)
    for key, value in {'max_conditions': None, 'genes_subset': None, 'knockout_genes': [],
                       'rich_medium_uptake': -1000.0, 'medium_completion_exclude': ['pnto__R', 'fol', 'hco3']}.items():
        equal(card['protocol']['params'][key], recipe['protocol'].get(key, value), arm + ' protocol default ' + key)
    equal(card['protocol']['intervention'], report['intervention'], arm + ' card intervention')
    equal(card['protocol']['study_fingerprint'], manifest['content_fingerprint'], arm + ' card fingerprint')
    equal(card['protocol']['evaluation_role'], 'development', arm + ' evaluation role')
    equal(card['warnings'], recipe['caveats'], arm + ' declared card caveats')
    equal(card['dataset_provenance']['n_genes_with_fitness'], str(len(fitness)), arm + ' dataset gene coverage')
    equal(card['dataset_provenance']['n_experiments'], str(len(fitness.columns)), arm + ' dataset experiment coverage')
    require(Path(card['model']['file']).resolve() == (directory / 'model.xml.gz').resolve(), arm + ' model provenance points elsewhere')
    for key, value in [('n_reactions', len(model.reactions)), ('n_metabolites', len(model.metabolites)), ('n_genes', len(model.genes)), ('model_id', model.id), ('sha256', sha(directory / 'model.xml.gz'))]:
        equal(card['model'][key], value, arm + ' model provenance ' + key)
    mapping = report['mapping']['model_to_browser']
    expected_pairs, seen = [], set()
    for gene in model.genes:
        locus = mapping.get(gene.id)
        if locus in fitness.index and locus not in seen:
            seen.add(locus)
            expected_pairs.append((gene.id, locus))
    equal(run['model_genes'].tolist(), [pair[0] for pair in expected_pairs], arm + ' native model-gene coverage')
    equal(run['browser_genes'].tolist(), [pair[1] for pair in expected_pairs], arm + ' native Browser-gene coverage')
    require(np.isfinite(run['sim_growth']).all() and np.isfinite(run['wt_growth']).all(), arm + ' saved numerical simulation failure')
    require((run['sim_growth'] >= -1e-6).all() and (run['wt_growth'] >= -1e-6).all(), arm + ' negative growth')
    grows = np.isfinite(run['wt_growth']) & (run['wt_growth'] >= run['st'])
    require((run['sim_growth'][:, ~grows] == 0).all(), arm + ' nonzero KO predictions where WT does not grow')
    table = pd.read_table(directory / 'conditions.tsv', keep_default_na=False)
    per_condition = pd.read_table(directory / 'per_condition_metrics.tsv')
    equal((table.condition + ' | ' + table.media).tolist(), run['conditions'].tolist(), arm + ' condition table axis')
    equal((per_condition.condition + ' | ' + per_condition.media).tolist(), run['conditions'].tolist(), arm + ' condition metric axis')
    if conditions is not None:
        expected_conditions = [key for key, record in conditions.items() if record['bigg_ids']]
        equal(run['conditions'].tolist(), expected_conditions, arm + ' complete mapped condition coverage')
    for j, key in enumerate(run['conditions']):
        name, media = key.split(' | ')
        selected = experiments[(experiments.expGroup == 'carbon source') & (experiments.condition_1 == name) &
                               (experiments.media == media) & experiments.condition_2.isin(['', 'Dimethyl Sulfoxide'])].expName
        selected = [value for value in selected if value in fitness.columns]
        require(selected, arm + ' condition has no exported fitness experiments: ' + key)
        observed = fitness.loc[run['browser_genes'], selected].mean(axis=1).to_numpy(dtype=float)
        require(np.allclose(observed, run['fitness'][:, j], atol=1e-10, rtol=0, equal_nan=True), arm + ' raw replicate-average mismatch: ' + key)
        equal(table.iloc[j]['wt_growth'], run['wt_growth'][j], arm + ' condition WT value')
        if conditions is not None:
            mapping_record = conditions[key]
            equal(table.iloc[j]['bigg_ids'], mapping_record['bigg_ids'], arm + ' condition metabolite mapping')
            equal(table.iloc[j]['mapping_confidence'], mapping_record['confidence'], arm + ' condition mapping confidence')
            equal(int(table.iloc[j]['n_experiments']), len(mapping_record['experiments']), arm + ' condition replicate coverage')
            missing = [f'EX_{met}_e' for met in mapping_record['bigg_ids'].split(';') if
                       f'EX_{met}_e' not in model.reactions and f'EX_{met}_e' not in report['physical']['medium_completion_added']]
            equal(table.iloc[j]['missing_exchanges'], ';'.join(missing), arm + ' absent carbon exchange inventory')
        equal(bool(table.iloc[j]['wt_grows']), bool(grows[j]), arm + ' condition WT class')
        equal(int(table.iloc[j]['n_ko_no_growth']), int((run['sim_growth'][:, j] < run['st']).sum()), arm + ' condition KO count')
        equal(int(table.iloc[j]['n_genes_fitness_below_threshold']), int((run['fitness'][:, j] < run['ft']).sum()), arm + ' condition fitness count')
        row = per_condition.iloc[j]
        equal(float(row['wt_growth']), run['wt_growth'][j], arm + ' per-condition WT value')
        finite = np.isfinite(run['sim_growth'][:, j]) & np.isfinite(run['fitness'][:, j])
        equal(int(row['n_genes']), int(finite.sum()), arm + ' per-condition finite pairs')
        if grows[j] and finite.any():
            for field, value in confusion(run['sim_growth'][:, j], run['fitness'][:, j], run['st'], run['ft']).items():
                equal(row[field], value, arm + ' per-condition ' + field)
            points = metric_points(run['sim_growth'][:, j], run['fitness'][:, j], run['st'], run['ft'])
            for field in ('mcc', 'aucpr_bernstein', 'auroc_standard'):
                equal(row[field], points[field], arm + ' per-condition ' + field)
        else:
            require(all(field not in row or pd.isna(row[field]) for field in ('tp', 'tn', 'fp', 'fn', 'mcc', 'aucpr_bernstein', 'auroc_standard')), arm + ' non-growing/unscored condition has metrics')
    native = native_score(run)
    equal(card['results']['gene_level_conditions_where_wt_grows'], native, arm + ' native score and intervals')
    equal(card['results']['condition_level'], {'n_conditions_mapped': len(run['conditions']), 'n_conditions_wt_grows': int(grows.sum())}, arm + ' native condition coverage')
    counts = card['results']['counts']
    if conditions is not None:
        equal(counts['conditions_total'], len(conditions), arm + ' all eligible condition coverage')
    for key, value in {'model_genes': len(model.genes), 'model_genes_mapped': len(mapping), 'genes_with_fitness': len(expected_pairs),
                       'genes_after_adjustment': len(expected_pairs), 'conditions_mapped': len(run['conditions']),
                       'conditions_wt_grows': int(grows.sum()), 'medium_completion_exchanges_added': len(report['physical']['medium_completion_added'])}.items():
        equal(counts[key], value, arm + ' coverage ' + key)
    physical_key = report['physical']['condition']
    require(physical_key in run['conditions'], arm + ' physical condition absent from benchmark')
    index = run['conditions'].tolist().index(physical_key)
    equal(run['wt_growth'][index], report['physical']['wild_type']['growth'], arm + ' physical/benchmark WT cross-check', atol=1e-5)
    for row in report['physical']['targeted_gene_predictions']:
        if row['locus'] in run['browser_genes'] and row['status'] == 'optimal' and grows[index]:
            gene_index = run['browser_genes'].tolist().index(row['locus'])
            equal(run['sim_growth'][gene_index, index], row['growth'], arm + ' physical/benchmark targeted KO', atol=1e-5)
    return run, native


def verify(root, directory):
    require(directory.is_dir(), 'Run directory does not exist')
    require(not (directory / 'failure.json').exists(), 'Run has a saved failure marker')
    manifest = verify_manifest(root, directory)
    recipe = manifest['plan']['recipe']
    equal([config['id'] for config in recipe['configurations']], ARM_IDS, 'Six-arm recipe inventory')
    # Independent verification and source-snapshot folders may coexist with
    # arms. A directory carrying any arm-defining artifact must be declared.
    arm_files = ('model.xml.gz', 'gene_map.json', 'intervention.json', 'physical.json', 'matrices.npz', 'card.json')
    equal(sorted(p.name for p in directory.iterdir() if p.is_dir() and any((p / name).exists() for name in arm_files)),
          sorted(ARM_IDS), 'Saved arm directory inventory')
    require(recipe['organism'] == 'Putida' and recipe['evaluation_role'] == 'development', 'Unexpected study scope')
    equal([config['run_benchmark'] for config in recipe['configurations']], [True] * 5 + [False], 'Benchmark arm inventory')
    candidate = read_json(root / recipe['candidate_path'])
    fitness, experiments, genes = read_fitness(root)
    conditions, media = condition_inventory(root, fitness, experiments)
    reports, coefficients = verify_models(root, directory, recipe, candidate, set(fitness.index))
    physical = verify_physical(directory, recipe, reports, set(fitness.index), coefficients, conditions, media)
    summary = read_json(directory / 'summary.json')
    equal(summary['study'], recipe, 'Summary recipe')
    equal(summary['study_fingerprint'], manifest['content_fingerprint'], 'Summary fingerprint')
    equal(summary['physical_arms'], physical, 'Summary physical arms')
    runs, native_reports = {}, {}
    equal([record['arm'] for record in summary['runs']], ARM_IDS[:5], 'Summary scored arm inventory')
    for arm, record in zip(ARM_IDS[:5], summary['runs']):
        run, native = verify_benchmark(directory / arm, reports[arm], recipe, manifest, fitness, experiments, conditions)
        runs[arm] = run
        equal({k: v for k, v in record.items() if k != 'arm'}, run['card']['results'], arm + ' summary/card results')
        native_reports[arm] = {'n_mapped_genes_with_exported_fitness': len(run['browser_genes']),
                              'n_conditions_mapped': len(run['conditions']), 'native_scored': native,
                              'mcc': native.get('mcc', {}).get('point'),
                              'mcc_defined': native.get('mcc', {}).get('point') is not None and np.isfinite(native.get('mcc', {}).get('point', np.nan))}
        print('Verified native scores, coverage, and saved model:', arm, flush=True)
    sequence_dir = directory / ARM_IDS[-1]
    require(not any((sequence_dir / name).exists() for name in ('card.json', 'card.md', 'matrices.npz', 'conditions.tsv', 'per_condition_metrics.tsv')), 'Sequence benchmark was run despite declared omission')
    comparisons, changes = {}, []
    equal([record['arm'] for record in summary['comparisons']], ARM_IDS[1:5], 'Comparison arm inventory')
    for record in summary['comparisons']:
        arm = record['arm']
        comparison = paired_comparison(runs['baseline'], runs[arm], resamples=recipe['comparison']['resamples'], seed=recipe['comparison']['seed'])
        for key, value in comparison.items():
            equal(record[key], value, arm + ' aligned comparison ' + key)
        require(record['evaluation_role'] == 'retrospective_development', arm + ' comparison exposure role')
        comparisons[arm] = comparison
        changes.extend(changed_predictions(runs['baseline'], runs[arm], arm))
    equal(read_json(directory / 'changed_predictions.json'), changes, 'Changed prediction records')
    equal(summary['n_changed_gene_condition_predictions'], len(changes), 'Changed prediction count')
    return {'study_fingerprint': manifest['content_fingerprint'], 'verified_arms': ARM_IDS,
            'fitness_coverage': {'annotated_genes': len(genes), 'exported_fitness_rows': len(fitness)},
            'native': native_reports, 'aligned_comparisons': comparisons, 'changed_prediction_count': len(changes),
            'model_inventories': {arm: {k: report[k] for k in ('n_reactions', 'n_metabolites', 'n_genes')} for arm, report in reports.items()},
            'physical_checks_scope': 'Stored solver statuses, finite values, threshold labels, pool row arithmetic, gene/reaction inventories, and cross-artifact agreement only; physical optima and energy probes were not independently re-solved.',
            'baseline_scope': 'All interventions checked against the saved baseline and hashed source candidate. Baseline preparation itself was not independently reconstructed.',
            'scientific_scope': 'Internal reproducibility on exposed development data; no physiological or independent predictive validation.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, default=ROOT / 'results/quinone_repair_2026_09_06/runs/main')
    parser.add_argument('--out', type=Path, help='New report path; default is independent_verification.json in the run directory')
    args = parser.parse_args()
    directory = args.run_dir.resolve()
    output = (args.out or directory / 'independent_verification.json').resolve()
    require(not output.exists(), f'Preserve the existing report; choose another --out: {output}')
    started = time.time()
    inputs = {str(path.relative_to(directory)): sha(path) for path in sorted(directory.rglob('*')) if path.is_file() and path.resolve() != output} if directory.is_dir() else {}
    report = {'created_utc': datetime.now(timezone.utc).isoformat(), 'passed': False,
              'verification_script': {'path': str(Path(__file__).resolve()), 'sha256': sha(__file__)},
              'runtime': {'python': platform.python_version(), 'platform': platform.platform(),
                          **{name: importlib.metadata.version(name) for name in ('numpy', 'pandas', 'cobra')}},
              'run_directory': str(directory), 'result_input_sha256': inputs}
    try:
        report.update(verify(ROOT.resolve(), directory))
        for relative, digest in inputs.items():
            equal(sha(directory / relative), digest, 'Result changed during verification: ' + relative)
        equal(sorted(str(path.relative_to(directory)) for path in directory.rglob('*') if path.is_file() and path.resolve() != output),
              sorted(inputs), 'Result artifact inventory changed during verification')
        report['passed'] = True
    except Exception as error:
        report['failure'] = {'type': type(error).__name__, 'message': str(error)}
    report['seconds'] = time.time() - started
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x') as stream:
        json.dump(safe(report), stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'passed': report['passed'], 'report': str(output), 'failure': report.get('failure')}))
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
