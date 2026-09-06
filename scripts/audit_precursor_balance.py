"""Independently audit saved precursor fluxes using XML and exact rational rows.

No optimizer, model library, runner helper, or phenotype parser is imported.
The audit covers the stated precursor equations, not global LP optimality.
Run after summary.json exists. Both output files must be fresh.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / 'results/quinone_precursor_2026_09_06/runs/main'
NS = {'s': 'http://www.sbml.org/sbml/level3/version1/core'}
PRECURSOR = ['4hbz_c', '4hbz_p', '4hbz_e', '4hbzcoa_c']
TAIL = ['3ophb_c', '2oph_c', '2ohph_c', '2omph_c', '2ombzl_c',
        '2ommbl_c', '2omhmbl_c', 'q8_c', 'q8h2_c']
COUMARATE = ['T4hcinnm_c', 'T4hcinnm_p', 'T4hcinnm_e', 'coucoa_c', '4hbald_c']
POOLS = {'total': PRECURSOR + TAIL + COUMARATE,
         'tail': TAIL, 'quinone': ['q8_c', 'q8h2_c'],
         'external_periplasm': ['4hbz_e', '4hbz_p'],
         'external': ['4hbz_e'], 'coumarate': COUMARATE}


def require(test, message):
    if not test:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, f'Duplicate JSON key: {key}')
            result[key] = value
        return result
    return json.loads(path.read_text(), object_pairs_hook=unique)


def exact_rows(path):
    model = ET.fromstring(gzip.decompress(path.read_bytes())).find('s:model', NS)
    require(model is not None, 'Unexpected SBML namespace')
    parameters = {r.attrib['id']: float(r.attrib['value'])
                  for r in model.findall('s:listOfParameters/s:parameter', NS)}
    reactions, bounds = {}, {}
    for reaction in model.findall('s:listOfReactions/s:reaction', NS):
        require(reaction.attrib['id'].startswith('R_'), 'Unexpected reaction encoding')
        rid = reaction.attrib['id'][2:]
        require(rid not in reactions, 'Duplicate reaction')
        coefficients = {}
        for side, sign in [('listOfReactants', -1), ('listOfProducts', 1)]:
            for species in reaction.findall(f's:{side}/s:speciesReference', NS):
                mid = species.attrib['species'][2:]
                value = float(species.attrib.get('stoichiometry', '1'))
                require(math.isfinite(value), 'Nonfinite stoichiometry')
                coefficients[mid] = coefficients.get(mid, Fraction(0)) + sign * Fraction.from_float(value)
        reactions[rid] = coefficients
        attrs = {key.split('}')[-1]: value for key, value in reaction.attrib.items()}
        bounds[rid] = [parameters[attrs[key]] for key in ('lowerFluxBound', 'upperFluxBound')]
    rows = {}
    for pool, mids in POOLS.items():
        rows[pool] = {rid: value for rid, coefficients in reactions.items()
                      if (value := sum((coefficients.get(mid, Fraction(0)) for mid in mids), Fraction(0)))}
    epsilon = -reactions['Growth']['q8h2_c']
    require(epsilon > 0, 'No positive quinone demand')
    require(rows['total'] == {'CHRPL': 1, 'EX_4hbz_e': -1, 'EX_T4hcinnm_e': -1,
                              'Growth': -epsilon, '4HBHYOX': -1, 'sink_2ohph_c': -1}, 'Changed total row')
    require(rows['tail'] == {'HBZOPT': 1, 'sink_2ohph_c': -1, 'Growth': -epsilon}, 'Changed tail row')
    require(rows['quinone'] == {'DMQMT': 1, 'Growth': -epsilon}, 'Changed quinone row')
    require(rows['external_periplasm'] == {'UHBZ1t_pp': -1, 'EX_4hbz_e': -1}, 'Changed transport row')
    require(rows['external'] == {'4HBZtex': -1, 'EX_4hbz_e': -1}, 'Changed outer transport row')
    for rid in ('CHRPL', 'Growth', '4HBHYOX', 'sink_2ohph_c'):
        require(bounds[rid][0] >= 0, f'Negative lower bound invalidates inequality: {rid}')
    species_ids = {r.attrib['id'][2:] for r in model.findall('s:listOfSpecies/s:species', NS)}
    return rows, bounds, epsilon, set(reactions), species_ids


def completion_rows(contexts, reaction_ids, species_ids):
    """Independently express the frozen protocol's exchange + inward MEDt recipe.

    These standard vitamin/ion components are outside every asserted pool. The
    source algorithm is read and hashed, not imported or executed by this audit.
    """
    expected = {'EX_4abz_e', 'EX_btn_e', 'EX_na1_e', 'EX_nac_e',
                'EX_ni2_e', 'EX_ribflv_e', 'EX_thm_e'}
    require(all(set(r['medium_completion_added']) == expected and
                len(r['medium_completion_added']) == len(expected) for r in contexts),
            'Unexpected completion inventory; independent source review required')
    added = {}
    for exchange in sorted(expected):
        extracellular, stem = exchange[3:], exchange[3:-2]
        cytoplasm, carrier = stem + '_c', 'MEDt_' + stem
        require(cytoplasm in species_ids, 'Completion lacks cytoplasmic target')
        require(exchange not in reaction_ids and carrier not in reaction_ids, 'Completion reaction collision')
        added[exchange] = {extracellular: -1}
        added[carrier] = {extracellular: -1, cytoplasm: 1}
    for coefficients in added.values():
        for mids in POOLS.values():
            require(sum(coefficients.get(mid, 0) for mid in mids) == 0, 'Completion changes a pool row')
    return added


def analyze(run):
    files = [run / name for name in ('summary.json', 'contexts.json', 'cases.jsonl.gz',
                                     'manifest.json', 'input_verification_after_run.json')]
    require(all(path.is_file() for path in files), 'Completed run artifacts are required')
    summary, contexts = read_json(files[0]), read_json(files[1])
    plan, manifest = summary['study'], read_json(files[3])
    parent = ROOT / plan['base_model']
    completion_source = ROOT / 'gembench/protocols/carbon_fitness_generic.py'
    files += [parent, Path(__file__).resolve(), completion_source]
    hashes = {str(path.relative_to(ROOT)): sha(path) for path in files}
    require(sha(parent) == plan['base_model_sha256'], 'Parent model hash mismatch')
    require(summary['study_fingerprint'] == manifest['content_fingerprint'], 'Fingerprint mismatch')
    require(read_json(files[4])['valid'], 'Main run reported changed inputs')
    for item in manifest['files']:
        path = ROOT / item['path']
        require(sha(path) == item['sha256'], f'Frozen input changed: {path}')
    rows, bounds, epsilon, reaction_ids, species_ids = exact_rows(parent)
    eps, tol = float(epsilon), plan['residual_tolerance']
    by_context = {r['condition']: r for r in contexts}
    require(len(by_context) == len(contexts), 'Duplicate condition context')
    added = completion_rows(contexts, reaction_ids, species_ids)
    reaction_ids |= set(added)
    with gzip.open(files[2], 'rt') as handle:
        cases = [json.loads(line) for line in handle if line.strip()]
    indexed = {(r['specification']['condition'], r['specification']['id']): r for r in cases}
    require(len(indexed) == len(cases) == summary['n_cases'], 'Duplicate/missing case')
    require(summary['n_solver_runs'] == 2 * len(cases), 'Solver case count mismatch')
    residuals = {pool: 0.0 for pool in rows}
    numerical_agreement = 0.0
    minima, doses, controls, donors = [], [], [], []
    kind_counts = Counter()
    for case in cases:
        spec = case['specification']
        context = by_context[spec['condition']]
        kind_counts[spec['kind']] += 1
        require({s['solver'] for s in case['solves']} == {'glpk', 'highs'} and len(case['solves']) == 2,
                'Require one witness from each solver')
        numerical_agreement = max(numerical_agreement, abs(case['solves'][0]['objective'] - case['solves'][1]['objective']))
        for solve in case['solves']:
            v = solve['full_primal_fluxes']
            require(solve['status'].lower() == 'optimal' and set(v) == reaction_ids, 'Unexpected witness status/reactions')
            require(all(math.isfinite(x) for x in v.values()), 'Nonfinite witness')
            require(abs(solve['objective'] - v[spec.get('objective', 'Growth')]) <= tol, 'Objective mismatch')
            require(abs(solve['growth_flux'] - v['Growth']) <= tol, 'Growth mismatch')
            for pool, coefficients in rows.items():
                residual = abs(float(sum((c * Fraction.from_float(v[r]) for r, c in coefficients.items()), Fraction(0))))
                residuals[pool] = max(residuals[pool], residual)
                require(residual <= tol, f'Balance residual: {pool}/{spec["condition"]}/{spec["id"]}')
            for rid in spec.get('close_reactions', []):
                require(abs(v[rid]) <= tol, 'Closed reaction carries flux')
            if 'PP_5317' in spec.get('loci', []):
                require(abs(v['CHRPL']) <= tol, 'Deleted CHRPL carries flux')
            if 'PP_5318' in spec.get('loci', []):
                require(abs(v['HBZOPT']) <= tol, 'Deleted HBZOPT carries flux')
            if 'PP_1376' in spec.get('loci', []):
                require(abs(v['UHBZ1t_pp']) <= tol, 'Deleted inner transporter carries flux')
            if 'growth_bounds' in spec:
                require(spec['growth_bounds'][0] - tol <= v['Growth'] <= spec['growth_bounds'][1] + tol, 'Growth constraint violated')
            lower, upper = spec.get('exchange_bounds', context['exchange_bounds']['EX_4hbz_e'])
            require(lower - tol <= v['EX_4hbz_e'] <= upper + tol, 'Exchange bound violated')
            if solve['solver'] == 'glpk':
                require(all(value == v[rid] for rid, value in case['selected_glpk_fluxes'].items()), 'Selected witness differs')
        v = next(s['full_primal_fluxes'] for s in case['solves'] if s['solver'] == 'glpk')
        require(case['growth'] == v['Growth'], 'Case growth disagrees with saved witness')
        wt = indexed[(spec['condition'], 'wild_type')]['growth']
        if spec['kind'] == 'minimum_uptake':
            coumarate = context['exchange_bounds']['EX_T4hcinnm_e'][0] < 0
            uptake, demand = -v['EX_4hbz_e'], eps * v['Growth']
            require(abs(v['Growth'] - spec['target_growth']) <= tol, 'Minimum target mismatch')
            if not coumarate:
                require(uptake >= demand - tol, 'Necessary uptake bound violated')
            minima.append({'condition': spec['condition'], 'fraction': spec['fraction'], 'growth': v['Growth'],
                           'uptake': uptake, 'epsilon_times_growth': demand,
                           'uptake_minus_demand': uptake - demand, 'coumarate_uptake_allowed': coumarate})
        elif spec['kind'] == 'dose_growth':
            require(context['exchange_bounds']['EX_T4hcinnm_e'][0] >= 0, 'Dose context admits coumarate')
            cap = -spec['exchange_bounds'][0]
            require(eps * v['Growth'] <= cap + tol, 'Dose material bound violated')
            doses.append({'multiplier': spec['dose_multiplier'], 'uptake_cap': cap, 'growth': v['Growth'],
                          'growth_over_unsupplemented_wt': v['Growth'] / wt,
                          'material_upper_bound_growth': cap / eps, 'uptake': -v['EX_4hbz_e']})
        elif spec['kind'] == 'growth_control':
            controls.append({'id': spec['id'], 'growth': v['Growth'], 'uptake': -v['EX_4hbz_e'],
                             'numerically_zero_growth_at_1e8': abs(v['Growth']) <= 1e-8,
                             'grows_at_protocol_threshold': v['Growth'] >= plan['protocol']['growth_threshold'],
                             'HBZOPT': v['HBZOPT'], 'catabolism': v['4HBHYOX']})
        elif spec['kind'] == 'donor_secretion':
            require(context['exchange_bounds']['EX_T4hcinnm_e'][0] >= 0, 'Donor context admits coumarate')
            require(v['EX_4hbz_e'] <= v['CHRPL'] - eps * v['Growth'] + tol, 'Donor material bound violated')
            donors.append({'id': spec['id'], 'maximum_secretion_objective': v['EX_4hbz_e'],
                           'growth': v['Growth'], 'growth_over_unsupplemented_wt': v['Growth'] / wt,
                           'growth_floor': spec['growth_bounds'][0], 'CHRPL': v['CHRPL'],
                           'quinone_demand': eps * v['Growth'], 'catabolism': v['4HBHYOX'],
                           'sink': v['sink_2ohph_c'], 'coumarate_export': v['EX_T4hcinnm_e']})
    for row in summary['condition_results']:
        require(row['wild_type'] == indexed[(row['condition'], 'wild_type')]['growth'], 'Summary WT mismatch')
        require(row['pp5317'] == indexed[(row['condition'], 'pp5317')]['growth'], 'Summary knockout mismatch')
        for item in row['minimum_uptake']:
            case = indexed[(row['condition'], f'minimum_{item["fraction"]}')]
            require(abs(item['minimum_uptake'] - max(0, -case['objective'])) <= tol, 'Summary uptake mismatch')
    require(hashes == {str(p.relative_to(ROOT)): sha(p) for p in files}, 'Audit input changed during reading')
    ordinary = [r for r in minima if not r['coumarate_uptake_allowed']]
    bypass = [r for r in minima if r['coumarate_uptake_allowed']]
    return {'status': 'passed', 'scope': 'Exact source-row derivation and saved-flux interpretation; no optimization or phenotype values.',
            'parent_sha256': sha(parent), 'epsilon_exact': str(epsilon), 'epsilon_float': eps,
            'fraction_semantics': 'Fraction.from_float of serialized numeric coefficients, matching the binary model representation.',
            'exact_rows': {pool: {'metabolites': POOLS[pool], 'nonzero_coefficients': {r: str(c) for r, c in row.items()}}
                           for pool, row in rows.items()},
            'parent_relevant_bounds': {r: bounds[r] for r in sorted(set().union(*rows.values()))},
            'completion_audit': {'n_exchanges_added_per_context': 7, 'n_reactions_added_per_context': len(added),
                'independently_reconstructed_stoichiometry': added, 'all_six_pool_projections_zero': True,
                'source': str(completion_source.relative_to(ROOT)),
                'scope': 'Inherited protocol completion outside the precursor pools; neither new 4HBZ supply nor new precursor transport.'},
            'n_conditions': len(contexts), 'n_cases': len(cases), 'n_saved_solver_witnesses': 2 * len(cases),
            'case_kind_counts': dict(kind_counts), 'max_abs_row_residual_by_pool': residuals,
            'max_solver_objective_difference': numerical_agreement,
            'minima': minima, 'dose_response': doses, 'growth_controls': controls, 'donor_controls': donors,
            'minimum_summary': {'without_coumarate_cases': len(ordinary),
                'max_abs_uptake_minus_required_demand_without_coumarate': max(abs(r['uptake_minus_demand']) for r in ordinary),
                'uptake_range_without_coumarate': [min(r['uptake'] for r in ordinary), max(r['uptake'] for r in ordinary)],
                'with_coumarate_cases': len(bypass),
                'max_abs_uptake_with_coumarate': max(abs(r['uptake']) for r in bypass)},
            'numerical_caveat': 'Absolute growth or exchange values <=1e-8 are operational numerical zeros, not proven biological zeros. A precursor residual allowance of 1e-8 corresponds to 1e-8/epsilon growth units in the derived bound.',
            'growth_equivalent_of_precursor_residual_tolerance': tol / eps,
            'interpretation_limits': [
                'Primal feasibility plus agreeing solver status/objectives is not a saved dual certificate of global optimality; this audit checks equations and interpretation, not the optimizer implementation.',
                'Minimum uptake equality shows feasibility under frozen assumptions, not uptake rates measured in a mutant.',
                'Dose caps are fluxes in mmol/gDW/h, not concentrations; no biomass/time/population balance is supplied.',
                'Maximum donor secretion under a growth floor is a deliberately chosen objective, not predicted spontaneous release, a measured donor rate, or evidence of sharing.',
                'No population coupling, recipient competition, regulation, or extracellular loss model is present.',
                'No replacement UbiC gene assignment, rescore, new reaction, or accepted biological correction follows.'],
            'experiment_input_count': len(manifest['files']), 'experiment_fingerprint': manifest['content_fingerprint'],
            'input_sha256': hashes,
            'auditor_timing': 'Independent additive audit; this script is not an input to the preceding frozen experiment.',
            'checker_refinement': 'Attempt 1 stopped before flux checks because its conservative guard required no medium-completion reactions. Attempt 2 independently reconstructs the seven existing protocol exchange/carrier pairs and proves zero projection onto all six pools. The initial checker source and failure record are preserved. This was checker development, not a failed or altered primary run.'}


def report(record):
    m = record['minimum_summary']
    lines = ['# Independent precursor balance and interpretation audit', '',
        'The saved results satisfy independently derived precursor balances. This audit uses the frozen SBML XML and both saved solver flux vectors; it imports no model library or runner helper and performs no optimization or phenotype analysis.', '',
        'The exact stored coefficient is ε = ' + record['epsilon_exact'] + ' = ' + f"{record['epsilon_float']:.14g}" + ' mmol/gDW.', '',
        '`CHRPL − EX_4hbz_e − EX_T4hcinnm_e = ε Growth + 4HBHYOX + sink_2ohph_c`.', '',
        'Additional rows establish `DMQMT = ε Growth`, `HBZOPT = ε Growth + sink`, and `EX_4hbz_e = −UHBZ1t_pp = −4HBZtex`. With CHRPL disabled and no coumarate import, the nonnegative disposal terms require uptake ≥ ε times growth. Transport closure prevents net supply; HBZOPT closure prevents growth independently of the external supply.', '',
        f"Verified {record['n_saved_solver_witnesses']} saved witnesses across {record['n_cases']} cases and {record['n_conditions']} conditions. The maximum absolute residual among these six independently summed rows was {max(record['max_abs_row_residual_by_pool'].values()):.4g}; the maximum difference between saved solver objectives was {record['max_solver_objective_difference']:.4g}. The inherited protocol adds seven exchange/carrier pairs for vitamins and ions. Their independently reconstructed stoichiometry has zero contribution to all six pools; bounds and deletions also preserve the rows.", '',
        '## Minimum supply and dose caps', '',
        f"Across {m['without_coumarate_cases']} minimum-supply cases without coumarate import, uptake differed from ε times target growth by at most {m['max_abs_uptake_minus_required_demand_without_coumarate']:.4g} mmol/gDW/h. Required uptake ranged from {m['uptake_range_without_coumarate'][0]:.6g} to {m['uptake_range_without_coumarate'][1]:.6g}. The {m['with_coumarate_cases']} coumarate cases required at most {m['max_abs_uptake_with_coumarate']:.4g} absolute 4HBZ uptake. The full condition/fraction records are in the JSON.", '',
        'A cap of d × ε × parental growth imposes mutant growth ≤ d × parental growth. It does not guarantee equality. At larger caps, the supplement can also alter carbon/precursor availability; an exact plateau at unsupplemented growth is not required.', '',
        '| Dose multiplier | Uptake cap (mmol/gDW/h) | Growth (h⁻¹) | Fraction of unsupplemented WT |',
        '|---:|---:|---:|---:|']
    for row in record['dose_response']:
        lines.append(f"| {row['multiplier']:g} | {row['uptake_cap']:.8g} | {row['growth']:.10g} | {row['growth_over_unsupplemented_wt']:.8g} |")
    lines += ['', '## Controls', '', '| Control | Growth (h⁻¹) | 4HBZ uptake (mmol/gDW/h) |', '|---|---:|---:|']
    for row in record['growth_controls']:
        lines.append(f"| {row['id']} | {row['growth']:.10g} | {row['uptake']:.8g} |")
    lines += ['', 'Positive rescue under a declared supply is a feasibility result. Loss of rescue under the transport/downstream controls identifies the represented route; preserved rescue after catabolism closure separates precursor incorporation from the need to degrade the supplement.', '',
              '## Deliberately maximized donor secretion', '', '| Donor control | Maximum secretion (mmol/gDW/h) | Growth (h⁻¹) | Fraction of WT |', '|---|---:|---:|---:|']
    for row in record['donor_controls']:
        lines.append(f"| {row['id']} | {row['maximum_secretion_objective']:.10g} | {row['growth']:.10g} | {row['growth_over_unsupplemented_wt']:.8g} |")
    lines += ['', 'The balance requires secretion ≤ CHRPL flux − ε times growth. These calculations maximize secretion subject to a growth floor. They do not predict spontaneous secretion, demonstrate release by cells, or show that a donor population supplies enough material to recipients. The model includes no sharing dynamics, biomass fractions, accumulation time or extracellular losses.', '',
        '## Numerical and biological limits', '',
        f"Absolute values ≤1e-8 are treated as operational numerical zeros, not biological zeros. Because ε is small, a precursor-balance residual allowance of 1e-8 translates to {record['growth_equivalent_of_precursor_residual_tolerance']:.4g} h⁻¹ in the derived growth inequality; the actual saved row residuals are reported above. The audit checks saved equations and agreement, not a dual optimality certificate.", '',
        'Fluxes in mmol/gDW/h cannot be converted into a medium concentration without a biomass/time balance. Neither successful rescue nor a positive donor optimum establishes carryover, cross-feeding or an alternative native enzyme. All inherited quinone, PpnP, transport and zero ATP-maintenance assumptions remain. No biological correction or new benchmark score is accepted here.', '',
        f"The preceding experiment froze {record['experiment_input_count']} files with fingerprint `{record['experiment_fingerprint']}`. This additive auditor is separate from those experiment inputs. Its own source hash and all directly read artifact hashes are recorded in `precursor_balance.json`; the frozen experiment input hashes were also rechecked.", '',
        'The first checker attempt conservatively stopped on the nonempty medium-completion inventory before checking fluxes. Its source and failure record were retained. The revised checker explicitly reconstructs that inherited completion and verifies pool disjointness. This refinement does not represent a failed primary run or a change to the experiment.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, default=DEFAULT)
    parser.add_argument('--out', type=Path, help='New output basename, default <run>/analysis/precursor_balance')
    args = parser.parse_args()
    prefix = args.out or args.run / 'analysis/precursor_balance'
    paths = [prefix.with_suffix(suffix) for suffix in ('.json', '.md')]
    require(all(not path.exists() for path in paths), 'Refusing to overwrite an audit artifact')
    result = analyze(args.run)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with paths[0].open('x') as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')
    with paths[1].open('x') as handle:
        handle.write(report(result))
    print(json.dumps({'status': result['status'], 'cases': result['n_cases'],
                      'max_row_residual': max(result['max_abs_row_residual_by_pool'].values()),
                      'output': str(prefix)}))


if __name__ == '__main__':
    main()
