"""Post-outcome PP_5317 assay-context audit; no optimization or model edits.

Only Putida fitness values are read. Reconstruct grouping and arithmetic means
independently of the benchmark implementation, retain experiment metadata, and
separate experimental scores from the existing model predictions. Literature
records below are source-reviewed interpretations, not machine-derived facts.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import statistics
import sys

import cobra
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STUDY = Path('results/ppnp_repair_2026_09_06/runs/main')
LOCUS = 'PP_5317'
META = {'orgId', 'locusId', 'sysName', 'geneName', 'desc'}
QUALITY = ('nMapped', 'nGenic', 'nUsed', 'gMed', 'gMedt0', 'cor12', 'mad12', 'gccor', 'adjcor')
HETEROGENEITY = ('mutantLibrary', 'dateStarted', 'setName', 'timeZeroSet', 'concentration_1',
                 'units_1', 'condition_2', 'concentration_2', 'units_2', 'vessel', 'shaking')

SOURCES = [
    {'id': 'wetmore2015', 'title': 'Rapid Quantification of Mutant Fitness in Diverse Bacteria by Sequencing Randomly Bar-Coded Transposons',
     'url': 'https://journals.asm.org/doi/10.1128/mbio.00306-15',
     'sections': ['Mutant fitness profiling', 'Analysis of BarSeq data', 'Assessment of experiment quality'],
     'finding': 'Competitive pooled assay; normalized log2 relative barcode change. Gene scores combine multiple insertion strains. Typical experiments span 4–6 generations. Global quality metrics include gene coverage, half-gene consistency and correlation biases.',
     'boundary': 'Method paper, not sample-specific evidence for KT2440 PP_5317.'},
    {'id': 'price2018', 'title': 'Mutant phenotypes for thousands of bacterial genes of unknown function',
     'url': 'https://doi.org/10.1038/s41586-018-0124-0',
     'author_manuscript_url': 'https://genomics.lbl.gov/supplemental/bigfit/authorfinal.pdf',
     'retrieved_pdf_sha256': '5ba13a263fd4a4493dde8275034a79fc7f09dd2d89398777584d271e82960b1d',
     'retrieved_pdf_size_bytes': 4101000,
     'sections': ['Author manuscript lines 692–727', 'Author manuscript lines 802–827'],
     'finding': 'Typically 4–8 population doublings to saturation; normalized weighted insertion-strain scores, central 10–90% insertions, abundance filtering. Replicate experiments need not use the same concentration.',
     'boundary': 'The later Putida metadata spans multiple projects. These typical methods do not supply exact generations for each exported sample.'},
    {'id': 'thompson2020', 'title': 'Fatty Acid and Alcohol Metabolism in Pseudomonas putida: Functional Analysis Using Random Barcode Transposon Sequencing',
     'url': 'https://journals.asm.org/doi/10.1128/aem.01665-20', 'sections': ['RB–Tn-Seq', 'Plate-based growth assays'],
     'finding': 'JBEI-1 recovered in LB with kanamycin to OD600 0.5; time-zero aliquots; one carbon-free MOPS wash; 1:50 inoculation into 10 mM carbon; 10 mL tubes, 30 °C, 200 rpm. Biological duplicates. Separate individual growth assays use two washes and 1:100.',
     'boundary': 'Protocol consistency with local set15/16 is not a unique experiment-ID mapping; endpoint generations and PP_5317 monoculture viability are not supplied.'},
    {'id': 'incha2020', 'title': 'Leveraging host metabolism for bisdemethoxycurcumin production in Pseudomonas putida',
     'url': 'https://doi.org/10.1016/j.mec.2019.e00119',
     'full_text_url': 'https://escholarship.org/content/qt4rh1g98j/qt4rh1g98j.pdf',
     'sections': ['4.5 RB-TnSeq experiments and analysis; printed page 7'],
     'finding': 'JBEI-1 recovered in LB/kanamycin to OD600 0.5, washed in MOPS, diluted 1:50 into 10 mM aromatic substrates or glucose; 600 µL deep wells at 30 °C/700 rpm; two 600 µL samples combined before BarSeq.',
     'boundary': 'Substrate/format agrees with set12, but no explicit expName link, number of washes, endpoint duration or generation count is given. Engineered production strains belong to separate assays.'},
    {'id': 'bernstein2023', 'title': 'Evaluating E. coli genome-scale metabolic model accuracy with high-throughput mutant fitness data',
     'url': 'https://link.springer.com/article/10.15252/msb.202311566',
     'sections': ['Vitamin/cofactor false-negative predictions and discussion'],
     'finding': 'Primary modelling analysis discusses cofactor carryover and cross-feeding as explanations for mismatches between pooled fitness and isolated-cell model requirements; it warns against adding unsupported biosynthesis solely to remove such mismatches.',
     'boundary': 'Different organism. Cited only for assay/model interpretation; no other-organism phenotype values are extracted and it is not evidence of PP_5317 rescue in KT2440.'},
]


def read_table(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def file_record(path):
    return {'path': str(path.relative_to(ROOT)), 'sha256': sha(path), 'size_bytes': path.stat().st_size}


def number(value):
    if value in ('', 'NA', 'NaN', 'nan'):
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError('Unexpected nonfinite numeric input')
    return result


def media_blocks(path):
    blocks = {}
    current = None
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        if fields[0] == 'Media':
            current = {'name': fields[1], 'properties': {}, 'components': []}
            blocks[fields[1]] = current
        elif fields[0] in ('Description', 'Minimal', 'X'):
            current['properties'][fields[0]] = fields[1]
        elif fields[0] != 'Controlled vocabulary':
            current['components'].append(dict(zip(('name', 'concentration', 'units'), fields)))
    return blocks


def quality_comparison(row):
    """Descriptive historic thresholds; not an official current pass flag."""
    return {'gMed_ge_50': float(row['gMed']) >= 50,
            'mad12_le_0_5': float(row['mad12']) <= .5,
            'cor12_ge_0_1': float(row['cor12']) >= .1,
            'abs_gccor_le_0_2': abs(float(row['gccor'])) <= .2,
            'abs_adjcor_le_0_25': abs(float(row['adjcor'])) <= .25}


def reaction_record(reaction):
    return {'id': reaction.id, 'equation': reaction.reaction,
            'metabolites': {m.id: float(v) for m, v in sorted(reaction.metabolites.items(), key=lambda item: item[0].id)},
            'bounds': list(reaction.bounds), 'gene_reaction_rule': reaction.gene_reaction_rule}


def extract(study):
    folder = ROOT / 'data/fitness_browser/Putida'
    arm = study / 'both_forward'
    input_paths = [folder / name for name in ('genes.tsv', 'experiments.tsv', 'fit_logratios.tsv', 'specific_phenotypes.tsv')]
    input_paths += [ROOT / 'data/reference' / name for name in ('feba_media_extract.tsv', 'fitness_browser_media_bigg.tsv', 'fitness_browser_carbon_sources_bigg.tsv')]
    input_paths += [arm / name for name in ('matrices.npz', 'gene_map.json', 'model.xml.gz', 'card.json', 'conditions.tsv', 'physical.json')]
    input_paths += [study / 'manifest.json', Path(__file__).resolve()]
    inputs = [file_record(path) for path in input_paths]
    genes = [row for row in read_table(folder / 'genes.tsv') if row['locusId'] == LOCUS]
    assert len(genes) == 1
    # Parse only the requested raw gene's numeric phenotype values.
    fitness_rows = [row for row in read_table(folder / 'fit_logratios.tsv') if row['locusId'] == LOCUS]
    assert len(fitness_rows) == 1
    raw = fitness_rows[0]
    columns = {name.split(' ')[0]: name for name in raw if name not in META}
    assert len(columns) == len(raw) - len(META)
    experiments = read_table(folder / 'experiments.tsv')
    media_map = {row['media']: row for row in read_table(ROOT / 'data/reference/fitness_browser_media_bigg.tsv')}
    groups = {}
    for row in experiments:
        if (row['expGroup'] == 'carbon source' and row['condition_2'] in ('', 'Dimethyl Sulfoxide') and
                row['media'] in media_map and row['expName'] in columns):
            key = row['condition_1'] + ' | ' + row['media']
            groups.setdefault(key, []).append(row)
    matrix = np.load(arm / 'matrices.npz', allow_pickle=True)
    indices = [i for i, gene in enumerate(matrix['browser_genes']) if str(gene) == LOCUS]
    assert len(indices) == 1
    i = indices[0]
    mapping = json.loads((arm / 'gene_map.json').read_text())
    aliases = [gene for gene, locus in mapping['model_to_browser'].items() if locus == LOCUS]
    assert aliases == ['NP_747418_1'] and str(matrix['model_genes'][i]) == aliases[0]
    params = json.loads((arm / 'card.json').read_text())['protocol']['params']
    growth_threshold, fitness_threshold = params['growth_threshold'], params['fitness_threshold']
    conditions = []
    for j, key in enumerate(matrix['conditions']):
        key = str(key)
        rows = groups[key]
        samples = []
        for row in rows:
            value = number(raw[columns[row['expName']]])
            samples.append({'expName': row['expName'], 'raw_export_column': columns[row['expName']],
                            'gene_log2_fitness': value, 'gene_t_statistic': None,
                            'global_quality_metrics': {field: number(row[field]) for field in QUALITY},
                            'historic_quality_reference_checks': quality_comparison(row),
                            'original_experiment_metadata': row})
        values = [sample['gene_log2_fitness'] for sample in samples if sample['gene_log2_fitness'] is not None]
        mean = statistics.mean(values) if values else None
        saved_fitness = float(matrix['fitness'][i, j])
        assert mean is not None and abs(mean - saved_fitness) <= 1e-12
        wt, ko = float(matrix['wt_growth'][j]), float(matrix['sim_growth'][i, j])
        assert math.isfinite(wt) and math.isfinite(ko)
        eligible = wt >= growth_threshold
        model_grows = ko >= growth_threshold
        observed_defect = mean < fitness_threshold
        role = 'excluded_model_wt_no_growth' if not eligible else 'disagreement' if model_grows == observed_defect else 'natural_rescue_agreement'
        assert not eligible or not observed_defect
        heterogeneity = {field: sorted({row[field] for row in rows}) for field in HETEROGENEITY
                         if len({row[field] for row in rows}) > 1}
        conditions.append({'condition': key, 'role': role, 'n_exported_samples': len(rows), 'n_finite_samples': len(values),
                           'mean_log2_fitness': mean, 'saved_mean_log2_fitness': saved_fitness,
                           'minimum_sample_fitness': min(values), 'maximum_sample_fitness': max(values),
                           'sample_standard_deviation': statistics.stdev(values) if len(values) > 1 else None,
                           'all_samples_above_defect_threshold': all(v >= fitness_threshold for v in values),
                           'wild_type_growth': wt, 'pp5317_model_growth': ko, 'model_wt_grows': eligible,
                           'model_knockout_grows': model_grows, 'experimental_defect_at_threshold': observed_defect,
                           'heterogeneous_metadata_fields': heterogeneity, 'samples': samples})
    by_role = {}
    for role in sorted({record['role'] for record in conditions}):
        selected = [record for record in conditions if record['role'] == role]
        samples = [sample for record in selected for sample in record['samples']]
        values = [sample['gene_log2_fitness'] for sample in samples]
        by_role[role] = {'n_conditions': len(selected), 'n_samples': len(samples),
                         'sample_fitness_range': [min(values), max(values)],
                         'condition_mean_range': [min(r['mean_log2_fitness'] for r in selected), max(r['mean_log2_fitness'] for r in selected)],
                         'library_labels': dict(Counter(sample['original_experiment_metadata']['mutantLibrary'] for sample in samples)),
                         'n_with_dmso': sum(sample['original_experiment_metadata']['condition_2'] == 'Dimethyl Sulfoxide' for sample in samples),
                         'all_historic_quality_reference_checks_pass': all(all(s['historic_quality_reference_checks'].values()) for s in samples),
                         'quality_metric_ranges': {field: [min(s['global_quality_metrics'][field] for s in samples),
                                                           max(s['global_quality_metrics'][field] for s in samples)] for field in QUALITY}}
    assert by_role['disagreement']['n_conditions'] == 31 and by_role['disagreement']['n_samples'] == 67
    assert by_role['natural_rescue_agreement']['n_conditions'] == 3 and by_role['natural_rescue_agreement']['n_samples'] == 6
    model = cobra.io.read_sbml_model(str(arm / 'model.xml.gz'))
    assert sorted(r.id for r in model.genes.get_by_id(aliases[0]).reactions) == ['CHRPL']
    blocks = media_blocks(ROOT / 'data/reference/feba_media_extract.tsv')
    relevant_media = [blocks[name] for name in ('MOPS minimal media_noCarbon', 'RCH2_defined_noCarbon',
                                               "Wolfe's mineral mix", "Wolfe's vitamin mix", 'MOPS Rich Defined media_noCarbon', 'LB')]
    physical = json.loads((arm / 'physical.json').read_text())
    result = {
        'scope': 'Post-outcome descriptive audit of exposed Putida PP_5317 data only; no optimization, interventions, new independent phenotype validation or other-organism numeric phenotype inspection.',
        'created_at': datetime.now(timezone.utc).isoformat(), 'gene': genes[0],
        'fitness_export_gene_metadata': {key: raw[key] for key in META},
        'model_gene_alias': aliases[0], 'model_gene_reactions': ['CHRPL'],
        'thresholds': {'model_growth_greater_than_or_equal': growth_threshold, 'experiment_defect_strictly_less_than': fitness_threshold},
        'aggregation': 'Unweighted arithmetic mean of finite experiment-level exported log2 gene fitness values; grouping is carbon-source name × medium. These exports already aggregate insertion strains and are not raw barcode counts.',
        'n_eligible_metadata_groups': len(groups), 'n_mapped_matrix_conditions': len(conditions),
        'summary_by_role': by_role, 'conditions': conditions,
        'quality_boundary': 'Metrics are experiment-wide, not PP_5317-specific coverage or uncertainty. Historical Wetmore2015 checks are descriptive; abs correlations used conservatively. No current official per-experiment pass flag or local PP_5317 t table is available. Do not infer statistical equivalence to zero.',
        'gene_specific_data_availability': {
            'local_fit_t_file_exists': (folder / 'fit_t.tsv').exists(),
            'specific_phenotype_rows': [row for row in read_table(folder / 'specific_phenotypes.tsv') if row['locusId'] == LOCUS],
            'raw_barcode_counts': None, 'central_insertion_count_and_positions': None,
            'per_experiment_population_doublings': None, 'pp5317_individual_monoculture_growth': None,
            't_table_retrieval_attempt': {'url': 'https://fit.genomics.lbl.gov/cgi-bin/createFitData.cgi?orgId=Putida&t=1',
                                          'date': '2026-09-06', 'outcome': 'HTTP 403 Forbidden on retrieval; existing local data untouched'}},
        'source_medium_blocks': relevant_media,
        'mapped_medium_records': [media_map[name] for name in sorted({r['condition'].split(' | ')[1] for r in conditions})],
        'declared_nitrate_experiments_outside_carbon_benchmark': [row for row in experiments if 'nitrate' in row['condition_1'].lower()],
        'existing_model_transport_and_precursor_reactions': [reaction_record(model.reactions.get_by_id(rid)) for rid in ('CHRPL', '4HBZtex', 'UHBZ1t_pp', 'EX_4hbz_e', 'HBZOPT')],
        'saved_glucose_physical_exchange_bounds': {rid: physical['final_exchange_bounds'].get(rid) for rid in ('EX_4hbz_e', 'EX_T4hcinnm_e', 'EX_no3_e', 'EX_o2_e')},
        'assay_primary_sources': SOURCES,
        'interpretation': {
            'supported': 'PP_5317 insertion mutants do not exhibit a severe competitive deficit at the chosen -2 threshold in any of the 67 samples behind the 31 disagreement groups.',
            'not_established': 'A complete PP_5317 deletion grows indefinitely in isolated defined-medium culture, precursor secretion or cross-feeding occurs, carryover suffices, or PP_5317 is dispensable for quinone synthesis.',
            'media': 'Neither selected minimal-medium recipe includes 4HBZ or nitrate. RCH2 includes 4-aminobenzoate in Wolfe vitamins; it is a distinct metabolite.4HBZ in MOPS Rich Defined does not establish use of that medium for inoculum or carryover.',
            'diagnostic_boundary': 'Existing reversible transport permits a hypothetical 4HBZ uptake/secretion calculation without adding transport. Feasibility would not establish actual secretion, abundance, community transfer or mutant viability.'},
        'provenance': {'inputs': inputs, 'python': sys.version, 'platform': platform.platform(),
                       'numpy_version': np.__version__, 'cobra_version': cobra.__version__},
    }
    assert all(sha(ROOT / item['path']) == item['sha256'] for item in inputs)
    return result


def render(result):
    disagreement = result['summary_by_role']['disagreement']
    rescue = [r for r in result['conditions'] if r['role'] == 'natural_rescue_agreement']
    lines = [
        '# PP_5317 fitness measurements and experimental context', '',
        'This is a post-outcome audit of exposed Putida development data. It adds no model reactions, runs no optimization and inspects no other organism’s numeric phenotypes. The companion JSON preserves all 43 condition groups, 92 experiment-level scores, original metadata and input hashes.', '',
        '**The 31 remaining disagreements are not caused by averaging away a strongly defective replicate.** Their 67 exported PP_5317 scores range from −1.053 to +0.553; every score is above the fixed −2 defect threshold. Condition means range from −0.842667 to +0.3625. These data establish an absence of a severe pooled competitive deficit at this threshold. They do not establish monoculture viability of a complete deletion.', '',
        '## What the assay measures', '',
        'RB-TnSeq compares mutant barcode abundance before and after competitive growth in a pooled library. Gene scores combine insertion strains and are normalized around the typical gene. A score near zero is therefore relative performance within the pool, rather than an absolute growth rate or direct viability measurement. Wetmore’s method paper describes typically 4–6 generations. [Wetmore et al., 2015](https://journals.asm.org/doi/10.1128/mbio.00306-15).', '',
        'Price’s later methods describe typically 4–8 population doublings to saturation, central 10–90% insertions and abundance filtering; replicates need not have identical concentrations. Exact generations for these Putida samples are absent from the local metadata. The benchmark’s growth threshold 0.001 and fitness threshold −2 compare different quantities and do not calibrate one into the other. [Price et al., 2018](https://doi.org/10.1038/s41586-018-0124-0).', '',
        'For the Putida fatty-acid/alcohol work, JBEI-1 was recovered in LB+kanamycin to OD600 0.5, sampled for time zero, washed **once** in carbon-free MOPS and diluted 1:50 into 10 mM carbon source; cultures used 10 mL tubes at 30 °C and 200 rpm. The separate individual-growth protocol used **two** washes and 1:100. Those assays must not be conflated. This is consistent with set15/set16 metadata, but is not a unique expName linkage. [Thompson et al., 2020](https://journals.asm.org/doi/10.1128/aem.01665-20).', '',
        'The aromatic-substrate paper instead describes JBEI-1 LB/kanamycin recovery, a MOPS wash and 1:50 transfer into 10 mM substrates, 600 µL deep wells at 30 °C/700 rpm, and combining two 600 µL samples before BarSeq. This fits set12’s substrate/format pattern; no explicit expName link, exact wash count, endpoint duration or generation count is supplied. Engineered production strains in that paper belong to separate assays. [Incha et al., 2020, §4.5](https://doi.org/10.1016/j.mec.2019.e00119).', '',
        'A primary GEM/fitness comparison identifies cofactor carryover and cross-feeding as possible explanations for discrepancies of this kind, and warns that unsupported biosynthetic additions can create errors elsewhere. That work concerns a different organism; it supplies a methodological caution, **not evidence that PP_5317 is cross-fed in KT2440**. [Bernstein et al., 2023](https://link.springer.com/article/10.15252/msb.202311566).', '',
        '## Mapping, replication and uncertainty', '',
        'The raw locus and sysName are both PP_5317, described as “Probable chorismate pyruvate-lyase”. The frozen model uniquely maps it to NP_747418_1, assigned only to CHRPL. The saved benchmark axis and every condition mean match an independent reconstruction from the raw export to within 1e-12. These exported gene-level scores already combine insertion strains; they are not raw barcode counts.', '',
        'The 43 mapped groups contain 34 model-WT-growing conditions. PP_5317 deletion blocks 31 and grows in the three conditions below. Nine groups with model WT growth below 0.001 are outside this gene-level comparison; their 19 scores are retained in JSON for coverage accounting. Fourteen additional eligible metadata groups are unmapped.', '',
        'The 67 disagreement samples include 52 Putida_ML5_JBEI and 15 Putida_ML5 labels; 20 include DMSO. Every sample is annotated aerobic and 30 °C. Some averages combine library labels, DMSO presence, dates, vessels or concentrations. For example, RCH2 glucose combines 20/40 mM and tube/plate experiments; RCH2 acetate combines 5/20 mM. MOPS benzoate mixes 5/10 mM and both library labels. Hence “replicate average” is broader than a matched repeat under every experimental detail. This heterogeneity does not explain away the −2 classification because every contributing score remains above it.', '',
        'All 67 disagreement samples pass a descriptive comparison to the historical experiment-wide quality criteria. Their median-gene counts range 68–368, half-gene correlation 0.1083–0.2979 and half-gene median absolute difference 0.1490–0.3351. These metrics do not give PP_5317-specific read depth or statistical uncertainty. Local t-scores, insertion positions/counts and raw barcode counts are unavailable; the public t-table endpoint returned HTTP 403 on 6 September 2026. The specific-phenotype export has no PP_5317 row, which is not evidence that its fitness is statistically equivalent to zero.', '',
        '## Three natural model rescues', '',
        '| Condition | Exported PP_5317 fitness values | Mean | Saved model knockout growth |',
        '|---|---|---:|---:|',
    ]
    for row in rescue:
        label = row['condition'].replace(' | ', ' / ')
        values = ', '.join(f"{sample['gene_log2_fitness']:+.3f}" for sample in row['samples'])
        lines.append(f"| {label} | {values} | {row['mean_log2_fitness']:+.4f} | {row['pp5317_model_growth']:.7f} |")
    lines += ['',
        'These six measurements use DMSO and are all above −2. The existing model can take up 4-hydroxybenzoate or form it during coumarate catabolism, bypassing CHRPL. “Rescue” here describes the saved model prediction under an already supplied substrate, not a dedicated experimental PP_5317 complementation or dose-response assay.', '',
        '## Medium, uptake and nitrate confounds', '',
        'The retained FEBA source recipes for MOPS minimal and RCH2_defined_noCarbon contain **no 4-hydroxybenzoate and no nitrate**. RCH2 contains Wolfe vitamins including **4-aminobenzoate**, a chemically distinct compound. MOPS Rich Defined lists 4-hydroxybenzoate at 0.01 mM, but none of the 43 mapped groups uses that medium, and the primary recovery protocols above use LB. This audit found no evidence that Rich Defined was an inoculum or that residual LB supplied enough precursor.', '',
        'The only nitrate-labelled Putida rows are set27IT016/set27IT017, classified as nitrogen-source experiments in MOPS minimal media_Glucose_noNitrogen with 10 mM sodium nitrate. They are outside this carbon benchmark. A reference to the soil isolate RCH2 or to its nitrate-reduction studies must not be imported as a condition of the KT2440 assay. All selected metadata is aerobic; dissolved oxygen histories are not available.', '',
        'Existing 4HBZtex and UHBZ1t_pp are reversible in the saved model. Its glucose physical context closes 4HBZ and nitrate uptake. The carbon-condition mapping separately opens 4HBZ or coumarate uptake at −10 mmol/gDW/h when supplied. That model bound is an assumption and is not the same physical quantity as an experimental concentration in mM. A hypothetical minimal-uptake or donor-secretion calculation could test model sufficiency without adding a transporter; it could not demonstrate measured secretion, cross-feeding or carryover.', '',
        '## Evidence that would distinguish the explanations', '',
        'Obtain PP_5317 barcode-level insertion coverage and t-statistics, then compare a sequence-verified clean deletion with its complemented strain in separately cultured, washed, defined medium across serial transfers. Test 4HBZ addition and, separately, cell-free donor-conditioned medium while measuring 4HBZ and quinone pools. Quantify inoculum and endpoint population doublings. These controls would distinguish incomplete disruption, transient stores, extracellular precursor supply and an unrepresented biosynthetic route; none is established by the current scores.', '',
        'The 31 disagreements remain valid benchmark errors under its fixed scoring rule. They should not be removed from the report, and they do not by themselves justify changing CHRPL’s GPR, adding a precursor source, or claiming a new biological pathway.', '',
        f"Reproduce with `python scripts/extract_quinone_precursor_fitness_context.py --out <fresh-directory>`. The JSON records {len(result['provenance']['inputs'])} input hashes, including the extraction script. No old inputs are overwritten.", '']
    return '\n'.join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, default=DEFAULT_STUDY)
    parser.add_argument('--out', type=Path, default=Path('results/quinone_precursor_2026_09_06/evidence'))
    args = parser.parse_args(argv)
    study = args.study if args.study.is_absolute() else ROOT / args.study
    out = args.out if args.out.is_absolute() else ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    for name in ('fitness_context.json', 'fitness_context.md'):
        if (out / name).exists():
            raise FileExistsError(f'Refusing to overwrite {out / name}')
    result = extract(study)
    (out / 'fitness_context.json').write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    (out / 'fitness_context.md').write_text(render(result))
    print(json.dumps(result['summary_by_role'], indent=2))


if __name__ == '__main__':
    main()
