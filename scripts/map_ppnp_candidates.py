"""Extract source-only PpnP reaction candidates and explicit model-design limits.

No model intervention, optimization, or numeric fitness access. The recovered
CarveMe universe must match the earlier gap-fill's recorded file hash. Source
reversibility is preserved as evidence, while the primary hypothesis is forward
only because the primary publication's reverse-activity statements conflict.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import cobra
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.map_quinone_pathway import met_record, reaction_record  # noqa: E402

BASE = 'results/quinone_repair_2026_09_06/runs/main/curated_path_template/model.xml.gz'
MAP = 'results/quinone_repair_2026_09_06/runs/main/curated_path_template/gene_map.json'
SOURCE = 'results/ppnp_repair_2026_09_06/evidence/source/universe_bacteria.xml.gz'
CURATED = 'models/bigg/iJN1463.xml'
GENES = 'data/fitness_browser/Putida/genes.tsv'
FEATURES = 'data/ncbi_feature_tables/GCF_000007565.2_Putida_feature_table.txt.gz'
SOURCE_IDS = ['PUNP1', 'PUNP5']
EXPECTED = {
    'PUNP1': {'adn_c': -1.0, 'pi_c': -1.0, 'ade_c': 1.0, 'r1p_c': 1.0},
    'PUNP5': {'ins_c': -1.0, 'pi_c': -1.0, 'hxan_c': 1.0, 'r1p_c': 1.0},
}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def stoichiometry(reaction):
    return {met.id: float(value) for met, value in reaction.metabolites.items()}


def adjacency(model, mids):
    return {mid: [{'reaction': reaction.id, 'coefficient': float(reaction.metabolites[model.metabolites.get_by_id(mid)]),
                   'bounds': list(reaction.bounds)}
                  for reaction in sorted(model.metabolites.get_by_id(mid).reactions, key=lambda r: r.id)] for mid in mids}


def signed_sum(model, fluxes):
    result = defaultdict(Fraction)
    for rid, flux in fluxes.items():
        reaction = model.reactions.get_by_id(rid)
        assert reaction.lower_bound <= flux <= reaction.upper_bound, rid
        for met, value in reaction.metabolites.items():
            result[met.id] += Fraction.from_float(float(value)) * Fraction(flux)
    return {mid: float(value) for mid, value in sorted(result.items()) if value}


def mapped_record(reaction, mapping):
    return {**reaction_record(reaction), 'model_to_locus': {gene.id: mapping.get(gene.id) for gene in reaction.genes}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'results/ppnp_repair_2026_09_06/evidence/reaction_candidates.json')
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('Preserve candidate evidence; choose a new --out path')
    input_paths = [BASE, MAP, SOURCE, CURATED, GENES, FEATURES, 'models/gapfilled/Putida_gapfill.json',
                   'data/genpept/Putida_genpept_map.tsv', 'scripts/map_quinone_pathway.py',
                   'scripts/map_ppnp_candidates.py', 'docs/evidence/10-putida-ribosyl-disposal-sources.md',
                   'results/ribosyl_disposal_2026_09_06/evidence/source_evidence.json',
                   'results/quinone_repair_2026_09_06/evidence/postrun_dependency_trace.json']
    hashes = {path: sha(ROOT / path) for path in input_paths}
    gapfill = json.loads((ROOT / 'models/gapfilled/Putida_gapfill.json').read_text())
    assert hashes[SOURCE] == gapfill['universe_sha256'], 'Recovered universe differs from earlier frozen source'
    source = cobra.io.read_sbml_model(str(ROOT / SOURCE))
    target = cobra.io.read_sbml_model(str(ROOT / BASE))
    curated = cobra.io.read_sbml_model(str(ROOT / CURATED))
    mapping = json.loads((ROOT / MAP).read_text())['model_to_browser']
    genes = pd.read_table(ROOT / GENES, dtype=str, keep_default_na=False)
    features = pd.read_table(ROOT / FEATURES, dtype=str, keep_default_na=False)
    records, all_mets, duplicates, equivalents = [], {}, {}, {}
    for rid in SOURCE_IDS:
        reaction = source.reactions.get_by_id(rid)
        assert stoichiometry(reaction) == EXPECTED[rid]
        assert not reaction.check_mass_balance() and not reaction.boundary
        assert reaction.bounds == (-1000.0, 1000.0) and not reaction.gene_reaction_rule
        assert rid not in target.reactions
        records.append(reaction_record(reaction))
        all_mets.update({met.id: met_record(met) for met in reaction.metabolites})
        duplicates[rid] = [reaction_record(other) for other in source.reactions
                           if other.id != rid and stoichiometry(other) == EXPECTED[rid]]
        equivalents[rid] = {label: [reaction_record(other) for other in model.reactions
                                    if stoichiometry(other) == EXPECTED[rid] or
                                    stoichiometry(other) == {mid: -v for mid, v in EXPECTED[rid].items()}]
                            for label, model in [('saved_target', target), ('local_iJN1463', curated)]}
        assert not equivalents[rid]['saved_target'] and not equivalents[rid]['local_iJN1463']
    expected_target = {mid: met_record(target.metabolites.get_by_id(mid)) for mid in all_mets}
    differences = {}
    for mid, original in all_mets.items():
        existing = expected_target[mid]
        assert existing['formula'] == original['formula'] and existing['compartment'] == original['compartment']
        changed = {field: {'source': original[field], 'target': existing[field]} for field in original if original[field] != existing[field]}
        if changed:
            differences[mid] = changed
    assert {mid: {'source': record['charge'], 'target': expected_target[mid]['charge']} for mid, record in all_mets.items()
            if record['charge'] != expected_target[mid]['charge']} == {'pi_c': {'source': -2, 'target': 0}, 'r1p_c': {'source': -2, 'target': 0}}
    local_hits = sorted(gid for gid, locus in mapping.items() if locus == 'PP_4248')
    assert not local_hits and 'PP_4248' not in target.genes
    current_gene = features[(features['# feature'] == 'gene') & features.attributes.str.contains('old_locus_tag=PP_4248', regex=False)]
    assert len(current_gene) == 1
    current_tag = current_gene.iloc[0].locus_tag
    current_cds = features[(features['# feature'] == 'CDS') & (features.locus_tag == current_tag)]
    assert len(current_cds) == 1
    current_id = current_cds.iloc[0].product_accession.replace('.', '_')
    assert current_id not in target.genes
    key_loci = ['PP_4248', 'PP_2458', 'PP_1777', 'PP_5288', 'PP_4976', 'PP_0591', 'PP_3254', 'PP_2460', 'PP_4218']
    locus_records = {}
    for locus in key_loci:
        hits = sorted(gid for gid, name in mapping.items() if name == locus)
        assert len(hits) <= 1, f'Ambiguous existing model identity for {locus}'
        locus_records[locus] = {'existing_model_gene_ids': hits,
                                'annotation_records': genes[genes.sysName == locus].to_dict('records'),
                                'associated_reactions': {gid: sorted(r.id for r in target.genes.get_by_id(gid).reactions) for gid in hits}}
    bypass = {'PNP': -1, 'RNMK': 1, 'NMNN': 1}
    bypass_net = signed_sum(target, bypass)
    assert bypass_net == {'adp_c': 1.0, 'atp_c': -1.0, 'h2o_c': -1.0, 'h_c': 1.0, 'pi_c': 1.0, 'r1p_c': -1.0, 'r5p_c': 1.0}
    relevant = ['AHCi', 'AHCYSNS', 'RHCYS', 'ADA', 'RBK', 'PPM', 'PGMT', 'PNP', 'PYNP2r', 'NP1', 'RNMK', 'NMNN',
                'NICRNS', 'NNATr', 'NADS1', 'NTD2', 'NTD7', 'NTD11', 'INSH', 'RIBabcpp', 'RIBtex', 'EX_rib__D_e']
    proposal = [{'id': rid, 'source_id': rid, 'substrate': substrate, 'metabolites': EXPECTED[rid]}
                for rid, substrate in zip(SOURCE_IDS, ['adenosine', 'inosine'])]
    artifact = {
        'schema_version': 1, 'candidate_id': 'ppnp_phosphorolysis_v1',
        'status': 'provisional_source_only_not_applied',
        'source_model': {'path': SOURCE, 'sha256': hashes[SOURCE], 'model_id': source.id},
        'source_origin': {'url': 'https://raw.githubusercontent.com/cdanielmachado/carveme/master/carveme/data/generated/universe_bacteria.xml.gz',
                          'repository_page': 'https://github.com/cdanielmachado/carveme/blob/master/carveme/data/generated/universe_bacteria.xml.gz',
                          'recovery_note': 'Original local external universe was absent. Downloaded bytes exactly match the pre-existing Putida gapfill universe_sha256; no changed-source substitution.',
                          'original_record': 'models/gapfilled/Putida_gapfill.json',
                          'scope': 'Universal reaction chemistry, not KT2440-specific gene or activity evidence.'},
        'baseline_model': {'path': BASE, 'sha256': hashes[BASE], 'model_id': target.id,
                           'quinone_biomass_coefficient': target.reactions.Growth.get_coefficient('q8h2_c')},
        'reaction_records': records, 'all_involved_metabolites': all_mets,
        'target_metabolite_map': {mid: mid for mid in all_mets},
        'expected_target_metabolites': expected_target, 'proposed_reactions': proposal,
        'gene_aliases': {'PP_4248': 'PP_4248'},
        'proposed_gene_rule': 'PP_4248',
        'proposed_direction_policy': {'primary_forward_only_bounds': [0.0, 1000.0],
                                      'separate_reversibility_sensitivity_bounds': [-1000.0, 1000.0],
                                      'basis': 'Original Sevin Fig S8 and Table 9 agree on forward adenosine/inosine phosphorolysis but conflict on reverse synthesis. Source reversibility is not sufficient evidence to transfer it into the primary arm.',
                                      'primary_source_doi': '10.1038/nmeth.4103',
                                      'supplement_url': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.4103/MediaObjects/41592_2017_BFnmeth4103_MOESM62_ESM.pdf'},
        'metadata_policy': {'preserve_existing_target_metadata': True, 'new_metabolites': [],
                            'compartment_aliases_required': {},
                            'expected_charge_differences': {'pi_c': {'source': -2, 'target': 0}, 'r1p_c': {'source': -2, 'target': 0}},
                            'all_record_differences': differences,
                            'caveat': 'Source reactions are atom/charge balanced. Target zero charges are inherited placeholders, and their cancellation does not validate charge metadata.'},
        'equivalent_reaction_audit': equivalents,
        'source_duplicates_not_selected': duplicates,
        'canonical_id_selection': 'PUNP1/PUNP5 carry database identifiers; exact _1 duplicates add no chemistry and are excluded before outcomes.',
        'candidate_identity': {'locus': 'PP_4248', 'existing_model_gene_ids': local_hits,
                               'target_model_gene_id': 'PP_4248', 'current_refseq_model_style_id_absent': current_id,
                               'current_gene_records': current_gene.to_dict('records'), 'current_cds_records': current_cds.to_dict('records'),
                               'caveat': 'These prove the unambiguous local alias choice, not enzyme function or sequence equality to Q88F51; the independent sequence audit is separate.'},
        'control_gene_records': locus_records,
        'existing_relevant_reactions': [mapped_record(target.reactions.get_by_id(rid), mapping) for rid in relevant],
        'relevant_complete_adjacency': adjacency(target, ['ahcys_c', 'rhcys_c', 'adn_c', 'ins_c', 'rib__D_c', 'r1p_c', 'r5p_c', 'rnam_c', 'nmn_c', 'nicrns_c', 'uri_c']),
        'ppm_assignment_discrepancy': {'saved_target': reaction_record(target.reactions.PPM), 'local_iJN1463': reaction_record(curated.reactions.PPM),
                                     'interpretation': 'Target PPM is PP_1777 only, curated PPM is PP_1777 or PP_5288; no GPR change is proposed. PP_1777 local annotation says phosphomannomutase. Functional ribose-P mutase assignment remains provisional.'},
        'ppm_bypass': {'signed_reaction_combination': bypass, 'summed_stoichiometry': bypass_net,
                       'interpretation': 'Permitted reaction directions yield an ATP-consuming r1p-to-r5p path. This symbolic route does not prove an optimized feasible state, but PPM deletion cannot be assumed to eliminate ribosyl disposal.'},
        'phosphate_control_limit': 'Closing external phosphate limits general biomass synthesis and is not a specific PpnP control. The exact phosphorolysis coefficient and a separately declared substrate-only enzyme module/assay with and without phosphate discriminate specificity; do not remove phosphate from the reaction equation.',
        'design_review': 'results/ppnp_repair_2026_09_06/evidence/design_review.md',
        'provenance': {'created_utc': datetime.now(timezone.utc).isoformat(),
                       'base_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                       'input_sha256': hashes, 'numeric_fitness_values_accessed': False,
                       'model_intervention_performed': False, 'optimization_performed': False,
                       'role': 'Post-outcome source-driven hypothesis design on exposed development organism; not an independent evaluation.'},
    }
    for path, digest in hashes.items():
        assert sha(ROOT / path) == digest, f'Input changed during extraction: {path}'
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(artifact, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(args.out)


if __name__ == '__main__':
    main()
