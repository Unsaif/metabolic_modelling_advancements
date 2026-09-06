"""Prespecified Putida biomass sensitivity; uses only the existing development data."""
from __future__ import annotations
import argparse
from dataclasses import asdict
import gzip
import json
from pathlib import Path
import sys
import time
import numpy as np
import cobra
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gembench.cards import BenchmarkCard, ModelProvenance, carbon_fitness_leakage, now, sha256_of
from gembench.comparison import compare_runs, load_run
from gembench.fitness_browser import load_organism, carbon_source_conditions
from gembench.patches import apply_universe_patches, apply_model_patches, apply_gpr_patches, load_patch_files
from gembench.protocols import carbon_fitness_generic as P
from scripts.run_carbon_fitness_generic import load_model, gene_map_for, score, per_condition
from scripts.audit_saved_results import json_safe

PLAN = ROOT / 'data/studies/quinone_biomass_sensitivity_v1.json'


def prepare(fb):
    model, path, source = load_model('Putida', 'gapfilled')
    apply_universe_patches(model, 'Putida', json.loads((ROOT/'data/reference/universe_patches_v0.1.json').read_text()))
    gm = gene_map_for('Putida', model, set(fb.genes['sysName']), 'gapfilled')
    apply_model_patches(model, 'Putida', json.loads((ROOT/'data/reference/model_patches_v0.4.json').read_text()), gm, verbose=False)
    gm = gene_map_for('Putida', model, set(fb.genes['sysName']), 'gapfilled')
    gp = load_patch_files([ROOT/f'data/reference/gpr_patches_v0.{i}.json' for i in [2,3,4]])
    apply_gpr_patches(model, 'Putida', gm, gp, verbose=False)
    gm = gene_map_for('Putida', model, set(fb.genes['sysName']), 'gapfilled')
    return model, gm, path, source


def write_json(path, value):
    path.write_text(json.dumps(json_safe(value), indent=2, allow_nan=False)+'\n')


def apply_biomass_intervention(model, config, template_coefficient, curated_coefficient):
    """Alter only the biomass substrate named in the frozen sensitivity arm."""
    biomass = model.reactions.get_by_id('Growth')
    for name in ('mql8_c', 'q8h2_c'):
        metabolite = model.metabolites.get_by_id(name)
        if biomass.metabolites.get(metabolite, 0.0) != 0:
            raise ValueError('Expected v0.4 to have no biomass quinone demand')
    choices = {'no_quinone_demand': (None, 0.0),
               'ubiquinol_template_amount': ('q8h2_c', template_coefficient),
               'ubiquinol_curated_amount': ('q8h2_c', curated_coefficient),
               'menaquinol_restored_control': ('mql8_c', template_coefficient)}
    metabolite_id, coefficient = choices[config['id']]
    if config['metabolite'] != metabolite_id:
        raise ValueError('Intervention identifier and metabolite disagree')
    if metabolite_id is not None:
        if not np.isfinite(coefficient) or coefficient >= 0:
            raise ValueError('A biomass substrate must have a finite negative coefficient')
        biomass.add_metabolites({model.metabolites.get_by_id(metabolite_id): coefficient})
    return coefficient


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=ROOT/'results/quinone_biomass_2026_09_06')
    ap.add_argument('--freeze-only', action='store_true')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    manifest_path = args.out/'prespecified_manifest.json'
    paths = [PLAN, Path(__file__), ROOT/'models/gapfilled/Putida.xml.gz', ROOT/'models/bigg/iJN1463.xml',
             ROOT/'data/reference/universe_patches_v0.1.json', ROOT/'data/reference/model_patches_v0.4.json',
             ROOT/'data/genpept/Putida_genpept_map.tsv', ROOT/'requirements-audit.txt']
    paths += [ROOT/f'data/reference/gpr_patches_v0.{i}.json' for i in [2,3,4]]
    paths += [p for p in sorted((ROOT/'gembench').rglob('*.py')) if p.name != 'study.py']
    paths += sorted((ROOT/'data/fitness_browser/Putida').glob('*'))
    paths += sorted((ROOT/'data/reference').glob('*.tsv'))
    paths += [ROOT/'scripts/run_carbon_fitness_generic.py', ROOT/'scripts/audit_saved_results.py']
    current = {str(p.relative_to(ROOT)): sha256_of(p) for p in paths if p.is_file()}
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        # Later unrelated additions to gembench may be ignored, but every frozen input must match.
        mismatched = [p for p,h in manifest['input_sha256'].items() if current.get(p) != h]
        if mismatched: raise ValueError(f'Frozen inputs changed: {mismatched}')
    else:
        manifest = {'created_before_simulation': now(), 'plan': json.loads(PLAN.read_text()), 'input_sha256': current}
        write_json(manifest_path, manifest)
    if args.freeze_only: return
    if (args.out/'summary.json').exists(): raise FileExistsError('Study already completed; use a fresh output path')
    with gzip.open(ROOT/'models/gapfilled/Putida.xml.gz','rt') as stream: template=cobra.io.read_sbml_model(stream)
    curated=cobra.io.read_sbml_model(str(ROOT/'models/bigg/iJN1463.xml'))
    template_coefficient=template.reactions.get_by_id('Growth').get_coefficient('mql8_c')
    curated_coefficient=curated.reactions.get_by_id('BIOMASS_KT2440_WT3').get_coefficient('q8h2_c')
    assert template_coefficient < 0 and curated_coefficient < 0
    fb=load_organism('Putida'); conditions=carbon_source_conditions(fb)
    configs=manifest['plan']['configurations']
    records=[]
    for config in configs:
        d=args.out/config['id']
        if d.exists(): raise FileExistsError(f'{d} exists; incomplete runs require a fresh output path')
        model,gm,path,source=prepare(fb)
        coefficient=apply_biomass_intervention(model,config,template_coefficient,curated_coefficient)
        model.id += '__'+config['id']
        params=P.GenericParams(**manifest['plan']['protocol'])
        started=time.time(); result=P.run(model,fb,conditions,gm,params)
        grows=result.wt_growth >= params.growth_threshold
        results={'condition_level':{'n_conditions_mapped':len(result.conditions),'n_conditions_wt_grows':int(grows.sum())},
                 'gene_level_conditions_where_wt_grows':score(result,grows)}
        card=BenchmarkCard(benchmark='carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Putida', created=now(),
            dataset_provenance=fb.provenance,
            model=ModelProvenance(model_id=model.id,file=str(Path(path).relative_to(ROOT)),source=source,sha256=sha256_of(path)),
            protocol={'params':asdict(params),'evaluation_role':'retrospective_development',
                      'intervention':{**config,'applied_coefficient':coefficient},'input_sha256':manifest['input_sha256']},
            leakage=carbon_fitness_leakage('Putida','gapfilled',patched=True,medium_completion=True),results=results)
        d.mkdir(); card.write(str(d/'card.json'),str(d/'card.md'))
        np.savez_compressed(d/'matrices.npz',sim_growth=result.sim_growth,wt_growth=result.wt_growth,fitness=result.fitness,
            model_genes=np.array(result.model_genes),browser_genes=np.array(result.browser_genes),conditions=np.array([c.key for c in result.conditions]))
        result.condition_table().to_csv(d/'conditions.tsv',sep='\t',index=False)
        per_condition(result).to_csv(d/'per_condition_metrics.tsv',sep='\t',index=False)
        records.append({'arm':config['id'],'coefficient':coefficient,'seconds':time.time()-started,**results})
        print('Completed',config['id'],results['condition_level'],flush=True)
    comparisons=[]; changes=[]
    base=load_run(args.out/configs[0]['id'])
    for config in configs[1:]:
        other=load_run(args.out/config['id'])
        comparisons.append({'arm':config['id'],**compare_runs(base,other)})
        # Gene identities are equal by construction; confirm before attribution.
        assert np.array_equal(base['browser_genes'],other['browser_genes'])
        assert np.array_equal(base['conditions'],other['conditions'])
        good=(base['wt_growth']>=0.001)&(other['wt_growth']>=0.001)
        diff=(base['sim_growth']>=0.001)!=(other['sim_growth']>=0.001)
        for i,j in zip(*np.where(diff & good[None,:])):
            changes.append({'arm':config['id'],'gene':base['browser_genes'][i], 'condition':base['conditions'][j],
              'before':float(base['sim_growth'][i,j]),'after':float(other['sim_growth'][i,j]),'fitness':float(base['fitness'][i,j])})
    write_json(args.out/'changed_predictions.json',changes)
    write_json(args.out/'summary.json',{'study':manifest['plan'],'runs':records,'comparisons':comparisons,
      'n_changed_gene_condition_predictions':len(changes),'interpretation':'Retrospective sensitivity analysis; coefficient values not fitted to fitness.'})


if __name__=='__main__': main()
