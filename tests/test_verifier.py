import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def test_nan_card_cannot_silently_pass_independent_verification(tmp_path):
    src = next((ROOT / 'results/carbon_fitness_multi/MR1').glob('*__gapfilled__nodroprich'))
    dest = tmp_path / 'MR1' / 'corrupt'
    shutil.copytree(src, dest)
    path = dest / 'card.json'
    card = json.loads(path.read_text())
    card['results']['gene_level_conditions_where_wt_grows']['mcc']['point'] = float('nan')
    path.write_text(json.dumps(card))
    result = subprocess.run([sys.executable, str(ROOT / 'scripts/verify_carbon_fitness_multi.py'),
                             '--results-dir', str(tmp_path)], capture_output=True, text=True)
    assert result.returncode != 0
    assert 'MISMATCH' in result.stdout


def test_no_cards_is_not_success(tmp_path):
    result = subprocess.run([sys.executable, str(ROOT / 'scripts/verify_carbon_fitness_multi.py'),
                             '--results-dir', str(tmp_path)], capture_output=True, text=True)
    assert result.returncode != 0
    assert 'nothing was verified' in result.stderr
