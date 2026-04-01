import json
import time
from pathlib import Path

from rdkit import Chem

from project.molecule import easy, medium, hard, smiles_to_name

DIFFICULTY_MAP = {"easy": easy, "medium": medium, "hard": hard}
POOL_SIZE = 2  # molecules pre-baked per difficulty

def load_data(path):
    try:
        with open(path, "r") as f:
            content = f.read().strip()
            if not content:
                return []
            return json.loads(content)
    except (FileNotFoundError, json.JSONDecodeError):
        return []

def append_to_json(path, new_item):
    data = load_data(path)
    data.append(new_item)

    with open(path, "w") as f:
        json.dump(data, f, indent=2)

def _bake_one(difficulty):
    """Generate one molecule dict, return it or None on failure."""
    BASE_DIR = Path(__file__).resolve().parent
    root = BASE_DIR / f"molecules/{difficulty}.json"
    generator = DIFFICULTY_MAP[difficulty]
    current = load_data(root)
    seen = [m["smiles"] for m in current]
    for attempt in range(POOL_SIZE):
        try:
            mol = generator()
            smiles = Chem.MolToSmiles(mol)
            if smiles in seen:
                continue
            seen.append(smiles)
            name = smiles_to_name(smiles)
            if name is None:
                continue
            data = {"name": name, "smiles": smiles}
            append_to_json(root, data)
        except Exception as e:
            print(f"  [{difficulty}] attempt {attempt+1} error: {e}")
    return None

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--difficulty")

    args = parser.parse_args()
    difficulty = args.difficulty

    while True:
        _bake_one(difficulty)
        time.sleep(0.5)