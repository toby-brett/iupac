
import base64
import io
import json
import random
import threading
import queue
from collections import deque

from PIL.ImageOps import contain
from flask import Flask, jsonify, request, render_template
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import Draw

from molecule import easy, medium, hard, smiles_to_name
from normalize import normalize, pick_hint

app = Flask(__name__,
            static_folder="static",
            template_folder="templates")
CORS(app)

DIFFICULTY_MAP = {"easy": easy, "medium": medium, "hard": hard}
MAX_ATTEMPTS = 5
POOL_SIZE = 2  # molecules pre-baked per difficulty

# One queue per difficulty — holds ready-to-serve dicts
_pools = {
    "easy":   queue.Queue(maxsize=POOL_SIZE),
    "medium": queue.Queue(maxsize=POOL_SIZE),
    "hard":   queue.Queue(maxsize=POOL_SIZE),
}

seen = deque(maxlen=10)

def _bake_one(difficulty):
    """Load one molecule dict, return it or None on failure."""
    with open(f"molecules/{difficulty}.json", "r") as f:
        while True:
            try:
                molecules = json.load(f)
                molecule = random.choice(molecules)
                smiles = molecule["smiles"]
                name = molecule["name"]
                if smiles in seen:
                    continue
                seen.append(smiles)
                if name is None:
                    continue
                mol = Chem.MolFromSmiles(smiles)
                img = Draw.MolToImage(mol, size=(400, 300))
                buf = io.BytesIO()
                img.save(buf, format="PNG")
                img_b64 = base64.b64encode(buf.getvalue()).decode()
                return {"image": img_b64, "name": name, "smiles": smiles}
            except Exception as e:
                print(f"[{difficulty}] error: {e}")
                break
    return None

import time # Add this at the top

def _worker(difficulty):
    q = _pools[difficulty]
    while True:
        item = _bake_one(difficulty)
        if item is not None:
            q.put(item)
            time.sleep(1) # Wait 1 second between successful generations
        else:
            time.sleep(5) # If it fails, wait 5 seconds before trying again

# Start one worker thread per difficulty
for diff in DIFFICULTY_MAP:
    t = threading.Thread(target=_worker, args=(diff,), daemon=True)
    t.start()

@app.route("/api/molecule", methods=["GET"])
def get_molecule():
    difficulty = request.args.get("difficulty", "medium").lower()
    if difficulty not in DIFFICULTY_MAP:
        return jsonify({"error": "difficulty must be easy, medium, or hard"}), 400

    q = _pools[difficulty]
    try:
        # Try to grab a pre-baked one instantly
        item = q.get(timeout=30)
        return jsonify(item)
    except queue.Empty:
        return jsonify({"error": "Timed out generating molecule. Try again."}), 500

@app.route("/api/try1", methods=["POST"])
def help_message():
    body = request.get_json(force=True)
    answer = body.get("answer", "")
    correct = body.get("correct", "")

    norm_answer = normalize(answer)
    norm_correct = normalize(correct)

   # message = get_message(norm_answer, norm_correct)
    message = pick_hint(norm_answer, norm_correct)
    return jsonify({
        "correct": norm_answer == norm_correct,
        "normalized_answer": norm_answer,
        "normalized_correct": norm_correct,
        "message": message
    })

@app.route("/")
def index():
    return render_template("index_real.html")

if __name__ == "__main__":
    app.run(debug=True, port=5000)