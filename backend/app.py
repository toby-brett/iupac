import base64
import io
import threading
import queue

from flask import Flask, jsonify, request, send_from_directory
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import Draw

from molecule import easy, medium, hard, smiles_to_name
from normalize import normalize, get_message

app = Flask(__name__)
CORS(app)

DIFFICULTY_MAP = {"easy": easy, "medium": medium, "hard": hard}
MAX_ATTEMPTS = 10
POOL_SIZE = 5  # molecules pre-baked per difficulty

# One queue per difficulty — holds ready-to-serve dicts
_pools = {
    "easy":   queue.Queue(maxsize=POOL_SIZE),
    "medium": queue.Queue(maxsize=POOL_SIZE),
    "hard":   queue.Queue(maxsize=POOL_SIZE),
}

def _bake_one(difficulty):
    """Generate one molecule dict, return it or None on failure."""
    generator = DIFFICULTY_MAP[difficulty]
    for attempt in range(MAX_ATTEMPTS):
        try:
            mol = generator()
            smiles = Chem.MolToSmiles(mol)
            name = smiles_to_name(smiles)
            if name is None:
                continue
            img = Draw.MolToImage(mol, size=(400, 300))
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            img_b64 = base64.b64encode(buf.getvalue()).decode()
            return {"image": img_b64, "name": name, "smiles": smiles}
        except Exception as e:
            print(f"  [{difficulty}] attempt {attempt+1} error: {e}")
    return None

def _worker(difficulty):
    """Background thread: keep the pool full."""
    q = _pools[difficulty]
    while True:
        # Block until there's space in the queue
        # (put blocks when full, so this naturally throttles)
        item = _bake_one(difficulty)
        if item is not None:
            q.put(item)   # blocks if queue is full — that's fine

# Start one worker thread per difficulty
for diff in DIFFICULTY_MAP:
    t = threading.Thread(target=_worker, args=(diff,), daemon=True)
    t.start()

@app.route("/molecule", methods=["GET"])
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

@app.route("/try1", methods=["POST"])
def help_message():
    body = request.get_json(force=True)
    answer = body.get("answer", "")
    correct = body.get("correct", "")

    norm_answer = normalize(answer)
    norm_correct = normalize(correct)

    message = get_message(norm_answer, norm_correct)
    print(norm_answer, norm_correct)
    return jsonify({
        "correct": norm_answer == norm_correct,
        "normalized_answer": norm_answer,
        "normalized_correct": norm_correct,
        "message": message
    })

@app.route("/")
def index():
    return send_from_directory("..", "index.html")

if __name__ == "__main__":
    app.run(debug=True, port=5000)