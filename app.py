import base64
import io

from flask import Flask, jsonify, request
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import Draw

from molecule import easy, medium, hard, smiles_to_name
from normalize import normalize, get_message

app = Flask(__name__)
CORS(app)

# api endpoints
DIFFICULTY_MAP = {"easy": easy, "medium": medium, "hard": hard}
MAX_ATTEMPTS = 10

@app.route("/molecule", methods=["GET"])
def get_molecule():
    difficulty = request.args.get("difficulty", "medium").lower()
    if difficulty not in DIFFICULTY_MAP:
        return jsonify({"error": "difficulty must be easy, medium, or hard"}), 400

    generator = DIFFICULTY_MAP[difficulty]

    for attempt in range(MAX_ATTEMPTS):
        try:
            print(f"Attempt {attempt+1}...")
            mol = generator()
            smiles = Chem.MolToSmiles(mol)
            print(f"  SMILES: {smiles}")
            name = smiles_to_name(smiles)
            print(f"  Name: {name}")
            if name is None:
                continue

            img = Draw.MolToImage(mol, size=(400, 300))
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            img_b64 = base64.b64encode(buf.getvalue()).decode()

            return jsonify({"image": img_b64, "name": name, "smiles": smiles})
        except Exception as e:
            print(f"  ERROR: {e}")
            continue

    return jsonify({"error": "Could not generate a valid molecule. Try again."}), 500

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

from flask import send_from_directory

@app.route("/")
def index():
    return send_from_directory(".", "index.html")

if __name__ == "__main__":
    app.run(debug=True, port=5000)