# IUPAC Trainer

A web app that quizzes you on IUPAC nomenclature. It generates random organic molecules, draws them, and you have to name them.

I made this to revise for my a level chemistry, as I didn't want to have to make my own flashcards, I wanted infinite generation. 

---

## What it does 

- generates molecules procedurally using RDKit (alkyl chains, benzene, esters, amides, and acid anhydrides as the "bases" with additional functional groups added ontop)
- Fetches the IUPAC name from PubChem, falling back to Cactus if PubChem doesn't have it (want to add my own procedural naming in future)
- Generates hints using regex pattern matching, to identify the roots of the answer, and what the user meant. i.e propan-1,3-ol vs propan-1,3-diol will generate the hint "how many alcohol groups are there?" (Because the user forgot the di-)
- Three difficulty levels, easy, medium and hard - harder levels use more functional groups and longer chains.
- Pre-saves molecules so there is no wait on website load

## Images

![screenshot](https://github.com/user-attachments/assets/19b3c0df-7779-4341-adb1-92d7357af26e)
![screenshot](https://github.com/user-attachments/assets/8ee1fe50-fcf3-4605-814d-c9deb2de49a0)
![screenshot](https://github.com/user-attachments/assets/01d21427-c2ea-45ee-90d1-01e2d9bb4990)

--- 
# Stack

- Python / Flask (backend)
- RDKit (molecule generation + rendering)
- PubChem and Cactus API (IUPAC name lookup)
- Vanilla JS + HTML (vibecoded frontend 😔)

---
# Setup 
Clone and run:
```bash
git clone https://github.com/toby-brett/iupac.git
cd iupac
pip install -r requirements.txt
cd project
python app.py
```
Then open `http://localhost:5000`

## Known issues
- Name lookup depends on PubChem availability
- Some valid molecules PubChem doesn't recognise and Cactus can't name - need to implement smiles -> IUPAC name myself

## To add
- Daily challenge + leaderboard
- Get a local naming pipeline (ML or alogrythmic?)
- Stereochemistry questions
- More molecules (very limited right now)
