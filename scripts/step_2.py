
'''
import pandas as pd
import json
import os
import time
from tqdm import tqdm

# === PATH BASE ===
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(BASE_DIR)

DATA_DIR = os.path.join(ROOT_DIR, "data")
OUTPUT_MAIN = os.path.join(ROOT_DIR, "output")
OUTPUT_DIR = os.path.join(OUTPUT_MAIN, "tpm_normalizzati")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# === FILE INPUT ===
TPM_PATH = os.path.join(DATA_DIR, "tpm_matrix_completo_batch_parallel.csv")
SAMPLES_PATH = os.path.join(DATA_DIR, "samples_raggruppati.json")
# === Timer globale ===
total_start = time.time()

print(">> Caricamento matrice TPM completa...")
tpm = pd.read_csv(TPM_PATH, index_col=0)
print(f">> TPM caricata: {tpm.shape[0]} geni × {tpm.shape[1]} sample\n")

print(">> Caricamento file samples_raggruppati.json...")
with open(SAMPLES_PATH, "r") as f:
    samples_by_tissue_age = json.load(f)
print(f">> Trovati {len(samples_by_tissue_age)} tessuti nel file\n")

# === Creazione cartella di output ===
os.makedirs(OUTPUT_DIR, exist_ok=True)

# === Elaborazione per tessuto ===
for tissue, age_groups in tqdm(samples_by_tissue_age.items(), desc="Tessuti"):
    start_time = time.time()
    print(f"\n>>> Inizio normalizzazione per tessuto: {tissue}")

    # Raccogli tutti i sample ID validi
    sample_ids = []
    for age_group in age_groups:
        sample_ids.extend(age_groups[age_group])
    sample_ids = [sid for sid in sample_ids if sid in tpm.columns]

    if not sample_ids:
        print(f"[!] Nessun sample valido per il tessuto: {tissue}")
        continue

    print(f">> Totale sample validi: {len(sample_ids)}")

    # Estrai submatrice TPM per il tessuto
    sub_tpm = tpm[sample_ids].copy()

    # Normalizzazione min-max per gene
    min_vals = sub_tpm.min(axis=1)
    max_vals = sub_tpm.max(axis=1)
    ranges = max_vals - min_vals
    ranges[ranges == 0] = 1
    norm_tpm = (sub_tpm.sub(min_vals, axis=0)).div(ranges, axis=0)

    print(">> Normalizzazione completata")

    # Mapping sample → fascia d'età
    sample_to_age = {}
    for age_group, ids in age_groups.items():
        for sid in ids:
            if sid in norm_tpm.columns:
                sample_to_age[sid] = age_group

    df_norm = norm_tpm.T  # sample × gene
    df_norm["age_group"] = df_norm.index.map(sample_to_age)
    df_norm = df_norm.dropna(subset=["age_group"])

    print(">> Raggruppamento per fascia d'età...")
    df_avg = df_norm.groupby("age_group").mean().T  # gene × età

    output_path = os.path.join(OUTPUT_DIR, f"{tissue}_normalized.csv")
    df_avg.to_csv(output_path)
    elapsed = time.time() - start_time
    print(f">>> Salvato: {output_path}")
    print(f">>> Tempo impiegato per {tissue}: {elapsed:.2f} secondi")

# === Riepilogo finale ===
total_elapsed = time.time() - total_start
print(f"\n=== FINE: normalizzazione completata per tutti i tessuti ===")
print(f"Tempo totale: {total_elapsed / 60:.2f} minuti")
'''


'''
# ================================
#   STEP 2 - NORMALIZZAZIONE TPM
# ================================

import pandas as pd
import json
import os
import time
from tqdm import tqdm

# === PATH BASE ===
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(BASE_DIR)

DATA_DIR = os.path.join(ROOT_DIR, "data")
OUTPUT_DIR = os.path.join(ROOT_DIR, "output", "tpm_normalizzati")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# === INPUT DATA ===
TPM_PATH = os.path.join(DATA_DIR, "tpm_matrix_completo_batch_parallel.csv")
SAMPLES_PATH = os.path.join(DATA_DIR, "samples_raggruppati.json")

# === TIMER ===
global_start = time.time()

print(">> Caricamento matrice TPM...")
tpm = pd.read_csv(TPM_PATH, index_col=0)
print(f">> TPM caricata: {tpm.shape[0]} geni × {tpm.shape[1]} sample\n")

print(">> Caricamento samples_raggruppati.json...")
with open(SAMPLES_PATH, "r") as f:
    samples_by_tissue_age = json.load(f)
print(f">> Trovati {len(samples_by_tissue_age)} tessuti\n")

# === ELABORAZIONE TESSUTI ===
for tissue, age_groups in tqdm(samples_by_tissue_age.items(), desc="Tessuti"):
    start_tissue = time.time()
    print(f"\n>>> Normalizzazione tessuto: {tissue}")

    sample_ids = []
    for age_group, ids in age_groups.items():
        sample_ids.extend(ids)

    sample_ids = [sid for sid in sample_ids if sid in tpm.columns]

    if not sample_ids:
        print(f"[!] Nessun sample valido per {tissue}")
        continue

    print(f">> Sample validi: {len(sample_ids)}")

    sub_tpm = tpm[sample_ids].copy()

    min_vals = sub_tpm.min(axis=1)
    max_vals = sub_tpm.max(axis=1)
    ranges = max_vals - min_vals
    ranges[ranges == 0] = 1

    norm = (sub_tpm.sub(min_vals, axis=0)).div(ranges, axis=0)

    sample_to_age = {}
    for age_group, ids in age_groups.items():
        for sid in ids:
            if sid in norm.columns:
                sample_to_age[sid] = age_group

    df_norm = norm.T
    df_norm["age_group"] = df_norm.index.map(sample_to_age)
    df_norm = df_norm.dropna(subset=["age_group"])

    df_avg = df_norm.groupby("age_group").mean().T

    out_path = os.path.join(OUTPUT_DIR, f"{tissue}_normalized.csv")
    df_avg.to_csv(out_path)

    print(f">>> Salvato: {out_path}")
    print(f">>> Tempo per {tissue}: {time.time()-start_tissue:.2f}s")

# === RIEPILOGO ===
total = time.time() - global_start
print("\n=== STEP 2 COMPLETATO ===")
print(f"Tempo totale: {total/60:.2f} minuti")
'''
# ================================
#   STEP 2 - NORMALIZZAZIONE TPM
#   (versione gui_leggero)
#   Supporta STAMP_TISSUES per normalizzare solo i tessuti richiesti.
# ================================

import pandas as pd
import json
import os
import time

# === PATH BASE ===
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(BASE_DIR)

DATA_DIR = os.path.join(ROOT_DIR, "data")
OUTPUT_DIR = os.path.join(ROOT_DIR, "output", "tpm_normalizzati")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# === INPUT DATA ===
TPM_PATH = os.path.join(DATA_DIR, "tpm_matrix_completo_batch_parallel.csv")
SAMPLES_PATH = os.path.join(DATA_DIR, "samples_raggruppati.json")

# === FILTRO TESSUTI (opzionale via env var) ===
# Se STAMP_TISSUES è impostata, processa solo i tessuti elencati (separati da virgola).
# Altrimenti processa tutti.
tissues_env = os.environ.get("STAMP_TISSUES", "")
tissues_filter = set(tissues_env.split(",")) if tissues_env.strip() else None

if tissues_filter:
    print(f">> Modalità selettiva: verrà normalizzato solo: {', '.join(tissues_filter)}\n")
else:
    print(">> Modalità completa: normalizzazione di tutti i tessuti.\n")

# === TIMER ===
global_start = time.time()

print(">> Caricamento matrice TPM...")
tpm = pd.read_csv(TPM_PATH, index_col=0)
print(">> TPM caricata correttamente.\n")

print(">> Caricamento samples_raggruppati.json...")
with open(SAMPLES_PATH, "r") as f:
    samples_by_tissue_age = json.load(f)
print(f">> Trovati {len(samples_by_tissue_age)} tessuti nel JSON.\n")

# === ELABORAZIONE TESSUTI ===
counter = 0
skipped = 0
total = len(samples_by_tissue_age)

for tissue, age_groups in samples_by_tissue_age.items():
    # Salta il tessuto se non è nella lista richiesta
    if tissues_filter and tissue not in tissues_filter:
        skipped += 1
        continue

    counter += 1
    print(f">> Normalizzazione tessuto ({counter}): {tissue}")

    # Raccogli sample validi
    sample_ids = []
    for age_group, ids in age_groups.items():
        sample_ids.extend(ids)

    sample_ids = [sid for sid in sample_ids if sid in tpm.columns]

    if not sample_ids:
        print("   [!] Nessun sample valido, saltato.\n")
        continue

    # Submatrice
    sub_tpm = tpm[sample_ids].copy()

    # Min-Max normalization
    min_vals = sub_tpm.min(axis=1)
    max_vals = sub_tpm.max(axis=1)
    ranges = max_vals - min_vals
    ranges[ranges == 0] = 1

    norm = (sub_tpm.sub(min_vals, axis=0)).div(ranges, axis=0)

    # Mapping sample → age
    sample_to_age = {}
    for age_group, ids in age_groups.items():
        for sid in ids:
            if sid in norm.columns:
                sample_to_age[sid] = age_group

    df_norm = norm.T
    df_norm["age_group"] = df_norm.index.map(sample_to_age)
    df_norm = df_norm.dropna(subset=["age_group"])

    # Media per fascia d'età
    df_avg = df_norm.groupby("age_group").mean().T

    # Output
    out_path = os.path.join(OUTPUT_DIR, f"{tissue}_normalized.csv")
    df_avg.to_csv(out_path)

    print("   -> Salvato.\n")

# === RIEPILOGO ===
total_time = time.time() - global_start
print("=== STEP 2 COMPLETATO ===")
print(f"Tessuti normalizzati: {counter} | Saltati (già esistenti): {skipped}")
print(f"Tempo totale: {total_time/60:.2f} minuti")
