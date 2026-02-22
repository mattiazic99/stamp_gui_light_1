'''# === step_3a_sets_stamp.py ===
import pandas as pd
import os
import time
from tqdm import tqdm


# === PATH BASE ===
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(BASE_DIR)

OUTPUT_DIR = os.path.join(ROOT_DIR, "output", "sets_stamp")
NORMALIZED_DIR = os.path.join(ROOT_DIR, "output", "tpm_normalizzati")

os.makedirs(OUTPUT_DIR, exist_ok=True)
THRESHOLD = 0.5  # threshold to consider a gene 'expressed'
AGE_GROUPS = ["30-39", "40-49", "50-59", "60-69", "70-79"]

# === Setup ===
os.makedirs(OUTPUT_DIR, exist_ok=True)
files = sorted([f for f in os.listdir(NORMALIZED_DIR) if f.endswith("_normalized.csv")])
num_files = len(files)

print(f">> Trovati {num_files} file normalizzati da elaborare...\n")

# === Timer globale ===
global_start = time.time()
elapsed_times = []

def is_switching_gene(binary_vector):
    changes = [binary_vector[i] != binary_vector[i+1] for i in range(len(binary_vector)-1)]
    return sum(changes) == 1

for i, file in enumerate(tqdm(files, desc="Tessuti")):
    start_time = time.time()
    tissue = file.replace("_normalized.csv", "")
    path = os.path.join(NORMALIZED_DIR, file)

    print(f"\n>>> Elaborazione: {tissue}")
    df = pd.read_csv(path, index_col=0)

    if df.empty:
        print(f"[!] File vuoto: {file}")
        continue

    # Binarizzazione con soglia
    df_bin = df.applymap(lambda x: 1 if x >= THRESHOLD else 0)

    # Mantieni solo le colonne delle 5 fasce età previste
    df_bin = df_bin[[col for col in AGE_GROUPS if col in df_bin.columns]]

    # Seleziona solo i geni switching (1 solo cambio)
    valid_genes = df_bin.apply(lambda row: is_switching_gene(row.tolist()), axis=1)
    switching_genes = df_bin[valid_genes]
    print(f">> Geni switching trovati: {switching_genes.shape[0]}")

    # Scrittura file: 5 righe esatte, una per ogni fascia (anche se vuota)
    output_path = os.path.join(OUTPUT_DIR, f"{tissue}_sets_stamp.txt")
    with open(output_path, "w") as f_out:
        for age_group in AGE_GROUPS:
            if age_group in switching_genes.columns:
                genes = switching_genes.index[switching_genes[age_group] == 1].tolist()
                f_out.write(" ".join(genes) + "\n")
            else:
                f_out.write("\n")  # fascia mancante → riga vuota

    elapsed = time.time() - start_time
    elapsed_times.append(elapsed)
    print(f">>> Salvato: {output_path}")
    print(f">>> Tempo impiegato per {tissue}: {elapsed:.2f} secondi")

    if (i + 1) % 5 == 0 and i + 1 < num_files:
        avg_time = sum(elapsed_times) / len(elapsed_times)
        remaining = avg_time * (num_files - i - 1)
        print(f">>> Tempo stimato rimanente: {remaining / 60:.2f} minuti")

# === Riepilogo finale ===
total_time = time.time() - global_start
print(f"\n=== FINE: calcolo geni switching STAMP completato per {num_files} tessuti ===")
print(f"Tempo totale: {total_time / 60:.2f} minuti")
'''
'''
# ===================================
#   STEP 3A - GENERAZIONE SETS STAMP
# ===================================

import pandas as pd
import os
import time
from tqdm import tqdm

# === PATH BASE ===
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(BASE_DIR)

NORMALIZED_DIR = os.path.join(ROOT_DIR, "output", "tpm_normalizzati")
OUTPUT_DIR = os.path.join(ROOT_DIR, "output", "sets_stamp")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Variabile DEFINITA al runtime
THRESHOLD = None

# Ordine fasce età STAMP
AGE_GROUPS = ["30-39", "40-49", "50-59", "60-69", "70-79"]


def is_switching_gene(row):
    vec = row.tolist()
    changes = sum(vec[i] != vec[i+1] for i in range(len(vec)-1))
    return changes == 1


def run_step3(threshold):
    print(f">> Threshold utilizzata: {threshold}")

    files = sorted([f for f in os.listdir(NORMALIZED_DIR) if f.endswith("_normalized.csv")])
    print(f">> File trovati: {len(files)}")

    start_global = time.time()

    for file in tqdm(files, desc="Tessuti"):
        tissue = file.replace("_normalized.csv", "")
        path = os.path.join(NORMALIZED_DIR, file)

        print(f"\n>>> Elaborazione tissue: {tissue}")
        start = time.time()

        df = pd.read_csv(path, index_col=0)
        if df.empty:
            print(f"[!] File vuoto: {file}")
            continue

        df_bin = df.applymap(lambda x: 1 if x >= threshold else 0)
        df_bin = df_bin[[c for c in AGE_GROUPS if c in df_bin.columns]]

        mask = df_bin.apply(is_switching_gene, axis=1)
        df_switch = df_bin[mask]

        print(f">> Geni switching: {df_switch.shape[0]}")

        out_path = os.path.join(OUTPUT_DIR, f"{tissue}_sets_stamp.txt")
        with open(out_path, "w") as f:
            for age in AGE_GROUPS:
                if age in df_switch.columns:
                    active_genes = df_switch.index[df_switch[age] == 1].tolist()
                    f.write(" ".join(active_genes) + "\n")
                else:
                    f.write("\n")

        print(f">>> Salvato: {out_path}")
        print(f">>> Tempo: {time.time()-start:.2f}s")

    print("\n=== STEP 3A COMPLETATO ===")
    print(f"Tempo totale: {(time.time()-start_global)/60:.2f} minuti")


if __name__ == "__main__":
    threshold_env = os.environ.get("STAMP_THRESHOLD")
    if threshold_env is None:
        raise ValueError("Errore: variabile STAMP_THRESHOLD non trovata.")
    THRESHOLD = float(threshold_env)

    run_step3(THRESHOLD)

    '''

# ===================================
#   STEP 3A - GENERAZIONE SETS STAMP
# ===================================

import pandas as pd
import os
import time

# === PATH BASE ===
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(BASE_DIR)

NORMALIZED_DIR = os.path.join(ROOT_DIR, "output", "tpm_normalizzati")
OUTPUT_DIR = os.path.join(ROOT_DIR, "output", "sets_stamp")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Valore definito via variabile d'ambiente
THRESHOLD = None

# Ordine fasce d’età STAMP
AGE_GROUPS = ["30-39", "40-49", "50-59", "60-69", "70-79"]


def is_switching_gene(row):
    vec = row.tolist()
    return sum(vec[i] != vec[i+1] for i in range(len(vec)-1)) == 1


def run_step3(threshold):
    print(f">> Threshold used: {threshold}")

    files = sorted(
        f for f in os.listdir(NORMALIZED_DIR)
        if f.endswith("_normalized.csv")
    )
    print(f">> Files found: {len(files)}\n")

    start_global = time.time()
    counter = 0
    total = len(files)

    for file in files:
        counter += 1
        tissue = file.replace("_normalized.csv", "")
        print(f">> Processing ({counter}/{total}): {tissue}")

        path = os.path.join(NORMALIZED_DIR, file)
        df = pd.read_csv(path, index_col=0)

        if df.empty:
            print("   [!] Empty file, skipped.\n")
            continue

        # Binarizzazione
        df_bin = df.map(lambda x: 1 if x >= threshold else 0)
        df_bin = df_bin[[c for c in AGE_GROUPS if c in df_bin.columns]]

        # Selezione geni switching
        mask = df_bin.apply(is_switching_gene, axis=1)
        df_switch = df_bin[mask]

        # Scrittura file
        out_path = os.path.join(OUTPUT_DIR, f"{tissue}_sets_stamp.txt")
        with open(out_path, "w") as f:
            for age in AGE_GROUPS:
                if age in df_switch.columns:
                    genes = df_switch.index[df_switch[age] == 1].tolist()
                    f.write(" ".join(genes) + "\n")
                else:
                    f.write("\n")

        print("   -> Saved.\n")

    print("=== Processing completed ===")
    print(f"Total time: {(time.time() - start_global)/60:.2f} minutes")
    print("")
    

if __name__ == "__main__":
    threshold_env = os.environ.get("STAMP_THRESHOLD")
    if threshold_env is None:
        raise ValueError("Errore: variabile STAMP_THRESHOLD non trovata.")
    THRESHOLD = float(threshold_env)

    run_step3(THRESHOLD)



