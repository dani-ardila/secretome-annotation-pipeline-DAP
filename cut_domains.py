# cut_simple_fixed.py
import os

def main():
    print("=== Cortador de subsecuencias (FASTA) ===")
    acc = input("Identificador (p.ej., ACC123) [opcional]: ").strip() or "SEQ"
    seq = input("Secuencia de aminoácidos (una sola línea, sin espacios): ").strip().upper()
    start = int(input("Posición de inicio (1-based): ").strip())
    end = int(input("Posición de fin (1-based, inclusive): ").strip())

    # --- Validaciones mínimas ---
    if start < 1 or end < 1 or start > end:
        raise ValueError("Rango inválido: asegúrate de que 1 <= inicio <= fin.")
    if end > len(seq):
        raise ValueError(f"El fin ({end}) excede la longitud de la secuencia ({len(seq)}).")

    # --- Corte 1-based -> índices Python 0-based ---
    subseq = seq[start-1:end]
    sublen = len(subseq)

    # --- Construir encabezado FASTA ---
    header = f">{acc}_dom_{start}_{end}|len={sublen}"

    # --- Ruta de salida fija ---
    out_dir = r"C:\Users\Daniela\OneDrive - Universidad Nacional de Colombia\Escritorio\Tareas\Trabajo C.elegans\domains"
    os.makedirs(out_dir, exist_ok=True)  # por si la carpeta no existe
    out_path = os.path.join(out_dir, "subsequence.fasta")

    # --- Guardar en FASTA ---
    with open(out_path, "a") as f:  # "a" = append, así no sobrescribe las secuencias anteriores
        f.write(header + "\n")
        for i in range(0, sublen, 70):  # envuelve en líneas de 70 caracteres
            f.write(subseq[i:i+70] + "\n")

    print(f"Listo ✅  Guardado en: {out_path}")
    print(f"Encabezado: {header}")

if __name__ == "__main__":
    main()
