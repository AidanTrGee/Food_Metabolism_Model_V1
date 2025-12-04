from __future__ import annotations
import csv, os
import matplotlib
matplotlib.use("Agg")  # headless backend
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

HERE = os.path.dirname(__file__)
ROOT = os.path.normpath(os.path.join(HERE, "../../.."))
CSV_PATH = os.path.join(ROOT, "out", "sim_output.csv")
PLOT_DIR = os.path.join(ROOT, "out", "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

def load_csv(path):
    with open(path, "r", newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        rows = list(r)
        times = [float(rec["time_min"]) for rec in rows]
        return times, rows, r.fieldnames

def series(times, rows, col):
    vals = []
    for rec in rows:
        v = rec.get(col, "")
        vals.append(float(v) if v != "" else 0.0)
    return vals

def main():
    if not os.path.exists(CSV_PATH):
        raise SystemExit(f"Expected CSV at {CSV_PATH}. Run the simulation first.")

    times, rows, cols = load_csv(CSV_PATH)

    targets = {
        "Arterial Blood Glucose (mg/dL)": "ArterialBlood.glucose_mg_dl",
        "Venous Blood Glucose (mg/dL)": "VenousBlood.glucose_mg_dl",
        "Portal Vein Glucose (mg/dL)": "PortalVein.glucose_mg_dl",
        "Liver Glucose (mg/dL)": "Liver.glucose_mg_dl",
        "Intestinal Lumen Glucose (mg/dL)": "IntestinalLumen.glucose_mg_dl",
        "Muscle Glucose (mg/dL)": "Muscle.glucose_mg_dl",
        "Adipose Lipid (mg/dL)": "Adipose.lipid_mg_dl",
        "Venous Blood Lipid (mg/dL)": "VenousBlood.lipid_mg_dl",
        "Lymph Lipid (mg/dL)": "Lymph.lipid_mg_dl",
        "Renal Glucose Excreted (g)": "signal.renal_glucose_excreted_total_g",
        "Insulin Signal": "signal.insulin",
        "Glucagon Signal": "signal.glucagon",
    }

    made = []
    pdf_path = os.path.join(PLOT_DIR, "all_plots.pdf")
    with PdfPages(pdf_path) as pdf:
        for title, col in targets.items():
            if cols is None or col not in cols:
                print(f"Skipping {title}: missing column '{col}' in CSV.")
                continue
            y = series(times, rows, col)
            plt.figure()
            plt.plot(times, y, linewidth=2.0)
            plt.xlabel("Time (min)")
            plt.ylabel(title)
            plt.title(title)
            plt.tight_layout()
            fname = title.lower().replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_") + ".png"
            out_png = os.path.join(PLOT_DIR, fname)
            plt.savefig(out_png, dpi=160)
            pdf.savefig()
            plt.close()
            made.append(out_png)

    if not made:
        print("No plots generated; verify the CSV contents and column names.")
        return

    print("Saved:")
    for p in made:
        print(" -", os.path.relpath(p, ROOT))
    print(" -", os.path.relpath(pdf_path, ROOT))

if __name__ == "__main__":
    main()
