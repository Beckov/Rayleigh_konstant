"""
Výpočet Rayleighovej konštanty η z hysteréznych slučiek pri nízkych poliach.

Model (Rayleighov zákon):
    B = μ0 [ (μi + η·Hm)·H ± (η/2)·(Hm² − H²) ]
    + … zostupná (horná) vetva, − … vzostupná (dolná) vetva

Každá slučka sa preloží metódou najmenších štvorcov (obe vetvy naraz).
Výsledné η je priemer cez slučky v zvolenom rozsahu Hm (predvolene Hm ≤ Hc/3).
Kontrola: lineárny fit μa(Hm) = μi + η·Hm.

Vstup: CSV z meracieho softvéru (tabulátorom oddelené bloky „Hₘ = … A/m“
so stĺpcami t, H(t), v_ind(t), B(t), μ_diff(t)), alebo súbor s dvoma stĺpcami H a B.

Použitie:
    python rayleigh_konstanta.py                 -> otvorí okno na výber súborov
    python rayleigh_konstanta.py subor1.csv ...  -> spracuje zadané súbory

Potrebné knižnice:  pip install numpy matplotlib
"""

import csv
import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt

# ======================= NASTAVENIA =======================
HC_PYTAT_SA = True     # True = pri každom súbore sa opýta na Hc (Enter = použije HC_RUCNE / automatiku)
HC_RUCNE = {}          # Hc v A/m pre jednotlivé vzorky podľa časti názvu súboru, napr.
                       #   {"Sm000": 99.5, "Sm002": 106.9, "Li05_Fe25": 757.3}
                       # alebo jedno číslo pre všetky súbory, napr. HC_RUCNE = 100
                       # Ak pre súbor nie je zadané nič, Hc sa určí z najväčšej slučky v súbore.
HM_MIN = None          # dolná hranica Hm v A/m; None = od najmenšej slučky
HM_MAX = None          # horná hranica Hm v A/m; None = Hc/3
VYLUCIT_HM = []        # Hm (A/m) slučiek, ktoré nechceš rátať, napr. [10, 20]
ULOZIT_GRAFY = True    # uloží PNG s grafmi vedľa vstupného súboru
ULOZIT_TABULKU = True  # uloží výsledky do *_rayleigh.csv vedľa vstupného súboru
ZOBRAZIT_GRAFY = True  # otvorí okno s grafmi
# ===========================================================

MU0 = 4e-7 * np.pi


def nacitaj_slucky(cesta):
    """Vráti zoznam (Hm_nominalne, H, B) zo súboru."""
    with open(cesta, encoding="utf-8-sig", errors="replace") as f:
        text = f.read()
    prvy = "\n".join(text.splitlines()[:8])
    delim = "\t" if "\t" in prvy else (";" if ";" in prvy else ",")
    rows = list(csv.reader(text.splitlines(), delimiter=delim))

    def cislo(s):
        s = (s or "").strip()
        if delim != ",":
            s = s.replace(",", ".")
        try:
            return float(s)
        except ValueError:
            return np.nan

    def vyber(hc, bc):
        H, B = [], []
        for r in rows:
            if len(r) <= max(hc, bc):
                if H:
                    break
                continue
            h, b = cislo(r[hc]), cislo(r[bc])
            if np.isfinite(h) and np.isfinite(b):
                H.append(h)
                B.append(b)
            elif H:
                break
        return np.array(H), np.array(B)

    hlavicka = rows[0] if rows else []
    bloky = []
    for c, bunka in enumerate(hlavicka):
        m = re.search(r"H.*?=\s*([0-9]+(?:[.,][0-9]+)?)", bunka or "")
        if m:
            bloky.append((c, float(m.group(1).replace(",", "."))))

    slucky = []
    if bloky:
        for i, (c, hm) in enumerate(bloky):
            koniec = bloky[i + 1][0] if i + 1 < len(bloky) else c + 5
            hc, bc = c + 1, c + 3
            for r in rows[1:5]:
                for k in range(c, min(koniec, len(r))):
                    t = r[k].strip()
                    if re.match(r"^H(\(|$|\s)", t):
                        hc = k
                    if re.match(r"^B(\(|$|\s)", t):
                        bc = k
            H, B = vyber(hc, bc)
            if len(H) > 20:
                slucky.append((hm, H, B))
    else:
        hc, bc = 0, 1
        for r in rows[:4]:
            for k, t in enumerate(r):
                t = t.strip()
                if re.match(r"^H(\(|$|\s|_)", t):
                    hc = k
                if re.match(r"^B(\(|$|\s|_)", t):
                    bc = k
        H, B = vyber(hc, bc)
        if len(H) > 20:
            slucky.append(((H.max() - H.min()) / 2, H, B))

    if not slucky:
        raise ValueError("V súbore som nenašiel stĺpce H a B s číslami.")
    return slucky


def priechody(X, Y):
    """Priemerná |X| v bodoch, kde Y mení znamienko (cyklicky)."""
    Xn, Yn = np.roll(X, -1), np.roll(Y, -1)
    idx = np.where(((Y <= 0) & (Yn > 0)) | ((Y >= 0) & (Yn < 0)))[0]
    if len(idx) == 0:
        return np.nan
    t = Y[idx] / (Y[idx] - Yn[idx])
    return np.mean(np.abs(X[idx] + t * (Xn[idx] - X[idx])))


def fit_rayleigh(H, B):
    """Fit dvoma parabolami. Vráti slovník s výsledkami pre jednu slučku."""
    n = len(H)
    Hm = (H.max() - H.min()) / 2
    H0 = (H.max() + H.min()) / 2
    k = max(1, round(n / 400))
    dH = np.roll(H, -k) - np.roll(H, k)
    s = np.where(dH < 0, 1.0, -1.0)            # +1 zostupná (horná) vetva
    u = (H - H0) / Hm                           # normované H kvôli presnosti
    X = np.column_stack([u, s * (1 - u ** 2), np.ones(n)])
    (p1, p2, d), *_ = np.linalg.lstsq(X, B, rcond=None)
    a, c = p1 / Hm, p2 / Hm ** 2
    rez = B - X @ np.array([p1, p2, d])
    eta = 2 * c / MU0
    return {
        "Hm": Hm, "H0": H0, "a": a, "c": c, "d": d,
        "eta": eta,
        "mui": a / MU0 - eta * Hm,
        "mua": (B.max() - B.min()) / 2 / (MU0 * Hm),
        "R2": 1 - np.sum(rez ** 2) / np.sum((B - B.mean()) ** 2),
        "Hc": priechody(H, B),
        "Br": priechody(B, H),
    }


def spracuj(cesta):
    meno = os.path.splitext(os.path.basename(cesta))[0]
    print(f"\n=== {meno} ===")
    slucky = nacitaj_slucky(cesta)
    vysl = []
    for hm_nom, H, B in slucky:
        r = fit_rayleigh(H, B)
        r["nom"], r["H"], r["B"] = hm_nom, H, B
        vysl.append(r)
    vysl.sort(key=lambda r: r["Hm"])

    hc = None
    if isinstance(HC_RUCNE, (int, float)) and HC_RUCNE > 0:
        hc = float(HC_RUCNE)
    elif isinstance(HC_RUCNE, dict):
        for kluc, hodnota in HC_RUCNE.items():
            if kluc.lower() in meno.lower():
                hc = float(hodnota)
                break
    hc_auto = vysl[-1]["Hc"]
    if HC_PYTAT_SA:
        navrh = hc if hc else hc_auto
        odp = input(f"Zadaj Hc pre {meno} v A/m [Enter = {navrh:.1f}]: ").strip().replace(",", ".")
        if odp:
            try:
                hc = float(odp)
            except ValueError:
                print("Neplatné číslo, použije sa navrhnutá hodnota.")
    if hc:
        print(f"Hc = {hc:.1f} A/m (zadané ručne)")
    else:
        hc = hc_auto
        print(f"Hc = {hc:.1f} A/m (z najväčšej slučky, Hm = {vysl[-1]['Hm']:.0f} A/m)")

    hm_min = HM_MIN if HM_MIN is not None else vysl[0]["Hm"]
    hm_max = HM_MAX if HM_MAX is not None else hc / 3
    pod_hc = [r for r in vysl if r["Hm"] < hc]

    def zahrnuta(r):
        vylucena = any(abs(r["nom"] - v) < 0.5 or abs(r["Hm"] - v) < 0.5 for v in VYLUCIT_HM)
        return hm_min - 1e-9 <= r["Hm"] <= hm_max + 1e-9 and not vylucena

    vyb = [r for r in pod_hc if zahrnuta(r)]
    if not vyb:
        print("V zvolenom rozsahu nie je žiadna slučka. Uprav HM_MIN / HM_MAX.")
        return None

    eta = np.array([r["eta"] for r in vyb])
    mui = np.array([r["mui"] for r in vyb])
    hm = np.array([r["Hm"] for r in vyb])
    mua = np.array([r["mua"] for r in vyb])
    sklon, icpt = np.polyfit(hm, mua, 1) if len(vyb) >= 2 else (np.nan, np.nan)
    sd = eta.std(ddof=1) if len(eta) > 1 else np.nan

    print(f"\n{'Hm (A/m)':>9} {'η (m/A)':>10} {'μi':>8} {'μa':>8} {'Br (mT)':>9} {'R²':>8}  zahrnutá")
    for r in pod_hc:
        print(f"{r['Hm']:9.1f} {r['eta']:10.4g} {r['mui']:8.2f} {r['mua']:8.2f} "
              f"{r['Br'] * 1e3:9.4g} {r['R2']:8.5f}  {'áno' if zahrnuta(r) else '-'}")
    print(f"\nRayleighova konštanta  η  = {eta.mean():.4g} ± {sd:.2g} m/A "
          f"(n = {len(vyb)}, Hm {hm.min():.1f}–{hm.max():.1f} A/m)")
    print(f"Počiatočná permeabilita μi = {mui.mean():.4g} ± {mui.std(ddof=1) if len(mui) > 1 else np.nan:.2g}")
    print(f"Kontrola μa(Hm):        η  = {sklon:.4g} m/A, μi = {icpt:.4g}")
    if len(eta) > 2 and sd / abs(eta.mean()) > 0.1:
        print("Pozor: η sa v rozsahu mení o viac ako 10 %. Zváž zúženie rozsahu "
              "alebo vylúčenie najnižších slučiek (VYLUCIT_HM).")

    zaklad = os.path.splitext(cesta)[0]
    if ULOZIT_TABULKU:
        with open(zaklad + "_rayleigh.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow([f"# eta = {eta.mean():.6g} +- {sd:.3g} m/A, mui = {mui.mean():.6g}, "
                        f"Hc = {hc:.4g} A/m, n = {len(vyb)}"])
            w.writerow(["Hm", "eta", "mui", "mua", "Br_mT", "R2", "zahrnuta"])
            for r in pod_hc:
                w.writerow([f"{r['Hm']:.6g}", f"{r['eta']:.6g}", f"{r['mui']:.6g}", f"{r['mua']:.6g}",
                            f"{r['Br'] * 1e3:.6g}", f"{r['R2']:.6g}", int(zahrnuta(r))])
        print(f"Tabuľka uložená: {zaklad}_rayleigh.csv")

    # ---------------- grafy ----------------
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    fig.suptitle(f"{meno}:  η = {eta.mean():.4g} ± {sd:.2g} m/A")
    farba = ["tab:blue" if zahrnuta(r) else "lightgray" for r in pod_hc]
    ax[0].scatter([r["Hm"] for r in pod_hc], [r["eta"] for r in pod_hc], c=farba, zorder=3)
    ax[0].hlines(eta.mean(), hm.min(), hm.max(), colors="tab:orange", linestyles="--")
    ax[0].set(xlabel="$H_m$ (A/m)", ylabel="η (m/A)", title="η pre každú slučku")

    r = vyb[len(vyb) // 2]
    Hm_, hh = r["Hm"], np.linspace(-r["Hm"], r["Hm"], 200)
    ax[1].plot(r["H"] - r["H0"], r["B"] * 1e3, ".", ms=2, label="meranie")
    ax[1].plot(hh, (r["a"] * hh + r["c"] * (Hm_ ** 2 - hh ** 2) + r["d"]) * 1e3, "tab:orange", label="fit")
    ax[1].plot(hh, (r["a"] * hh - r["c"] * (Hm_ ** 2 - hh ** 2) + r["d"]) * 1e3, "tab:orange")
    ax[1].axhline(0, color="k", lw=0.8)
    ax[1].axvline(0, color="k", lw=0.8)
    ax[1].set(xlabel="H (A/m)", ylabel="B (mT)", title=f"Slučka $H_m$ = {Hm_:.1f} A/m, R² = {r['R2']:.5f}")
    ax[1].legend()

    ax[2].scatter([r["Hm"] for r in pod_hc], [r["mua"] for r in pod_hc], c=farba, zorder=3)
    if np.isfinite(sklon):
        ax[2].plot(hm, icpt + sklon * hm, "--", color="tab:orange", label=f"sklon η = {sklon:.4g} m/A")
        ax[2].legend()
    ax[2].set(xlabel="$H_m$ (A/m)", ylabel="$μ_a$ (-)", title="Kontrola: $μ_a = μ_i + η H_m$")
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    if ULOZIT_GRAFY:
        fig.savefig(zaklad + "_rayleigh.png", dpi=150)
        print(f"Graf uložený: {zaklad}_rayleigh.png")
    return {"vzorka": meno, "eta": eta.mean(), "eta_sd": sd, "mui": mui.mean(),
            "eta_mua": sklon, "Hc": hc, "n": len(vyb)}


def main():
    subory = sys.argv[1:]
    if not subory:
        try:
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk()
            root.withdraw()
            subory = list(filedialog.askopenfilenames(
                title="Vyber CSV súbory so slučkami",
                filetypes=[("CSV / TXT", "*.csv *.tsv *.txt *.dat"), ("Všetky súbory", "*.*")]))
            root.destroy()
        except Exception:
            pass
    if not subory:
        print("Nezadal si žiadny súbor. Použitie: python rayleigh_konstanta.py subor.csv")
        return

    suhrn = []
    for s in subory:
        try:
            v = spracuj(s)
            if v:
                suhrn.append(v)
        except Exception as e:
            print(f"{s}: chyba – {e}")

    if len(suhrn) > 1:
        print("\n=== Súhrn ===")
        print(f"{'vzorka':<60} {'η (m/A)':>10} {'± sd':>8} {'μi':>8} {'η z μa':>10} {'n':>3}")
        for v in suhrn:
            print(f"{v['vzorka'][:60]:<60} {v['eta']:10.4g} {v['eta_sd']:8.2g} "
                  f"{v['mui']:8.2f} {v['eta_mua']:10.4g} {v['n']:3d}")
    if ZOBRAZIT_GRAFY and suhrn:
        plt.show()


if __name__ == "__main__":
    main()
