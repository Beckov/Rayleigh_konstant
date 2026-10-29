# Rayleighova konštanta z hysteréznych slučiek

Webová aplikácia a Python skript na výpočet **Rayleighovej konštanty η** a **počiatočnej permeability μᵢ**
z hysteréznych slučiek nameraných pri nízkych poliach (pod koercivitou).

**Aplikácia online:** `https://<používateľ>.github.io/rayleigh-konstanta/`

## Model

Pri nízkych poliach (Rayleighova oblasť) sa hysterézna slučka skladá z dvoch parabol:

```
B = μ0 [ (μi + η·Hm)·H ± (η/2)·(Hm² − H²) ]
```

| Veličina | Jednotka | Význam |
|---|---|---|
| B | T | magnetická indukcia |
| H | A/m | okamžitá intenzita magnetického poľa |
| Hm | A/m | amplitúda intenzity poľa danej slučky |
| μ0 | H/m | permeabilita vákua, 4π·10⁻⁷ |
| μi | – | počiatočná relatívna permeabilita |
| η | m/A | Rayleighova konštanta |
| ± | | + zostupná (horná) vetva, − vzostupná (dolná) vetva |

## Výpočet

1. Každá slučka sa podľa smeru zmeny H(t) rozdelí na zostupnú a vzostupnú vetvu.
2. Obe vetvy sa naraz preložia metódou najmenších štvorcov modelom
   `B = a·H ± c·(Hm² − H²) + d`.
3. Z koeficientov: `η = 2c/μ0`, `μi = a/μ0 − η·Hm`.
4. Výsledné η a μi sú priemery cez slučky vo zvolenom rozsahu Hm (predvolene Hm ≤ Hc/3),
   ± je výberová smerodajná odchýlka.
5. Kontrola: lineárny fit amplitúdovej permeability `μa = Bm/(μ0·Hm) = μi + η·Hm`.

Rayleighov zákon platí iba hlboko pod koercivitou, kde je η konštantné. Na grafe η(Hm)
treba vybrať oblasť plató. Najnižšie slučky môžu mať η mierne vyššie, pretože fázový posun
alebo vírivé prúdy slučku otvárajú elipticky.

## Vstupné dáta

- CSV z meracieho softvéru: tabulátorom oddelené bloky `Hₘ = … A/m` so stĺpcami
  `t, H(t), v_ind(t), B(t), μ_diff(t)` (jedna slučka = jeden blok).
- Alebo jednoduchý súbor s dvoma stĺpcami `H` a `B` (jedna slučka).

## Webová aplikácia (`index.html`)

Otvor stránku, nahraj jeden alebo viac CSV súborov (kliknutím alebo pretiahnutím).
Dáta sa spracúvajú len v prehliadači, nikam sa neodosielajú.

- **Hc** sa navrhne z najväčšej slučky v súbore, dá sa prepísať ručne.
- **Hm od / do** nastavuje rozsah slučiek; tlačidlo *Hm ≤ Hc/3* ho nastaví automaticky.
- Kliknutím na bod v grafe alebo na zaškrtávacie pole v tabuľke slučku vylúčiš alebo vrátiš.
- Tabuľku výsledkov skopíruješ s desatinnou čiarkou (Excel) alebo bodkou (Veusz).

## Python skript (`rayleigh_konstanta.py`)

```
pip install numpy matplotlib
python rayleigh_konstanta.py                  # otvorí okno na výber súborov
python rayleigh_konstanta.py subor1.csv ...   # spracuje zadané súbory
```

Nastavenia sú na začiatku skriptu (`HC_PYTAT_SA`, `HC_RUCNE`, `HM_MIN`, `HM_MAX`, `VYLUCIT_HM`).
Ku každému súboru uloží `*_rayleigh.csv` s výsledkami a `*_rayleigh.png` s grafmi.
