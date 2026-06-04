# Data Dictionary — California Proposition 99 Panel

Source: Abadie, Diamond, and Hainmueller (2010), distributed with the Arkhangelsky et al. (2021)
replication package (OpenICPSR DOI: 10.3886/E146381V1).

File: `input/california_prop99.csv`  
Delimiter: semicolon (`;`)  
Rows: 1,209 (39 states × 31 years, 1970–2000)

| Variable        | Type    | Units              | Description                                                                 |
|-----------------|---------|--------------------|-----------------------------------------------------------------------------|
| `State`         | string  | —                  | U.S. state name (39 states; California is the treated unit)                 |
| `Year`          | integer | year               | Calendar year, 1970–2000                                                    |
| `PacksPerCapita`| float   | packs per person   | Annual per-capita cigarette sales (outcome variable)                        |
| `treated`       | integer | 0 / 1              | Binary treatment indicator: 1 for California from 1989 onward, 0 otherwise |

Notes:
- California is the sole treated unit (Proposition 99 took effect January 1, 1989).
- Pre-treatment period: 1970–1988 (T₀ = 19 years).
- Post-treatment period: 1989–2000 (T₁ = 12 years).
- 38 donor (control) states.
