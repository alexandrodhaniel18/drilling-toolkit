# Drilling Toolkit: Trajectory

Modul Python untuk menghitung trajectory sumur dari data survei dengan metode
minimum curvature: TVD, Northing, Easting, dogleg severity (DLS), vertical
section, dan closure, lengkap dengan grafik 2D dan 3D.

## Persamaan
- Sudut dogleg: cos β = cos(I2 - I1) - sin I1 sin I2 (1 - cos(A2 - A1))
- Ratio factor: RF = (2/β) tan(β/2)
- ΔN, ΔE, ΔTVD = (ΔMD/2)(komponen arah stasiun 1 + stasiun 2) x RF
- DLS = β / ΔMD x 30 m

Referensi: Sawaryn & Thorogood (2005), SPE 84246; Taylor & Mason (1972), SPE 3362;
Bourgoyne dkk. (1986), Applied Drilling Engineering.

## Instalasi dan pemakaian
    python -m venv .venv
    .venv\Scripts\Activate.ps1
    pip install -r requirements.txt
    python main.py

## Konvensi
Inklinasi dari vertikal (deg), azimuth dari utara searah jarum jam (deg),
MD dan TVD dalam meter, DLS dalam deg/30 m. Azimuth diasumsikan sudah
dikoreksi ke acuan koordinat yang sama.

## Validasi
`pytest` menjalankan 8 tes: sumur vertikal, interval lurus, busur 90° (hasil
harus sama dengan jari-jari), sifat fisik (MD >= TVD, jarak lurus <= MD),
closure, vertical section, dan pembacaan CSV.

## Contoh hasil
![Plan dan section view](outputs/well_2d.png)
![Trajectory 3D](outputs/well_3d.png)

## Keterbatasan
Model deterministik: tanpa ketidakpastian survei (ellipse of uncertainty).
Data contoh bersifat sintetis.