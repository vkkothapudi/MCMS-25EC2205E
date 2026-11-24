# MCMS-25EC2205E: Modern Communication Systems

Course pages for **25EC2205E Modern Communication Systems**: lecture notes, derivations, interactive simulations, solved numericals, programs, assessment banks and lab material, one page per module.

| Page | Module | Course outcome |
|---|---|---|
| [`index.html`](index.html) | **M1 Analog Modulation**: AM, AM demodulation, FM, FM demodulation, PM, broadcast contexts, GNU Radio labs. Its LABS tab also contains M2, M3, practicals P1–P13 (Practical 3 has two parts: 3A FM and 3B PM) and WiSpec skilling S1–S24. | CO1 |
| [`m2_digital_modulation.html`](m2_digital_modulation.html) | **M2 Digital Modulation**: why digital, BPSK, QPSK/M-PSK, M-QAM, BER analysis, Gray coding, 5G NR / Wi-Fi modulation, Python BER lab | CO2 |
| [`m3_pulse_shaping_isi.html`](m3_pulse_shaping_isi.html) | **M3 Pulse Shaping and ISI**: ISI, Nyquist criterion, raised cosine, matched filter, eye diagrams, equalization, Python eye-diagram lab | CO3 |

Each module page has a START tab with a live simulation, one lecture per syllabus topic, the syllabus lab, verified programs, a quiz and assignment bank, a syllabus map with references, and a LABS tab.

## Opening the pages

Open the HTML files in any modern browser. Keep the folder structure as it is: the pages load their figures from `img/` and `labs/`. An internet connection is needed for the equation renderer (MathJax) and fonts.

## Folders

- `img/`: figures, program plots and GNU Radio screenshots
- `labs/prac/`, `labs/wis/`: images for the lab practicals and WiSpec skilling sheets
- `labs/m2/`, `labs/m3/`: images for the M2 and M3 pages shown inside the M1 LABS tab
- `labs/prac/p3a/`, `labs/prac/p3b/`: programs for Practical 3A (FM: Bessel spectrum and Carson's rule) and Practical 3B (PM vs FM and PM demodulation), each as Python, MATLAB/Octave and a GNU Radio flowgraph
- `labs/grc/`: GNU Radio Companion flowgraphs (`.grc`)

### Practical numbering

Practical 3 was added in two parts, so the former P3–P12 are now P4–P13 on the M1 page. The M2 and M3 pages and the separate lab manual and record book still use the old P1–P12 numbering.
