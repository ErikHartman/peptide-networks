# peptide-networks

Datan är uppdelad i _ninf_ (non infected), _inf_ (infected), och _WF_ (akut steril sårvätska).
Äger inte datan så sprid inte osv osv osv

## Idéer

### skapa nätverk

Detta är vår "viktigaste" uppgift: hur skapar vi nätverk över peptidomet på ett bra sätt. Om vi vill gruppera peptider
beroende på region kan vi använda Levenshtein-sträckan. Vi kan antingen analysera alla peptider i ett prov,
eller peptider som kommer från ett visst protein. Om vi däremot vill gruppera peptider beroende på funktion skulle vi behöva
en matris där lika aminosyror (e.g. alanin och glycin) ges ett lågt värde, men olika ett högt värde. Jag tror att
dessa två metoder är bra startpunkter. Då vi skapar nätverken kommer vi behöva undersöka vilka thresholds som är rimliga
att applicera för att skapa en kant mellan två noder. Hur vi gör detta bör vi diskutera.

### nätverkskarakteristik

Det första vi förmodligen kommer titta på när vi väl skapat de olika nätverken är dess struktur och karakteristik.
Kan vi säga ngt i allmänhet om peptidomnätverken. Vad fångar vi etc?? Här kan vi använda mått inom grafteori.
Därefter vill vi se om dessa mått skiljer sig mellan nätverken. Vad signifierar de olika måtten?

### enzym-prediktion

Enzym-prediktion är svårt då alla peptider tas i åtanke pga. "brus" och ospecifitet i datan då exopeptidas skymmer de mer intressanta enzymen.
Genom att identifiera den längsta peptiden i varje community kan vi skapa ett subset där vi vet att exopeptidasen har minimal verkan.
Algoritmen blir alltså: network --> community prediction --> get longest peptides --> predict proteases

---
## Installing python package
To run pytest from other folder one needs to run:

pip install -e .

in source folder.

## Folder structure

- /data
- /peptide_networks
  - /experiment 1
  - /experiment 2
- /test\
  - /experiment 1
  - /experiment 2
- /results
  - /experiment 1
  - /experiment 2

---

## Workflow

Set up pytest

### Mattias

1. Proteasprediktion
	Lite frygt 30 % av alla peptider är direkta derivat av andra peptider i populationen.
	Ca 70 % av alla peptider kan härledas till andra peptider. 
### Erik

1. skapar edge lists
   - BLOSUM62
   - Testa andra matriser
2. skapa nätverk från edge list och mata ut lite karakteristik
3. plotta nätverk
