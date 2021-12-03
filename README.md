# peptide-networks

Datan är uppdelad i _ninf_ (non infected), _inf_ (infected), och _WF_ (akut steril sårvätska).
Äger inte datan så sprid inte osv osv osv

## Idéer
### enzym-prediktion
Enzym-prediktion är svårt då alla peptider tas i åtanke pga. "brus" och ospecifitet i datan då exopeptidas skymmer de mer intressanta enzymen.
Genom att identifiera den längsta peptiden i varje community kan vi skapa ett subset där vi vet att exopeptidasen har minimal verkan.
Algoritmen blir alltså: network --> community prediction --> get longest peptides --> predict proteases


