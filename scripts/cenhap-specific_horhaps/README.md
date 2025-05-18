# cenhap-specific_horhaps

Usage: `./cenhap-specific_horhaps.sh stv.bed asseblies.fa 300`
- stv.bed = output of stv tool
- assemblies.fa must contain all assablies of a given chromosome with cenhap assignments (last character in fasta header line after "_". e.g. >HG02572.pat.cen10h1_2). They must be same assemblies used to produse stv.bed
- 300 - size of random sample of given stv from given cenhap assemblies. More samples => longer alignmnent
