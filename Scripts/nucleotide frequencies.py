from collections import Counter

fa = "psph_5mC_cov30_pm10.fa"
seqs = []
with open(fa) as f:
    for line in f:
        if line.startswith(">"):
            s = next(f).strip().upper()
            seqs.append(s)

L = len(seqs[0]); c = L//2
pos = {
    "-3": Counter(s[c-3] for s in seqs),
    "-2": Counter(s[c-2] for s in seqs),
    "-1": Counter(s[c-1] for s in seqs),
    "+1": Counter(s[c+1] for s in seqs),
    "+2": Counter(s[c+2] for s in seqs),
    "+3": Counter(s[c+3] for s in seqs),
}
for k in ["-3","-2","-1","+1","+2","+3"]:
    print(k, pos[k])
