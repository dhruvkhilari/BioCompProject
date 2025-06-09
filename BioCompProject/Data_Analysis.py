from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# loading the genbank
gbk_path = r"C:\Users\sauri\OneDrive\Desktop\BioComp Club\Data.gbk"
record = SeqIO.read(gbk_path, "genbank")
seq = str(record.seq).upper()
length = len(seq)

#calc cpg density
cpg_sites = [i + 1 for i in range(length - 1) if seq[i:i+2] == "CG"]
window, step = 200, 50
win_mid, density = [], []
for s in range(0, length - window + 1, step):
    win_mid.append(s + window // 2)
    density.append(seq[s:s + window].count("CG"))

#features
features = []
for feat in record.features:
    if feat.type in ("CDS", "gene", "misc_feature"):
        name = feat.qualifiers.get("gene",
               feat.qualifiers.get("product",
               feat.qualifiers.get("label", [""])))[0]
        if name:
            start = int(feat.location.start) + 1
            end = int(feat.location.end)
            features.append((start, end, name))

# periodicity score
window_p, step_p = 147, 10
positions, scores = [], []
for i in range(0, length - window_p + 1, step_p):
    wseq = seq[i:i+window_p]
    score = sum(1 for k in range(0, window_p - 1, 10)
                if wseq[k:k+2] in ("AA", "TT", "TA", "AT"))
    positions.append(i + window_p // 2 + 1)
    scores.append(score)
positions = np.array(positions)
scores = np.array(scores)
mean_score = scores.mean()

# formatting and gui
plt.style.use("ggplot")
plt.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
})

# creating the plot
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(4, 1, height_ratios=[2, 1, 0.6, 1.8], hspace=0.5)

# Panel A: CpG Density
ax1 = fig.add_subplot(gs[0])
ax1.plot(win_mid, density, color="#1f77b4", lw=1.5)
ax1.fill_between(win_mid, density, color="#1f77b4", alpha=0.2)
ax1.set_ylabel("CpG / 200 bp")
ax1.set_title("A) CpG Density Across Plasmid")
ax1.grid(axis="y", linestyle="--", alpha=0.4)

# Panel B: CpG Site Events
ax2 = fig.add_subplot(gs[1], sharex=ax1)
ax2.eventplot(cpg_sites, orientation="horizontal",
              colors="#2ca02c", linelengths=0.8)
ax2.set_yticks([])
ax2.set_title("B) Individual CpG Site Positions")
ax2.grid(axis="y", linestyle="--", alpha=0.3)

# Panel C: Annotated Features
ax3 = fig.add_subplot(gs[2], sharex=ax1)
ax3.set_ylim(0, 1)
for start, end, name in features:
    ax3.broken_barh([(start, end - start)], (0, 1),
                    facecolors="#ff7f0e", edgecolors="#d62728", alpha=0.5)
    ax3.text((start + end) / 2, 0.5, name,
             rotation=45, ha="right", va="center",
             fontsize=7, color="#4b0082", fontweight="bold")
ax3.set_yticks([])
ax3.set_title("C) Gene Locations")
ax3.grid(axis="x", linestyle="--", alpha=0.3)

# Panel D: Nucleosome Periodicity Score
ax4 = fig.add_subplot(gs[3], sharex=ax1)
ax4.plot(positions, scores, color="#6a5acd", lw=1.5)
ax4.fill_between(positions, scores, color="#6a5acd", alpha=0.3)
ax4.axhline(mean_score, color="red", linestyle="--", lw=1,
            label=f"Mean = {mean_score:.2f}")
ax4.set_ylabel("Affinity Score")
ax4.set_xlabel("Position on Plasmid (bp)")
ax4.set_title("D) Predicted Nucleosome Affinity")
ax4.legend(loc="upper right")
ax4.grid(axis="y", linestyle="--", alpha=0.4)
plt.setp(ax4.get_xticklabels(), rotation=45, ha="right")

#showing the plot
plt.tight_layout()
plt.show()