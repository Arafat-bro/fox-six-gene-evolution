import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

genes = {
    "Fol4287":     [("SIX10", 5000, 5516, "-"), ("SIX12", 6743, 7171, "-"), ("SIX7", 9482, 10141, "+")],
    "race3":       [("SIX7", 5000, 5659, "-"), ("SIX12", 7970, 8398, "-"), ("SIX10", 9625, 10141, "+")],
    "cepae":       [("SIX12", 5000, 5359, "+"), ("SIX10", 8600, 8962, "+")],
    "lini_block1": [("SIX7", 5000, 5641, "-"), ("SIX12", 6547, 6975, "+"), ("SIX10", 8874, 9390, "+")],
}
colors = {"SIX7": "#1f77b4", "SIX10": "#d62728", "SIX12": "#2ca02c"}
row_y = {"Fol4287": 3, "race3": 2, "cepae": 1, "lini_block1": 0}

hits = {
    "race3": [(2, 15142, 15142, 4, 99.941, 15143)],
    "cepae": [(8103,9203,5863,4767,96.370,1102), (4734,5785,7371,6312,95.476,1061),
              (2701,3857,10379,9231,92.876,1165), (7277,7492,12437,12656,90.909,220),
              (7026,7192,6224,6058,97.006,167), (7882,8110,12655,12421,87.764,237),
              (4127,4352,7935,8161,85.903,227)],
    "lini_block1": [(6350,7488,7371,6223,95.490,1153), (4735,5936,10390,9180,93.355,1219),
                     (8532,9544,5860,4848,96.252,1014), (7697,7863,6224,6058,97.605,167),
                     (5947,6125,7730,7559,89.385,179)],
}

fig, ax = plt.subplots(figsize=(13, 7))

for row in ["race3", "cepae", "lini_block1"]:
    y_top, y_bot = row_y["Fol4287"], row_y[row]
    for qs, qe, ss, se, pident, length in hits[row]:
        reverse = ss > se
        # normalize subject coords so the polygon never self-crosses
        s_left, s_right = (se, ss) if reverse else (ss, se)
        alpha = min(0.15 + (pident - 75) / 100, 0.85)
        color = "#4477aa" if not reverse else "#cc6677"  # blue = same strand, red-gray = inverted vs Fol4287
        verts = [(qs, y_bot + 0.15), (qe, y_bot + 0.15), (s_right, y_top - 0.15), (s_left, y_top - 0.15), (qs, y_bot + 0.15)]
        ax.add_patch(mpatches.Polygon(verts, closed=True, facecolor=color, alpha=alpha, edgecolor="none"))

for name, glist in genes.items():
    y = row_y[name]
    ax.axhline(y, color="black", lw=0.5, zorder=1)
    for gname, start, end, strand in glist:
        x0, x1 = (start, end) if strand == "+" else (end, start)
        ax.annotate("", xy=(x1, y), xytext=(x0, y),
                    arrowprops=dict(arrowstyle="-|>", color=colors[gname], lw=7, mutation_scale=22))
        offset = 0.35 if name == "Fol4287" else -0.35
        va = "bottom" if name == "Fol4287" else "top"
        ax.text((start+end)/2, y + offset, gname, ha="center", va=va, fontsize=9, weight="bold")
    ax.text(-600, y, name, ha="right", va="center", fontsize=11, weight="bold")

ax.set_xlim(-1500, 16000)
ax.set_ylim(-0.8, 3.8)
ax.set_yticks([])
ax.set_xlabel("position in window (bp)")
ax.set_title("SIX7-SIX10-SIX12 region: gene arrangement and flanking-sequence conservation vs Fol4287")
handles = [
    mpatches.Patch(color="#4477aa", alpha=0.7, label="match, same orientation as Fol4287"),
    mpatches.Patch(color="#cc6677", alpha=0.7, label="match, inverted relative to Fol4287"),
]
ax.legend(handles=handles, loc="upper right", fontsize=8)
plt.tight_layout()
plt.savefig("results/synteny/six7_10_12_flanking_comparison.png", dpi=200)
print("saved")
