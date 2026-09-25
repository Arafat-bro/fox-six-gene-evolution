"""
Quick coordinate-based schematic of SIX7/SIX10/SIX12 gene arrangement
across four genomes carrying >=2 of these genes. Not a full clinker
synteny diagram (no flanking-gene homology shading) - this is a
to-scale positional/orientation schematic built directly from BLAST
merged-hit coordinates, intended as an honest, time-boxed substitute.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# genome_label, contig, [(gene, start, end, strand)]
data = [
    ("Fol4287 (lycopersici)\nQESU01000255.1", [
        ("SIX10", 5729, 6245, "-"),
        ("SIX7", 10211, 10870, "+"),
    ]),
    ("D11 (lycopersici)\nCM012196.1", [
        ("SIX7", 363265, 363924, "-"),
        ("SIX10", 367890, 368406, "+"),
    ]),
    ("f. sp. lini\nCM023981.1 (linked copies only)", [
        ("SIX7_copy2", 304130, 304771, "-"),
        ("SIX10_copy1", 308004, 308520, "+"),
    ]),
    ("f. sp. cepae\nMRCU01000021.1", [
        ("SIX12", 244625, 244984, "+"),
        ("SIX10", 248225, 248587, "+"),
    ]),
]

colors = {"SIX7": "#2b6cb0", "SIX7_copy2": "#2b6cb0",
          "SIX10": "#c05621", "SIX10_copy1": "#c05621",
          "SIX12": "#2f855a"}

fig, axes = plt.subplots(len(data), 1, figsize=(9, 6.5), sharex=False)

for ax, (label, genes) in zip(axes, data):
    starts = [g[1] for g in genes]
    ends = [g[2] for g in genes]
    lo, hi = min(starts) - 500, max(ends) + 500
    ax.set_xlim(lo, hi)
    ax.set_ylim(0, 1)
    ax.axhline(0.5, color="grey", lw=1, zorder=1)
    for gene, start, end, strand in genes:
        color = colors.get(gene, "#718096")
        width = end - start
        x = start if strand == "+" else end
        dx = width if strand == "+" else -width
        ax.add_patch(mpatches.FancyArrow(
            x, 0.5, dx, 0, width=0.18, length_includes_head=True,
            head_width=0.32, head_length=min(400, width * 0.3),
            color=color, zorder=2))
        ax.text((start + end) / 2, 0.85, gene, ha="center", fontsize=8)
    ax.set_yticks([])
    ax.set_title(label, fontsize=9, loc="left")
    ax.set_xlabel("genomic coordinate (bp)", fontsize=7)

plt.tight_layout()
plt.savefig("/home/claude/six7_10_12_schematic.png", dpi=180)
print("saved")
