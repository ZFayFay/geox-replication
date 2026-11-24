"""Draw figures/cover.png, the cover image of the blog post. Run from the repo root."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BG, INK, MUTED = "#10151f", "#f2f2f0", "#9aa3b2"


def cover(path, title, rows, note):
    fig = plt.figure(figsize=(8, 4.5), dpi=200, facecolor=BG)
    fig.text(0.06, 0.86, title, color=INK, fontsize=19, fontweight="bold", va="top")
    ax = fig.add_axes([0.06, 0.16, 0.88, 0.52], facecolor=BG)
    top = max(v for _, v, _, _ in rows)
    for i, (label, value, text, color) in enumerate(rows):
        y = (len(rows) - 1 - i) * 1.5
        ax.barh(y, value, height=0.6, color=color)
        ax.text(
            value + top * 0.015, y, text, va="center", color=INK, fontsize=17, fontweight="bold"
        )
        ax.text(0, y + 0.38, label, va="bottom", color=MUTED, fontsize=11)
    ax.set_xlim(0, top * 1.22)
    ax.set_ylim(-0.45, (len(rows) - 1) * 1.5 + 0.9)
    ax.axis("off")
    fig.text(0.06, 0.06, note, color=MUTED, fontsize=9)
    fig.savefig(path, dpi=200, facecolor=BG)


cover(
    "figures/cover.png",
    "Experiments where nothing happened,\nflagged as a lift",
    [
        ("Target at a 10% level", 10, "10%", "#6b7280"),
        ("My rebuild of GeoX, Iowa sales", 22, "22%", "#4c8fe0"),
        ("Google's library, clean synthetic data", 28, "28%", "#e0884c"),
    ],
    "Testing Google's Meridian GeoX on ten years of public Iowa sales",
)
