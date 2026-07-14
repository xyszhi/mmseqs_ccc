"""
Sankey diagram for KaMinPar partition results across selected k values.
Selected k values: k10, k30, k60, k100
"""

import os
import collections
import plotly.graph_objects as go

BASE_DIR = "/Users/zhixy/CLionProjects/mmseqs_ccc/kaminpar"
SELECTED_K = ["k10", "k20", "k30", "k40", "k50", "k60", "k70", "k80", "k90", "k100", "k110", "k120", "k130", "k140"]
HEADER_LINES = 7  # lines to skip in partition.txt
TOP_N = 200  # keep top N blocks per k, merge rest as "other"

# ── 1. Load partition data ──────────────────────────────────────────────────
def load_partition(k_name):
    path = os.path.join(BASE_DIR, k_name, "partition.txt")
    assignments = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i < HEADER_LINES:
                continue
            assignments.append(int(line.strip()))
    return assignments

print("Loading partitions...")
partitions = {k: load_partition(k) for k in SELECTED_K}
n_nodes = len(partitions[SELECTED_K[0]])
print(f"  {n_nodes} nodes loaded for each k.")

# ── 2. For each k, determine top-N blocks; map others to "other" ────────────
def get_top_blocks(assignments, top_n):
    counter = collections.Counter(assignments)
    top = [block for block, _ in counter.most_common(top_n)]
    return set(top)

top_blocks = {k: get_top_blocks(partitions[k], TOP_N) for k in SELECTED_K}

def mapped_label(k, block):
    if block in top_blocks[k]:
        return f"{k}_B{block}"
    return f"{k}_other"

# ── 3. Build node list and index ────────────────────────────────────────────
all_node_labels = []
node_index = {}

for k in SELECTED_K:
    counter = collections.Counter(partitions[k])
    # top blocks sorted by size descending
    top_sorted = sorted(top_blocks[k], key=lambda b: -counter[b])
    for b in top_sorted:
        lbl = f"{k}_B{b}"
        if lbl not in node_index:
            node_index[lbl] = len(all_node_labels)
            all_node_labels.append(lbl)
    # other node
    other_lbl = f"{k}_other"
    if other_lbl not in node_index:
        node_index[other_lbl] = len(all_node_labels)
        all_node_labels.append(other_lbl)

# ── 4. Pre-compute block sizes ──────────────────────────────────────────────
counter_all = {k: collections.Counter(partitions[k]) for k in SELECTED_K}

# ── 5. Build flows between adjacent k values ────────────────────────────────
print("Computing flows...")
source_list = []
target_list = []
value_list  = []

for i in range(len(SELECTED_K) - 1):
    k_src = SELECTED_K[i]
    k_tgt = SELECTED_K[i + 1]
    flow = collections.Counter()
    for node_id in range(n_nodes):
        src_lbl = mapped_label(k_src, partitions[k_src][node_id])
        tgt_lbl = mapped_label(k_tgt, partitions[k_tgt][node_id])
        flow[(src_lbl, tgt_lbl)] += 1
    for (s, t), v in flow.items():
        if v < 20000:
            continue
        # also skip if source or target block is too small
        src_block = s.split("_B")[-1] if "_B" in s else None
        tgt_block = t.split("_B")[-1] if "_B" in t else None
        if src_block and counter_all[k_src][int(src_block)] < 20000:
            continue
        if tgt_block and counter_all[k_tgt][int(tgt_block)] < 20000:
            continue
        source_list.append(node_index[s])
        target_list.append(node_index[t])
        value_list.append(v)

# ── 5. Assign colors ────────────────────────────────────────────────────────
import colorsys

def k_color_palette(k_name, n):
    """Generate n distinct colors for a given k layer."""
    hue_base = {"k10": 0.0, "k20": 0.07, "k30": 0.14, "k40": 0.21, "k50": 0.28, "k60": 0.35, "k70": 0.45, "k80": 0.52, "k90": 0.59, "k100": 0.65, "k110": 0.70, "k120": 0.75, "k130": 0.82, "k140": 0.88}
    h = hue_base.get(k_name, 0.0)
    colors = []
    for i in range(n):
        lightness = 0.45 + 0.25 * (i / max(n - 1, 1))
        r, g, b = colorsys.hls_to_rgb(h, lightness, 0.6)
        colors.append(f"rgba({int(r*255)},{int(g*255)},{int(b*255)},0.8)")
    return colors

node_colors = ["rgba(180,180,180,0.5)"] * len(all_node_labels)
for k in SELECTED_K:
    counter = collections.Counter(partitions[k])
    top_sorted = sorted(top_blocks[k], key=lambda b: -counter[b])
    palette = k_color_palette(k, len(top_sorted))
    for idx, b in enumerate(top_sorted):
        lbl = f"{k}_B{b}"
        node_colors[node_index[lbl]] = palette[idx]

# ── 6. Build display labels ─────────────────────────────────────────────────
display_labels = []
for lbl in all_node_labels:
    k = lbl.split("_")[0]
    if lbl.endswith("_other"):
        total = sum(v for b, v in counter_all[k].items() if b not in top_blocks[k])
        display_labels.append(f"other\n({total:,})")
    else:
        block = int(lbl.split("_B")[1])
        size = counter_all[k][block]
        display_labels.append(f"B{block}\n({size:,})")

# ── 7. Remove isolated nodes (not connected to any link) ────────────────────
connected = set(source_list) | set(target_list)
# Build remapping: old index -> new index
old_to_new = {}
new_labels = []
new_colors = []
for old_idx, lbl in enumerate(all_node_labels):
    if old_idx in connected:
        old_to_new[old_idx] = len(new_labels)
        new_labels.append(display_labels[old_idx])
        new_colors.append(node_colors[old_idx])

source_list = [old_to_new[s] for s in source_list]
target_list = [old_to_new[t] for t in target_list]
display_labels = new_labels
node_colors = new_colors

# ── 8. Plot ─────────────────────────────────────────────────────────────────
print("Plotting...")
fig = go.Figure(go.Sankey(
    arrangement="perpendicular",
    node=dict(
        pad=12,
        thickness=18,
        line=dict(color="white", width=0.5),
        label=display_labels,
        color=node_colors,
    ),
    link=dict(
        source=source_list,
        target=target_list,
        value=value_list,
        color=[node_colors[s].replace(",0.8)", ",0.3)").replace(",0.5)", ",0.2)") for s in source_list],
    )
))

# Add k-value annotations as column headers
for i, k in enumerate(SELECTED_K):
    fig.add_annotation(
        x=i / (len(SELECTED_K) - 1),
        y=1.05,
        xref="paper", yref="paper",
        text=f"<b>{k}</b>",
        showarrow=False,
        font=dict(size=14),
    )

fig.update_layout(
    title_text="KaMinPar Partition Sankey: k10 → k20 → ... → k140 (all k values)",
    title_font_size=16,
    font_size=10,
    width=5000,
    height=2400,
)

out_path = "/Users/zhixy/CLionProjects/mmseqs_ccc/sankey_all_k_thresh20000_clean.png"
fig.write_image(out_path, scale=2)
print(f"Saved: {out_path}")
