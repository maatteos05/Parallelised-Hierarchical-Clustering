import argparse
import pandas as pd
from pathlib import Path

class TreeNode:
    def __init__(self, cluster_id, left=None, right=None, dist=0.0, leaves=None):
        self.cluster_id = cluster_id
        self.left = left
        self.right = right
        self.dist = dist
        self.leaves = leaves if leaves is not None else {cluster_id}

def build_tree_from_dendrogram(df, N):
    # Map from cluster ID to Node
    nodes = {i: TreeNode(i) for i in range(N)}
    
    # List of merge events
    merge_events = []
    
    current_id = N
    for index, row in df.iterrows():
        cl1 = int(row['cl1'])
        cl2 = int(row['cl2'])
        dist = float(row['dist'])
        
        left_node = nodes[cl1]
        right_node = nodes[cl2]
        
        merged_leaves = left_node.leaves | right_node.leaves
        new_node = TreeNode(current_id, left_node, right_node, dist, merged_leaves)
        
        nodes[current_id] = new_node
        merge_events.append((current_id, cl1, cl2, dist, new_node))
        
        current_id += 1
        
    root = nodes[current_id - 1] if N > 1 else nodes[0]
    return root, merge_events

def get_leaf_order(node):
    if node.left is None and node.right is None:
        return [node.cluster_id]
    return get_leaf_order(node.left) + get_leaf_order(node.right)

def draw_ascii_dendrogram(root, merge_events, N):
    """Draws a vertical ASCII-art dendrogram tree with level/distance labels on the left."""
    if N < 2:
        return "Not enough points to draw a dendrogram.\n"

    # 1. Determine leaf order to avoid crossing lines
    leaf_order = get_leaf_order(root)
    
    # 2. Configure sizes and spacing dynamically based on size
    col_spacing = 4 if N > 15 else 6
    left_margin = 32
    
    # Column mapping for each leaf
    col_map = {leaf_id: i * col_spacing + left_margin for i, leaf_id in enumerate(leaf_order)}
    
    # Row assignments:
    # y = 0: top margin / root stem
    # Each merge level L (from 1 to N-1) gets a dedicated horizontal row.
    # We leave 1 empty row between merges for vertical stems.
    # So y_merge(L) = L * 2
    # y_leaves = (N) * 2
    # y_labels = y_leaves + 1, y_leaves + 2
    grid_height = N * 2 + 3
    grid_width = left_margin + N * col_spacing + 5
    
    y_leaves = grid_height - 3
    
    # Initialize the character grid
    grid = [[" " for _ in range(grid_width)] for _ in range(grid_height)]
    
    # x and y coordinates map for each cluster ID
    x_map = {}
    y_map = {}
    
    # Initialize leaves on the grid
    for leaf in leaf_order:
        x = col_map[leaf]
        grid[y_leaves][x] = "○"
        # Center the leaf label
        lbl = str(leaf)
        for offset, char in enumerate(lbl):
            lx = x - len(lbl) // 2 + offset
            if 0 <= lx < grid_width:
                grid[y_leaves + 1][lx] = char
        x_map[leaf] = x
        y_map[leaf] = y_leaves
        
    # 3. Draw each merge step-by-step from Level 1 (bottom) to Level N-1 (top)
    for level_idx, (node_id, cl1, cl2, dist, node) in enumerate(merge_events):
        level_num = level_idx + 1
        
        # Row index for this merge level
        # Level 1 is at the bottom (just above leaves), Level N-1 is at the top.
        y_merge = y_leaves - level_num * 2
        
        x_left = x_map[cl1]
        x_right = x_map[cl2]
        
        # Keep them left-to-right
        if x_left > x_right:
            x_left, x_right = x_right, x_left
            
        x_mid = (x_left + x_right) // 2
        
        # Save node coordinates
        x_map[node_id] = x_mid
        y_map[node_id] = y_merge
        
        # Draw level labels on the left margin
        lbl_str = f"LEVEL {level_num:<2d} (Dist: {dist:.6f})"
        for idx, char in enumerate(lbl_str):
            grid[y_merge][idx] = char
        # Draw arrow pointing to the merge bar
        arrow = " ──────> "
        for idx, char in enumerate(arrow):
            grid[y_merge][len(lbl_str) + idx] = char
            
        # Draw horizontal merge bar
        grid[y_merge][x_left] = "┌"
        grid[y_merge][x_right] = "┐"
        for x in range(x_left + 1, x_right):
            grid[y_merge][x] = "─"
        grid[y_merge][x_mid] = "┴"
        
        # Draw vertical drop stems from this merge bar down to the children
        y_left_child = y_map[cl1]
        for y in range(y_merge + 1, y_left_child):
            grid[y][x_left] = "│"
            
        y_right_child = y_map[cl2]
        for y in range(y_merge + 1, y_right_child):
            grid[y][x_right] = "│"
            
        # If this is the root, draw the top stem going up
        if level_num == N - 1:
            for y in range(0, y_merge):
                grid[y][x_mid] = "│"
                
    # Convert grid to string
    lines = []
    for row in grid:
        lines.append("".join(row).rstrip())
    return "\n".join(lines) + "\n"

def main():
    parser = argparse.ArgumentParser(description="Visualize a dendrogram CSV as a text tree structure.")
    parser.add_argument("input_csv", type=str, help="Path to the dendrogram CSV file")
    parser.add_argument("output_txt", type=str, help="Path to the output visualization TXT file")
    args = parser.parse_args()
    
    csv_path = Path(args.input_csv)
    if not csv_path.exists():
        print(f"Error: Input file {csv_path} does not exist.")
        return
        
    # Read the dendrogram CSV
    df = pd.read_csv(csv_path)
    
    # Calculate N: N merges is N-1 rows in the CSV, so N = len(df) + 1
    N = len(df) + 1
    
    root, merge_events = build_tree_from_dendrogram(df, N)
    
    # Draw the beautiful ASCII dendrogram
    dendrogram_art = draw_ascii_dendrogram(root, merge_events, N)
    
    # Place all visualization text files in the 'data/visual' directory
    repo_root = Path(__file__).resolve().parent.parent
    out_name = Path(args.output_txt).name
    out_path = repo_root / "data" / "visual" / out_name
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("========================================================================\n")
        f.write(f"DENDROGRAM PLOT FOR {csv_path.name} (N={N} points)\n")
        f.write("========================================================================\n\n")
        
        f.write(dendrogram_art)
        
    print(f"Dendrogram visualization successfully written to {out_path}")

if __name__ == "__main__":
    main()
