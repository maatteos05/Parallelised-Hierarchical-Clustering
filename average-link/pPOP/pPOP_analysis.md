# Deep analysis: our nested pPOP (average-link) vs. Dash et al. 2007

This maps our final implementation onto the paper's algorithm (Table 2) and
complexity model (Table 4, priority-queue type), states exactly where we follow
the paper and where the average-link choice forces a documented deviation, and
explains the benchmark numbers in those terms.

Final benchmark (8-core MacBook, speedup vs. sequential baseline):

```
pPOP:  N=500  N=1000  N=2000  N=5000
 t=1   1.16x  1.93x   4.01x   7.65x
 t=2   1.08x  2.02x   4.78x   9.22x
 t=4   0.99x  2.03x   5.02x  10.24x
 t=8   0.81x  2.21x   6.48x  17.42x
```

---

## 1. Algorithm-level correspondence (paper Table 2 vs. ours)

The paper's nested pPOP iteration has these steps. Ours matches step for step:

| Paper Table 2 step | Our implementation | Match |
|---|---|---|
| 1-2. Each P creates priority queues for its chunk of cells Ccell(P), in parallel | `build_heaps_thread` dispatched on the pool; threads claim cells via atomic `next_cell` | yes, parallel + dynamic |
| 3. Repeat while overall_min_dist < delta | inner `while (n_active>1)` with break on `best_dist >= delta` | yes |
| 4-5. Each P finds closest pair for each cell in its chunk | `find_min_thread`, atomic cell counter, dynamic scheduling | yes, parallel + dynamic |
| 6-7. Designated P finds overall closest pair + its cell | sequential reduction over `results[]` | yes |
| 8. Merge CL1, CL2 | create `k`, set `sz[k]`, weighted `rep_x/rep_y` | yes |
| 9-10. Each P updates priority queues of Cclus(C) | `update_heaps_thread` on the pool over `affected_clusters` | yes, parallel |
| 11-12. Designated P determines affected neighbouring cells | folded into `affected_clusters` build (scan k's cells) | yes (sequential, cheap) |
| 13-15. Each P updates priority queues of affected cells | same `update_heaps_thread` dispatch | yes, parallel |

Partitioning scheme matches POP: a c x c grid over the 2-D bounding box, with
each cell expanded by delta/2 on every side so adjacent cells overlap by exactly
delta. This is the POP invariant that guarantees correctness: any two clusters
closer than delta share at least one cell, so the merge that would join them is
never missed. We verified this invariant empirically -- every dendrogram matches
scipy's linkage(method='average') to 1e-5 across N=50..800 and t=1,4,8.

The nested control matches Section 2: delta starts at 0 and grows x1.1 per outer
iteration (D_INCR=0.1); c shrinks per iteration (C_DECR=0.1), the 90-10 rule.
The final iteration (small c, large delta) plays the role of the paper's Phase 2,
merging whatever remains.

---

## 2. The one fundamental deviation: average-link vs. centroid

The paper's pPOP is built on the centroid method. We use average-link (UPGMA).
This single choice produces every difference between our constant factors and
the paper's, and it is worth stating precisely.

Centroid: a merged cluster has a representative POINT (its centroid). The
distance from any other cluster to the merged cluster depends only on that one
point. Consequence in the paper (Table 2 step 10): after a merge you only have
to update the priority queues of clusters in the container cell -- O(n/c) work --
because a cluster in a different cell is, by the delta-overlap construction,
farther than delta and thus irrelevant this iteration.

Average-link: there is no representative point. The distance D[k][m] between the
merged cluster k and another cluster m is the Lance-Williams combination
D[k][m] = (sz_i*D[i][m] + sz_j*D[j][m]) / (sz_i+sz_j), which references the old
pairwise distances to BOTH parents. To keep these correct for later iterations
(when the grid is rebuilt and k may be compared against clusters that are in
other cells now), D[k][m] must be maintained globally, for every active m --
O(n) per merge, not O(n/c).

So we keep the paper's O(N/(c*p)) advantage on find-min (Table 4 step 3), but on
the update step (Table 4 step 4) we are a factor of c heavier than the paper:
our global Lance-Williams update is O(n/p) where their container-only centroid
update is O((n/c)*log n / p). We parallelized the update so it is still O(n/p)
rather than O(n) -- but we cannot recover the O(c) reduction, because that
reduction is a property of the centroid method, not of the parallelization.

We followed the paper's parallelization strategy faithfully; the average-link 
distance measure costs us the O(c) reduction factor on exactly one of the two 
critical steps.

---

## 3. Step-by-step complexity: paper vs. ours (priority-queue type, Table 4)

Per merge, n = remaining clusters, c = cells, p = processors:

| Step | Paper pPOP (centroid) | Ours (average-link) | Same? |
|---|---|---|---|
| Create PQs (per outer iter) | O(N^2/(c*p)) | O(N^2/(c*p)) | yes |
| find closest pair | O(n/p) | O(n/p) | yes |
| merge | O(1) | O(1) | yes |
| update PQs / distances | O((n/c)*log n / p) | O((n/p)*log n) global | NO -- factor c |

So our overall per-iteration cost carries the paper's O(c) reduction on
create-PQ and find-min, and loses it only on the update. The paper itself says
(Table 4 discussion) that Step 1 (create PQs) is where most wall-clock time
goes -- "the time taken by Step 1 is substantial" -- and the RF for Step 1 is
O(c). We keep that one. That is why we still see strong speedup despite the
heavier update.

---

## 4. What actually moved the benchmark (and why it matches the theory)

Four implementation facts, in the order we found them, each tied to the model:

1. Global parallel Lance-Williams update. Necessary for average-link
   correctness (Section 2). Parallelized so it is O(n/p), matching the paper's
   step-4 parallelization even though not its O(c) factor.

2. Thread pool (created once). The paper's model assumes p persistent
   processors with negligible synchronization ("assuming no synchronization
   delays"). Per-merge std::thread creation violated that assumption badly;
   the pool restores it. This is an implementation requirement to make the
   paper's complexity model actually apply on a real machine.

3. O(1) heap clear at rebuild (swap-to-empty). This was the single biggest
   real-world cost and is invisible to the paper's asymptotic model: popping
   a priority queue empty element-by-element is O(size*log size), and we did it
   for every active cluster every outer iteration. Swapping with a fresh empty
   heap is O(1). This is a constant-factor fix the paper never has to mention
   because its model counts the rebuild as one O(N^2/(c*p)) construction, not
   the teardown. N=2000 dropped ~7100ms -> ~2745ms from this one line.

4. Cell-count decay floor -- the load-balance fix. This is the paper's own
   warning made concrete. The paper assumes "equal load distribution, each
   processor responsible for c/p cells." Our naive x0.9 decay with integer
   truncation collapsed c from 14 to 1 in three iterations, after which a single
   cell held ~95% of all merges and exactly one thread did all the work -- the
   degenerate load imbalance the paper flags on p.10-12. Flooring c so that
   total_cells stays ~2p during the heavy-merge phase (while capping by the
   cluster count so we don't force more cells than clusters) restores the
   "c/p cells per processor" assumption. This is what turned flat scaling
   (t1=t8=14000ms at N=5000) into 2848->1251ms.

Note that fixes 3 and 4 are NOT visible in the paper's complexity tables: one is
a constant-factor teardown cost, the other is the gap between the model's
"assume equal load distribution" and what a small-N nested run actually does.
Both are legitimate, citable implementation findings.

---

## 5. Reading the speedup numbers through the paper's lens

- Super-linear vs. sequential (17.4x on 8 threads at N=5000). Not anomalous.
  Two independent factors multiply: (a) the partitioning reduces total work by
  the O(c) reduction factor on create-PQ and find-min -- this shows up already
  at t=1 as 7.65x, i.e. single-thread pPOP is 7.65x the sequential baseline
  purely from doing less work; (b) the 8 threads then divide that reduced work,
  adding the further 7.65->17.4 factor. Speedup-vs-sequential conflates
  algorithmic gain with parallel gain; that product is expected and matches the
  paper's claim that the O(c) reduction is "independent of the number of
  processors."

- Scaling improves with N (t8: 0.81x -> 2.21x -> 6.48x -> 17.42x as N grows
  500->5000). Exactly the paper's central experimental finding: speedup is
  small for small datasets and grows with N, because larger N keeps more
  clusters per cell during more of the run, so the "equal load distribution"
  assumption holds longer before the endgame collapses to few cells.

- Small-N still weak / sub-1 at t=8, N=500 (0.81x). Also the paper's finding:
  at small N the per-iteration parallel work is tiny relative to synchronization
  and the single-cell endgame dominates, so adding processors can reduce
  speedup. The paper observed the same drop at 3K-5K on their hardware and
  attributed it to processors idling on cells with few clusters.

---

## 6. Honest limitations to state in the report

- We could only benchmark to N=5000 on a laptop: the stored matrix is
  O((2N)^2) doubles ~ 800 MB at N=5000 and grows quadratically, so N=30K (the
  regime where the paper reports near-linear speedup) needs ~28 GB and a
  cluster. The paper hit the same wall -- they could only run the stored-matrix
  baseline to 5K for the same memory reason.

- Our baseline is a lean flat-matrix scan (naive), not the paper's pTRAD
  (priority-queue per cluster). pTRAD carries heap overhead comparable to pPOP,
  so the paper's relative speedup (pPOP/pTRAD) starts >1 even at small N. Our
  naive baseline is lighter, so the crossover where pPOP overtakes it sits at
  larger N. We report speedup vs. sequential and vs. naive separately to keep
  this distinction clear.

- The O(c) update reduction is permanently unavailable to average-link. A
  centroid implementation would update only the container cell and would be
  faster by that factor, at the cost of using a different (centroid) linkage.
  This is a measure-vs-performance trade, not a bug.