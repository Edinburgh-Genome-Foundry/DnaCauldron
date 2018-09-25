"""This is some very experimental code tht shouldn't be here, but hey."""

from itertools import count
from collections import Counter


import numpy as np
import networkx as nx
import tatapov

from .tools import reverse_complement

def weighted_pick(weights, n_picks):
    return np.searchsorted(np.cumsum(weights),
                           np.random.rand(n_picks) * np.sum(weights))

def predict_good_clones_proportion(parts_overhangs, annealing_data):
    overhangs = [
        o
        for pos, o1, o2 in parts_overhangs
        for o in [o1, o2]
        if o not in ["LEFT", "RIGHT"]
    ]
    probas = tatapov.data_subset(annealing_data, overhangs, add_reverse=True)
    precomputed_random = {
        oh:  (p for p in weighted_pick(probas[oh], 15000))
        for oh in probas
    }
    parts_overhangs_dict = {
        p: {
            "left": "LEFT" if left == "LEFT" else reverse_complement(left),
            "right": right
        }
        for p, left, right in parts_overhangs
    }
    free_overhangs = {
        p + "_" + side: {
            "part": p,
            "sequence": parts_overhangs_dict[p][side],
            "other_side": parts_overhangs_dict[p][other_side]
        }
        for p in parts_overhangs_dict
        for (side, other_side) in [("left", "right"), ("right", "left")]
    }
    all_free_overhangs = [
        o for o in sorted(free_overhangs)
        if free_overhangs[o]["sequence"] not in ["LEFT", "RIGHT"]
    ]
    all_free_overhangs_sequences = [
        free_overhangs[o]["sequence"]
        for o in all_free_overhangs
    ]
    probas = probas.loc[all_free_overhangs_sequences]
    
    expected = tuple(pos for pos, o1, o2 in parts_overhangs)

    def simulate_construct(start, end, direction):
        
        construct = [start]
        last_free_overhang = parts_overhangs_dict[start][direction]
        next_part = None
        while next_part not in (start, end):
            next_overhang = all_free_overhangs[
                next(precomputed_random[last_free_overhang])]
            next_part = free_overhangs[next_overhang]["part"]
            next_sequence = free_overhangs[next_overhang]["sequence"]
            last_free_overhang = free_overhangs[next_overhang]["other_side"]
            construct.append(next_part)
        return tuple(construct)

    constructs = [
        c for c in [
            simulate_construct("backbone_left", "backbone_right", "right")
            for i in range(5000)
        ]
        if c[-1] == "backbone_right"
    ]
    constructs += [
        c for c in [
            simulate_construct("backbone_right", "backbone_left", "left")
            for i in range(5000)
        ]
        if c[0][-1] == "backbone_left"
    ]
    cnt = Counter(constructs)
    for cst, c in cnt.items():
        cnt[cst] = 1.0 * c / (len(cst)**3)
        
    return (cnt[expected] + cnt[expected[::-1]]) / sum(list(cnt.values()))