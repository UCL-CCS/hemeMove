"""Since the line-by-line ordering of snapshots can vary depending on
decomposition, this enforces a canonical ordering.

The ordering is based on the grid positions, z fastest varying.

"""

import numpy as N

def computeGridIdx(snap):
    return (snap.grid[:,0] * snap.bb_len[1] +
            snap.grid[:,1]) * snap.bb_len[2] + \
            snap.grid[:,2]

def order(snap):
    idx = computeGridIdx(snap)
    idxIdx = idx.argsort()
    
    new = snap.copy()

    for i, ind in enumerate(idxIdx):
        new[i] = snap[ind]
        continue
    
    return new

# This creates an array of indices which give an ordering to a 3D vector array.
def computeIncreasingGridIndex(snap):
    # We first move each position to be in a range [0, max-min]
    # Then can order by i * (max-min)^2 + j * (max-min) + k. Sweet!
    # THIS IS A HACK BUT I JUST WANT IT TO WORK. Better than multiplying by 100 would be
    # to do it some other way, i.e. define an ordering, pass that into some sort function.
    min = N.min(snap.position) * 100.0
    max = N.max(snap.position) * 100.0

    return ((snap.position[:,0])  * (max-min) +
            (snap.position[:,1])) * (max-min) + \
            (snap.position[:,2])

# This takes an (i,j,k) indexed array and arranges it in a natural order (for 
# each i (increasing) go through each j (increasing) and for each (i,j) go 
# through each k (increasing).
def continuousOrder(snap):
    idx = computeIncreasingGridIndex(snap)
    idxIdx = idx.argsort()
    
    new = snap.copy()

    for i, ind in enumerate(idxIdx):
        new[i] = snap[ind]
        continue
    
    return new

