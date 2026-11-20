#!/usr/bin/env python3
"""Exact checks for a Hamiltonian ray through all of Z^d, d >= 2."""

import argparse
from itertools import product
from math import prod


def box_after_cycle(dimension, cycle):
    """Cube of side 2*cycle+1, with a one-unit shift in odd cycles of odd d."""
    shift = -(dimension % 2) * (cycle % 2)
    low = [-cycle] * (dimension - 1) + [-cycle + shift]
    high = [cycle] * (dimension - 1) + [cycle + shift]
    return tuple(low), tuple(high)


def cycle_slabs(dimension, cycle):
    """Yield disjoint slabs in walk order, each with its starting corner."""
    low, high = map(list, box_after_cycle(dimension, cycle - 1))
    parity_sign = -1 if cycle % 2 == 0 else 1
    signs = [1] + [parity_sign] * (dimension - 2)
    signs.append(parity_sign if dimension % 2 else 1)
    yield from expand_round(low, high, signs)
    assert (tuple(low), tuple(high)) == box_after_cycle(dimension, cycle)


def expand_round(low, high, signs):
    """Perform one two-pass expansion of an odd-sided rectangle."""
    dimension = len(low)
    for second_pass in (False, True):
        for axis in range(dimension):
            if signs[axis] > 0:
                high[axis] += 1
                level = high[axis]
            else:
                low[axis] -= 1
                level = low[axis]
            slab_low = low[:]
            slab_high = high[:]
            slab_low[axis] = slab_high[axis] = level
            face_axes = [k for k in range(dimension) if k != axis]
            if not second_pass and axis > 0:
                slow = axis - 1
            elif second_pass and axis < dimension - 1:
                slow = dimension - 1
            else:
                slow = None
            if slow is None:
                order = face_axes
                flips = face_axes
                assert all((high[k] - low[k] + 1) % 2 for k in face_axes)
            else:
                order = [slow] + [k for k in face_axes if k != slow]
                flips = [slow]
                assert (high[slow] - low[slow] + 1) % 2 == 0
            size = prod(slab_high[k] - slab_low[k] + 1 for k in order)
            yield slab_low, slab_high, order, signs[:], size
            for k in flips:
                signs[k] *= -1


def integer_root(value, power):
    low, high = 0, value + 1
    while low + 1 < high:
        middle = (low + high) // 2
        if middle**power <= value:
            low = middle
        else:
            high = middle
    return low


def point_at(dimension, tick):
    assert dimension >= 2 and tick >= 0
    if tick == 0:
        return (0,) * dimension
    cycle = (integer_root(tick, dimension) + 1) // 2
    index = tick - (2*cycle - 1)**dimension
    for low, high, order, signs, size in cycle_slabs(dimension, cycle):
        if index >= size:
            index -= size
            continue
        return slab_point(low, high, order, signs, size, index)
    raise AssertionError("tick outside the cycle")


def slab_point(low, high, order, signs, size, index):
    point = low[:]
    stride = size
    for axis in order:
        stride //= high[axis] - low[axis] + 1
        layer, index = divmod(index, stride)
        point[axis] = high[axis] - layer if signs[axis] > 0 else low[axis] + layer
        if layer % 2:
            index = stride - 1 - index
    return tuple(point)


def in_box(point, low, high):
    return all(a <= x <= b for x, a, b in zip(point, low, high))


def time_of(point):
    dimension = len(point)
    radius = max(map(abs, point))
    if radius == 0:
        return 0
    cycle = next(
        m for m in range(max(1, radius - 1), radius + 2)
        if in_box(point, *box_after_cycle(dimension, m))
    )
    time = (2*cycle - 1)**dimension
    for low, high, order, signs, size in cycle_slabs(dimension, cycle):
        if not in_box(point, low, high):
            time += size
            continue
        stride = 1
        index = 0
        for axis in reversed(order):
            layer = high[axis] - point[axis] if signs[axis] > 0 else point[axis] - low[axis]
            if layer % 2:
                index = stride - 1 - index
            index += layer * stride
            stride *= high[axis] - low[axis] + 1
        return time + index
    raise AssertionError("point outside the cycle")


def hamiltonian_bandwidth_at_most(side, bound):
    """Exhaust all square start orbits using the deadline of each visited vertex."""
    count = side**2
    neighbors = []
    for x, y in product(range(side), repeat=2):
        neighbors.append(sum(
            1 << (a*side + b)
            for a, b in ((x-1,y), (x+1,y), (x,y-1), (x,y+1))
            if 0 <= a < side and 0 <= b < side
        ))
    full = (1 << count) - 1
    path = []
    calls = 0

    def search(used):
        nonlocal calls
        calls += 1
        if used == full:
            return True
        candidates = neighbors[path[-1]] & ~used
        deadline = len(path) - bound
        if deadline >= 0:
            due = neighbors[path[deadline]] & ~used
            if due & (due-1):
                return False
            if due:
                candidates &= due
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            path.append(bit.bit_length()-1)
            if search(used | bit):
                return True
            path.pop()
        return False

    middle = (side-1)//2
    for a in range(middle+1):
        for b in range(a, middle+1):
            start = a*side+b
            path[:] = [start]
            if search(1 << start):
                return True, calls
    return False, calls


def check_ray(dimension, last_cycle):
    total = (2*last_cycle+1)**dimension
    visited = set()
    previous = None
    for tick in range(total):
        point = point_at(dimension, tick)
        assert point not in visited
        assert time_of(point) == tick
        if previous is not None:
            assert sum(abs(a-b) for a,b in zip(previous, point)) == 1
        visited.add(point)
        previous = point
        root = integer_root(tick + 1, dimension)
        if root % 2 == 1 and root**dimension == tick + 1:
            low, high = box_after_cycle(dimension, (root-1)//2)
            assert all(in_box(x, low, high) for x in visited)
    print(f"d={dimension}: {total} vertices, bijection and unit steps verified")


def growing_count(stage):
    return (4*stage+1)**2 * prod(4*j+1 for j in range(1, stage))


def growing_point(tick):
    if tick == 0:
        return (0, 0)
    stage = 1
    while growing_count(stage) <= tick:
        stage += 1
    index = tick - growing_count(stage-1)
    radii = [2*(stage-1)]*2 + [2*j for j in range(stage-2, -1, -1)]
    low, high = [-r for r in radii], radii[:]
    signs = [1]*len(radii)
    for _ in range(2):
        for lo, hi, order, corner, size in expand_round(low, high, signs):
            if index < size:
                return slab_point(lo, hi, order, corner, size, index)
            index -= size
    raise AssertionError("tick outside the stage")


def check_growing(last_stage):
    count = growing_count(last_stage)
    previous = None
    visited = set()
    for tick in range(count):
        point = growing_point(tick)
        point += (0,) * (last_stage+1-len(point))
        assert point not in visited
        assert previous is None or sum(abs(a-b) for a,b in zip(previous, point)) == 1
        visited.add(point)
        previous = point
    radii = [2*last_stage]*2 + [2*j for j in range(last_stage-1, 0, -1)]
    assert all(all(abs(x) <= r for x,r in zip(point, radii)) for point in visited)
    print(f"growing dimension: {count} vertices, no repeats and unit steps verified")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-square', type=int, default=8)
    args = parser.parse_args()
    for dimension, cycles in ((2, 20), (3, 6), (4, 4), (5, 2), (6, 2)):
        check_ray(dimension, cycles)
    # Hamiltonian path labels respect the checkerboard bipartition.
    # A failure at the previous odd value proves the snake's upper bound sharp.
    check_growing(3)
    for side in range(2, args.max_square+1):
        exists, calls = hamiltonian_bandwidth_at_most(side, 2*side-3)
        assert not exists
        print(f"{side} by {side}: exact Hamiltonian bandwidth {2*side-1}, {calls} search nodes")


if __name__ == "__main__":
    main()
