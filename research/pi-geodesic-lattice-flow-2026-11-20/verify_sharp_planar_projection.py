#!/usr/bin/env python3
"""Exact checks for splicing a planar spiral into a Hamiltonian ray."""

from itertools import count, islice, pairwise, product
from math import isqrt


def vertex(planar, tail):
    return planar, tuple(sorted((j, value) for j, value in tail.items() if value))


def targets():
    """Fair enumeration of full targets (plane, canonical nonzero tail)."""
    seen = set()
    for bound in count(1):
        digits = range(-bound, bound + 1)
        for xy in product(digits, repeat=2):
            for dense in product(digits, repeat=bound):
                if any(dense):
                    point = vertex(xy, dict(enumerate(dense, 1)))
                    if point not in seen:
                        seen.add(point)
                        yield point


def planar_walk(start, end):
    """Walk on the first axis and then the second, including both endpoints."""
    x, y = start
    yield x, y
    while x != end[0]:
        x += 1 if end[0] > x else -1
        yield x, y
    while y != end[1]:
        y += 1 if end[1] > y else -1
        yield x, y


def portal_path(radius, target, entry, exit_axis):
    p, q = (radius, 0), (radius, 1)
    a, h = target
    approach = list(planar_walk(p, a))
    tail = {}
    yield vertex(p, tail)
    tail[entry] = 1
    yield vertex(p, tail)
    tail[entry] = 2
    yield vertex(p, tail)
    for xy in approach[1:]:
        yield vertex(xy, tail)
    for j, value in h:
        for k in range(1, abs(value)+1):
            tail[j] = k if value > 0 else -k
            yield vertex(a, tail)
    tail[entry] = 1
    yield vertex(a, tail)
    del tail[entry]
    yield vertex(a, tail)
    tail[exit_axis] = 1
    yield vertex(a, tail)
    for j, value in reversed(h):
        for k in range(abs(value)-1, -1, -1):
            tail[j] = k if value > 0 else -k
            yield vertex(a, tail)
        del tail[j]
    for xy in reversed(approach[:-1]):
        yield vertex(xy, tail)
    yield vertex(q, tail)
    yield vertex(q, {})


def upper_budget(target):
    xy, h = target
    return 2*sum(map(abs, xy)) + 2*sum(abs(v) for _, v in h) + 6


def overhead(radius, target):
    (x, y), h = target
    return 2*(radius-x+abs(y)) + 2*sum(abs(v) for _, v in h) + 6


def first_radius(min_radius, budget, allowance):
    if allowance(min_radius) >= budget:
        return min_radius
    low, high = min_radius, 2*min_radius
    while allowance(high) < budget:
        low, high = high, 2*high
    while high-low > 1:
        mid = (low+high)//2
        if allowance(mid) < budget:
            low = mid
        else:
            high = mid
    return high


def spiral(radius):
    x = y = 0
    yield (x, y), 0
    for r in range(1, radius+1):
        for dx, dy, steps in (
            (1,0,1), (0,1,2*r-1), (-1,0,2*r), (0,-1,2*r), (1,0,2*r)
        ):
            for _ in range(steps):
                x, y = x+dx, y+dy
                yield (x, y), r


def unit_distance(left, right):
    left_xy, left_tail = left
    right_xy, right_tail = right
    difference = sum(abs(a-b) for a,b in zip(left_xy, right_xy))
    tail = dict(left_tail)
    for j,v in right_tail:
        tail[j] = tail.get(j,0)-v
    return difference + sum(map(abs, tail.values()))


def check_full_prefix(num_targets):
    jobs = {}
    nonzero_seen = set()
    radius = elapsed_overhead = skipped = largest_axis = 0
    for target in islice(targets(), num_targets):
        if target in nonzero_seen:
            skipped += 1
            continue
        xy, h = target
        lower = max(radius+1, max(map(abs, xy))+1)
        radius = first_radius(lower, elapsed_overhead+upper_budget(target), lambda r: 100*r)
        largest_axis = max(largest_axis, h[-1][0])
        path = list(portal_path(radius, target, largest_axis+1, largest_axis+2))
        largest_axis += 2
        assert target in path
        mid = set(path[1:-1])
        assert len(mid) == len(path)-2
        assert all(p[1] for p in mid) and not (mid & nonzero_seen)
        assert len(path)-2 == overhead(radius, target)
        for (x,y), _ in mid:
            assert max(abs(x),abs(y)) < radius or (x,y) in ((radius,0), (radius,1))
        nonzero_seen.update(mid)
        jobs[radius] = path
        elapsed_overhead += len(path)-2
        assert elapsed_overhead <= 102*radius
    assert len(jobs)+skipped == num_targets

    started = ((0,0), ())
    seen = {started}
    plane = {(0,0)}
    shell_seen = {0: 1}
    completed = {0: 0}
    time = 0
    cost_at_end = 0
    base = spiral(radius)
    prev,_ = next(base)
    for q,r in base:
        path = jobs.get(r) if prev == (r,0) and q == (r,1) else None
        add = path[1:] if path is not None else [(q,())]
        current = (prev,())
        for point in add:
            assert point not in seen and unit_distance(current,point) == 1
            seen.add(point)
            time += 1
            current = point
            xy = point[0]
            if xy not in plane:
                plane.add(xy)
                ring = max(map(abs,xy))
                shell_seen[ring] = shell_seen.get(ring,0)+1
                if shell_seen[ring] == 8*ring:
                    assert ring-1 in completed
                    completed[ring] = time
        if path is not None:
            cost_at_end += len(path)-2
        if q == (r,-r):
            assert completed[r] == time == 4*r*(r+1)+cost_at_end
        prev = q
    assert nonzero_seen.issubset(seen)
    assert set(islice(targets(), num_targets)).issubset(seen)
    print(f"{num_targets} full targets verified, {skipped} covered before their turn.")
    print(f"radius {radius}: {len(seen)} distinct vertices, exact T={time}.")


def check_sharp_budget():
    # Distinct first-plane coordinates and the same tail e_1 ensure these
    # targets cannot occur in another detour's nonzero-gate part.
    scheduled = []
    total = previous = 0
    for index in range(6):
        target = ((-index, -1), ((1,1),))
        r = first_radius(max(previous+1,index+1,2), total+upper_budget(target), isqrt)
        total += overhead(r,target)
        assert total <= 2*r+isqrt(r)
        scheduled.append((r, total))
        previous = r
    r,s = scheduled[0]
    t = ((0,-1), ((1,1),))
    path = list(portal_path(r,t,2,3))
    assert len(path)-2 == s
    assert all(unit_distance(a,b) == 1 for a,b in pairwise(path))
    assert len(set(path)) == len(path)
    print("sqrt(R) allowance: six exact schedule certificates.")
    for r,s in scheduled:
        print(f"  radius={r}   excess={s}   budget={2*r+isqrt(r)}")


if __name__ == "__main__":
    check_full_prefix(80)
    check_sharp_budget()
