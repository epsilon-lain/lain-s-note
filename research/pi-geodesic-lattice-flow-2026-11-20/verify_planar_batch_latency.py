#!/usr/bin/env python3
"""Exact finite checks for batching a fixed reference ray into a planar spiral."""

from bisect import bisect_right

from verify_sharp_planar_projection import planar_walk, spiral, unit_distance, vertex
from verify_single_line_hamiltonian import growing_point


def reference_points(count):
    points = []
    for index in range(count):
        x, y, *tail = growing_point(index)
        points.append(vertex((x,y), dict(enumerate(tail,1))))
    return points


def link(start, end, marker, enter_outer=False):
    """Join with a protected fresh axis. Yield all vertices after start."""
    start_xy, start_tail = start
    end_xy, end_tail = end
    tail = dict(start_tail)
    tail[marker] = 1
    yield vertex(start_xy, tail)

    if enter_outer:
        # Change y before x so the only new point of this ring is the portal.
        planar = [(y,x) for x,y in planar_walk(start_xy[::-1],end_xy[::-1])]
    else:
        planar = list(planar_walk(start_xy,end_xy))
    for xy in planar[1:]:
        yield vertex(xy, tail)

    end_values = dict(end_tail)
    for axis in sorted(dict(start_tail).keys() | end_values.keys()):
        old = tail.get(axis,0)
        new = end_values.get(axis,0)
        step = 1 if new > old else -1
        for value in range(old+step,new+step,step):
            tail[axis] = value
            yield vertex(end_xy, tail)
    yield end


def coverage_check(radii):
    cap = radii[-2]
    targets = reference_points(cap)
    jobs = dict(zip(radii[1:], radii[:-1]))
    max_axis = max((j for _,tail in targets for j,_ in tail),default=0)
    initial = ((0,0),())
    when = {initial: 0}
    first_projected = {(0,0)}
    counts = {0:1}
    coverage = [0]
    delays = 0
    time = 0
    current = initial
    base = spiral(radii[-1])
    prev,_ = next(base)

    def record(point):
        nonlocal current,time
        assert unit_distance(current,point) == 1 and point not in when
        current = point
        time += 1
        when[point] = time
        xy = point[0]
        if xy not in first_projected:
            first_projected.add(xy)
            k = max(map(abs,xy))
            counts[k] = counts.get(k,0)+1
            if counts[k] == 8*k:
                assert len(coverage) == k
                coverage.append(time)

    for q,r in base:
        assert current == (prev,())
        if r in jobs and prev == (r,0) and q == (r,1):
            goal = jobs[r]
            outstanding = [point for point in targets[:goal] if point not in when]
            assert all(tail and max(map(abs,xy)) < r for xy,tail in outstanding)
            dest = (q,())
            start_time = time
            if outstanding:
                for end in outstanding + [dest]:
                    max_axis += 1
                    for point in link(current,end,max_axis,enter_outer=(end==dest)):
                        assert (max(map(abs,point[0])) < r
                                or point[0] in ((r,0),(r,1)))
                        record(point)
            else:
                record(dest)
            extra = time-start_time-1
            assert extra <= 2*r+4*goal+2
            delays += extra
            assert all(target in when for target in targets[:goal])
            print(f"batch at r={r}: prefix={goal}, new={len(outstanding)}, extra={extra}")
        else:
            record((q,()))
        if q == (r,-r):
            assert time == coverage[r] == 4*r*(r+1)+delays
        prev = q

    footprint = set()
    radius = latest = tested = 0
    for n,point in enumerate(targets,1):
        latest = max(latest,when[point])
        xy,_ = point
        footprint.add(xy)
        radius = max(radius,max(map(abs,xy)))
        r = bisect_right(coverage,latest)-1
        if radius < r < radii[-1]:
            extra = coverage[r+1]-4*(r+1)*(r+2)
            assert extra >= 2*(r-radius)-1+n-len(footprint)
            tested += 1
    print(f"verified {len(when)} distinct unit-step vertices and {tested} lower-bound instances")


def log_schedule(stop):
    radii = [2,4]
    while len(radii) < stop:
        # g(N) = bit_length(N+2)^4 in the general schedule.
        multiplier = (radii[-2]+3).bit_length()
        radii.append(max(2,multiplier)*radii[-1])
    return radii


def verify_schedules():
    log_radii = log_schedule(9)
    for j in range(2,len(log_radii)):
        smallest = log_radii[j-2]+1
        g = (smallest+2).bit_length()**4
        assert log_radii[j]**2 <= 16*smallest**2*g
    power_radii = [2,4]
    for _ in range(6):
        power_radii.append(power_radii[-1]**2)
    for old,r in zip(power_radii,power_radii[2:]):
        assert r*r == old**8
    print("Verified logarithmic and power schedule inequalities.")


if __name__ == "__main__":
    coverage_check([2,84,100,128])
    verify_schedules()
