import sys
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as LA
import json
from os.path import basename
import argparse
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import itertools




def plot_unitcell(data):
    A = np.array(data['index_cell_vectors'])
    aLinv = np.array(data['primitive_cell_vectors'])

    for f in get_faces(A):
        ax.plot(*np.array(f, dtype=np.float64).T, color='k')

    for f in get_faces(aLinv):
        ax.plot(*np.array(f, dtype=np.float64).T, color='green')


def get_faces(M):
    faces = []

    for i in range(3):
        vec_a = M[:, i]
        vec_b = M[:, (i+1) % 3]
        vec_c = M[:, (i+2) % 3]
        f = [[0, 0, 0], vec_a, vec_a + vec_b, vec_b, [0, 0, 0]]
        f = [np.array(x, dtype=np.int64) for x in f]
        faces.append(f)
        faces.append([x + vec_c for x in f])
    return faces


def unwrap(dx, A):
    '''
    Searches for the unitcell-equivalent point that makes 'dx' the shortest vector.
    '''
    candidate_DX = dx
    for idx in itertools.product((-1, 0, 1), (-1, 0, 1), (-1, 0, 1)):
        tmp = dx + A@idx
        if LA.norm(tmp) < LA.norm(candidate_DX):
            candidate_DX = tmp
    return candidate_DX

def are_parallel(dx0, dx1, atol=1e-10):
    return LA.norm(np.cross(dx0, dx1)) < atol


def link2startstop(x):
    # expects 'x' to be in link format, i.e.
    # x = {'position', 'boundary': [ [r0, m0], [r1, m1])]}

    if len(x["boundary"]) > 2:
        print("Malformed boundary: link should have two or zero ends")
    r0, m0 = x["boundary"][0]
    r1, m1 = x["boundary"][1]
    assert (m0 + m1 == 0), "Malformed link: multipliers do not sum to 0"
    if m0 == -1:
        # swap
        m0, m1 = m1, m0
        r0, r1 = r1, r0

    # it should now be possible to deduce the correct wrapping of p0, p1,
    # such that all are in same cell
    return np.array(r0), np.array(r1)


def plot_directed_link(x, A):
    r0, r1 = link2startstop(x)

    dx0, dx1 = (r0 - x['pos']), (r1 - x['pos'])

    if not are_parallel(dx0, dx1):
        dx0 = unwrap(dx0, A)
        dx1 = unwrap(dx1, A)

    ax.quiver(*r0, *(-dx0), color='k')
    ax.quiver(*r1, *(-dx1), color='k', arrow_length_ratio=0)


def plot_undirected_link(x, A):
    r0, r1 = link2startstop(x)

    dx0, dx1 = (r0 - x['pos']), (r1 - x['pos'])

    if not are_parallel(dx0, dx1):
        dx0 = unwrap(dx0, A)
        dx1 = unwrap(dx1, A)

    ax.quiver(*r0, *(-dx0), color='k', arrow_length_ratio=0)
    ax.quiver(*r1, *(-dx1), color='k', arrow_length_ratio=0)


def plot_idx(x, i):
    ax.text(*x['pos'], "%d" % i)

def plot_pos(x):
    ax.text(*x['pos'], "%d %d %d" % tuple(l for l in x['pos']))

def plot_points(data, args):
    point_data = data['points']

    if point_data is None:
        print("No points in latfile.")
        return

    xyz = []
    for i, x in enumerate(point_data):
        xyz.append(x['pos'])

        if args.show_idx:
            plot_idx(x, i)
        if args.show_pos:
            plot_pos(x)

    ax.scatter(*np.array(xyz).T, color='r', marker='o')


def plot_links(data, args):
    A = np.array(data['index_cell_vectors'])
    link_data = data['links']

    if link_data is None:
        print("No links in latfile.")
        return

    for i, x in enumerate(link_data):
        if args.undirected:
            plot_undirected_link(x,A)
        else:
            plot_directed_link(x, A)

        if args.show_idx:
            plot_idx(x, i)
        if args.show_pos:
            plot_pos(x)


def find_link(linkpos, link_data):
    for link in link_data:
        if link['pos'] == linkpos:
            return link
    raise LookupError(f"No link at {linkpos}")


def plot_plaqs(data, args):
    link_data = data["links"]
    plaq_data = data["plaqs"]
    
    if plaq_data is None:
        print("No plaquettes in latfile.")
        return
    A = np.array(data['index_cell_vectors'])

    for i, x in enumerate(plaq_data):
        # Plot plaquette centroid
        pos = np.array(x['pos'])

        triangles = []
        for link_pos, m in x['boundary']:
            link = find_link(link_pos, link_data)
            r0, r1 = link2startstop(link)
            # decide if we crossed a boundary
            dx0 = unwrap(r0 - pos, A)
            dx1 = unwrap(r1 - pos, A)

            triangles.append([pos, pos+dx0, pos+dx1])

        poly = Poly3DCollection(triangles, alpha=0.5)
        poly.set_facecolor('cyan')  # Set surface color
        ax.add_collection3d(poly)

        if args.show_idx:
            plot_idx(x, i)
        if args.show_pos:
            plot_pos(x)


def plot_vols(data, args):
    vol_data = data['vols']
    if vol_data is None:
        print("No volumes in latfile.")
        return
    xyz = []
    for i, x in enumerate(vol_data):
        xyz.append(x['pos'])

        if args.show_idx:
            plot_idx(x, i)
        if args.show_pos:
            plot_pos(x)

    ax.scatter(*np.array(xyz).T, color='b', marker='o')


func_to_run = {
        'points': plot_points,
        'links': plot_links,
        'plaqs': plot_plaqs,
        'vols': plot_vols
        }

ap = argparse.ArgumentParser()
ap.add_argument("file", help="The .lat.json file specifying the lattice",
                type=str)
ap.add_argument("objects", nargs='+', choices=func_to_run.keys())
ap.add_argument("--show_idx", action='store_true')
ap.add_argument("--show_pos", action='store_true')
ap.add_argument("--undirected", action='store_true')
ap.add_argument("--save")

args = ap.parse_args()

data = None
with open(args.file, 'r') as f:
    data = json.load(f)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

for arg in args.objects:
    func_to_run[arg](data, args)

def parse_filename(fname):
    retval = {}
    for tok in basename(fname).split(';'):
        t = tok.split('=')
        if len(t) != 2:
            continue
        retval[t[0]] = t[1]
    return retval

opts = parse_filename(args.file)
try:
    ax.set_title(f'{float(opts['p'])*100}% dilution, removed {opts['nn']}-neighbours')
except KeyError:
    ax.set_title(basename(args.file))

if args.save:
    print("Saving to "+args.save)
    fig.savefig(args.save)
else:
    plt.show()
