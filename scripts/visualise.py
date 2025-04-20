import sys
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as LA
import json
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import itertools

if len(sys.argv) < 3:
    print(f"USAGE: {sys.argv[0]} <dump.json> [points|links|plaqs|vols]")

assert (sys.argv[2] in ["points", "links", "plaqs", "vols"])

data = None
with open(sys.argv[1], 'r') as f:
    data = json.load(f)


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# plot the unit cell
A = np.array(data['index_cell_vectors'])
aLinv = np.array(data['primitive_cell_vectors'])


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


for f in get_faces(A):
    ax.plot(*np.array(f, dtype=np.float64).T, color='k')


for f in get_faces(aLinv):
    ax.plot(*np.array(f, dtype=np.float64).T, color='green')


def plot_points(data):
    point_data = data['points']
    xyz = []
    for x in point_data:
        xyz.append(x['pos'])

    ax.scatter(*np.array(xyz).T, color='r', marker='o')




def plot_link(x):

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
    dx0 = (np.array(r0) - x['pos'])
    dx1 = (np.array(r1) - x['pos'])

    n = np.cross(dx0, dx1)
    if LA.norm(n) > 1e-10:
        candidate_X1 = dx1
        candidate_X0 = dx0
        for idx in itertools.product((-1, 0, 1), (-1, 0, 1), (-1, 0, 1)):
            tmp = dx1 + A@idx
            if LA.norm(tmp) < LA.norm(candidate_X1):
                candidate_X1 = tmp

            tmp = dx0 + A@idx
            if LA.norm(tmp) < LA.norm(candidate_X0):
                candidate_X0 = tmp

        dx0 = candidate_X0
        dx1 = candidate_X1

    ax.quiver(*r0, *(-dx0), color='k')
    ax.quiver(*r1, *(-dx1), color='k', arrow_length_ratio=0)

def plot_links(data):
    link_data = data['links']
    for x in link_data:
        plot_link(x)


def find_link(linkpos, link_data):
    for l in link_data:
        if l['pos'] == linkpos:
            return l
    raise LookupError(f"No link at {linkpos}")

def plot_plaqs(data):
    link_data = data["links"]
    plaq_data = data["plaqs"]
    
    # Color mapping for orientation
    color_map = {1: 'blue', -1: 'red'}
    
    for plaq in plaq_data:
        # Plot plaquette centroid
        pos = plaq['pos']






func_to_run = {
        'points': plot_points,
        'links': plot_links,
        'plaqs': plot_plaqs
        }

for arg in sys.argv[2:]:
    func_to_run[arg](data)


plt.show()
