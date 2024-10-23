import sys
import matplotlib.pyplot as plt
import numpy as np
import json

if len(sys.argv) < 3:
    print(f"USAGE: {sys.argv[0]} <dump.json> [points|links|plaqs|vols]")

assert (sys.argv[2] in ["points", "links", "plaqs", "vols"])

data = None
with open(sys.argv[1], 'r') as f:
    data = json.load(f)

print(data)

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



def plot_points(point_data):
    xyz = []
    for x in point_data:
        xyz.append(x['pos'])


    ax.scatter(*np.array(xyz).T, color='r', marker='o')

def plot_links(link_data): 
    xyz = []
    for x in link_data:
        xyz.append(x['pos'])

        if len(x["boundary"]) > 2:
            print("Malformed boundary: link should have two or zero ends")
        p0 = x["boundary"][0]
        p1 = x["boundary"][1]
        assert (p0["mult"] * p1["mult"] == -1)
        if p0["mult"] == -1:
            tmp = p0
            p0 = p1
            p1 = tmp

        # it should now be possible to deduce the correct wrapping of p0, p1,
        # such that all are in same cell
        ax.quiver(*p0['pos'], *(np.array(x['pos']) - p0['pos']), color='k')
        ax.quiver(*x['pos'], *(np.array(p1['pos'])-x['pos']), color='k', arrow_length_ratio=0)


func_to_run = {
        'points': plot_points,
        'links': plot_links
        }

for arg in sys.argv[2:]:
    func_to_run[arg](data[arg])


plt.show()
