import sys
import numpy as np
import numpy.linalg as LA
import json
from os.path import basename
import argparse
from fury import window, actor, ui, pick
import itertools


actors = []

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




def gen_points(data, args):
    point_data = data['points']
    xyz = []
    for i, pdata in enumerate(point_data):
        x = pdata['pos']
        xyz.append(x)
    xyz = np.array(xyz, dtype=np.float64)
    return [actor.sphere(xyz,
                         colors=np.repeat([[0.8, 0., 0.]], xyz.shape[0], axis=0),
                         radii=0.5)
            ]




def gen_link(x, A):
    r0, r1 = link2startstop(x)
    dx0, dx1 = (r0 - x['pos']), (r1 - x['pos'])
    if not are_parallel(dx0, dx1):
        dx0 = unwrap(dx0, A)
        dx1 = unwrap(dx1, A)

    return np.array([dx0 + x['pos'], dx1 + x['pos']], dtype=np.float64)

def gen_links(data, args):
    A = np.array(data['index_cell_vectors'])
    link_data = data['links']

    if args.thick_links:
        center_list = []
        dir_list = []
        for i, cell in enumerate(link_data):
            center_list.append(cell['pos'])
            r0, r1 = link2startstop(cell)
            dir_list.append(unwrap(r0 - cell['pos'], A))
        center_list = np.array(center_list, dtype=np.float64)
        dir_list = np.array(dir_list, dtype=np.float64)
        height_list = np.linalg.norm(dir_list, axis=1)*2
        return [actor.cylinder(center_list,
                               directions=dir_list,
                               heights=height_list,
                           colors=np.repeat([[1., 1., 1]], len(center_list),
                                            axis=0),
                           )]
    else:
        linklist = []
        for i, cell in enumerate(link_data):
            linklist.append(gen_link(cell, A))
        return [actor.line(linklist,
                           colors=np.repeat([[1., 1., 1]], len(linklist),
                                            axis=0)
                           )]





def find_link(linkpos, link_data):
    for link in link_data:
        if link['pos'] == linkpos:
            return link
    raise LookupError(f"No link at {linkpos}")


def gen_plaqs(data, args):
    link_data = data["links"]
    plaq_data = data["plaqs"]
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
            gen_idx(x, i)
        if args.show_pos:
            gen_pos(x)


def gen_vols(data, args):
    vol_data = data['vols']
    xyz = []
    for i, x in enumerate(vol_data):
        xyz.append(x['pos'])

        if args.show_idx:
            gen_idx(x, i)
        if args.show_pos:
            gen_pos(x)

    ax.scatter(*np.array(xyz).T, color='b', marker='o')


func_to_run = {
        'points': gen_points,
        'links': gen_links,
        'plaqs': gen_plaqs,
        'vols': gen_vols
        }

ap = argparse.ArgumentParser()
ap.add_argument("file", help="The .lat.json file specifying the lattice",
                type=str)
ap.add_argument("objects", nargs='+', choices=func_to_run.keys())
ap.add_argument("--show_idx", action='store_true')
ap.add_argument("--show_pos", action='store_true')
ap.add_argument("--undirected", action='store_true')
ap.add_argument("--window_size", type=int, nargs=2, default=(2400,1200))
ap.add_argument("--save")
ap.add_argument("--thick_links", action="store_true")

args = ap.parse_args()

data = None
with open(args.file, 'r') as f:
    data = json.load(f)


def parse_filename(fname):
    retval = {}
    for tok in basename(fname).split(';'):
        t = tok.split('=')
        if len(t) != 2:
            continue
        retval[t[0]] = t[1]
    return retval


###############################################################################
# Instantiating the actors

actors = []
for arg in args.objects:
    actors += func_to_run[arg](data, args)


###############################################################################
# Define the info panel

info_panel = ui.Panel2D(size=(300, 100), align="right")
info_panel.center = (150, 200)


opts = parse_filename(args.file)

title = basename(args.file)
try:
    title = f'{float(opts['p'])*100}% dilution, removed {opts['nn']}-neighbours'
except KeyError:
    pass

text_block = ui.TextBlock2D(text=title, font_size=20)
info_panel.add_element(text_block, (0.3, 0.3))

###############################################################################
# Builds the scene and initialises
scene = window.Scene()


for a in actors:
    scene.add(a)


# Create the show manager
showm = window.ShowManager(scene=scene, size=(800, 600), reset_camera=True)
scene.add(info_panel)



###############################################################################
# Picking

# Chore: Build a table to associate vertices with objects





label_actor = actor.vector_text(text="Test")
pickm = pick.PickingManager()

def left_click_callback(obj, event):
    # Get the event position on display and pick
    event_pos = pickm.event_position(showm.iren)
    picked_info = pickm.pick(event_pos, showm.scene)
    print(picked_info['actor'])

for a in actors:
    a.AddObserver("LeftButtonPressEvent", left_click_callback, 1)

# window.show(scene, size=args.window_size)

if args.save:
    print("Saving to "+args.save)
    window.record(scene=scene, out_path=args.save, size=args.window_size)
else:
    showm.start()


