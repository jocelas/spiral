from helper_functions import *
import numpy as np
import matplotlib.pyplot as plt
from helper_functions import *
from plotting_helix import *
import plotly.graph_objects as go

def vis_tangent_vectors(segments, L, theta, tau):
    radius, polar_angle, deltaz = find_helix_parameters(L, theta, tau)
    normal_vectors = []
    for n in range(segments + 1):
        x = - radius * polar_angle * np.sin(n * polar_angle)
        y = radius * polar_angle * np.cos(n * polar_angle)
        z = deltaz
        normal_vectors.append(unit(np.array([x,y,z])))
    return np.asarray(normal_vectors)

def xy_perp(v):
    ## find the vector orthogonal to v lying in the xy plane
    vx, vy, vz = v
    # Handle special case where v has no xy component
    if vx == 0 and vy == 0:
        # any xy vector is fine, choose (1, 0, 0)
        return np.array([1, 0, 0], dtype=float)
    w = np.array([-vy, vx, 0], dtype=float)
    return unit(w)

def vis_local_xy_plane(vectors):
    basis = []
    for vector in vectors:
        vector = unit(vector.copy())
        u1 = xy_perp(vector)
        u2 = unit(np.cross(vector, u1))
        basis.append([u1, u2])
    return np.array(basis)




def vis_circle(angle, point, basis, pipe_diameter):
    r = pipe_diameter / 2
    return point + r * (np.cos(angle)*basis[0] + np.sin(angle)*basis[1])

def vis_whole_circle(point, basis, pipe_diameter, sides = 10):
    points = []
    angles = np.linspace(0,2*np.pi, sides)
    for angle in angles:
        points.append(vis_circle(angle, point, basis, pipe_diameter))
    return np.asarray(points)

def connect_lines_mesh3d(line1, line2, color='lightblue', opacity=0.6, name='surface'):
    """
    Create a Mesh3d surface filling the area between two 3D lines.
    
    Parameters
    ----------
    line1, line2 : (N,3) array_like
        Arrays of 3D points with the same length.
    color : str
        Color of the surface.
    opacity : float
        Opacity of the surface.
    name : str
        Trace name for Plotly legend.
    
    Returns
    -------
    go.Mesh3d
        Plotly Mesh3d trace connecting the two lines.
    """
    line1 = np.asarray(line1)
    line2 = np.asarray(line2)
    assert line1.shape == line2.shape, "Both lines must have same number of points"

    n = len(line1)
    x = np.concatenate([line1[:,0], line2[:,0]])
    y = np.concatenate([line1[:,1], line2[:,1]])
    z = np.concatenate([line1[:,2], line2[:,2]])

    # Triangles between the two lines
    triangles = []
    for i in range(n - 1):
        triangles += [
            (i, i + 1, n + i),
            (i + 1, n + i + 1, n + i)
        ]
    i, j, k = np.array(triangles).T

    return go.Mesh3d(
        x=x, y=y, z=z,
        i=i, j=j, k=k,
        color=color, opacity=opacity,
        name=name
    )



def plot_helix(L, theta, tau, segments, pipe_diameter, sides = 10, double_helix = True, return_mode = False):

    fig = go.Figure()
    origin = np.zeros(3)

    points_straight = generate_straight_helix(L, theta, tau, segments)
    vectors = vis_tangent_vectors(segments, L, theta, tau)

    vectors = np.gradient(points_straight, axis=0)
    basis = vis_local_xy_plane(vectors)

    circles = []
    for i in range(len(basis)):
        circles.append(vis_whole_circle(points_straight[i], basis[i], pipe_diameter, sides = sides))
        
    circles = np.asarray(circles)
    lines = np.transpose(circles, axes=[1,0,2])

    for i in range(len(lines) -1):
        fig.add_trace(connect_lines_mesh3d(lines[i], lines[i+1], color='blue', opacity = .6))


    for circle in lines:
        x,y,z = circle.T
        fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode = 'lines', line = dict(color = 'blue')))
        i = i+1


    if double_helix:
            second_points = make_second_helix(points_straight)
            vectors = vis_tangent_vectors(segments, L, theta, tau)


            vectors = np.gradient(second_points, axis=0)
            basis = vis_local_xy_plane(vectors)

            circles = []
            for i in range(len(basis)):
                circles.append(vis_whole_circle(second_points[i], basis[i], pipe_diameter, sides = sides))
                
            circles = np.asarray(circles)
            lines = np.transpose(circles, axes=[1,0,2])

            for i in range(len(lines) -1):
                fig.add_trace(connect_lines_mesh3d(lines[i], lines[i+1], color='red', opacity = .6))



            for circle in lines:
                x,y,z = circle.T
                fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode = 'lines', line = dict(color = 'red')))
                i = i+1



    fig.update_layout(
    scene=dict(
        aspectmode='data'   # ensures equal scaling for x, y, z
    )
    )    

    fig.update_layout(
    dragmode="turntable",        # intuitive rotation
    scene_camera=dict(eye=dict(x=0, y=1.8, z=1.2)),  # pleasant starting view
    )

    if return_mode: return fig

    fig.show()