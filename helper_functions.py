import numpy as np

def homogeneous_to_cartesian(v):
    """Convert a homogeneous coordinate to a Cartesian coordinate.

    Args:
        v (array-like): A 4D vector in homogeneous coordinates.

    Returns:
        np.ndarray: A 3D vector in Cartesian coordinates.
    """
    v = np.asarray(v)
    if v.shape != (4,):
        raise ValueError("Input must be a 4D vector.")
    if v[3] == 0:
        raise ValueError("The last component of the homogeneous coordinate cannot be zero.")
    return v[:3] / v[3]

toCart = homogeneous_to_cartesian



def generate_next_point(point, Transformation):
    if len(point) == 3:
        point = np.append(point, 1)
    elif len(point) != 4:
        raise ValueError("point must be a 3D or 4D vector.")
    v = Transformation @ point
    return homogeneous_to_cartesian(v)

def generate_point_on_helix(n_points :int, point, Transformation):
    if n_points < 1:
        raise ValueError("n_points must be at least 1.")
    if len(point) == 4:
        point = homogeneous_to_cartesian(point)
    elif len(point) != 3:
        raise ValueError("point must be a 3D or 4D vector.")
    points = [point]
    for _ in range(n_points-1):
        point = generate_next_point(point, Transformation)
        points.append(point)
    return np.array(points)

def make_line(point, vector, length = 100): 
    """Generate two points defining a line in 3D space.

    Args:
        point (array-like): A point on the line (3D vector).
        vector (array-like): The direction vector of the line (3D vector).
        length (float, optional): The length of the line segment. Defaults to 100.

    Returns:
        np.ndarray: An array containing two points defining the line segment.
    """
    point = np.asarray(point)
    vector = np.asarray(vector)
    if point.shape != (3,) or vector.shape != (3,):
        raise ValueError("Both point and vector must be 3D vectors.")
    vector = vector / np.linalg.norm(vector) * (length / 2)
    return np.array([point - vector, point + vector])





def unit(v, eps=1e-12):
    n = np.linalg.norm(v)
    if n < eps: 
        raise ValueError("Zero-length vector encountered.")
    return v / n

def _bisector(vertex, p1, p2, eps=1e-9):
    # Internal angle bisector direction at 'vertex' for angle p1-vertex-p2
    u1 = unit(p1 - vertex)
    u2 = unit(p2 - vertex)
    d  = u1 + u2
    if np.linalg.norm(d) < eps:
        # Degenerate (straight angle) – fall back to external bisector
        d = u1 - u2
    return unit(d)

def closest_line_between_lines(P1, u, P2, v, eps=1e-12):
    """
    Lines: L1: P1 + t u,  L2: P2 + s v  (u,v assumed unit)
    Returns midpoint M of the shortest segment and direction vector w (from L1 to L2).
    """
    a = 1.0             # u·u
    b = float(u @ v)    # u·v
    c = 1.0             # v·v
    w0 = P1 - P2
    d = float(u @ w0)
    e = float(v @ w0)
    denom = a*c - b*b   # = 1 - (u·v)^2

    if denom < eps:
        # Lines are parallel (or nearly). Shortest vector is w0 minus its projection on u.
        w = w0 - (d) * u
        if np.linalg.norm(w) < eps:
            # Same line – arbitrary: return that line
            return P1.copy(), u
        Pq = P1
        Pr = P2 + (e) * v  # any point on L2; not strictly needed
    else:
        t = (b*e - c*d) / denom
        s = (a*e - b*d) / denom
        Pq = P1 + t*u
        Pr = P2 + s*v
        w  = Pr - Pq

    M = 0.5*(Pq + Pr)          # a point on the connecting line
    if np.linalg.norm(w) < eps:
        # Lines intersect – return intersection point and any perpendicular direction
        return M, unit(np.cross(u, np.array([1,0,0])))
    dir_vec = unit(w)          # direction of the connecting line
    return M, dir_vec

def connecting_line_for_bisectors(a, b, c, d):
    """
    a,b,c,d: (3,) numpy arrays
    Returns (point, direction) for the line that connects (at shortest distance)
    the angle-bisector at B of angle ABC and the angle-bisector at C of angle BCD.
    """
    # Bisector at B for angle ABC
    uB = _bisector(b, a, c)     # line L1: x = b + t*uB
    # Bisector at C for angle BCD
    uC = _bisector(c, b, d)     # line L2: x = c + s*uC

    M, dir_vec = closest_line_between_lines(b, uB, c, uC)
    return M, dir_vec


def find_helix_axis(point, transformation):
    points = generate_point_on_helix(4, point, transformation)
    a, b, c, d = points
    point, direction = connecting_line_for_bisectors(a, b, c, d)
    return point, direction




# -------- Example usage --------
# a = np.array([0.0, 0.0, 0.0])
# b = np.array([1.0, 0.0, 0.0])
# c = np.array([1.0, 1.0, 0.0])
# d = np.array([2.0, 1.0, 1.0])
# point, direction = connecting_line_for_bisectors(a, b, c, d)
# print("Point on connecting line:", point)
# print("Direction (unit):", direction)


def rotation_to_z(v):
    v = np.asarray(v, dtype=float)
    v = v / np.linalg.norm(v)          # normalize
    z = np.array([0, 0, 1], float)

    if np.allclose(v, z):              # already aligned
        return np.eye(3)
    if np.allclose(v, -z):             # opposite direction
        # rotate 180° about any perpendicular axis
        return np.array([[-1, 0, 0],
                         [ 0,-1, 0],
                         [ 0, 0, 1]])

    axis = np.cross(v, z)
    axis /= np.linalg.norm(axis)
    angle = np.arccos(np.clip(np.dot(v, z), -1.0, 1.0))

    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])

    R = np.eye(3) + np.sin(angle)*K + (1 - np.cos(angle))*(K @ K)
    return R


def translation_to_origin(point):
    point = np.asarray(point, dtype=float)
    if point.shape != (3,):
        raise ValueError("point must be a 3D vector.")
    T = np.eye(4)
    T[:3, 3] = -point
    return T

def make_helix_axis_z_axis(point, direction, return_transformation=True):
    point = np.asarray(point, dtype=float)
    direction = np.asarray(direction, dtype=float)
    if point.shape == (3,):
        ponint = np.append(point, 1)
    elif point.shape != (4,):
        raise ValueError("point must be a 3D or 4D vector.")
    if direction.shape != (3,):
        raise ValueError("direction must be a 3D vector.")
    direction = np.append(direction, 1)
    direction = direction / np.linalg.norm(direction)  # normalize
    T = translation_to_origin(point)
    R = rotation_to_z(direction)
    RT = np.eye(4)
    RT[:3, :3] = R
    point_in_z = RT @ T @ point
    direction_in_z = RT @ T @ direction
    assert np.allclose(direction_in_z, np.array([0, 0, 1, 1])), "Direction not aligned with z-axis"
    assert np.allclose(point_in_z[2], 0), "Point not in xy-plane"
    if return_transformation:
        return point_in_z, direction_in_z, R
    return point_in_z, direction_in_z

def transform_points(points, axis_origin, axis_direction):
    points -= axis_origin
    R = rotation_to_z(axis_direction)
    points = points @ R.T
    return points



def turn_spiral_about_z_axis(points):
    """
    Rotate a list or Nx3 array of 3D points by 180° around the z-axis.
    """
    points = np.asarray(points, dtype=float)
    R = np.array([[-1,  0, 0],
                  [ 0, -1, 0],
                  [ 0,  0, 1]])  # 180° about z

    return points @ R.T  # apply rotation to all points
