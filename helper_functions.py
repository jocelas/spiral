import numpy as np
import warnings

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


def transformation_matrix(L, theta, tau):
    '''To construct the helix, we will assume that the first point B lies at the origin. 
The first segment goes in the x direction, so C is at $(L,0,0,1)$.
Further, the joint faces of the first segment are normal to the $xy$-plane, meaning
$$N_B = (-\cos(\theta/2), \sin(\theta/2), 0, 1)\quad N_C = (\cos(\theta/2), \sin(\theta/2), 0, 1)$$
which are the normal vectors to the joint faces.

Note that this means that the axis of the helix is not the $z$-axis, and it is in fact not even parallel to it unless the twist angle $\tau$ is 0.

To determine the third point A we do the following:

1. translate C to the origin ($\vec t = (-L, 0, 0), R = \mathbf 1$)
2. rotate by $-\theta/2$ along the $z$-axis. Now $\vec N_B$ is aligned with the $x$-axis.
3. rotate by $\tau$ around the $x$-axis.
4. rotate by $-\theta/2$ around the $z$-axis again. Now the point D would lie where C was before, and B has moved to A.

Doing this transformation to B will give you A, and doing the inverse to C will give you D'''
    tau_rad = np.deg2rad(tau)
    theta_rad = np.deg2rad(theta)

    ## translate b to origin
    T1 = np.array([[1,0,0,-L],
                    [0,1,0,0],
                    [0,0,1,0],
                    [0,0,0,1]])

    ## rotate by -theta/2 around z axis
    R1 = np.array([[np.cos(-theta_rad/2), -np.sin(-theta_rad/2), 0, 0],
                    [np.sin(-theta_rad/2), np.cos(-theta_rad/2), 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]])

    ## rotate by tau around x axis
    R2 = np.array([[1, 0, 0, 0],
                    [0, np.cos(tau_rad), -np.sin(tau_rad), 0],
                    [0, np.sin(tau_rad), np.cos(tau_rad), 0],
                    [0, 0, 0, 1]])

    ## rotate by -tau around x axis
    R2_inv = np.array([[1, 0, 0, 0],
                        [0, np.cos(-tau_rad), -np.sin(-tau_rad), 0],
                        [0, np.sin(-tau_rad), np.cos(-tau_rad), 0],
                        [0, 0, 0, 1]])  

    ## rotate by theta/2 around z axis
    R3 = np.array([[np.cos(theta_rad/2), -np.sin(theta_rad/2), 0, 0],
                    [np.sin(theta_rad/2), np.cos(theta_rad/2), 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]])

    ## translate back
    T2 = np.array([[1,0,0,L],
                    [0,1,0,0],        
                    [0,0,1,0],
                    [0,0,0,1]])


    M = R1@R2@R1@T1
    inverse_M = np.linalg.inv(M)
    inverse_M_from_scratch = T2@R3@R2_inv@R3

    if np.allclose(inverse_M, inverse_M_from_scratch):
       # print("OK. Inverse matrix calculation is consistent")
       pass
    else:  
        print("Inverse matrix calculation is NOT consistent")
        raise ValueError("Inverse matrix calculation is NOT consistent")
    
    return M, inverse_M



def generate_next_point(point, Transformation):

    if len(point) == 3:
        point = np.append(point, 1)
    elif len(point) != 4:
        raise ValueError("point must be a 3D or 4D vector.")
    v = Transformation @ point
    return homogeneous_to_cartesian(v)

def generate_points_on_helix(segments :int, Transformation):
    # generates N segments based on the transformation
    if segments < 1:
        raise ValueError("n_points must be at least 1.")
    point = np.array([0,0,0]) # Transformation matrix works such that the starting point is on the origin!
    points = [point]
    for _ in range(segments):
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




### These functions are important for finding the axis of the Helix

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


def find_axis_dir_ev(M_in):
    ## This finds the direction of the axis, as the direction is an eigenvector of the rotation Part of the trafo matrix
    M = M_in.copy()[:3,:3]
    vals, vecs = np.linalg.eig(M)

    axis_vec = np.real(vecs[:, np.isclose(vals, 1.)])

    #test if it actually an eigenvector
    assert np.allclose(M@axis_vec, axis_vec), "Computation of eigenvector failed"

    unit(axis_vec)

    return axis_vec[:,0]



def find_helix_axis(transformation):
    points = generate_points_on_helix(3, transformation)
    a, b, c, d = points
    point, direction = connecting_line_for_bisectors(a, b, c, d)
    ev_direction = find_axis_dir_ev(transformation)

    assert np.allclose(unit(direction), unit(ev_direction), atol = 1e-5) or np.allclose(unit(direction), -unit(ev_direction), atol = 1e-5), f"Direction from bisectors does not match eigenvector direction \n {np.abs(np.dot(unit(direction), unit(ev_direction)))}"   

    return point, direction

# end of axis finding


def rotation_to_z(v):
    ## take a vector and return the transformation that rotates it to the z axis
    # done with ChatGPT
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

    # Check if det R is 1:
    assert np.isclose(np.linalg.det(R), 1), "The rotation to z was not correctly caluclated."

    return R


def translation_to_origin(point):
    # find homogeneous matrix that translates a point to the origin.
    point = np.asarray(point, dtype=float)
    if point.shape != (3,):
        raise ValueError("point must be a 3D vector.")
    T = np.eye(4)
    T[:3, 3] = -point
    return T

def transform_helix_to_z_axis(points, axis_origin, axis_direction):
    points -= axis_origin
    R = rotation_to_z(axis_direction)
    points = points @ R.T # rotate all the points such that the axis is z axis
    points -= np.array([0,0,points[0,2]]) # put the first segment on the xy plane
    return points


def make_second_helix(points):
    # makes the second helix for a double helix
    points = np.asarray(points, dtype=float)
    R = np.array([[-1,  0, 0],
                  [ 0, -1, 0],
                  [ 0,  0, 1.]])  # 180° about z

    return points @ R.T  # apply rotation to all points

def point_line_distance(P, A, v):
    #Distance between point P and line passing through A with direction v.
    P, A, v = map(lambda x: np.asarray(x, float), (P, A, v))
    return np.linalg.norm(np.cross(P - A, v)) / np.linalg.norm(v)

def angle(v1, v2):
    # Angle between two vectors v1 and v2
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    dot_product = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    angle_rad = np.arccos(dot_product)
    return np.rad2deg(angle_rad)

def find_number_of_segments_for_length(points_in, length):
    
    delta_points = np.diff(points_in.copy(), axis=0)
    assert np.allclose(delta_points[:,2], delta_points[0,2])
    deltaz = delta_points[0,2]
    number = int(np.floor(length/deltaz))

    return length, deltaz, number, deltaz*number


def full_helix_calculation(L, theta, tau, length = None, segments = None):

    if length is None and segments is None:
        warnings.warn("Either length or segments should be provided. Defaulting to segments = 25.")
        segments = 25

    if length is not None and length  <= 0:
        warnings.warn("Lenght is invalid, defaulting to automatic length calculation")
        length = None

    M, inverse_M = transformation_matrix(L, theta, tau)

    points = generate_points_on_helix(4, M)
    axis_origin, axis_direction = find_helix_axis(M)

    helix_radius = point_line_distance(points[0], axis_origin, axis_direction)
    if length is not None and segments is None:
        points_straight = transform_helix_to_z_axis(points.copy(), axis_origin, axis_direction)
        _, deltaz, segments, actual_length  = find_number_of_segments_for_length(points_straight, length)

    points = generate_points_on_helix(segments, M)
    points_straight = transform_helix_to_z_axis(points.copy(), axis_origin, axis_direction)

    if length is not None and segments is None:
        return points_straight, helix_radius, segments, actual_length, deltaz
    
    if length is None and segments is not None:
        actual_length = np.abs(points_straight[-1,2] - points_straight[0,2])

        deltaz = np.abs(points_straight[1,2] - points_straight[0,2])
        assert np.isclose(actual_length/(segments), deltaz), f"Calculated pitch (deltaz) does not match actual pitch from points.\n Deltaz = {deltaz}, length/segments = {actual_length/segments}" 
        return points_straight, helix_radius, segments, actual_length, deltaz

    if length is not None and segments is not None:
        actual_length = np.abs(points_straight[-1,2] - points_straight[0,2])
        _, deltaz, segments_from_length, actual_length_from_length  = find_number_of_segments_for_length(points_straight, length)
        return points_straight, helix_radius, segments, actual_length, deltaz, segments_from_length, actual_length_from_length

    raise NotImplementedError("Something has gone wrong in the full helix calculation!")

def generate_straight_helix(L, theta, tau, segments):
    # generates straight matrix points
    M, _ = transformation_matrix(L, theta, tau)
    points = generate_points_on_helix(segments, M)
    axis_origin, axis_direction = find_helix_axis(M)
    points_straight = transform_helix_to_z_axis(points, axis_origin, axis_direction)
    return points_straight



def center_length_from_outer(outer_length, cut_angle, diameter):
    cut_angle_rad = np.deg2rad(cut_angle)
    return outer_length - diameter * np.tan(cut_angle_rad)


def polar_angle_between(p1, p2):
    # returns polar angle in **RADIANS**
    x1, y1, _ = p1
    x2, y2, _ = p2
    theta1 = np.arctan2(y1, x1)
    theta2 = np.arctan2(y2, x2)
    dtheta = theta2 - theta1
    # normalize to [-pi, pi)
    dtheta = (dtheta + np.pi) % (2*np.pi) - np.pi
    return dtheta

def find_number_of_turns(polar_angle, number):
    return polar_angle * number / 2 / np.pi


def find_helix_parameters(L, theta, tau):

    points_straight = generate_straight_helix(L, theta, tau, 4)

    radius = point_line_distance(points_straight[0], np.array([0.,0,0]), np.array([0,0,1.]))
    radius_test = point_line_distance(points_straight[1], np.array([0.,0,0]), np.array([0,0,1.]))

    assert np.isclose(radius, radius_test), "radius calculation wrong"

    polar_angle = polar_angle_between(points_straight[1], points_straight[0])

    deltaz = points_straight[1,2] - points_straight[0,2]
    assert np.isclose (deltaz, points_straight[-1,2] - points_straight[-2,2]), "deltaz calculation wrong"

    return radius, polar_angle, deltaz



def analytical_angle_to_xy_plane(radius, polar_angle, deltaz):
    # Calculates angle to xy plane of analytical spiral
    pitch = deltaz * 2 * np.pi / polar_angle # height increase for 1 entire turn
    tanalpha = pitch / (2 * np.pi * radius)
    alpha = np.arctan(tanalpha)
    return np.rad2deg(alpha)

def segment_angle_to_xy_plane(L, deltaz):
    # Calculates angle of 2 points on helix with z as the helix axis to the xy plane (plane orthogonal to z axis)
    sinalpha = deltaz / L
    alpha = np.arcsin(sinalpha)
    return np.rad2deg(alpha)



# r, theta, d = find_helix_parameters(30,20,15)

# ang = tangent_angle(r, theta, d)

# M, _ = transformation_matrix(30,20,15)
# points = generate_points_on_helix(5, M)
# axis_origin, axis_direction = find_helix_axis(points[0], M)
# points_straight = transform_points(points, axis_origin, axis_direction)

# angle(points_straight[1]-points_straight[0], [0,0,1]), ang

# distance = 2 * np.pi * d * np.cos(np.deg2rad(ang))

# ang, distance, 2 * np.pi * d



def is_helix_loose_enough(L, theta, tau, pipe_diameter, eps = 1):
    radius, polar_angle, deltaz = find_helix_parameters(L, theta, tau)
    nfull = 2*np.pi / polar_angle
    if nfull * deltaz > 2 * (pipe_diameter + eps):
        return True
    else:
        return False

def is_helix_wide_enough(L, theta, tau, pipe_diameter, eps = 1):
    radius, polar_anlge, deltaz = find_helix_parameters(L, theta, tau)
    if 2 * radius > pipe_diameter + eps:
        return True
    else:
        return False
    
def is_helix_allowed(L, theta, tau, pipe_diameter, eps_wide = 3, eps_loose = 3):
    return is_helix_loose_enough(L, theta, tau, pipe_diameter, eps = eps_loose) * is_helix_wide_enough(L, theta, tau, pipe_diameter, eps = eps_wide)


def find_allowed_twist_range(L, theta, pipe_diameter, 
                             eps_loose=3, eps_wide=3,
                             tau_min=0.0, tau_max=180.0, tol=1e-4):
    """
    Find the tau interval [tau_lo, tau_hi] for which is_helix_allowed() is True.
    Uses binary search on each edge.
    """
    # Done with ChatGPT
    def f(tau):
        return is_helix_allowed(L, theta, tau, pipe_diameter, eps_loose=eps_loose, eps_wide=eps_wide)

    # First, locate roughly where f becomes True
    N = 50  # coarse scan
    taus = np.linspace(tau_min, tau_max, N)
    vals = [f(t) for t in taus]

    if not any(vals):
        return None, None  # no valid tau range at all

    # Find coarse True region indices
    idx_true = np.where(vals)[0]
    i_start, i_end = idx_true[0], idx_true[-1]
    lo_left, hi_left = taus[max(i_start - 1, 0)], taus[i_start]
    lo_right, hi_right = taus[i_end], taus[min(i_end + 1, N - 1)]

    # --- refine lower edge ---
    while hi_left - lo_left > tol:
        mid = 0.5 * (lo_left + hi_left)
        if f(mid):
            hi_left = mid
        else:
            lo_left = mid
    tau_lo = hi_left

    # --- refine upper edge ---
    while hi_right - lo_right > tol:
        mid = 0.5 * (lo_right + hi_right)
        if f(mid):
            lo_right = mid
        else:
            hi_right = mid
    tau_hi = lo_right

    return tau_lo, tau_hi


def list_possible_helices(L, theta, pipe_diameter, desired_length, tau_increment = 1, bottom_offset = 2, top_offset = 0):
    taumin, taumax = find_allowed_twist_range(L, theta, pipe_diameter, eps_loose=1, eps_wide=1)
    taumin, taumax = np.ceil(taumin) + bottom_offset, np.floor(taumax) - top_offset

    list_of_taus = np.arange(taumin, taumax, tau_increment)

    outputs = []

    for tau in list_of_taus:
        radius, polar_angle, deltaz = find_helix_parameters(L, theta, tau)

        Nsegments_for_length = np.floor(desired_length / deltaz)
        actual_length = deltaz * Nsegments_for_length

        donut_hole_diameter = np.abs(radius * 2 - pipe_diameter)

        Nturns = find_number_of_turns(polar_angle, Nsegments_for_length)

        turns_per_segment = find_number_of_turns(polar_angle, 1)

        outer_diameter = pipe_diameter + radius * 2

        outputs.append([tau, Nsegments_for_length, actual_length, Nturns, outer_diameter, donut_hole_diameter, deltaz, turns_per_segment])

    return np.array(outputs)