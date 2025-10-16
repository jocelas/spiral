/***** helperfunctions.js *****/
// Basic vector operations (3D vectors represented as [x, y, z])
function add(u, v)   { return [u[0] + v[0], u[1] + v[1], u[2] + v[2]]; }
function sub(u, v)   { return [u[0] - v[0], u[1] - v[1], u[2] - v[2]]; }
function mul(v, s)   { return [v[0] * s, v[1] * s, v[2] * s]; }
function dot(u, v)   { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; }
function cross(u, v) {
  return [
    u[1]*v[2] - u[2]*v[1],
    u[2]*v[0] - u[0]*v[2],
    u[0]*v[1] - u[1]*v[0]
  ];
}
function norm(v)     { return Math.sqrt(dot(v, v)); }
function unit(v) {
  const n = norm(v);
  if (n < 1e-12) throw new Error("Zero-length vector");
  return [v[0]/n, v[1]/n, v[2]/n];
}

// Compute centerline length from outer length, given cut angle and pipe diameter
function centerLengthFromOuter(outer_length, cut_angle_deg, pipe_diameter_inch) {
  const cut_angle_rad = cut_angle_deg * Math.PI / 180;
  // Convert pipe diameter to same unit (mm) as outer_length before calculation
  const pipe_diameter_mm = pipe_diameter_inch * 25.4;
  return outer_length - pipe_diameter_mm * Math.tan(cut_angle_rad);
}

// Construct the homogeneous transformation matrix for one helix segment.
// L: segment center length (straight), theta: total cut angle (deg), tau: twist angle (deg).
function transformationMatrix(L, theta_deg, tau_deg) {
  const tau = tau_deg * Math.PI / 180;
  const theta = theta_deg * Math.PI / 180;
  // Define 4x4 matrices as nested arrays
  const T1 = [  // Translation of B to origin: translate by (-L, 0, 0)
    [1, 0, 0, -L],
    [0, 1, 0,  0],
    [0, 0, 1,  0],
    [0, 0, 0,  1]
  ];
  const cos_half = Math.cos(-theta/2), sin_half = Math.sin(-theta/2);
  const R1 = [  // Rotate by -θ/2 about Z-axis
    [ cos_half, -sin_half, 0, 0],
    [ sin_half,  cos_half, 0, 0],
    [ 0,        0,        1, 0],
    [ 0,        0,        0, 1]
  ];
  const cosT = Math.cos(tau), sinT = Math.sin(tau);
  const R2 = [  // Rotate by τ around X-axis
    [1,    0,     0,    0],
    [0,  cosT, -sinT,   0],
    [0,  sinT,  cosT,   0],
    [0,    0,     0,    1]
  ];
  // Compute M = R1 * R2 * R1 * T1 (matrix multiplication)
  function multiplyMatrices(A, B) {
    const result = Array.from({ length: 4 }, () => [0, 0, 0, 0]);
    for (let i = 0; i < 4; i++) {
      for (let j = 0; j < 4; j++) {
        let sum = 0;
        for (let k = 0; k < 4; k++) {
          sum += A[i][k] * B[k][j];
        }
        result[i][j] = sum;
      }
    }
    return result;
  }
  // Multiply in sequence: R1 * R2
  const R1_R2 = multiplyMatrices(R1, R2);
  // Then (R1*R2) * R1
  const R1R2_R1 = multiplyMatrices(R1_R2, R1);
  // Finally * T1
  const M = multiplyMatrices(R1R2_R1, T1);
  return M;
}

// Apply a 4x4 transformation matrix to a 3D point (homogeneous coordinate)
function applyTransformation(point, M) {
  // point is [x, y, z, 1] or [x,y,z] (we handle both)
  const x = point[0], y = point[1], z = point[2];
  const X = M[0][0]*x + M[0][1]*y + M[0][2]*z + M[0][3]*1;
  const Y = M[1][0]*x + M[1][1]*y + M[1][2]*z + M[1][3]*1;
  const Z = M[2][0]*x + M[2][1]*y + M[2][2]*z + M[2][3]*1;
  const W = M[3][0]*x + M[3][1]*y + M[3][2]*z + M[3][3]*1;
  // Return Cartesian coordinates (divide by W, though W should be 1 for rigid transforms)
  return [X/W, Y/W, Z/W];
}

// Generate points along the helix by repeatedly applying transformation matrix M.
// Returns an array of [x,y,z] points of length (n+1) (including the starting point).
function generatePointsOnHelix(n, M) {
  const points = [];
  // Start at B (origin of first segment) at (0,0,0)
  let currentPoint = [0, 0, 0];  // we'll treat this as [x,y,z] with implicit homogeneous 1
  points.push([...currentPoint]);  // initial point
  for (let i = 0; i < n; i++) {
    const newPoint = applyTransformation(currentPoint, M);
    points.push(newPoint);
    currentPoint = newPoint;
  }
  return points;
}

// Internal helper: compute bisector direction at 'vertex' for angle (p1-vertex-p2)
function bisector(vertex, p1, p2) {
  const u1 = unit(sub(p1, vertex));
  const u2 = unit(sub(p2, vertex));
  let d = add(u1, u2);
  if (norm(d) < 1e-9) {
    // Degenerate case (straight line): use external bisector (difference)
    d = sub(u1, u2);
  }
  return unit(d);
}

// Compute the shortest connecting line between two lines in space.
// Each line is given by a point (P1, P2) and direction (u, v). Assumes u, v are unit vectors.
function closestLineBetweenLines(P1, u, P2, v) {
  const w0 = sub(P1, P2);
  const a = 1.0;             // dot(u,u)
  const b = dot(u, v);
  const c = 1.0;             // dot(v,v)
  const d = dot(u, w0);
  const e = dot(v, w0);
  const denom = a*c - b*b;
  let Pq, Pr, w;
  if (Math.abs(denom) < 1e-12) {
    // Lines are nearly parallel
    // Take vector connecting P1 to P2, project it onto u for offset
    const w_parallel = mul(u, d);
    w = sub(w0, w_parallel);
    if (norm(w) < 1e-12) {
      // Lines are coincident (overlap) – choose an arbitrary perpendicular direction
      return { point: P1, direction: unit(cross(u, [1, 0, 0])) };
    }
    Pq = P1;
    Pr = add(P2, mul(v, e));
  } else {
    // Compute closest points on each line
    const t = (b*e - c*d) / denom;
    const s = (a*e - b*d) / denom;
    Pq = add(P1, mul(u, t));
    Pr = add(P2, mul(v, s));
    w = sub(Pr, Pq);
  }
  const midpoint = mul(add(Pq, Pr), 0.5);
  if (norm(w) < 1e-12) {
    // Lines intersect – use any perpendicular direction (ensure not zero vector)
    return { point: midpoint, direction: unit(cross(u, [1, 0, 0])) };
  }
  return { point: midpoint, direction: unit(w) };
}

// Determine the helix axis for the transformation matrix M.
// Returns an object {origin: [x,y,z], direction: [dx,dy,dz]} describing the axis line.
function findHelixAxis(M) {
  // Generate four points along the helix (3 segments)
  const pts = generatePointsOnHelix(3, M);  // yields [a, b, c, d]
  const [a, b, c, d] = pts;
  // Compute bisector directions at b (for angle a-b-c) and at c (for angle b-c-d)
  const uB = bisector(b, a, c);
  const uC = bisector(c, b, d);
  // Find closest connecting line between line through b (dir uB) and line through c (dir uC)
  const axisLine = closestLineBetweenLines(b, uB, c, uC);
  // axisLine.point and axisLine.direction are the axis origin and direction.
  return { origin: axisLine.point, direction: unit(axisLine.direction) };
}

// Compute a rotation matrix (3x3) that rotates a given unit vector 'v' to align with [0,0,1] (the z-axis).
function rotationToZ(v) {
  // Ensure v is a unit vector
  const nv = unit(v);
  const z = [0, 0, 1];
  // If already aligned with z-axis (or opposite), handle directly
  if (Math.abs(nv[0] - z[0]) < 1e-9 && Math.abs(nv[1] - z[1]) < 1e-9 && Math.abs(nv[2] - z[2]) < 1e-9) {
    return [ [1,0,0], [0,1,0], [0,0,1] ];  // identity
  }
  if (Math.abs(nv[0] + z[0]) < 1e-9 && Math.abs(nv[1] + z[1]) < 1e-9 && Math.abs(nv[2] + z[2]) < 1e-9) {
    // v is opposite of z: rotate 180° about any axis in the XY plane (choose X-axis)
    return [ [1,0,0], [0,-1,0], [0,0,-1] ];
  }
  // Otherwise, find rotation axis and angle using Rodrigues' formula
  const axis = unit(cross(nv, z));                // axis to rotate around
  const cosA = dot(nv, z);                        // cos(theta) where theta is angle between nv and z
  const angle = Math.acos(Math.max(-1, Math.min(1, cosA)));  // clamp cosA to [-1,1] and get angle
  const K = [
    [0,        -axis[2],  axis[1]],
    [ axis[2],  0,       -axis[0]],
    [-axis[1],  axis[0],  0      ]
  ];
  // Rodrigues' rotation formula: R = I + sinθ*K + (1-cosθ)*K^2
  const I = [ [1,0,0], [0,1,0], [0,0,1] ];
  const K2 = [  // K squared
    [
      K[0][0]*K[0][0] + K[0][1]*K[1][0] + K[0][2]*K[2][0],
      K[0][0]*K[0][1] + K[0][1]*K[1][1] + K[0][2]*K[2][1],
      K[0][0]*K[0][2] + K[0][1]*K[1][2] + K[0][2]*K[2][2]
    ],
    [
      K[1][0]*K[0][0] + K[1][1]*K[1][0] + K[1][2]*K[2][0],
      K[1][0]*K[0][1] + K[1][1]*K[1][1] + K[1][2]*K[2][1],
      K[1][0]*K[0][2] + K[1][1]*K[1][2] + K[1][2]*K[2][2]
    ],
    [
      K[2][0]*K[0][0] + K[2][1]*K[1][0] + K[2][2]*K[2][0],
      K[2][0]*K[0][1] + K[2][1]*K[1][1] + K[2][2]*K[2][1],
      K[2][0]*K[0][2] + K[2][1]*K[1][2] + K[2][2]*K[2][2]
    ]
  ];
  const sinAngle = Math.sin(angle), cosAngle = Math.cos(angle);
  // R = I + sinθ*K + (1-cosθ)*K^2
  const R = [ [0,0,0], [0,0,0], [0,0,0] ];
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      R[i][j] = I[i][j] + sinAngle * K[i][j] + (1 - cosAngle) * K2[i][j];
    }
  }
  return R;
}

// Transform helix points to align the helix axis with the z-axis.
// points: array of [x,y,z]; axisOrigin: a point on the helix axis; axisDir: axis direction vector.
function transformHelixToZAxis(points, axisOrigin, axisDir) {
  // Translate points so that axisOrigin goes to (0,0,0)
  const translated = points.map(p => sub(p, axisOrigin));
  // Rotate all points so that axisDir aligns with [0,0,1]
  const R = rotationToZ(axisDir);
  function applyRotation(vec) {
    // Multiply 3x3 matrix R by vector vec
    const [x, y, z] = vec;
    return [
      R[0][0]*x + R[0][1]*y + R[0][2]*z,
      R[1][0]*x + R[1][1]*y + R[1][2]*z,
      R[2][0]*x + R[2][1]*y + R[2][2]*z
    ];
  }
  const rotated = translated.map(applyRotation);
  // Shift vertically so that the first point has z = 0 (put the base of helix on XY-plane)
  if (rotated.length > 0) {
    const baseZ = rotated[0][2];
    for (let i = 0; i < rotated.length; i++) {
      rotated[i][2] -= baseZ;
    }
  }
  return rotated;
}

// Generate a helix aligned with the z-axis (vertical), given dimensions and segment count.
// Returns an array of [x,y,z] points (segment endpoints along the helix).
function generateStraightHelix(L, theta, tau, segments) {
  const M = transformationMatrix(L, theta, tau);
  const points = generatePointsOnHelix(segments, M);
  const axis = findHelixAxis(M);
  const straightPoints = transformHelixToZAxis(points, axis.origin, axis.direction);
  return straightPoints;
}

// Compute helix radius, polar angle, and vertical rise per segment (deltaz) for given parameters.
function findHelixParameters(L, theta, tau) {
  // Use 4 segments to approximate parameters
  const pts = generateStraightHelix(L, theta, tau, 4);
  // Radius: distance from z-axis (0,0) to first point
  const radius = Math.hypot(pts[0][0], pts[0][1]);
  // Check second point's radius (should be same if axis alignment is correct)
  const radius_test = Math.hypot(pts[1][0], pts[1][1]);
  // Polar angle between first two points around z-axis (in radians)
  let dtheta = Math.atan2(pts[1][1], pts[1][0]) - Math.atan2(pts[0][1], pts[0][0]);
  // Normalize angle to [-π, π)
  dtheta = ((dtheta + Math.PI) % (2*Math.PI)) - Math.PI;
  // Vertical rise per segment
  const deltaz = pts[1][2] - pts[0][2];
  // (We assume radius ≈ radius_test; small differences may occur due to numerical tolerance)
  return { radius, polar_angle: dtheta, deltaz };
}

// Check if helix has enough vertical spacing between coils (loose enough)
function isHelixLooseEnough(L, theta, tau, pipeDiam, eps = 1) {
  const { polar_angle, deltaz } = findHelixParameters(L, theta, tau);
  const turns_per_segment = polar_angle / (2 * Math.PI);
  // nfull = number of segments per full 360° turn = 1/turns_per_segment (if polar_angle is fraction of 2π)
  const nfull = (Math.abs(turns_per_segment) < 1e-9) ? Infinity : 1 / turns_per_segment;
  // Condition: vertical rise for one full turn > 2*(pipeDiam + eps)
  return nfull * Math.abs(deltaz) > 2 * (pipeDiam + eps);
}

// Check if helix is wide enough (radius sufficiently larger than pipe radius)
function isHelixWideEnough(L, theta, tau, pipeDiam, eps = 1) {
  const { radius } = findHelixParameters(L, theta, tau);
  return 2 * radius > pipeDiam + eps;
}

// Check if a helix configuration meets both loose and wide criteria
function isHelixAllowed(L, theta, tau, pipeDiam, eps_wide = 3, eps_loose = 3) {
  return isHelixLooseEnough(L, theta, tau, pipeDiam, eps_loose) &&
         isHelixWideEnough(L, theta, tau, pipeDiam, eps_wide);
}

// Find the range [tau_lo, tau_hi] of twist angles (tau in degrees) that produce a valid helix.
// Scans from 0 to 180 degrees and uses binary search to refine the edges.
function findAllowedTwistRange(L, theta, pipeDiam, eps_loose = 3, eps_wide = 3, tau_min = 0, tau_max = 180) {
  const N = 50;
  const taus = [];
  const vals = [];
  // Coarse scan
  for (let i = 0; i < N; i++) {
    const tau = tau_min + (tau_max - tau_min) * (i / (N - 1));
    taus.push(tau);
    vals.push(isHelixAllowed(L, theta, tau, pipeDiam, eps_wide, eps_loose));
  }
  if (!vals.includes(true)) {
    return [null, null];  // no valid tau range
  }
  const idx_true = vals.map((v,i) => v ? i : -1).filter(i => i >= 0);
  const i_start = idx_true[0];
  const i_end = idx_true[idx_true.length - 1];
  // Bracket the true region with one false on each side if possible
  let lo_left = taus[Math.max(i_start - 1, 0)];
  let hi_left = taus[i_start];
  let lo_right = taus[i_end];
  let hi_right = taus[Math.min(i_end + 1, N - 1)];
  // Refine lower bound
  while (hi_left - lo_left > 1e-4) {
    const mid = 0.5 * (lo_left + hi_left);
    if (isHelixAllowed(L, theta, mid, pipeDiam, eps_wide, eps_loose)) {
      hi_left = mid;
    } else {
      lo_left = mid;
    }
  }
  const tau_lo = hi_left;
  // Refine upper bound
  while (hi_right - lo_right > 1e-4) {
    const mid = 0.5 * (lo_right + hi_right);
    if (isHelixAllowed(L, theta, mid, pipeDiam, eps_wide, eps_loose)) {
      lo_right = mid;
    } else {
      hi_right = mid;
    }
  }
  const tau_hi = lo_right;
  return [tau_lo, tau_hi];
}

// List all possible helix configurations (twist angles) for given dimensions.
// Returns an array of objects for each feasible tau in [tau_min, tau_max] (incremented by tau_increment).
function listPossibleHelices(L, theta, pipeDiam, desiredLength, tau_increment = 1, bottom_offset = 1, top_offset = 0) {
  const [tau_min, tau_max] = findAllowedTwistRange(L, theta, pipeDiam, 1, 1);
  if (tau_min === null || tau_max === null) {
    return [];  // no valid configurations
  }
  // Apply offsets and rounding to bounds
  let startTau = Math.ceil(tau_min) + bottom_offset;
  let endTau = Math.floor(tau_max) - top_offset;
  const results = [];
  for (let tau = startTau; tau < endTau; tau += tau_increment) {
    const { radius, polar_angle, deltaz } = findHelixParameters(L, theta, tau);
    if (deltaz === 0) continue;  // avoid division by zero if it ever occurs
    const Nsegments = Math.floor(Math.abs(desiredLength) / Math.abs(deltaz));
    const actual_length = deltaz * Nsegments;
    const donut_hole_diameter = Math.abs(2 * radius - pipeDiam);
    const turns_per_segment = polar_angle / (2 * Math.PI);
    const Nturns = turns_per_segment * Nsegments;
    const outer_diameter = pipeDiam + 2 * radius;
    results.push({
      tau: tau,
      Nsegments: Nsegments,
      actual_length: actual_length,
      Nturns: Nturns,
      outer_diameter: outer_diameter,
      donut_hole_diameter: donut_hole_diameter,
      deltaz: deltaz,
      turns_per_segment: turns_per_segment
    });
  }
  return results;
}
