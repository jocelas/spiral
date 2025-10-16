


function homogeneousToCartesian(v) {
  // Convert a homogeneous coordinate to a Cartesian coordinate.
  //
  // Args:
  //   v (Array): A 4D vector in homogeneous coordinates.
  //
  // Returns:
  //   Array: A 3D vector in Cartesian coordinates.

  if (!Array.isArray(v) || v.length !== 4) {
    throw new Error("Input must be a 4D vector.");
  }

  const w = v[3];
  if (w === 0) {
    throw new Error("The last component of the homogeneous coordinate cannot be zero.");
  }

  return [v[0] / w, v[1] / w, v[2] / w];
}

function transformationMatrix(L, theta, tau) {
  // Convert degrees to radians
  const tauRad = (tau * Math.PI) / 180;
  const thetaRad = (theta * Math.PI) / 180;

  // 4×4 translation matrix: translate B to origin
  const T1 = [
    [1, 0, 0, -L],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
  ];

  // Rotate by -θ/2 around z-axis
  const R1 = [
    [Math.cos(-thetaRad / 2), -Math.sin(-thetaRad / 2), 0, 0],
    [Math.sin(-thetaRad / 2), Math.cos(-thetaRad / 2), 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
  ];

  // Rotate by τ around x-axis
  const R2 = [
    [1, 0, 0, 0],
    [0, Math.cos(tauRad), -Math.sin(tauRad), 0],
    [0, Math.sin(tauRad), Math.cos(tauRad), 0],
    [0, 0, 0, 1],
  ];

  // Rotate by -τ around x-axis
  const R2_inv = [
    [1, 0, 0, 0],
    [0, Math.cos(-tauRad), -Math.sin(-tauRad), 0],
    [0, Math.sin(-tauRad), Math.cos(-tauRad), 0],
    [0, 0, 0, 1],
  ];

  // Rotate by θ/2 around z-axis
  const R3 = [
    [Math.cos(thetaRad / 2), -Math.sin(thetaRad / 2), 0, 0],
    [Math.sin(thetaRad / 2), Math.cos(thetaRad / 2), 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
  ];

  // Translate back
  const T2 = [
    [1, 0, 0, L],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
  ];

  // Helper for 4×4 matrix multiplication
  function matMul(A, B) {
    const result = Array.from({ length: 4 }, () => Array(4).fill(0));
    for (let i = 0; i < 4; i++) {
      for (let j = 0; j < 4; j++) {
        for (let k = 0; k < 4; k++) {
          result[i][j] += A[i][k] * B[k][j];
        }
      }
    }
    return result;
  }

  // Compute main and inverse transformations
  const M = matMul(matMul(R1, R2), matMul(R1, T1));
  const inverseMFromScratch = matMul(matMul(T2, R3), matMul(R2_inv, R3));

  // Return both matrices
  return { M, inverseM: inverseMFromScratch };
}

function generateNextPoint(point, Transformation) {
  if (!Array.isArray(point)) {
    throw new Error("Point must be an array.");
  }

  let p;
  if (point.length === 3) {
    p = [...point, 1];
  } else if (point.length === 4) {
    p = [...point];
  } else {
    throw new Error("Point must be a 3D or 4D vector.");
  }

  // 4×4 matrix × 4D vector multiplication
  const v = new Array(4).fill(0);
  for (let i = 0; i < 4; i++) {
    for (let j = 0; j < 4; j++) {
      v[i] += Transformation[i][j] * p[j];
    }
  }

  return homogeneousToCartesian(v);
}

function generatePointsOnHelix(segments, Transformation) {
  if (segments < 1) {
    throw new Error("segments must be at least 1.");
  }

  let point = [0, 0, 0];
  const points = [point];

  for (let i = 0; i < segments; i++) {
    point = generateNextPoint(point, Transformation);
    points.push(point);
  }

  return points;
}

function makeLine(point, vector, length = 100) {
  if (point.length !== 3 || vector.length !== 3) {
    throw new Error("Both point and vector must be 3D vectors.");
  }

  // normalize vector
  const norm = Math.hypot(...vector);
  if (norm === 0) throw new Error("Direction vector cannot be zero.");

  const scaled = vector.map(v => (v / norm) * (length / 2));

  const p1 = point.map((p, i) => p - scaled[i]);
  const p2 = point.map((p, i) => p + scaled[i]);

  return [p1, p2];
}


function unit(v, eps = 1e-12) {
  const n = Math.hypot(...v);
  if (n < eps) throw new Error("Zero-length vector encountered.");
  return v.map(x => x / n);
}

function _bisector(vertex, p1, p2, eps = 1e-9) {
  const u1 = unit(p1.map((x, i) => p1[i] - vertex[i]));
  const u2 = unit(p2.map((x, i) => p2[i] - vertex[i]));
  let d = u1.map((x, i) => x + u2[i]);

  const normD = Math.hypot(...d);
  if (normD < eps) {
    // Degenerate (straight angle) – fall back to external bisector
    d = u1.map((x, i) => x - u2[i]);
  }

  return unit(d);
}

function dot(u, v) {
  return u.reduce((sum, x, i) => sum + x * v[i], 0);
}

function cross(u, v) {
  return [
    u[1] * v[2] - u[2] * v[1],
    u[2] * v[0] - u[0] * v[2],
    u[0] * v[1] - u[1] * v[0],
  ];
}

function closestLineBetweenLines(P1, u, P2, v, eps = 1e-12) {
  const a = 1.0;
  const b = dot(u, v);
  const c = 1.0;
  const w0 = P1.map((x, i) => x - P2[i]);
  const d = dot(u, w0);
  const e = dot(v, w0);
  const denom = a * c - b * b; // = 1 - (u·v)^2

  let Pq, Pr, w;

  if (denom < eps) {
    // Lines are parallel (or nearly)
    w = w0.map((x, i) => x - d * u[i]);
    const normW = Math.hypot(...w);
    if (normW < eps) {
      // Same line
      return [P1.slice(), u.slice()];
    }
    Pq = P1.slice();
    Pr = P2.map((x, i) => x + e * v[i]);
  } else {
    const t = (b * e - c * d) / denom;
    const s = (a * e - b * d) / denom;
    Pq = P1.map((x, i) => x + t * u[i]);
    Pr = P2.map((x, i) => x + s * v[i]);
    w = Pr.map((x, i) => x - Pq[i]);
  }

  const M = Pq.map((x, i) => 0.5 * (x + Pr[i]));
  const normW = Math.hypot(...w);

  if (normW < eps) {
    // Lines intersect – return intersection point and any perpendicular direction
    return [M, unit(cross(u, [1, 0, 0]))];
  }

  const dirVec = unit(w);
  return [M, dirVec];
}

function connectingLineForBisectors(a, b, c, d) {
  const uB = _bisector(b, a, c);
  const uC = _bisector(c, b, d);

  const [M, dirVec] = closestLineBetweenLines(b, uB, c, uC);
  return [M, dirVec];
}

function findAxisDirEv(M_in) {
  const M = M_in.slice(0, 3).map(row => row.slice(0, 3));

  // Axis of rotation for a proper rotation matrix:
  const axis = [
    M[2][1] - M[1][2],
    M[0][2] - M[2][0],
    M[1][0] - M[0][1],
  ];
  const n = Math.hypot(...axis);
  if (n < 1e-12) throw new Error("Rotation axis could not be determined.");
  return axis.map(x => x / n);
}

function findHelixAxis(transformation) {
  const points = generatePointsOnHelix(3, transformation);
  const [a, b, c, d] = points;

  const [point, direction] = connectingLineForBisectors(a, b, c, d);
  const evDirection = findAxisDirEv(transformation);

  const dirNorm = unit(direction);
  const evNorm = unit(evDirection);

  const dotVal = Math.abs(dot(dirNorm, evNorm));

  if (!(Math.abs(dotVal - 1) < 1e-5)) {
    throw new Error(
      `Direction from bisectors does not match eigenvector direction\n ${dotVal}`
    );
  }

  return [point, direction];
}


function rotationToZ(v) {
  let vec = v.map(Number);
  const norm = Math.hypot(...vec);
  vec = vec.map(x => x / norm);
  const z = [0, 0, 1];

  if (vec.every((x, i) => Math.abs(x - z[i]) < 1e-9)) {
    // already aligned
    return [
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
    ];
  }

  if (vec.every((x, i) => Math.abs(x + z[i]) < 1e-9)) {
    // opposite direction — rotate 180° around any perpendicular axis
    return [
      [-1, 0, 0],
      [0, -1, 0],
      [0, 0, 1],
    ];
  }

  let axis = cross(vec, z);
  axis = unit(axis);
  const angle = Math.acos(Math.min(Math.max(dot(vec, z), -1.0), 1.0));

  const [ax, ay, az] = axis;
  const K = [
    [0, -az, ay],
    [az, 0, -ax],
    [-ay, ax, 0],
  ];

  // K² = K @ K
  const K2 = Array.from({ length: 3 }, () => Array(3).fill(0));
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      for (let k = 0; k < 3; k++) {
        K2[i][j] += K[i][k] * K[k][j];
      }
    }
  }

  // R = I + sin(angle)*K + (1 - cos(angle))*K²
  const R = Array.from({ length: 3 }, (_, i) =>
    Array.from({ length: 3 }, (_, j) =>
      (i === j ? 1 : 0) +
      Math.sin(angle) * K[i][j] +
      (1 - Math.cos(angle)) * K2[i][j]
    )
  );

  // Determinant check (rough numeric)
  const det =
    R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1]) -
    R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0]) +
    R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]);
  if (Math.abs(det - 1) > 1e-6)
    throw new Error("Rotation to z was not correctly calculated.");

  return R;
}

function translationToOrigin(point) {
  if (point.length !== 3) throw new Error("point must be a 3D vector.");
  const T = [
    [1, 0, 0, -point[0]],
    [0, 1, 0, -point[1]],
    [0, 0, 1, -point[2]],
    [0, 0, 0, 1],
  ];
  return T;
}

function transformHelixToZAxis(points, axisOrigin, axisDirection) {
  // subtract origin
  const shifted = points.map(p => p.map((x, i) => x - axisOrigin[i]));

  // rotate so axisDirection aligns with z-axis
  const R = rotationToZ(axisDirection);

  // multiply each point by Rᵀ
  const rotated = shifted.map(p => {
    return [
      dot([R[0][0], R[1][0], R[2][0]], p),
      dot([R[0][1], R[1][1], R[2][1]], p),
      dot([R[0][2], R[1][2], R[2][2]], p),
    ];
  });

  // shift to place first segment on xy-plane
  const zShift = rotated[0][2];
  return rotated.map(p => [p[0], p[1], p[2] - zShift]);
}

function makeSecondHelix(points) {
  const R = [
    [-1, 0, 0],
    [0, -1, 0],
    [0, 0, 1],
  ];

  // points @ Rᵀ
  return points.map(p => [
    dot([R[0][0], R[1][0], R[2][0]], p),
    dot([R[0][1], R[1][1], R[2][1]], p),
    dot([R[0][2], R[1][2], R[2][2]], p),
  ]);
}

function pointLineDistance(P, A, v) {
  // Distance between point P and line passing through A with direction v
  const crossProd = cross(
    P.map((x, i) => x - A[i]),
    v
  );
  return Math.hypot(...crossProd) / Math.hypot(...v);
}

function angle(v1, v2) {
  const v1_u = unit(v1);
  const v2_u = unit(v2);
  const dotProduct = Math.min(Math.max(dot(v1_u, v2_u), -1.0), 1.0);
  const angleRad = Math.acos(dotProduct);
  return (angleRad * 180) / Math.PI;
}

function findNumberOfSegmentsForLength(pointsIn, length) {
  // assumes pointsIn is an array of [x, y, z] arrays
  const deltaPoints = [];
  for (let i = 0; i < pointsIn.length - 1; i++) {
    const a = pointsIn[i];
    const b = pointsIn[i + 1];
    deltaPoints.push(b.map((x, j) => x - a[j]));
  }

  // check z increments are constant
  const deltazValues = deltaPoints.map(dp => dp[2]);
  const allEqual = deltazValues.every(z => Math.abs(z - deltazValues[0]) < 1e-9);
  if (!allEqual) throw new Error("Z differences are not consistent.");

  const deltaz = deltazValues[0];
  const number = Math.floor(length / deltaz);

  return [length, deltaz, number, deltaz * number];
}



function generateStraightHelix(L, theta, tau, segments) {
  const { M } = transformationMatrix(L, theta, tau);
  const points = generatePointsOnHelix(segments, M);
  const [axisOrigin, axisDirection] = findHelixAxis(M);
  const pointsStraight = transformHelixToZAxis(points, axisOrigin, axisDirection);
  return pointsStraight;
}

function centerLengthFromOuter(outerLength, cutAngle, diameter) {
  const cutAngleRad = (cutAngle * Math.PI) / 180;
  return outerLength - diameter * Math.tan(cutAngleRad);
}

function polarAngleBetween(p1, p2) {
  const [x1, y1] = p1;
  const [x2, y2] = p2;
  const theta1 = Math.atan2(y1, x1);
  const theta2 = Math.atan2(y2, x2);
  let dtheta = theta2 - theta1;
  // normalize to [-π, π)
  dtheta = ((dtheta + Math.PI) % (2 * Math.PI)) - Math.PI;
  return dtheta; // in radians
}

function findNumberOfTurns(polarAngle, number) {
  return (polarAngle * number) / (2 * Math.PI);
}

function findHelixParameters(L, theta, tau) {
  const pointsStraight = generateStraightHelix(L, theta, tau, 4);

  const radius = pointLineDistance(pointsStraight[0], [0, 0, 0], [0, 0, 1]);
  const radiusTest = pointLineDistance(pointsStraight[1], [0, 0, 0], [0, 0, 1]);

  if (Math.abs(radius - radiusTest) > 1)
    throw new Error("radius calculation wrong");

  const polarAngle = polarAngleBetween(pointsStraight[1], pointsStraight[0]);

  const deltaz = pointsStraight[1][2] - pointsStraight[0][2];
  const deltazTest =
    pointsStraight[pointsStraight.length - 1][2] -
    pointsStraight[pointsStraight.length - 2][2];

  if (Math.abs(deltaz - deltazTest) > 1)
    // throw new Error("deltaz calculation wrong");

  return [radius, polarAngle, deltaz];
}

function analyticalAngleToXYPlane(radius, polarAngle, deltaz) {
  // Calculates angle to xy-plane of analytical spiral
  const pitch = (deltaz * 2 * Math.PI) / polarAngle; // height increase per full turn
  const tanAlpha = pitch / (2 * Math.PI * radius);
  const alpha = Math.atan(tanAlpha);
  return (alpha * 180) / Math.PI;
}

function segmentAngleToXYPlane(L, deltaz) {
  const sinAlpha = deltaz / L;
  const alpha = Math.asin(sinAlpha);
  return (alpha * 180) / Math.PI;
}

function isHelixLooseEnough(L, theta, tau, pipeDiameter, eps = 1) {
  const [radius, polarAngle, deltaz] = findHelixParameters(L, theta, tau);
  const nfull = (2 * Math.PI) / polarAngle;
  return nfull * deltaz > 2 * (pipeDiameter + eps);
}

function isHelixWideEnough(L, theta, tau, pipeDiameter, eps = 1) {
  const [radius, polarAngle, deltaz] = findHelixParameters(L, theta, tau);
  return 2 * radius > pipeDiameter + eps;
}

function isHelixAllowed(L, theta, tau, pipeDiameter, epsWide = 3, epsLoose = 3) {
  return (
    (isHelixLooseEnough(L, theta, tau, pipeDiameter, epsLoose) ? 1 : 0) *
    (isHelixWideEnough(L, theta, tau, pipeDiameter, epsWide) ? 1 : 0)
  );
}

function findAllowedTwistRange(
  L,
  theta,
  pipeDiameter,
  epsLoose = 3,
  epsWide = 3,
  tauMin = 0.0,
  tauMax = 180.0,
  tol = 1e-4
) {
  function f(tau) {
    return !!isHelixAllowed(L, theta, tau, pipeDiameter, epsWide, epsLoose);
  }

  // Coarse scan
  const N = 50;
  const taus = Array.from({ length: N }, (_, i) => tauMin + (i * (tauMax - tauMin)) / (N - 1));
  const vals = taus.map(f);

  if (!vals.some(Boolean)) return [null, null];

  const idxTrue = vals
    .map((val, i) => (val ? i : -1))
    .filter(i => i !== -1);
  const iStart = idxTrue[0];
  const iEnd = idxTrue[idxTrue.length - 1];

  let loLeft = taus[Math.max(iStart - 1, 0)];
  let hiLeft = taus[iStart];
  let loRight = taus[iEnd];
  let hiRight = taus[Math.min(iEnd + 1, N - 1)];

  // Refine lower edge
  while (hiLeft - loLeft > tol) {
    const mid = 0.5 * (loLeft + hiLeft);
    if (f(mid)) hiLeft = mid;
    else loLeft = mid;
  }
  const tauLo = hiLeft;

  // Refine upper edge
  while (hiRight - loRight > tol) {
    const mid = 0.5 * (loRight + hiRight);
    if (f(mid)) loRight = mid;
    else hiRight = mid;
  }
  const tauHi = loRight;

  return [tauLo, tauHi];
}

function listPossibleHelices(
  L,
  theta,
  pipeDiameter,
  desiredLength,
  tauIncrement = 1,
  bottomOffset = 1,
  topOffset = 0
) {
  let [tauMin, tauMax] = findAllowedTwistRange(L, theta, pipeDiameter, 1, 1);

  tauMin = Math.ceil(tauMin) + bottomOffset;
  tauMax = Math.floor(tauMax) - topOffset;

  const listOfTaus = [];
  for (let tau = tauMin; tau < tauMax; tau += tauIncrement) listOfTaus.push(tau);

  const outputs = [];

  for (const tau of listOfTaus) {
    const [radius, polarAngle, deltaz] = findHelixParameters(L, theta, tau);
    const nSegmentsForLength = Math.floor(desiredLength / deltaz);
    const actualLength = deltaz * nSegmentsForLength;

    const donutHoleDiameter = Math.abs(radius * 2 - pipeDiameter);
    const nTurns = findNumberOfTurns(polarAngle, nSegmentsForLength);
    const turnsPerSegment = findNumberOfTurns(polarAngle, 1);
    const outerDiameter = pipeDiameter + radius * 2;

    outputs.push([
      tau,
      nSegmentsForLength,
      actualLength,
      nTurns,
      outerDiameter,
      donutHoleDiameter,
      deltaz,
      turnsPerSegment,
    ]);
  }

  return outputs;
}
