
// helper_functions.js
// Ported from helper_functions.py to plain JavaScript (no external math libs).
// Vectors are plain arrays [x,y,z] or [x,y,z,w]. Matrices are 2D arrays.

// ---- Basic math helpers ----
const deg2rad = (deg) => deg * Math.PI / 180;
const nearlyEqual = (a, b, eps=1e-9) => Math.abs(a - b) <= eps;

function dot(a, b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
function cross(a, b) {
  return [
    a[1]*b[2] - a[2]*b[1],
    a[2]*b[0] - a[0]*b[2],
    a[0]*b[1] - a[1]*b[0],
  ];
}
function add(a, b) { return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]; }
function sub(a, b) { return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]; }
function scale(v, s) { return [v[0]*s, v[1]*s, v[2]*s]; }
function norm(v) { return Math.hypot(v[0], v[1], v[2]); }
function unit(v, eps=1e-12) {
  const n = norm(v);
  if (n < eps) throw new Error("Zero-length vector encountered.");
  return [v[0]/n, v[1]/n, v[2]/n];
}

// ---- Matrix helpers ----
function matMul(A, B) {
  // general matrix multiply
  const m = A.length, n = B[0].length, p = B.length;
  const C = Array.from({length: m}, () => Array(n).fill(0));
  for (let i=0;i<m;i++){
    for (let k=0;k<p;k++){
      const aik = A[i][k];
      for (let j=0;j<n;j++){
        C[i][j] += aik * B[k][j];
      }
    }
  }
  return C;
}
function matVec(M, v) {
  // v is column vector as array; returns array
  const r = M.length, c = M[0].length;
  if (v.length !== c) throw new Error("Dimension mismatch in matVec.");
  const out = Array(r).fill(0);
  for (let i=0;i<r;i++){
    let s = 0;
    for (let j=0;j<c;j++) s += M[i][j]*v[j];
    out[i] = s;
  }
  return out;
}
function eye(n) {
  const I = Array.from({length:n},(_,i)=>{
    const row = Array(n).fill(0); row[i]=1; return row;
  });
  return I;
}

// ---- Homogeneous ↔ Cartesian ----
function homogeneous_to_cartesian(v) {
  if (v.length !== 4) throw new Error("Input must be a 4D vector.");
  if (nearlyEqual(v[3], 0)) throw new Error("The last component cannot be zero.");
  return [v[0]/v[3], v[1]/v[3], v[2]/v[3]];
}
const toCart = homogeneous_to_cartesian;

// ---- Transformation construction ----
function transformation_matrix(L, theta, tau) {
  const tau_rad = deg2rad(tau);
  const theta_rad = deg2rad(theta);

  const T1 = [
    [1,0,0,-L],
    [0,1,0, 0],
    [0,0,1, 0],
    [0,0,0, 1],
  ];
  const c1 = Math.cos(-theta_rad/2), s1 = Math.sin(-theta_rad/2);
  const R1 = [
    [ c1,-s1,0,0],
    [ s1, c1,0,0],
    [  0,  0,1,0],
    [  0,  0,0,1],
  ];
  const ct = Math.cos(tau_rad), st = Math.sin(tau_rad);
  const R2 = [
    [1,0, 0,0],
    [0,ct,-st,0],
    [0,st, ct,0],
    [0,0, 0,1],
  ];
  const ct_i = Math.cos(-tau_rad), st_i = Math.sin(-tau_rad);
  const R2_inv = [
    [1,0, 0,0],
    [0,ct_i,-st_i,0],
    [0,st_i, ct_i,0],
    [0,0, 0,1],
  ];
  const c3 = Math.cos(theta_rad/2), s3 = Math.sin(theta_rad/2);
  const R3 = [
    [ c3,-s3,0,0],
    [ s3, c3,0,0],
    [  0,  0,1,0],
    [  0,  0,0,1],
  ];
  const T2 = [
    [1,0,0,L],
    [0,1,0,0],
    [0,0,1,0],
    [0,0,0,1],
  ];

  // M = R1 @ R2 @ R1 @ T1
  let M = matMul(R1, matMul(R2, matMul(R1, T1)));
  // inverse_M_from_scratch = T2@R3@R2_inv@R3
  let inverse_M = matMul(T2, matMul(R3, matMul(R2_inv, R3)));
  return {M, inverse_M};
}

function generate_next_point(point, Transformation) {
  let p4;
  if (point.length === 3) p4 = [...point, 1];
  else if (point.length === 4) p4 = point.slice();
  else throw new Error("point must be a 3D or 4D vector.");
  const v = matVec(Transformation, p4);
  return homogeneous_to_cartesian(v);
}

function generate_points_on_helix(n_points, point, M) {
  if (n_points < 1) throw new Error("n_points must be at least 1.");
  let p = (point.length === 4) ? homogeneous_to_cartesian(point) :
          (point.length === 3) ? point.slice() :
          (()=>{throw new Error("point must be 3D/4D");})();
  const pts = [p];
  for (let i=0;i<n_points-1;i++){
    p = generate_next_point(p, M);
    pts.push(p);
  }
  return pts; // array of [x,y,z]
}

// ---- Geometry helpers ----
function _bisector(vertex, p1, p2) {
  const u1 = unit(sub(p1, vertex));
  const u2 = unit(sub(p2, vertex));
  let d = add(u1, u2);
  if (norm(d) < 1e-9) d = sub(u1, u2);
  return unit(d);
}

// shortest line segment between two 3D lines (unit directions)
function closest_line_between_lines(P1, u, P2, v, eps=1e-12) {
  const a = 1.0;
  const b = dot(u,v);
  const c = 1.0;
  const w0 = sub(P1, P2);
  const d = dot(u, w0);
  const e = dot(v, w0);
  const denom = a*c - b*b;

  let Pq, Pr, w;
  if (denom < eps) {
    // parallel
    w = sub(w0, scale(u, d));
    if (norm(w) < eps) return {M: P1.slice(), dir: u.slice()};
    Pq = P1.slice();
    Pr = add(P2, scale(v, e));
  } else {
    const t = (b*e - c*d) / denom;
    const s = (a*e - b*d) / denom;
    Pq = add(P1, scale(u, t));
    Pr = add(P2, scale(v, s));
    w = sub(Pr, Pq);
  }
  const M = scale(add(Pq, Pr), 0.5);
  if (norm(w) < eps) {
    return {M, dir: unit(cross(u, [1,0,0]))};
  }
  return {M, dir: unit(w)};
}

function connecting_line_for_bisectors(a,b,c,d){
  const uB = _bisector(b, a, c); // line at B
  const uC = _bisector(c, b, d); // line at C
  const {M, dir} = closest_line_between_lines(b, uB, c, uC);
  return {point: M, direction: dir};
}

// NOTE: The original Python also verified the eigenvector of the 3x3 rotation
// matrix has eigenvalue 1. We omit that expensive eigen computation here and
// rely on the bisector construction, which yields the helix axis reliably.
function find_helix_axis(point, M4) {
  const pts4 = generate_points_on_helix(4, point, M4);
  const [a,b,c,d] = pts4;
  const {point: P, direction: dir} = connecting_line_for_bisectors(a,b,c,d);
  return {axis_origin: P, axis_direction: dir};
}

function rotation_to_z(vin) {
  const v = unit(vin);
  const z = [0,0,1];
  if (Math.abs(v[0]-0)<1e-12 && Math.abs(v[1]-0)<1e-12 && Math.abs(v[2]-1)<1e-12) {
    return eye(3);
  }
  if (Math.abs(v[0]-0)<1e-12 && Math.abs(v[1]-0)<1e-12 && Math.abs(v[2]+1)<1e-12) {
    return [[-1,0,0],[0,-1,0],[0,0,1]];
  }
  let axis = cross(v, z);
  axis = unit(axis);
  const angle = Math.acos(Math.max(-1, Math.min(1, dot(v, z))));
  const K = [
    [0,        -axis[2],  axis[1]],
    [axis[2],   0,       -axis[0]],
    [-axis[1],  axis[0],  0      ],
  ];
  const I = eye(3);
  // R = I + sin(a)K + (1-cos(a))K^2
  const K2 = matMul(K, K);
  const R = addMat(addMat(I, scaleMat(K, Math.sin(angle))), scaleMat(K2, (1-Math.cos(angle))));
  return R;
}
function addMat(A,B){
  return A.map((row,i)=>row.map((v,j)=>v+B[i][j]));
}
function scaleMat(A,s){
  return A.map(row=>row.map(v=>v*s));
}

function translation_to_origin(point){
  if (point.length!==3) throw new Error("point must be a 3D vector.");
  const T = eye(4);
  T[0][3] = -point[0];
  T[1][3] = -point[1];
  T[2][3] = -point[2];
  return T;
}

function transform_points(points, axis_origin, axis_direction){
  // points: array of [x,y,z]
  const shifted = points.map(p => sub(p, axis_origin));
  const R = rotation_to_z(axis_direction); // 3x3
  const rot = shifted.map(p => {
    const x = R[0][0]*p[0] + R[0][1]*p[1] + R[0][2]*p[2];
    const y = R[1][0]*p[0] + R[1][1]*p[1] + R[1][2]*p[2];
    const z = R[2][0]*p[0] + R[2][1]*p[1] + R[2][2]*p[2];
    return [x,y,z];
  });
  const z0 = rot[0][2];
  return rot.map(p => [p[0], p[1], p[2]-z0]);
}

function turn_spiral_about_z_axis(points){
  const R = [[-1,0,0],[0,-1,0],[0,0,1]];
  return points.map(p=>[
    R[0][0]*p[0]+R[0][1]*p[1]+R[0][2]*p[2],
    R[1][0]*p[0]+R[1][1]*p[1]+R[1][2]*p[2],
    R[2][0]*p[0]+R[2][1]*p[1]+R[2][2]*p[2],
  ]);
}

function point_line_distance(P, A, v){
  const num = norm(cross(sub(P,A), v));
  const den = norm(v);
  return num/den;
}

function angle(v1, v2){
  const u1 = unit(v1), u2 = unit(v2);
  const dp = Math.max(-1, Math.min(1, dot(u1,u2)));
  return (Math.acos(dp) * 180 / Math.PI);
}

function find_number_of_segments_for_length(points_in, length){
  const delta = [];
  for (let i=0;i<points_in.length-1;i++){
    delta.push(sub(points_in[i+1], points_in[i]));
  }
  const deltaz = delta[0][2];
  const number = Math.ceil(length / deltaz);
  return {length, deltaz, number, actual: deltaz*number};
}

function full_helix_calculation(L, theta, tau, {length=null, segments=null}={}){
  if (length===null && segments===null){
    console.warn("Either length or segments should be provided. Defaulting to segments = 25.");
    segments = 25;
  }
  if (length!==null && length<=0){
    console.warn("Length is invalid, defaulting to automatic length calculation");
    length = null;
  }

  const {M} = transformation_matrix(L, theta, tau);
  const pts4 = generate_points_on_helix(4, [0,0,0,1], M);
  const {axis_origin, axis_direction} = find_helix_axis(pts4[0], M);
  const helix_radius = point_line_distance(pts4[0], axis_origin, axis_direction);

  let segments_from_length=null, actual_length_from_length=null, deltaz=null;
  let ptsN, ptsStraight;

  if (length!==null && segments===null){
    ptsN = generate_points_on_helix(4, [0,0,0,1], M);
    const tmpStraight = transform_points(ptsN, axis_origin, axis_direction);
    const calc = find_number_of_segments_for_length(tmpStraight, length);
    ({number: segments_from_length, actual: actual_length_from_length, deltaz} = calc);
    segments = segments_from_length;
  }

  ptsN = generate_points_on_helix(segments, [0,0,0,1], M);
  ptsStraight = transform_points(ptsN, axis_origin, axis_direction);

  let actual_length;
  if (length===null && segments!==null){
    actual_length = Math.abs(ptsStraight[ptsStraight.length-1][2] - ptsStraight[0][2]);
    deltaz = Math.abs(ptsStraight[1][2] - ptsStraight[0][2]);
    return {points_straight: ptsStraight, helix_radius, segments, actual_length, deltaz};
  }

  if (length!==null && segments!==null){
    actual_length = Math.abs(ptsStraight[ptsStraight.length-1][2] - ptsStraight[0][2]);
    if (segments_from_length===null){
      const tmp = find_number_of_segments_for_length(ptsStraight, length);
      segments_from_length = tmp.number;
      actual_length_from_length = tmp.actual;
      deltaz = tmp.deltaz;
    }
    return {points_straight: ptsStraight, helix_radius, segments, actual_length, deltaz, segments_from_length, actual_length_from_length};
  }

  throw new Error("full_helix_calculation: unexpected branch");
}

function center_length_from_outer(outer_length, cut_angle, diameter){
  const cut_angle_rad = deg2rad(cut_angle);
  return outer_length - diameter * Math.tan(cut_angle_rad);
}

function polar_angle_between(p1, p2){
  const [x1,y1] = p1, [x2,y2] = p2;
  const t1 = Math.atan2(y1, x1);
  const t2 = Math.atan2(y2, x2);
  let d = t2 - t1;
  d = (d + Math.PI) % (2*Math.PI) - Math.PI;
  return d;
}
function find_number_of_turns(polar_angle, number){
  return polar_angle * number / (2*Math.PI);
}

function find_helix_parameters(L, theta, tau){
  const {M} = transformation_matrix(L, theta, tau);
  const initial_point = [0,0,0];
  const {axis_origin, axis_direction} = find_helix_axis(initial_point, M);
  const radius = point_line_distance(initial_point, axis_origin, axis_direction);
  const pts = generate_points_on_helix(2, initial_point, M);
  const ptsStraight = transform_points(pts, axis_origin, axis_direction);
  const polar_angle = polar_angle_between(ptsStraight[1], ptsStraight[0]);
  const deltaz = ptsStraight[1][2] - ptsStraight[0][2];
  return {radius, polar_angle, deltaz};
}

function tangent_angle(radius, polar_angle, deltaz){
  return Math.acos(deltaz / Math.sqrt((radius*polar_angle)**2 + deltaz**2)) * 180 / Math.PI;
}

function is_helix_loose_enough(L, theta, tau, pipe_diameter, eps=1){
  const {radius, polar_angle, deltaz} = find_helix_parameters(L, theta, tau);
  const tang_ang = tangent_angle(radius, polar_angle, deltaz);
  const {M} = transformation_matrix(L, theta, tau);
  const pts = generate_points_on_helix(2, [0,0,0], M);
  const {axis_origin, axis_direction} = find_helix_axis(pts[0], M);
  const discrete_angle = angle(sub(pts[1], pts[0]), axis_direction);
  const distance = 2 * Math.PI * deltaz * Math.cos(discrete_angle * Math.PI / 180);
  return !(distance < 2 * pipe_diameter + eps);
}

function is_helix_wide_enough(radius, pipe_diameter, eps=1){
  return radius > pipe_diameter / 2 + eps;
}

// Export for global usage
window.HelixHelpers = {
  transformation_matrix,
  generate_points_on_helix,
  full_helix_calculation,
  center_length_from_outer,
  polar_angle_between,
  find_number_of_turns,
  is_helix_loose_enough,
  is_helix_wide_enough,
  turn_spiral_about_z_axis,
};
