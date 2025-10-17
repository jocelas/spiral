
// helper_functions.js — translated from the updated helper_functions.py
// Data conventions: vectors are arrays [x,y,z] or [x,y,z,w], matrices are nested arrays.

(function(){
  const EPS = 1e-12;

  // ---------- Basic math helpers ----------
  const deg2rad = (deg) => deg * Math.PI / 180;
  const nearlyEqual = (a,b,eps=1e-9) => Math.abs(a-b) <= eps;

  function add(a,b){ return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]; }
  function sub(a,b){ return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]; }
  function scale(v,s){ return [v[0]*s, v[1]*s, v[2]*s]; }
  function dot(a,b){ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
  function cross(a,b){
    return [
      a[1]*b[2] - a[2]*b[1],
      a[2]*b[0] - a[0]*b[2],
      a[0]*b[1] - a[1]*b[0],
    ];
  }
  function norm(v){ return Math.hypot(v[0], v[1], v[2]); }
  function unit(v, eps=1e-12){
    const n = norm(v);
    if (n < eps) throw new Error("Zero-length vector encountered.");
    return [v[0]/n, v[1]/n, v[2]/n];
  }

  // ---------- Matrix helpers ----------
  function eye(n){
    return Array.from({length:n},(_,i)=>{
      const r = Array(n).fill(0); r[i]=1; return r;
    });
  }
  function matMul(A,B){
    const m=A.length, p=B.length, n=B[0].length;
    const C = Array.from({length:m},()=>Array(n).fill(0));
    for(let i=0;i<m;i++){
      for(let k=0;k<p;k++){
        const aik=A[i][k];
        for(let j=0;j<n;j++) C[i][j]+=aik*B[k][j];
      }
    }
    return C;
  }
  function matVec(M,v){
    const r=M.length, c=M[0].length;
    if (v.length!==c) throw new Error("matVec dimension mismatch");
    const out = Array(r).fill(0);
    for(let i=0;i<r;i++){
      let s=0;
      for(let j=0;j<c;j++) s+=M[i][j]*v[j];
      out[i]=s;
    }
    return out;
  }
  function det3(R){
    // determinant of 3x3
    return R[0][0]*(R[1][1]*R[2][2]-R[1][2]*R[2][1])
         - R[0][1]*(R[1][0]*R[2][2]-R[1][2]*R[2][0])
         + R[0][2]*(R[1][0]*R[2][1]-R[1][1]*R[2][0]);
  }

  // ---------- Homogeneous → Cartesian ----------
  function homogeneous_to_cartesian(v){
    if (v.length!==4) throw new Error("Input must be a 4D vector.");
    if (nearlyEqual(v[3],0)) throw new Error("The last component of the homogeneous coordinate cannot be zero.");
    return [v[0]/v[3], v[1]/v[3], v[2]/v[3]];
  }

  // ---------- Transformation matrix ----------
  function transformation_matrix(L, theta, tau){
    const tau_rad = Math.PI * tau / 180;
    const theta_rad = Math.PI * theta / 180;

    const T1 = [
      [1,0,0,-L],
      [0,1,0, 0],
      [0,0,1, 0],
      [0,0,0, 1],
    ];

    const c1=Math.cos(-theta_rad/2), s1=Math.sin(-theta_rad/2);
    const R1 = [
      [ c1,-s1,0,0],
      [ s1, c1,0,0],
      [  0,  0,1,0],
      [  0,  0,0,1],
    ];

    const ct=Math.cos(tau_rad), st=Math.sin(tau_rad);
    const R2=[
      [1,0, 0,0],
      [0,ct,-st,0],
      [0,st, ct,0],
      [0,0, 0,1],
    ];

    const c3=Math.cos(theta_rad/2), s3=Math.sin(theta_rad/2);
    const R3=[
      [ c3,-s3,0,0],
      [ s3, c3,0,0],
      [  0,  0,1,0],
      [  0,  0,0,1],
    ];
    const T2=[
      [1,0,0,L],
      [0,1,0,0],
      [0,0,1,0],
      [0,0,0,1],
    ];

    const M = matMul(R1, matMul(R2, matMul(R1, T1)));

    const ct_i=Math.cos(-tau_rad), st_i=Math.sin(-tau_rad);
    const R2_inv=[
      [1,0, 0,0],
      [0,ct_i,-st_i,0],
      [0,st_i, ct_i,0],
      [0,0, 0,1],
    ];
    const inverse_M = matMul(T2, matMul(R3, matMul(R2_inv, R3)));
    return {M, inverse_M};
  }

  // ---------- Helix point generation ----------
  function generate_next_point(point, Transformation){
    let p4;
    if (point.length===3) p4 = [...point,1];
    else if (point.length===4) p4 = point.slice();
    else throw new Error("point must be a 3D or 4D vector.");
    const v = matVec(Transformation, p4);
    return homogeneous_to_cartesian(v);
  }

  function generate_points_on_helix(segments, Transformation){
    if (segments < 1) throw new Error("n_points must be at least 1.");
    let point = [0,0,0]; // starts at origin
    const pts = [point];
    for (let i=0;i<segments;i++){
      point = generate_next_point(point, Transformation);
      pts.push(point);
    }
    return pts; // length = segments+1
  }

  function make_line(point, vector, length=100){
    if (point.length!==3 || vector.length!==3) throw new Error("Both point and vector must be 3D vectors.");
    const v = unit(vector);
    const half = length/2;
    const seg = scale(v, half);
    return [ sub(point, seg), add(point, seg) ];
  }

  // ---------- Axis finding helpers ----------
  function _bisector(vertex, p1, p2){
    const u1 = unit(sub(p1, vertex));
    const u2 = unit(sub(p2, vertex));
    let d = add(u1, u2);
    if (norm(d) < 1e-9) d = sub(u1, u2);
    return unit(d);
  }

  function closest_line_between_lines(P1, u, P2, v, eps=1e-12){
    const a = 1.0;
    const b = (u[0]*v[0]+u[1]*v[1]+u[2]*v[2]);
    const c = 1.0;
    const w0 = sub(P1, P2);
    const d = (u[0]*w0[0]+u[1]*w0[1]+u[2]*w0[2]);
    const e = (v[0]*w0[0]+v[1]*w0[1]+v[2]*w0[2]);
    const denom = a*c - b*b;

    let Pq, Pr, w;
    if (denom < eps){
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
    if (norm(w) < eps){
      return {M, dir: unit(cross(u, [1,0,0]))};
    }
    return {M, dir: unit(w)};
  }

  function connecting_line_for_bisectors(a,b,c,d){
    const uB = _bisector(b, a, c);
    const uC = _bisector(c, b, d);
    const {M, dir} = closest_line_between_lines(b, uB, c, uC);
    return {point: M, direction: dir};
  }

  function find_axis_dir_ev(M_in){
    const R = [
      [M_in[0][0], M_in[0][1], M_in[0][2]],
      [M_in[1][0], M_in[1][1], M_in[1][2]],
      [M_in[2][0], M_in[2][1], M_in[2][2]],
    ];
    const A = [
      [R[0][0]-1, R[0][1],   R[0][2]],
      [R[1][0],   R[1][1]-1, R[1][2]],
      [R[2][0],   R[2][1],   R[2][2]-1],
    ];
    function crossRows(i,j){
      const ri = A[i], rj = A[j];
      return [
        ri[1]*rj[2]-ri[2]*rj[1],
        ri[2]*rj[0]-ri[0]*rj[2],
        ri[0]*rj[1]-ri[1]*rj[0],
      ];
    }
    let v = crossRows(0,1);
    if (norm(v) < 1e-8) v = crossRows(0,2);
    if (norm(v) < 1e-8) v = crossRows(1,2);
    if (norm(v) < 1e-8) v = [0,0,1];
    v = unit(v);
    const Rv = [
      R[0][0]*v[0] + R[0][1]*v[1] + R[0][2]*v[2],
      R[1][0]*v[0] + R[1][1]*v[1] + R[1][2]*v[2],
      R[2][0]*v[0] + R[2][1]*v[1] + R[2][2]*v[2],
    ];
    const err = norm(sub(Rv, v));
    if (err > 1e-4) {
      const n = norm(Rv);
      if (n>1e-12) v = [Rv[0]/n, Rv[1]/n, Rv[2]/n];
    }
    return v;
  }

  function find_helix_axis(Transformation){
    const pts = generate_points_on_helix(3, Transformation); // 4 points
    const [a,b,c,d] = pts;
    const {point, direction} = connecting_line_for_bisectors(a,b,c,d);
    const ev_direction = find_axis_dir_ev(Transformation);
    // Optional alignment check omitted in JS
    return {axis_origin: point, axis_direction: direction};
  }

  // ---------- Rotations & transforms ----------
  function rotation_to_z(v){
    let vn = unit(v);
    const z = [0,0,1];
    if (Math.abs(vn[0])<1e-12 && Math.abs(vn[1])<1e-12 && Math.abs(vn[2]-1)<1e-12) return eye(3);
    if (Math.abs(vn[0])<1e-12 && Math.abs(vn[1])<1e-12 && Math.abs(vn[2]+1)<1e-12) {
      return [[-1,0,0],[0,-1,0],[0,0,1]];
    }
    let axis = cross(vn, z);
    axis = unit(axis);
    const angle = Math.acos(Math.max(-1, Math.min(1, dot(vn, z))));
    const K = [
      [0,        -axis[2],  axis[1]],
      [axis[2],   0,       -axis[0]],
      [-axis[1],  axis[0],  0      ],
    ];
    const I = eye(3);
    const K2 = matMul(K, K);
    const sinA = Math.sin(angle), oneMinusCos = 1 - Math.cos(angle);
    const R = addMat(addMat(I, scaleMat(K, sinA)), scaleMat(K2, oneMinusCos));
    const d = R[0][0]*(R[1][1]*R[2][2]-R[1][2]*R[2][1])
            - R[0][1]*(R[1][0]*R[2][2]-R[1][2]*R[2][0])
            + R[0][2]*(R[1][0]*R[2][1]-R[1][1]*R[2][0]);
    // if (Math.abs(d-1) >= 1e-6) console.warn("Rotation det != 1", d);
    return R;
  }
  function addMat(A,B){ return A.map((row,i)=>row.map((v,j)=>v+B[i][j])); }
  function scaleMat(A,s){ return A.map(row=>row.map(v=>v*s)); }

  function translation_to_origin(point){
    if (point.length!==3) throw new Error("point must be a 3D vector.");
    const T = eye(4);
    T[0][3] = -point[0];
    T[1][3] = -point[1];
    T[2][3] = -point[2];
    return T;
  }

  function transform_helix_to_z_axis(points, axis_origin, axis_direction){
    const shifted = points.map(p => sub(p, axis_origin));
    const R = rotation_to_z(axis_direction);
    const rot = shifted.map(p => {
      const x = R[0][0]*p[0] + R[0][1]*p[1] + R[0][2]*p[2];
      const y = R[1][0]*p[0] + R[1][1]*p[1] + R[1][2]*p[2];
      const z = R[2][0]*p[0] + R[2][1]*p[1] + R[2][2]*p[2];
      return [x,y,z];
    });
    const z0 = rot[0][2];
    return rot.map(p => [p[0], p[1], p[2]-z0]);
  }

  function make_second_helix(points){
    const R = [[-1,0,0],[0,-1,0],[0,0,1]]; // 180° about z
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
    return Math.acos(dp) * 180 / Math.PI;
  }

  function find_number_of_segments_for_length(points_in, length){
    const delta = [];
    for (let i=0;i<points_in.length-1;i++){
      delta.push(sub(points_in[i+1], points_in[i]));
    }
    const deltaz = delta[0][2];
    const number = Math.floor(length / deltaz);
    return {length, deltaz, number, actual: deltaz * number};
  }

  // ---------- High-level helix ops ----------
  function generate_straight_helix(L, theta, tau, segments){
    const {M} = transformation_matrix(L, theta, tau);
    const points = generate_points_on_helix(segments, M);
    const {axis_origin, axis_direction} = find_helix_axis(M);
    const points_straight = transform_helix_to_z_axis(points, axis_origin, axis_direction);
    return points_straight;
  }

  function center_length_from_outer(outer_length, cut_angle, diameter){
    const cut_angle_rad = Math.PI * cut_angle / 180;
    return outer_length - diameter * Math.tan(cut_angle_rad);
  }

  function polar_angle_between(p1, p2){
    const [x1,y1] = p1, [x2,y2] = p2;
    const t1 = Math.atan2(y1, x1);
    const t2 = Math.atan2(y2, x2);
    let d = t2 - t1;
    d = (d + Math.PI) % (2*Math.PI) - Math.PI;
    return d; // radians
  }

  function find_number_of_turns(polar_angle, number){
    return polar_angle * number / (2*Math.PI);
  }

  function find_helix_parameters(L, theta, tau){
    const pts = generate_straight_helix(L, theta, tau, 4);
    const radius = point_line_distance(pts[0], [0,0,0], [0,0,1]);
    const radius_test = point_line_distance(pts[1], [0,0,0], [0,0,1]);
    const polar_angle = polar_angle_between(pts[1], pts[0]);
    const deltaz = pts[1][2] - pts[0][2];
    return {radius, polar_angle, deltaz};
  }

  function analytical_angle_to_xy_plane(radius, polar_angle, deltaz){
    const pitch = deltaz * 2 * Math.PI / polar_angle;
    const tanalpha = pitch / (2 * Math.PI * radius);
    const alpha = Math.atan(tanalpha);
    return alpha * 180 / Math.PI;
  }

  function segment_angle_to_xy_plane(L, deltaz){
    const sinalpha = deltaz / L;
    const alpha = Math.asin(sinalpha);
    return alpha * 180 / Math.PI;
  }

  function is_helix_loose_enough(L, theta, tau, pipe_diameter, eps=1){
    const {polar_angle, deltaz} = find_helix_parameters(L, theta, tau);
    const nfull = 2*Math.PI / polar_angle;
    return (nfull * deltaz > 2 * (pipe_diameter + eps));
  }

  function is_helix_wide_enough(L, theta, tau, pipe_diameter, eps=1){
    const {radius} = find_helix_parameters(L, theta, tau);
    return (2 * radius > pipe_diameter + eps);
  }

  function is_helix_allowed(L, theta, tau, pipe_diameter, eps_wide=3, eps_loose=3){
    return is_helix_loose_enough(L, theta, tau, pipe_diameter, eps_loose)
        && is_helix_wide_enough(L, theta, tau, pipe_diameter, eps_wide);
  }

  // function find_allowed_twist_range(L, theta, pipe_diameter,
  //                                   eps_loose=3, eps_wide=3,
  //                                   tau_min=0.0, tau_max=180.0, tol=1e-4){
  //   function f(tau){
  //     return is_helix_allowed(L, theta, tau, pipe_diameter, eps_wide, eps_loose);
  //   }
  //   const N = 50;
  //   const taus = Array.from({length:N},(_,i)=> tau_min + i*(tau_max - tau_min)/(N-1));
  //   const vals = taus.map(t=>f(t));

  //   if (!vals.some(Boolean)) return {tau_lo: null, tau_hi: null};

  //   const idx_true = [];
  //   for (let i=0;i<N;i++) if (vals[i]) idx_true.push(i);
  //   let i_start = idx_true[0], i_end = idx_true[idx_true.length-1];

  //   let lo_left = taus[Math.max(i_start-1, 0)], hi_left = taus[i_start];
  //   let lo_right = taus[i_end], hi_right = taus[Math.min(i_end+1, N-1)];

  //   while (hi_left - lo_left > tol){
  //     const mid = 0.5*(lo_left + hi_left);
  //     if (f(mid)) hi_left = mid; else lo_left = mid;
  //   }
  //   const tau_lo = hi_left;

  //   while (hi_right - lo_right > tol){
  //     const mid = 0.5*(lo_right + hi_right);
  //     if (f(mid)) lo_right = mid; else hi_right = mid;
  //   }
  //   const tau_hi = lo_right;
  //   return {tau_lo, tau_hi};
  // }

  function find_allowed_twist_range(
  L, theta, pipe_diameter,
  eps_loose = 0, eps_wide = 0,
  tau_min = 0.0, tau_max = 180.0, tau_step = 0.1
) {
  // Simple sweep: check allowed helices from tau_min to tau_max
  const taus = [];
  for (let t = tau_min; t <= tau_max; t += tau_step) taus.push(t);

  const allowed = taus.map(t =>
    is_helix_allowed(L, theta, t, pipe_diameter, eps_wide, eps_loose)
  );

  let tau_lo = null, tau_hi = null;
  let inRange = false;

  for (let i = 0; i < taus.length; i++) {
    if (allowed[i] && !inRange) {
      // Found start of valid range
      tau_lo = taus[i];
      inRange = true;
    } else if (!allowed[i] && inRange) {
      // Found end of valid range
      tau_hi = taus[i - 1];
      break;
    }
  }

  // If still valid at the end, close at tau_max
  if (inRange && tau_hi === null) tau_hi = tau_max;

  return { tau_lo, tau_hi };
}


  function list_possible_helices(L, theta, pipe_diameter, desired_length, tau_increment=1, bottom_offset=1, top_offset=0){
    let {tau_lo, tau_hi} = find_allowed_twist_range(L, theta, pipe_diameter);
    if (tau_lo==null || tau_hi==null) return [];
    tau_lo = Math.ceil(tau_lo) + bottom_offset;
    tau_hi = Math.floor(tau_hi) - top_offset;

    const outputs = [];
    for (let tau = tau_lo; tau < tau_hi; tau += tau_increment){
      const {radius, polar_angle, deltaz} = find_helix_parameters(L, theta, tau);
      const Nsegments_for_length = Math.floor(desired_length / deltaz);
      const actual_length = deltaz * Nsegments_for_length;
      const donut_hole_diameter = Math.abs(radius*2 - pipe_diameter);
      const Nturns = find_number_of_turns(polar_angle, Nsegments_for_length);
      const turns_per_segment = find_number_of_turns(polar_angle, 1);
      const outer_diameter = pipe_diameter + radius * 2;
      outputs.push([tau, Nsegments_for_length, actual_length, Nturns, outer_diameter, donut_hole_diameter, deltaz, turns_per_segment]);
    }
    return outputs;
  }

  window.HelixHelpers = {
    transformation_matrix,
    generate_points_on_helix,
    generate_straight_helix,
    make_line,
    transform_helix_to_z_axis,
    make_second_helix,
    point_line_distance,
    angle,
    find_number_of_segments_for_length,
    center_length_from_outer,
    polar_angle_between,
    find_number_of_turns,
    find_helix_parameters,
    analytical_angle_to_xy_plane,
    segment_angle_to_xy_plane,
    is_helix_loose_enough,
    is_helix_wide_enough,
    is_helix_allowed,
    find_allowed_twist_range,
    list_possible_helices,
  };
})();
