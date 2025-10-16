/***** plothelix.js *****/
// Helper to compute a vector in the XY-plane perpendicular to a given vector v.
function xy_perp(v) {
  const [vx, vy, vz] = v;
  // If v has no horizontal component, return a default unit vector in XY-plane
  if (Math.abs(vx) < 1e-9 && Math.abs(vy) < 1e-9) {
    return [1, 0, 0];  // arbitrary horizontal direction
  }
  // A perpendicular in XY-plane can be (-vy, vx, 0)
  return unit([-vy, vx, 0]);
}

// Given an array of tangent vectors, compute local perpendicular basis (u1, u2) for each.
function visLocalXYPlane(vectors) {
  return vectors.map(vec => {
    const v_unit = unit(vec.slice());  // normalize copy of vec
    const u1 = xy_perp(v_unit);
    const u2 = unit(cross(v_unit, u1));
    return [u1, u2];
  });
}

// Compute all points of a circle at a given helix point, using basis vectors for orientation.
function visWholeCircle(centerPoint, basis, pipe_diameter, sides = 10) {
  const r = pipe_diameter / 2;
  const [u1, u2] = basis;
  const circlePoints = [];
  // Generate points around 0 to 2π (inclusive to close the loop)
  for (let i = 0; i < sides; i++) {
    const angle = 2 * Math.PI * (i / (sides - 1));
    const cosA = Math.cos(angle), sinA = Math.sin(angle);
    // point = center + r*(cosA*u1 + sinA*u2)
    const offset = [
      r * (cosA * u1[0] + sinA * u2[0]),
      r * (cosA * u1[1] + sinA * u2[1]),
      r * (cosA * u1[2] + sinA * u2[2])
    ];
    circlePoints.push(add(centerPoint, offset));
  }
  return circlePoints;
}

// Create a Mesh3D trace connecting two parallel rings of points (line1 and line2).
function connectLinesMesh3d(line1, line2, color = 'lightblue', opacity = 0.6, name = 'surface') {
  const n = line1.length;
  // Flatten x, y, z coordinates for both lines
  const x = [], y = [], z = [];
  line1.forEach(p => { x.push(p[0]); y.push(p[1]); z.push(p[2]); });
  line2.forEach(p => { x.push(p[0]); y.push(p[1]); z.push(p[2]); });
  // Triangles between the two rings:
  const I = [], J = [], K = [];
  for (let i = 0; i < n - 1; i++) {
    // Two triangles for each rectangular panel between line1[i]-line1[i+1] and line2[i]-line2[i+1]
    I.push(i, i+1);
    J.push(i+1, n + i + 1);
    K.push(n + i, n + i);
    // Triangle indices:
    //  (i, i+1, n+i) and (i+1, n+i+1, n+i)
    // Pushing in pairs as above corresponds to these two triangles.
  }
  return {
    type: 'mesh3d',
    x: x, y: y, z: z,
    i: I, j: J, k: K,
    opacity: opacity,
    color: color,
    name: name
  };
}

// Generate Plotly traces for a helix (and its optional second helix) and render the plot.
function plotHelix(L, theta, tau, segments, pipe_diameter, doubleHelix = true, containerId = 'plotContainer') {
  // Generate base helix points (centerline aligned to z-axis)
  const basePoints = generateStraightHelix(L, theta, tau, segments);
  // Compute tangent vectors along helix (using central differences for interior points)
  const tangents = [];
  for (let i = 0; i < basePoints.length; i++) {
    if (i === 0) {
      tangents.push(sub(basePoints[1], basePoints[0]));  // forward difference at start
    } else if (i === basePoints.length - 1) {
      tangents.push(sub(basePoints[i], basePoints[i-1]));  // backward difference at end
    } else {
      // average of forward and backward differences for smoother tangent
      const fwd = sub(basePoints[i+1], basePoints[i]);
      const bwd = sub(basePoints[i], basePoints[i-1]);
      tangents.push(mul(add(fwd, bwd), 0.5));
    }
  }
  // Local perpendicular bases at each point
  const bases = visLocalXYPlane(tangents);
  // Compute cross-section circles for each segment point
  const circles = basePoints.map((pt, idx) => visWholeCircle(pt, bases[idx], pipe_diameter));
  // Transpose circles to get "meridian" lines along the helix
  const linesAlong = [];
  const sides = circles[0].length;  // number of points per circle
  for (let j = 0; j < sides; j++) {
    const line = [];
    for (let i = 0; i < circles.length; i++) {
      line.push(circles[i][j]);
    }
    linesAlong.push(line);
  }
  // Prepare Plotly traces
  const traces = [];
  // Mesh surfaces for each adjacent pair of meridians (skip connecting last to first to avoid duplicate seam)
  for (let j = 0; j < linesAlong.length - 1; j++) {
    traces.push(connectLinesMesh3d(linesAlong[j], linesAlong[j+1], 'blue', 0.6, 'Helix Surface'));
  }
  // Outline lines along the helix (meridians)
  linesAlong.forEach(line => {
    const xs = line.map(p => p[0]);
    const ys = line.map(p => p[1]);
    const zs = line.map(p => p[2]);
    traces.push({
      type: 'scatter3d',
      mode: 'lines',
      x: xs, y: ys, z: zs,
      line: { color: 'blue', width: 2 },
      name: 'Helix Edge'
    });
  });
  // If double helix, generate mirrored helix and add its traces
  if (doubleHelix) {
    // Mirror the centerline points across the z-axis
    const secondPoints = basePoints.map(([x, y, z]) => [-x, -y, z]);
    // Compute tangents and bases for second helix
    const tangents2 = [];
    for (let i = 0; i < secondPoints.length; i++) {
      if (i === 0) tangents2.push(sub(secondPoints[1], secondPoints[0]));
      else if (i === secondPoints.length - 1) tangents2.push(sub(secondPoints[i], secondPoints[i-1]));
      else {
        const fwd = sub(secondPoints[i+1], secondPoints[i]);
        const bwd = sub(secondPoints[i], secondPoints[i-1]);
        tangents2.push(mul(add(fwd, bwd), 0.5));
      }
    }
    const bases2 = visLocalXYPlane(tangents2);
    const circles2 = secondPoints.map((pt, idx) => visWholeCircle(pt, bases2[idx], pipe_diameter));
    const linesAlong2 = [];
    const sides2 = circles2[0].length;
    for (let j = 0; j < sides2 - 1; j++) {  // note: skip last index if it's duplicate 2π point
      const line = [];
      for (let i = 0; i < circles2.length; i++) {
        line.push(circles2[i][j]);
      }
      linesAlong2.push(line);
    }
    // Mesh surfaces for second helix
    for (let j = 0; j < linesAlong2.length - 1; j++) {
      traces.push(connectLinesMesh3d(linesAlong2[j], linesAlong2[j+1], 'red', 0.6, 'Helix Surface 2'));
    }
    // Lines for second helix
    linesAlong2.forEach(line => {
      const xs = line.map(p => p[0]);
      const ys = line.map(p => p[1]);
      const zs = line.map(p => p[2]);
      traces.push({
        type: 'scatter3d',
        mode: 'lines',
        x: xs, y: ys, z: zs,
        line: { color: 'red', width: 2 },
        name: 'Helix Edge 2'
      });
    });
  }
  // Define layout for 3D plot
  const layout = {
    scene: {
      aspectmode: 'data',    // equal axis scaling
      camera: { eye: { x: 0, y: 1.8, z: 1.2 } }
    },
    dragmode: 'turntable',
    margin: { l: 0, r: 0, b: 0, t: 0 }  // full viewport usage
  };
  // Render the plot into the container
  Plotly.newPlot(containerId, traces, layout);
}
