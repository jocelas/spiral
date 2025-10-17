
// plottinghelix.js — translated from plotting_helix.py
// Depends on window.HelixHelpers and Plotly.js

(function(){
  const H = window.HelixHelpers;

  function vis_tangent_vectors(segments, L, theta, tau){
    const {radius, polar_angle, deltaz} = H.find_helix_parameters(L, theta, tau);
    const normal_vectors = [];
    for (let n = 0; n <= segments; n++){
      const x = - radius * polar_angle * Math.sin(n * polar_angle);
      const y =   radius * polar_angle * Math.cos(n * polar_angle);
      const z =   deltaz;
      const v = [x, y, z];
      normal_vectors.push(H._unit ? H._unit(v) : (function(vv){ const n=Math.hypot(...vv); return vv.map(x=>x/n);})(v));
    }
    return normal_vectors;
  }

  function xy_perp(v){
    const [vx, vy, vz] = v;
    if (vx === 0 && vy === 0){
      return [1, 0, 0];
    }
    const w = [-vy, vx, 0];
    const n = Math.hypot(...w);
    return w.map(x=>x/n);
  }

  function vis_local_xy_plane(vectors){
    const basis = [];
    for (const vec of vectors){
      const v = (()=>{const n=Math.hypot(...vec); return vec.map(x=>x/n);})();
      const u1 = xy_perp(v);
      const crossv = [
        v[1]*u1[2] - v[2]*u1[1],
        v[2]*u1[0] - v[0]*u1[2],
        v[0]*u1[1] - v[1]*u1[0],
      ];
      const n2 = Math.hypot(...crossv);
      const u2 = crossv.map(x=>x/n2);
      basis.push([u1, u2]);
    }
    return basis;
  }

  function vis_circle(angle, point, basis, pipe_diameter){
    const r = pipe_diameter / 2;
    const [u1, u2] = basis;
    return [
      point[0] + r * (Math.cos(angle)*u1[0] + Math.sin(angle)*u2[0]),
      point[1] + r * (Math.cos(angle)*u1[1] + Math.sin(angle)*u2[1]),
      point[2] + r * (Math.cos(angle)*u1[2] + Math.sin(angle)*u2[2]),
    ];
  }

  function vis_whole_circle(point, basis, pipe_diameter, sides=10){
    const pts = [];
    for (let i=0;i<sides;i++){
      const angle = 2*Math.PI*i/sides;
      pts.push(vis_circle(angle, point, basis, pipe_diameter));
    }
    return pts;
  }

  function connect_lines_mesh3d(line1, line2, color='lightblue', opacity=1.0, name='surface'){
    if (line1.length !== line2.length) throw new Error("Both lines must have same number of points");
    const n = line1.length;
    const x = [], y = [], z = [];
    for (const p of line1){x.push(p[0]); y.push(p[1]); z.push(p[2]);}
    for (const p of line2){x.push(p[0]); y.push(p[1]); z.push(p[2]);}
    const i = [], j = [], k = [];
    for (let idx=0; idx<n-1; idx++){
      i.push(idx); j.push(idx+1); k.push(n+idx);
      i.push(idx+1); j.push(n+idx+1); k.push(n+idx);
    }
    return {
      type: 'mesh3d',
      x, y, z, i, j, k,
      color, opacity, name,
    };
  }

  function plotHelix(L, theta, tau, segments, pipe_diameter, sides=10, double_helix=true, returnMode=false){
    const points_straight = H.generate_straight_helix(L, theta, tau, segments);
    let vectors = [];
    for (let i=0;i<points_straight.length-1;i++){
      const p1=points_straight[i], p2=points_straight[i+1];
      vectors.push([p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]);
    }
    vectors.push(vectors[vectors.length-1].slice());
    const basis = vis_local_xy_plane(vectors);

    const circles = [];
    for (let i=0;i<basis.length;i++){
      circles.push(vis_whole_circle(points_straight[i], basis[i], pipe_diameter, sides));
    }
    const lines = [];
    for (let j=0;j<sides;j++){
      const arr=[];
      for (let i=0;i<circles.length;i++) arr.push(circles[i][j]);
      lines.push(arr);
    }

    const data = [];
    for (let i=0;i<lines.length-1;i++){
      data.push(connect_lines_mesh3d(lines[i], lines[i+1], 'blue'));
      // This fix added by Jonas because it didn't work otherwise
      if (i == lines.length-2){
        data.push(connect_lines_mesh3d(lines[0], lines[i+1], 'blue', 'surface'));
      }
    }
    for (const circle of lines){
      const x = circle.map(p=>p[0]), y = circle.map(p=>p[1]), z = circle.map(p=>p[2]);
      data.push({type:'scatter3d', mode:'lines', x, y, z, line:{color:'blue'}});
    }

    if (double_helix){
      const second_points = H.make_second_helix(points_straight);
      const vectors2 = [];
      for (let i=0;i<second_points.length-1;i++){
        const p1=second_points[i], p2=second_points[i+1];
        vectors2.push([p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]);
      }
      vectors2.push(vectors2[vectors2.length-1].slice());
      const basis2 = vis_local_xy_plane(vectors2);

      const circles2 = [];
      for (let i=0;i<basis2.length;i++){
        circles2.push(vis_whole_circle(second_points[i], basis2[i], pipe_diameter, sides));
      }
      const lines2 = [];
      for (let j=0;j<sides;j++){
        const arr=[];
        for (let i=0;i<circles2.length;i++) arr.push(circles2[i][j]);
        lines2.push(arr);
      }

      for (let i=0;i<lines2.length-1;i++){
        data.push(connect_lines_mesh3d(lines2[i], lines2[i+1], 'red'));
        if (i == lines2.length-2){
        data.push(connect_lines_mesh3d(lines2[0], lines2[i+1], 'red'));
      }
      }
      for (const circle of lines2){
        const x = circle.map(p=>p[0]), y = circle.map(p=>p[1]), z = circle.map(p=>p[2]);
        data.push({type:'scatter3d', mode:'lines', x, y, z, line:{color:'red'}});
      }
    }

    const layout = {
      scene: {aspectmode:'data', camera:{eye:{x:0, y:1.8, z:1.2}}},
      dragmode: 'turntable',
      showlegend: false
    };

    if (returnMode) return {data, layout};
    Plotly.newPlot('plot', data, layout);
  }

  window.PlottingHelix = {
    vis_tangent_vectors,
    xy_perp,
    vis_local_xy_plane,
    vis_circle,
    vis_whole_circle,
    connect_lines_mesh3d,
    plotHelix
  };
})();
