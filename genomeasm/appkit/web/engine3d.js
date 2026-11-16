/* DNA Studio, 3D engine: a small molecular renderer written for this project (no library, works
   offline). Spheres are ray-cast per pixel ("impostors"), so every atom is a perfect lit sphere
   with correct depth; bonds are instanced cylinders. A plain 2D-canvas renderer takes over when
   WebGL 2 is not available. The scene uses a right-handed coordinate system and is never
   mirrored, so a right-handed helix is drawn right-handed. */
(function (global) {
  "use strict";
  var NS = (global.DNAStudio = global.DNAStudio || {});

  // ---- quaternions (x, y, z, w) ------------------------------------------------------------------
  function qAxis(x, y, z, angle) {
    var s = Math.sin(angle / 2);
    return [x * s, y * s, z * s, Math.cos(angle / 2)];
  }
  function qMul(a, b) {
    return [
      a[3] * b[0] + a[0] * b[3] + a[1] * b[2] - a[2] * b[1],
      a[3] * b[1] - a[0] * b[2] + a[1] * b[3] + a[2] * b[0],
      a[3] * b[2] + a[0] * b[1] - a[1] * b[0] + a[2] * b[3],
      a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2],
    ];
  }
  function qNorm(q) {
    var n = Math.hypot(q[0], q[1], q[2], q[3]) || 1;
    return [q[0] / n, q[1] / n, q[2] / n, q[3] / n];
  }
  function qToRows(q) {
    var x = q[0], y = q[1], z = q[2], w = q[3];
    return [
      1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w),
      2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w),
      2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y),
    ];
  }
  NS.math = { qAxis: qAxis, qMul: qMul, qNorm: qNorm, qToRows: qToRows };

  function hexToRgb(hex) {
    var v = parseInt(hex.slice(1), 16);
    return [((v >> 16) & 255) / 255, ((v >> 8) & 255) / 255, (v & 255) / 255];
  }
  NS.hexToRgb = hexToRgb;

  // ---- shaders -----------------------------------------------------------------------------------
  var LIGHTING = [
    "uniform vec3 uBg; uniform vec3 uFog; uniform vec2 uGloss;",
    "vec3 shade(vec3 base, vec3 n, vec3 pos, float emphasis) {",
    "  vec3 v = normalize(-pos);",
    "  vec3 l1 = normalize(vec3(-0.45, 0.65, 0.62));",
    "  vec3 l2 = normalize(vec3(0.70, -0.20, 0.55));",
    "  float d1 = max(dot(n, l1), 0.0);",
    "  float d2 = max(dot(n, l2), 0.0);",
    "  float spec = pow(max(dot(n, normalize(l1 + v)), 0.0), uGloss.x) * uGloss.y;",
    "  vec3 ambient = mix(vec3(0.20, 0.21, 0.26), vec3(0.42, 0.44, 0.48), n.y * 0.5 + 0.5);",
    "  vec3 col = base * (ambient + 0.74 * d1 + 0.20 * d2) + vec3(spec);",
    "  col *= mix(0.50, 1.0, smoothstep(0.0, 0.55, dot(n, v)));",
    "  float fog = clamp((-pos.z - uFog.x) / (uFog.y - uFog.x), 0.0, 1.0) * uFog.z;",
    "  col = mix(col, uBg, fog);",
    "  return mix(uBg, col, mix(0.16, 1.0, emphasis));",
    "}",
  ].join("\n");

  var SPHERE_VS = [
    "#version 300 es",
    "layout(location=0) in vec2 aCorner;",
    "layout(location=1) in vec4 aCenterRadius;",
    "layout(location=2) in vec4 aColor;",
    "uniform mat4 uView; uniform mat4 uProj;",
    "out vec3 vCenter; out float vRadius; out vec4 vColor; out vec3 vPos;",
    "void main() {",
    "  vec4 c = uView * vec4(aCenterRadius.xyz, 1.0);",
    "  vCenter = c.xyz; vRadius = aCenterRadius.w; vColor = aColor;",
    "  vPos = c.xyz + vec3(aCorner * aCenterRadius.w * 1.4, 0.0);",
    "  gl_Position = uProj * vec4(vPos, 1.0);",
    "}",
  ].join("\n");

  var SPHERE_FS = [
    "#version 300 es",
    "precision highp float;",
    "in vec3 vCenter; in float vRadius; in vec4 vColor; in vec3 vPos;",
    "uniform mat4 uProj;",
    LIGHTING,
    "out vec4 fragColor;",
    "void main() {",
    "  vec3 ray = normalize(vPos);",
    "  float b = dot(ray, vCenter);",
    "  float disc = b * b - dot(vCenter, vCenter) + vRadius * vRadius;",
    "  float edge = disc / (vRadius * vRadius);",
    "  float width = fwidth(edge);",
    "  if (disc < 0.0) discard;",
    "  vec3 hit = ray * (b - sqrt(disc));",
    "  vec3 n = normalize(hit - vCenter);",
    "  vec4 clip = uProj * vec4(hit, 1.0);",
    "  gl_FragDepth = clip.z / clip.w * 0.5 + 0.5;",
    "  fragColor = vec4(shade(vColor.rgb, n, hit, vColor.a), smoothstep(0.0, max(width, 1e-5), edge));",
    "}",
  ].join("\n");

  var CYL_VS = [
    "#version 300 es",
    "layout(location=0) in vec3 aPos;",
    "layout(location=1) in vec3 aStart;",
    "layout(location=2) in vec3 aEnd;",
    "layout(location=3) in float aRadius;",
    "layout(location=4) in vec4 aColor;",
    "uniform mat4 uView; uniform mat4 uProj;",
    "out vec3 vNormal; out vec3 vPos; out vec4 vColor;",
    "void main() {",
    "  vec3 axis = aEnd - aStart;",
    "  vec3 dir = normalize(axis);",
    "  vec3 helper = abs(dir.y) < 0.99 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);",
    "  vec3 u = normalize(cross(dir, helper));",
    "  vec3 w = cross(dir, u);",
    "  vec3 radial = u * aPos.x + w * aPos.z;",
    "  vec4 p = uView * vec4(aStart + axis * aPos.y + radial * aRadius, 1.0);",
    "  vNormal = mat3(uView) * radial; vPos = p.xyz; vColor = aColor;",
    "  gl_Position = uProj * p;",
    "}",
  ].join("\n");

  var CYL_FS = [
    "#version 300 es",
    "precision highp float;",
    "in vec3 vNormal; in vec3 vPos; in vec4 vColor;",
    LIGHTING,
    "out vec4 fragColor;",
    "void main() { fragColor = vec4(shade(vColor.rgb, normalize(vNormal), vPos, vColor.a), 1.0); }",
  ].join("\n");

  function compile(gl, type, source) {
    var shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
      throw new Error("shader: " + gl.getShaderInfoLog(shader));
    }
    return shader;
  }
  function program(gl, vs, fs) {
    var p = gl.createProgram();
    gl.attachShader(p, compile(gl, gl.VERTEX_SHADER, vs));
    gl.attachShader(p, compile(gl, gl.FRAGMENT_SHADER, fs));
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS)) {
      throw new Error("program: " + gl.getProgramInfoLog(p));
    }
    var uniforms = {};
    var count = gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS);
    for (var i = 0; i < count; i++) {
      var info = gl.getActiveUniform(p, i);
      uniforms[info.name] = gl.getUniformLocation(p, info.name);
    }
    return { id: p, u: uniforms };
  }

  var SEGMENTS = 16;
  var SPHERE_STRIDE = 8; // centre xyz, radius, rgb, emphasis
  var CYL_STRIDE = 11; // start xyz, end xyz, radius, rgb, emphasis

  /* Viewer3D(container, options)
       options.background: "#rrggbb" stage colour (used for fog and dimming)
       options.onPick(index | -1): a click on sphere number index (or on nothing)
       options.onFrame(viewer): called after each frame (for labels drawn on the overlay)
       options.forceCanvas: use the 2D fallback renderer (for testing) */
  function Viewer3D(container, options) {
    this.options = options || {};
    this.container = container;
    this.spheres = new Float32Array(0);
    this.cyls = new Float32Array(0);
    this.nSpheres = 0;
    this.nCyls = 0;
    this.rotation = [0, 0, 0, 1];
    this.center = [0, 0, 0];
    this.sceneRadius = 20;
    this.distance = 100;
    this.fov = (26 * Math.PI) / 180;
    this.spin = false;
    this.spinSpeed = 0.35; // radians per second about the vertical axis of the scene
    this.bg = hexToRgb(this.options.background || "#0b1626");
    this.active = false;
    this.visible = true;
    this.disposed = false;
    this.frameRequested = false;
    this.lastTime = 0;
    this.frames = 0;

    container.style.position = container.style.position || "relative";
    this.canvas = document.createElement("canvas");
    this.overlay = document.createElement("canvas");
    [this.canvas, this.overlay].forEach(function (c) {
      c.style.cssText = "position:absolute;inset:0;width:100%;height:100%;display:block;";
    });
    this.overlay.style.pointerEvents = "none";
    this.canvas.style.touchAction = "none";
    this.canvas.style.cursor = "grab";
    this.canvas.setAttribute("role", "img");
    this.canvas.setAttribute("aria-label", this.options.label || "Interactive 3D model of DNA");
    container.appendChild(this.canvas);
    container.appendChild(this.overlay);
    this.overlayCtx = this.overlay.getContext("2d");

    this.gl = null;
    if (!this.options.forceCanvas) {
      try {
        this.gl = this.canvas.getContext("webgl2", {
          antialias: true, alpha: false, depth: true, powerPreference: "high-performance",
        });
        if (this.gl) this._initGL();
      } catch (error) {
        if (global.console) console.warn("DNA Studio: WebGL 2 unavailable, using the 2D renderer.", error);
        this.gl = null;
      }
    }
    if (!this.gl) {
      // a canvas that was asked for a WebGL context cannot give a 2D context: replace it
      var fresh = this.canvas.cloneNode(false);
      container.replaceChild(fresh, this.canvas);
      this.canvas = fresh;
      this.ctx2d = this.canvas.getContext("2d");
    }
    this.renderer = this.gl ? "webgl2" : "canvas2d";
    this._bindEvents();
    this.resize();
  }

  Viewer3D.prototype._initGL = function () {
    var gl = this.gl;
    this.sphereProgram = program(gl, SPHERE_VS, SPHERE_FS);
    this.cylProgram = program(gl, CYL_VS, CYL_FS);

    this.sphereVao = gl.createVertexArray();
    gl.bindVertexArray(this.sphereVao);
    var quad = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, quad);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]), gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
    this.sphereBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.sphereBuffer);
    var s = SPHERE_STRIDE * 4;
    gl.enableVertexAttribArray(1);
    gl.vertexAttribPointer(1, 4, gl.FLOAT, false, s, 0);
    gl.vertexAttribDivisor(1, 1);
    gl.enableVertexAttribArray(2);
    gl.vertexAttribPointer(2, 4, gl.FLOAT, false, s, 16);
    gl.vertexAttribDivisor(2, 1);

    this.cylVao = gl.createVertexArray();
    gl.bindVertexArray(this.cylVao);
    var mesh = [];
    for (var i = 0; i <= SEGMENTS; i++) {
      var a = (i / SEGMENTS) * Math.PI * 2;
      mesh.push(Math.cos(a), 0, Math.sin(a), Math.cos(a), 1, Math.sin(a));
    }
    this.cylVertices = mesh.length / 3;
    var meshBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, meshBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(mesh), gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
    this.cylBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.cylBuffer);
    var c = CYL_STRIDE * 4;
    [[1, 3, 0], [2, 3, 12], [3, 1, 24], [4, 4, 28]].forEach(function (spec) {
      gl.enableVertexAttribArray(spec[0]);
      gl.vertexAttribPointer(spec[0], spec[1], gl.FLOAT, false, c, spec[2]);
      gl.vertexAttribDivisor(spec[0], 1);
    });
    gl.bindVertexArray(null);
    gl.enable(gl.DEPTH_TEST);
    gl.enable(gl.SAMPLE_ALPHA_TO_COVERAGE);
    var self = this;
    this.canvas.addEventListener("webglcontextlost", function (event) {
      event.preventDefault();
      self.contextLost = true;
    });
    this.canvas.addEventListener("webglcontextrestored", function () {
      self.contextLost = false;
      self._initGL();
      self._upload();
      self.invalidate();
    });
  };

  /* setData(spheres, cylinders): Float32Arrays with 8 and 11 numbers per item (see the strides). */
  Viewer3D.prototype.setData = function (spheres, cyls) {
    this.spheres = spheres;
    this.cyls = cyls;
    this.nSpheres = spheres.length / SPHERE_STRIDE;
    this.nCyls = cyls.length / CYL_STRIDE;
    this._upload();
    this.invalidate();
  };

  Viewer3D.prototype._upload = function () {
    var gl = this.gl;
    if (!gl) return;
    gl.bindBuffer(gl.ARRAY_BUFFER, this.sphereBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, this.spheres, gl.DYNAMIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.cylBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, this.cyls, gl.DYNAMIC_DRAW);
  };

  /* frame(center, radius): look at a sphere of the given radius; keeps the current rotation. */
  Viewer3D.prototype.frame = function (center, radius, zoom, extent) {
    this.center = center.slice();
    this.sceneRadius = radius;
    this.extent = extent || null; // [half width, half height, half depth] to keep in view
    this.homeZoom = zoom || 1;
    this._fit();
    this.invalidate();
  };

  Viewer3D.prototype._fit = function () {
    var aspect = this.width / Math.max(this.height, 1);
    var half = Math.tan(this.fov / 2);
    var limit = Math.min(half, half * aspect);
    var distance = this.sceneRadius / Math.sin(Math.atan(limit));
    if (this.extent) {
      distance = Math.max(this.extent[1] / half, this.extent[0] / (half * aspect)) + this.extent[2];
    }
    this.homeDistance = distance / (this.homeZoom || 1);
    if (!this.userZoomed) this.distance = this.homeDistance;
  };

  Viewer3D.prototype.setRotation = function (q) {
    this.rotation = qNorm(q);
    this.invalidate();
  };

  Viewer3D.prototype.reset = function (q) {
    this.userZoomed = false;
    this.distance = this.homeDistance;
    this.rotation = q ? qNorm(q) : [0, 0, 0, 1];
    this.invalidate();
  };

  Viewer3D.prototype.zoomBy = function (factor) {
    this.userZoomed = true;
    var d = this.distance * factor;
    this.distance = Math.min(this.homeDistance * 3, Math.max(this.homeDistance * 0.18, d));
    this.invalidate();
  };

  Viewer3D.prototype.setSpin = function (on) {
    this.spin = !!on;
    this.lastTime = 0;
    this.invalidate();
  };

  Viewer3D.prototype.setBackground = function (hex) {
    this.bg = hexToRgb(hex);
    this.invalidate();
  };

  Viewer3D.prototype.resize = function () {
    var rect = this.container.getBoundingClientRect();
    var dpr = Math.min(global.devicePixelRatio || 1, 2);
    this.width = Math.max(1, Math.round(rect.width));
    this.height = Math.max(1, Math.round(rect.height));
    this.dpr = dpr;
    var w = Math.round(this.width * dpr), h = Math.round(this.height * dpr);
    [this.canvas, this.overlay].forEach(function (c) {
      if (c.width !== w || c.height !== h) { c.width = w; c.height = h; }
    });
    this._fit();
    this.invalidate();
  };

  Viewer3D.prototype._matrices = function () {
    var r = qToRows(this.rotation);
    var c = this.center;
    var tx = -(r[0] * c[0] + r[1] * c[1] + r[2] * c[2]);
    var ty = -(r[3] * c[0] + r[4] * c[1] + r[5] * c[2]);
    var tz = -(r[6] * c[0] + r[7] * c[1] + r[8] * c[2]) - this.distance;
    var aspect = this.width / Math.max(this.height, 1);
    var f = 1 / Math.tan(this.fov / 2);
    var near = Math.max(0.5, this.distance - this.sceneRadius * 1.6);
    var far = this.distance + this.sceneRadius * 1.6;
    this.rows = r;
    this.translation = [tx, ty, tz];
    this.f = f;
    this.aspect = aspect;
    this.view = new Float32Array([
      r[0], r[3], r[6], 0, r[1], r[4], r[7], 0, r[2], r[5], r[8], 0, tx, ty, tz, 1,
    ]);
    this.proj = new Float32Array([
      f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) / (near - far), -1,
      0, 0, (2 * far * near) / (near - far), 0,
    ]);
    this.fog = [this.distance - this.sceneRadius, this.distance + this.sceneRadius * 1.15, 0.62];
  };

  /* toView([x, y, z]) -> position in camera space (camera at the origin, looking along -z). */
  Viewer3D.prototype.toView = function (p) {
    var r = this.rows, t = this.translation;
    return [
      r[0] * p[0] + r[1] * p[1] + r[2] * p[2] + t[0],
      r[3] * p[0] + r[4] * p[1] + r[5] * p[2] + t[1],
      r[6] * p[0] + r[7] * p[1] + r[8] * p[2] + t[2],
    ];
  };

  /* project([x, y, z]) -> {x, y, depth, scale} in CSS pixels of the canvas, or null if behind. */
  Viewer3D.prototype.project = function (p) {
    if (!this.rows) this._matrices();
    var v = this.toView(p);
    if (v[2] > -0.01) return null;
    var k = this.f / -v[2];
    return {
      x: (1 + (k / this.aspect) * v[0]) * this.width / 2,
      y: (1 - k * v[1]) * this.height / 2,
      depth: -v[2],
      scale: (k * this.height) / 2, // CSS pixels per angstrom at that depth
    };
  };

  /* pick(x, y): index of the nearest sphere under a canvas position (CSS pixels), or -1. */
  Viewer3D.prototype.pick = function (x, y) {
    if (!this.rows) this._matrices();
    var dx = ((2 * x) / this.width - 1) * this.aspect / this.f;
    var dy = (1 - (2 * y) / this.height) / this.f;
    var n = Math.hypot(dx, dy, 1);
    var ray = [dx / n, dy / n, -1 / n];
    var best = -1, bestT = Infinity, s = this.spheres;
    for (var i = 0; i < this.nSpheres; i++) {
      var o = i * SPHERE_STRIDE;
      if (this.pickMask && !this.pickMask[i]) continue; // decoration, not pickable
      var c = this.toView([s[o], s[o + 1], s[o + 2]]);
      var b = ray[0] * c[0] + ray[1] * c[1] + ray[2] * c[2];
      var disc = b * b - (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]) + s[o + 3] * s[o + 3];
      if (disc >= 0) {
        var t = b - Math.sqrt(disc);
        if (t > 0 && t < bestT) { bestT = t; best = i; }
      }
    }
    return best;
  };

  Viewer3D.prototype.invalidate = function () {
    if (this.frameRequested || this.disposed) return;
    this.frameRequested = true;
    var self = this;
    global.requestAnimationFrame(function (time) {
      self.frameRequested = false;
      self._frame(time);
    });
  };

  Viewer3D.prototype._frame = function (time) {
    if (this.disposed) return;
    if (this.spin && this.visible && !this.dragging) {
      var dt = this.lastTime ? Math.min(0.05, (time - this.lastTime) / 1000) : 0;
      this.rotation = qNorm(qMul(this.rotation, qAxis(0, 1, 0, dt * this.spinSpeed)));
    }
    this.lastTime = time;
    this.render();
    if (this.spin && this.visible) this.invalidate();
  };

  Viewer3D.prototype.render = function () {
    this._matrices();
    this.frames++;
    if (this.gl) this._renderGL(); else this._renderCanvas();
    var ctx = this.overlayCtx;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.clearRect(0, 0, this.width, this.height);
    if (this.options.onFrame) this.options.onFrame(this, ctx);
  };

  Viewer3D.prototype._renderGL = function () {
    var gl = this.gl;
    if (this.contextLost) return;
    gl.viewport(0, 0, this.canvas.width, this.canvas.height);
    gl.clearColor(this.bg[0], this.bg[1], this.bg[2], 1);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    var self = this;
    function common(p) {
      gl.useProgram(p.id);
      gl.uniformMatrix4fv(p.u.uView, false, self.view);
      gl.uniformMatrix4fv(p.u.uProj, false, self.proj);
      gl.uniform3fv(p.u.uBg, self.bg);
      gl.uniform3fv(p.u.uFog, self.fog);
      gl.uniform2f(p.u.uGloss, 42.0, 0.22);
    }
    if (this.nCyls) {
      common(this.cylProgram);
      gl.bindVertexArray(this.cylVao);
      gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, this.cylVertices, this.nCyls);
    }
    if (this.nSpheres) {
      common(this.sphereProgram);
      gl.bindVertexArray(this.sphereVao);
      gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, this.nSpheres);
    }
    gl.bindVertexArray(null);
  };

  // 2D fallback: painter's algorithm (far to near), spheres as shaded discs, bonds as thick lines.
  Viewer3D.prototype._renderCanvas = function () {
    var ctx = this.ctx2d, items = [], i, o, a, b;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    var bg = this.bg;
    function css(r, g, bl, k) {
      function ch(v, base) { return Math.round(255 * (base + (v - base) * k)); }
      return "rgb(" + ch(r, bg[0]) + "," + ch(g, bg[1]) + "," + ch(bl, bg[2]) + ")";
    }
    ctx.fillStyle = css(bg[0], bg[1], bg[2], 1);
    ctx.fillRect(0, 0, this.width, this.height);
    var s = this.spheres, c = this.cyls;
    for (i = 0; i < this.nSpheres; i++) {
      o = i * SPHERE_STRIDE;
      a = this.project([s[o], s[o + 1], s[o + 2]]);
      if (a) items.push({ z: a.depth, kind: 0, a: a, r: s[o + 3] * a.scale, o: o });
    }
    for (i = 0; i < this.nCyls; i++) {
      o = i * CYL_STRIDE;
      a = this.project([c[o], c[o + 1], c[o + 2]]);
      b = this.project([c[o + 3], c[o + 4], c[o + 5]]);
      if (a && b) items.push({ z: (a.depth + b.depth) / 2, kind: 1, a: a, b: b, r: c[o + 6] * a.scale, o: o });
    }
    items.sort(function (p, q) { return q.z - p.z; });
    var near = this.fog[0], far = this.fog[1];
    for (i = 0; i < items.length; i++) {
      var it = items[i];
      var fade = 1 - Math.min(1, Math.max(0, (it.z - near) / (far - near))) * this.fog[2];
      if (it.kind === 0) {
        var k = fade * (0.16 + 0.84 * Math.max(0, s[it.o + 7]));
        var g = ctx.createRadialGradient(it.a.x - it.r * 0.35, it.a.y - it.r * 0.4, it.r * 0.1, it.a.x, it.a.y, it.r);
        g.addColorStop(0, css(Math.min(1, s[it.o + 4] + 0.35), Math.min(1, s[it.o + 5] + 0.35), Math.min(1, s[it.o + 6] + 0.35), k));
        g.addColorStop(0.6, css(s[it.o + 4], s[it.o + 5], s[it.o + 6], k));
        g.addColorStop(1, css(s[it.o + 4] * 0.45, s[it.o + 5] * 0.45, s[it.o + 6] * 0.45, k));
        ctx.fillStyle = g;
        ctx.beginPath();
        ctx.arc(it.a.x, it.a.y, Math.max(it.r, 0.5), 0, Math.PI * 2);
        ctx.fill();
      } else {
        ctx.strokeStyle = css(c[it.o + 7] * 0.85, c[it.o + 8] * 0.85, c[it.o + 9] * 0.85, fade * (0.16 + 0.84 * c[it.o + 10]));
        ctx.lineWidth = Math.max(1, it.r * 2);
        ctx.lineCap = "round";
        ctx.beginPath();
        ctx.moveTo(it.a.x, it.a.y);
        ctx.lineTo(it.b.x, it.b.y);
        ctx.stroke();
      }
    }
  };

  Viewer3D.prototype._bindEvents = function () {
    var self = this, canvas = this.canvas, last = null, moved = 0, pointers = {};
    function position(event) {
      var rect = canvas.getBoundingClientRect();
      return [event.clientX - rect.left, event.clientY - rect.top];
    }
    this._onDown = function (event) {
      pointers[event.pointerId] = position(event);
      last = position(event);
      moved = 0;
      self.dragging = true;
      self.active = true;
      canvas.style.cursor = "grabbing";
      if (canvas.setPointerCapture) { try { canvas.setPointerCapture(event.pointerId); } catch (e) { /* ignore */ } }
    };
    this._onMove = function (event) {
      if (!self.dragging || !last) return;
      var p = position(event);
      var ids = Object.keys(pointers);
      if (ids.length >= 2) { // pinch to zoom
        var before = Math.hypot(pointers[ids[0]][0] - pointers[ids[1]][0], pointers[ids[0]][1] - pointers[ids[1]][1]);
        pointers[event.pointerId] = p;
        var after = Math.hypot(pointers[ids[0]][0] - pointers[ids[1]][0], pointers[ids[0]][1] - pointers[ids[1]][1]);
        if (before > 0 && after > 0) self.zoomBy(before / after);
        return;
      }
      var dx = p[0] - last[0], dy = p[1] - last[1];
      moved += Math.abs(dx) + Math.abs(dy);
      last = p;
      var k = 3.2 / Math.max(self.height, 1);
      self.rotation = qNorm(qMul(qMul(qAxis(0, 1, 0, dx * k), qAxis(1, 0, 0, dy * k)), self.rotation));
      self.invalidate();
    };
    this._onUp = function (event) {
      delete pointers[event.pointerId];
      if (!self.dragging) return;
      self.dragging = false;
      canvas.style.cursor = "grab";
      if (moved < 4 && self.options.onPick) {
        var p = position(event);
        self.options.onPick(self.pick(p[0], p[1]));
      }
      self.lastTime = 0;
      self.invalidate();
    };
    this._onWheel = function (event) {
      if (!self.active && !self.options.wheelAlways) return; // do not steal the page scroll
      event.preventDefault();
      self.zoomBy(Math.exp(event.deltaY * 0.0012));
    };
    this._onLeave = function () { self.active = false; };
    this._onDouble = function () { if (self.options.onReset) self.options.onReset(); else self.reset(); };
    canvas.addEventListener("pointerdown", this._onDown);
    canvas.addEventListener("pointermove", this._onMove);
    canvas.addEventListener("pointerup", this._onUp);
    canvas.addEventListener("pointercancel", this._onUp);
    canvas.addEventListener("wheel", this._onWheel, { passive: false });
    canvas.addEventListener("pointerleave", this._onLeave);
    canvas.addEventListener("dblclick", this._onDouble);
    if (global.ResizeObserver) {
      this.resizeObserver = new ResizeObserver(function () { self.resize(); });
      this.resizeObserver.observe(this.container);
    }
    if (global.IntersectionObserver) {
      this.intersectionObserver = new IntersectionObserver(function (entries) {
        self.visible = entries[0].isIntersecting;
        if (self.visible) { self.lastTime = 0; self.invalidate(); }
      });
      this.intersectionObserver.observe(this.container);
    }
    this._onVisibility = function () {
      if (!document.hidden) { self.lastTime = 0; self.invalidate(); }
    };
    document.addEventListener("visibilitychange", this._onVisibility);
  };

  /* dispose(): stop drawing and free the WebGL context (a page may hold only a few of them). */
  Viewer3D.prototype.dispose = function () {
    this.disposed = true;
    if (this.resizeObserver) this.resizeObserver.disconnect();
    if (this.intersectionObserver) this.intersectionObserver.disconnect();
    document.removeEventListener("visibilitychange", this._onVisibility);
    if (this.gl) {
      var ext = this.gl.getExtension("WEBGL_lose_context");
      if (ext) ext.loseContext();
    }
    if (this.canvas.parentNode) this.canvas.parentNode.removeChild(this.canvas);
    if (this.overlay.parentNode) this.overlay.parentNode.removeChild(this.overlay);
  };

  NS.Viewer3D = Viewer3D;
  NS.SPHERE_STRIDE = SPHERE_STRIDE;
  NS.CYL_STRIDE = CYL_STRIDE;
})(typeof window !== "undefined" ? window : globalThis);
