/* DNA Studio, molecule view: turns a scene made by genomeasm.molbio.model into spheres and
   cylinders, and adds the controls (representation, colours, hydrogen bonds, labels, selection).
   All geometry and all numbers in the info card come from the scene; nothing is invented here. */
(function (global) {
  "use strict";
  var NS = (global.DNAStudio = global.DNAStudio || {});

  var BASE_COLOUR = { A: "#00b383", T: "#e8662a", G: "#f2ae1c", C: "#2f8fdc", U: "#cc79a7" };
  var ELEMENT_COLOUR = { C: "#b9c2cc", N: "#4a86f0", O: "#ee5548", P: "#f2a02e", S: "#9a8fd0", B: "#b9c2cc" };
  var COMPONENT_COLOUR = ["#f2a02e", "#9a8fd0", "#3fa7e8"]; // phosphate, sugar, base
  var COMPONENT_NAME = ["phosphate group", "deoxyribose (sugar)", "nitrogenous base"];
  var STRAND_COLOUR = ["#56b4e9", "#f0a73a"];
  var NEUTRAL = "#8f9bab";
  var HBOND_COLOUR = "#ffe6a3";
  var VDW = { C: 1.7, N: 1.55, O: 1.52, P: 1.8 };
  var BASE_NAME = { A: "adenine", T: "thymine", G: "guanine", C: "cytosine" };
  var BASE_FAMILY = { A: "purine (two rings)", G: "purine (two rings)", C: "pyrimidine (one ring)", T: "pyrimidine (one ring)" };
  NS.BASE_COLOUR = BASE_COLOUR;

  function el(tag, cls, text) {
    var node = document.createElement(tag);
    if (cls) node.className = cls;
    if (text !== undefined) node.textContent = text;
    return node;
  }
  NS.el = el;

  /* MoleculeView(root, scene, options)
       options.style   "ballstick" | "spacefill" | "sticks"   (atomic scenes)
       options.colour  "element" | "base" | "component" | "strand"
       options.hbonds, options.labels, options.spin: booleans
       options.atomLabels: write the atom names next to the atoms (small scenes)
       options.toolbar: false hides the controls (used inside slides)
       options.rotation: starting orientation as a quaternion
       options.zoom: magnification of the default view */
  function MoleculeView(root, scene, options) {
    var o = options || {};
    this.root = root;
    this.scene = scene;
    this.atomic = scene.kind === "atomic";
    this.state = {
      style: o.style || "ballstick",
      colour: o.colour || (this.atomic ? "element" : "base"),
      hbonds: o.hbonds !== false,
      labels: o.labels !== false,
      spin: !!o.spin,
      selected: -1,
    };
    this.options = o;
    this.n = scene.atoms.el.length;
    this.xyz = scene.atoms.xyz;
    root.classList.add("ds-mol");
    this.stage = el("div", "ds-stage");
    this.stage.dataset.theme = o.theme || "dark";
    root.appendChild(this.stage);
    var self = this;
    this.viewer = new NS.Viewer3D(this.stage, {
      background: o.background || "#0b1626",
      label: o.label,
      forceCanvas: o.forceCanvas,
      onPick: function (index) { self._picked(index); },
      onFrame: function (viewer, ctx) { self._overlay(viewer, ctx); },
      onReset: function () { self.resetView(); },
    });
    this.home = o.rotation || [0, 0, 0, 1];
    this.viewer.setRotation(this.home);
    this.viewer.frame([0, 0, 0], scene.radius, o.zoom || 1, scene.extent);
    this.badge = el("div", "ds-badge");
    this.stage.appendChild(this.badge);
    this.card = el("div", "ds-card");
    this.card.hidden = true;
    this.stage.appendChild(this.card);
    this.hint = el("div", "ds-hint", o.hint || "Drag to rotate · click the model, then scroll to zoom · click a nucleotide for details");
    this.stage.appendChild(this.hint);
    if (o.toolbar !== false) this._toolbar();
    this.legend = el("div", "ds-legend");
    root.appendChild(this.legend);
    this.rebuild();
    this.viewer.setSpin(this.state.spin);
    if (this.viewer.renderer !== "webgl2") {
      this.badge.textContent = "Simplified drawing: this browser has no WebGL 2";
    }
  }

  MoleculeView.prototype._toolbar = function () {
    var self = this, bar = el("div", "ds-toolbar");
    function group(label, key, choices) {
      var box = el("div", "ds-group");
      box.appendChild(el("span", "ds-group-label", label));
      var buttons = el("div", "ds-seg");
      choices.forEach(function (choice) {
        var b = el("button", "", choice[1]);
        b.type = "button";
        b.dataset.value = choice[0];
        b.setAttribute("aria-pressed", String(self.state[key] === choice[0]));
        b.addEventListener("click", function () {
          self.state[key] = choice[0];
          Array.prototype.forEach.call(buttons.children, function (c) {
            c.setAttribute("aria-pressed", String(c === b));
          });
          self.rebuild();
        });
        buttons.appendChild(b);
      });
      box.appendChild(buttons);
      bar.appendChild(box);
    }
    function toggle(label, key, apply) {
      var b = el("button", "ds-toggle", label);
      b.type = "button";
      b.setAttribute("aria-pressed", String(!!self.state[key]));
      b.addEventListener("click", function () {
        self.state[key] = !self.state[key];
        b.setAttribute("aria-pressed", String(self.state[key]));
        apply();
      });
      return b;
    }
    if (this.atomic) {
      group("Model", "style", [["ballstick", "Ball and stick"], ["spacefill", "Space-filling"], ["sticks", "Sticks"]]);
      group("Colour by", "colour", [["element", "Element"], ["base", "Base"], ["component", "Part"], ["strand", "Strand"]]);
    } else {
      group("Colour by", "colour", [["base", "Base"], ["component", "Part"], ["strand", "Strand"]]);
    }
    var toggles = el("div", "ds-group");
    toggles.appendChild(toggle("Hydrogen bonds", "hbonds", function () { self.rebuild(); }));
    toggles.appendChild(toggle("Labels", "labels", function () { self.viewer.invalidate(); }));
    toggles.appendChild(toggle("Rotate", "spin", function () { self.viewer.setSpin(self.state.spin); }));
    var zoomIn = el("button", "ds-icon", "+"), zoomOut = el("button", "ds-icon", "−"), reset = el("button", "ds-toggle", "Reset view");
    zoomIn.type = zoomOut.type = reset.type = "button";
    zoomIn.title = "Zoom in"; zoomOut.title = "Zoom out";
    zoomIn.setAttribute("aria-label", "Zoom in"); zoomOut.setAttribute("aria-label", "Zoom out");
    zoomIn.addEventListener("click", function () { self.viewer.zoomBy(0.8); });
    zoomOut.addEventListener("click", function () { self.viewer.zoomBy(1.25); });
    reset.addEventListener("click", function () { self.resetView(); });
    toggles.appendChild(zoomOut); toggles.appendChild(zoomIn); toggles.appendChild(reset);
    bar.appendChild(toggles);
    this.root.insertBefore(bar, this.stage);
  };

  MoleculeView.prototype.resetView = function () {
    this.viewer.reset(this.home);
    this.select(-1);
  };

  MoleculeView.prototype._atomColour = function (i) {
    var a = this.scene.atoms, residue = this.scene.residues[a.res[i]], mode = this.state.colour;
    if (mode === "element") {
      if (a.el[i] === "B") return BASE_COLOUR[residue.base];
      return ELEMENT_COLOUR[a.el[i]] || NEUTRAL;
    }
    if (mode === "component") return COMPONENT_COLOUR[a.comp[i]];
    if (mode === "strand") return STRAND_COLOUR[residue.strand];
    return a.comp[i] === 2 ? BASE_COLOUR[residue.base] : NEUTRAL; // by base
  };

  MoleculeView.prototype._emphasis = function (residue) {
    var sel = this.state.selected;
    if (sel < 0) return 1;
    if (residue === sel) return 1;
    return this.scene.residues[sel].partner === residue ? 0.72 : 0.0;
  };

  /* rebuild(): scene + current settings -> the two instance arrays of the engine. */
  MoleculeView.prototype.rebuild = function () {
    var scene = this.scene, a = scene.atoms, n = this.n, xyz = this.xyz, st = this.state;
    var spheres = [], cyls = [], mask = [], sphereAtom = [], i, c, e;
    var radius, stick;
    for (i = 0; i < n; i++) {
      if (this.atomic) {
        var vdw = VDW[a.el[i]] || 1.6;
        radius = st.style === "spacefill" ? vdw : st.style === "sticks" ? 0.2 : vdw * 0.27;
      } else {
        radius = a.el[i] === "P" ? 1.55 : a.el[i] === "S" ? 1.25 : 0.95;
      }
      c = NS.hexToRgb(this._atomColour(i));
      e = this._emphasis(a.res[i]);
      spheres.push(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2], radius, c[0], c[1], c[2], e);
      mask.push(1);
      sphereAtom.push(i);
    }
    if (!(this.atomic && st.style === "spacefill")) {
      for (var k = 0; k < scene.bonds.length; k += 2) {
        var p = scene.bonds[k], q = scene.bonds[k + 1];
        if (this.atomic) stick = st.style === "sticks" ? 0.2 : 0.13;
        else stick = a.comp[p] === 2 || a.comp[q] === 2 ? 0.62 : 0.5;
        var cp = NS.hexToRgb(this._atomColour(p)), cq = NS.hexToRgb(this._atomColour(q));
        if (!this.atomic && (a.comp[p] === 2 || a.comp[q] === 2)) {
          cp = cq = NS.hexToRgb(this._atomColour(a.comp[p] === 2 ? p : q)); // the base rod has one colour
        }
        var mx = (xyz[3 * p] + xyz[3 * q]) / 2, my = (xyz[3 * p + 1] + xyz[3 * q + 1]) / 2, mz = (xyz[3 * p + 2] + xyz[3 * q + 2]) / 2;
        cyls.push(xyz[3 * p], xyz[3 * p + 1], xyz[3 * p + 2], mx, my, mz, stick, cp[0], cp[1], cp[2], this._emphasis(a.res[p]));
        cyls.push(mx, my, mz, xyz[3 * q], xyz[3 * q + 1], xyz[3 * q + 2], stick, cq[0], cq[1], cq[2], this._emphasis(a.res[q]));
      }
    }
    if (st.hbonds) {
      var hc = NS.hexToRgb(HBOND_COLOUR), dot = this.atomic ? 0.085 : 0.2, gap = this.atomic ? 0.34 : 0.62;
      scene.hbonds.forEach(function (h) {
        var dx = h.b[0] - h.a[0], dy = h.b[1] - h.a[1], dz = h.b[2] - h.a[2];
        var length = Math.hypot(dx, dy, dz), count = Math.max(2, Math.round(length / gap));
        var em = Math.max(this._emphasis(h.ra), this._emphasis(h.rb)) > 0.9 ? 1 : (st.selected < 0 ? 1 : 0);
        for (var j = 1; j < count; j++) {
          var t = j / count;
          spheres.push(h.a[0] + dx * t, h.a[1] + dy * t, h.a[2] + dz * t, dot, hc[0], hc[1], hc[2], em);
          mask.push(0);
          sphereAtom.push(-1);
        }
      }, this);
    }
    this.sphereAtom = sphereAtom;
    this.viewer.pickMask = new Uint8Array(mask);
    this.viewer.setData(new Float32Array(spheres), new Float32Array(cyls));
    this._legend();
  };

  MoleculeView.prototype._legend = function () {
    var items = [], mode = this.state.colour, self = this;
    if (mode === "element") {
      items = this.atomic
        ? [["C", "carbon"], ["N", "nitrogen"], ["O", "oxygen"], ["P", "phosphorus"]].map(function (x) { return [ELEMENT_COLOUR[x[0]], x[1]]; })
        : [[ELEMENT_COLOUR.P, "phosphate"], [ELEMENT_COLOUR.S, "sugar"]];
    } else if (mode === "component") {
      items = COMPONENT_COLOUR.map(function (c, i) { return [c, COMPONENT_NAME[i]]; });
    } else if (mode === "strand") {
      items = [[STRAND_COLOUR[0], "strand 1"], [STRAND_COLOUR[1], "strand 2 (antiparallel)"]];
    } else {
      items = ["A", "T", "G", "C"].map(function (b) { return [BASE_COLOUR[b], b + " " + BASE_NAME[b]]; });
      items.push([NEUTRAL, "sugar-phosphate backbone"]);
    }
    if (this.state.hbonds && this.scene.hbonds.length) items.push([HBOND_COLOUR, "hydrogen bond", true]);
    this.legend.textContent = "";
    items.forEach(function (item) {
      var chip = el("span", "ds-chip");
      var swatch = el("i", item[2] ? "ds-swatch ds-dotted" : "ds-swatch");
      swatch.style.background = item[2] ? "transparent" : item[0];
      swatch.style.borderColor = item[0];
      chip.appendChild(swatch);
      chip.appendChild(document.createTextNode(item[1]));
      self.legend.appendChild(chip);
    });
    if (this.options.note) this.legend.appendChild(el("span", "ds-note", this.options.note));
  };

  MoleculeView.prototype._picked = function (sphere) {
    var atom = sphere >= 0 ? this.sphereAtom[sphere] : -1;
    this.pickedAtom = atom;
    this.select(atom >= 0 ? this.scene.atoms.res[atom] : -1);
  };

  /* select(residue index or -1): highlight a nucleotide and its partner, fill the info card. */
  MoleculeView.prototype.select = function (residue) {
    if (this.options.selectable === false) return;
    this.state.selected = residue;
    this.rebuild();
    var card = this.card;
    card.textContent = "";
    card.hidden = residue < 0;
    this.hint.hidden = residue >= 0;
    if (residue < 0) { if (this.options.onSelect) this.options.onSelect(null); return; }
    var r = this.scene.residues[residue], a = this.scene.atoms, self = this;
    var title = el("div", "ds-card-title");
    var dot = el("i", "ds-swatch");
    dot.style.background = BASE_COLOUR[r.base];
    title.appendChild(dot);
    title.appendChild(document.createTextNode("d" + r.base + r.num + " · strand " + (r.strand + 1)));
    card.appendChild(title);
    card.appendChild(el("div", "ds-card-sub", r.name));
    card.appendChild(el("div", "ds-card-row", "Base: " + BASE_NAME[r.base] + ", a " + BASE_FAMILY[r.base]));
    if (r.partner !== null && r.partner !== undefined) {
      var partner = this.scene.residues[r.partner];
      var bonds = this.scene.hbonds.filter(function (h) {
        return (h.ra === residue && h.rb === r.partner) || (h.rb === residue && h.ra === r.partner);
      });
      card.appendChild(el("div", "ds-card-row ds-strong", "Pairs with d" + partner.base + partner.num + " through " + bonds.length + " hydrogen bonds"));
      if (this.atomic) {
        var list = el("ul", "ds-card-list");
        bonds.forEach(function (h) {
          var li = el("li", "", h.text + "  " + h.d.toFixed(2) + " Å");
          if (h.role) li.title = h.role;
          list.appendChild(li);
        });
        card.appendChild(list);
        card.appendChild(el("div", "ds-card-foot", "Distances measured between the atoms of PDB " + this.scene.pdb + "."));
      }
    }
    if (this.atomic) {
      var counts = [0, 0, 0];
      for (var i = 0; i < this.n; i++) if (a.res[i] === residue) counts[a.comp[i]]++;
      var parts = el("div", "ds-card-parts");
      counts.forEach(function (count, index) {
        if (!count) return;
        var chip = el("span", "ds-chip");
        var sw = el("i", "ds-swatch");
        sw.style.background = COMPONENT_COLOUR[index];
        chip.appendChild(sw);
        chip.appendChild(document.createTextNode(COMPONENT_NAME[index] + ": " + count + " atoms"));
        parts.appendChild(chip);
      });
      card.appendChild(parts);
    }
    var close = el("button", "ds-card-close", "×");
    close.type = "button";
    close.setAttribute("aria-label", "Clear the selection");
    close.addEventListener("click", function () { self.select(-1); });
    card.appendChild(close);
    if (this.options.onSelect) this.options.onSelect(r);
  };

  function textWithHalo(ctx, text, x, y, colour, halo) {
    ctx.lineWidth = 3.5;
    ctx.strokeStyle = halo;
    ctx.strokeText(text, x, y);
    ctx.fillStyle = colour;
    ctx.fillText(text, x, y);
  }

  MoleculeView.prototype._overlay = function (viewer, ctx) {
    var scene = this.scene, st = this.state, self = this, a = scene.atoms, i, p;
    var halo = "rgba(8,16,30,0.85)", fs = this.options.fontScale || 1;
    function px(big, small) { return Math.round((self.options.bigLabels ? big : small) * fs * 10) / 10 + "px 'Segoe UI', system-ui, sans-serif"; }
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    if (st.labels) {
      scene.labels.forEach(function (label) {
        p = viewer.project(label.xyz);
        if (!p) return;
        if (label.kind === "end") {
          ctx.font = "700 " + px(24, 17);
          textWithHalo(ctx, label.text, p.x, p.y, STRAND_COLOUR[label.strand] || "#ffffff", halo);
        } else {
          // a groove label is shown while it stands beside the helix, not in front of or behind it
          var middle = viewer.project([0, label.xyz[1], 0]);
          if (middle && Math.abs(p.x - middle.x) < 9.5 * middle.scale) return;
          ctx.font = "600 " + px(19, 13);
          textWithHalo(ctx, label.text, p.x + (middle && p.x < middle.x ? -34 : 34), p.y, "#dfe8f3", halo);
        }
      });
    }
    if (this.options.atomLabels && st.labels) {
      ctx.font = "600 " + px(17, 12);
      for (i = 0; i < this.n; i++) {
        if (st.selected >= 0 && a.res[i] !== st.selected) continue;
        p = viewer.project([this.xyz[3 * i], this.xyz[3 * i + 1], this.xyz[3 * i + 2]]);
        if (p) textWithHalo(ctx, a.name[i], p.x, p.y - 0.62 * p.scale - 7 * fs, "#f4f7fb", halo);
      }
    }
    if (this.options.distanceLabels && st.hbonds) {
      ctx.font = "700 " + px(19, 13);
      scene.hbonds.forEach(function (h) {
        p = viewer.project([(h.a[0] + h.b[0]) / 2, (h.a[1] + h.b[1]) / 2, (h.a[2] + h.b[2]) / 2]);
        if (p) textWithHalo(ctx, h.d.toFixed(2) + " Å", p.x, p.y - 12 * fs, HBOND_COLOUR, halo);
      });
    }
    if (st.selected >= 0 && this.pickedAtom >= 0 && this.atomic && !this.options.atomLabels) {
      i = this.pickedAtom;
      p = viewer.project([this.xyz[3 * i], this.xyz[3 * i + 1], this.xyz[3 * i + 2]]);
      if (p) {
        ctx.font = "600 12px 'Segoe UI', system-ui, sans-serif";
        textWithHalo(ctx, a.name[i] + " (" + COMPONENT_NAME[a.comp[i]] + ")", p.x, p.y - 0.5 * p.scale - 10, "#ffffff", halo);
      }
    }
  };

  MoleculeView.prototype.dispose = function () {
    this.viewer.dispose();
    this.root.textContent = "";
  };

  NS.MoleculeView = MoleculeView;
})(typeof window !== "undefined" ? window : globalThis);
