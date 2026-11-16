/* DNA Studio, process player: step-by-step animations of replication, transcription and
   translation. A scenario (genomeasm.molbio.processes) gives the sequences and, for every step,
   the state at its end; the state at any time is interpolated, so each frame depends on time
   only and the animation can be paused, stepped and scrubbed. The drawings are schematic. */
(function (global) {
  "use strict";
  var NS = (global.DNAStudio = global.DNAStudio || {});
  var W = 1200, H = 520;
  var BASE = { A: "#00b383", T: "#e8662a", G: "#f2ae1c", C: "#2f8fdc", U: "#cc79a7" };
  var INK = "#e8eef6", SOFT = "#9fb0c5", OLD = "#8fa3bd", NEW = "#5ad1ff", RNA = "#f29ccb", HB = "#ffe6a3";
  var CLASS_COLOUR = { nonpolar: "#93a4bb", polar: "#56b4e9", acidic: "#e8662a", basic: "#00b383", stop: "#d94f4f" };
  var FONT = "'Segoe UI', system-ui, sans-serif";

  function clamp(v) { return v < 0 ? 0 : v > 1 ? 1 : v; }
  function ease(t) { t = clamp(t); return t * t * (3 - 2 * t); }
  function lerp(a, b, t) { return a + (b - a) * t; }

  function rounded(ctx, x, y, w, h, r) {
    ctx.beginPath();
    ctx.moveTo(x + r, y);
    ctx.arcTo(x + w, y, x + w, y + h, r);
    ctx.arcTo(x + w, y + h, x, y + h, r);
    ctx.arcTo(x, y + h, x, y, r);
    ctx.arcTo(x, y, x + w, y, r);
    ctx.closePath();
  }

  /* One nucleotide: a coloured tile with its letter. RNA nucleotides get a light ring. */
  function tile(ctx, x, y, letter, o) {
    o = o || {};
    var s = o.size || 30, w = o.width || s, alpha = o.alpha === undefined ? 1 : o.alpha;
    if (alpha <= 0.01) return;
    ctx.save();
    ctx.globalAlpha *= alpha;
    rounded(ctx, x - w / 2, y - s / 2, w, s, o.rna ? s * 0.34 : 6);
    ctx.fillStyle = BASE[letter] || "#667";
    ctx.fill();
    if (o.rna) { ctx.lineWidth = 2.2; ctx.strokeStyle = "rgba(255,255,255,0.85)"; ctx.stroke(); }
    if (o.ring) { ctx.lineWidth = 3; ctx.strokeStyle = o.ring; ctx.stroke(); }
    ctx.fillStyle = "#ffffff";
    ctx.font = "700 " + Math.round(s * 0.56) + "px " + FONT;
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText(letter, x, y + 1);
    ctx.restore();
  }

  function text(ctx, value, x, y, o) {
    o = o || {};
    ctx.save();
    if (o.alpha !== undefined) ctx.globalAlpha *= o.alpha;
    ctx.font = (o.weight || 600) + " " + (o.size || 15) + "px " + FONT;
    ctx.textAlign = o.align || "center";
    ctx.textBaseline = o.baseline || "middle";
    if (o.pill) {
      var m = ctx.measureText(value), pad = 8, h = (o.size || 15) + 10;
      var left = o.align === "left" ? x - pad : o.align === "right" ? x - m.width - pad : x - m.width / 2 - pad;
      rounded(ctx, left, y - h / 2, m.width + 2 * pad, h, h / 2);
      ctx.fillStyle = o.pill;
      ctx.fill();
    }
    ctx.fillStyle = o.color || INK;
    ctx.fillText(value, x, y);
    ctx.restore();
  }

  function dots(ctx, x, y, count, alpha, vertical) {
    if (alpha <= 0.01) return;
    ctx.save();
    ctx.globalAlpha *= alpha;
    ctx.fillStyle = HB;
    for (var k = 0; k < count; k++) {
      var off = (k - (count - 1) / 2) * 7;
      ctx.beginPath();
      ctx.arc(vertical ? x : x + off, vertical ? y + off : y, 2.3, 0, Math.PI * 2);
      ctx.fill();
    }
    ctx.restore();
  }

  function blob(ctx, x, y, rx, ry, colour, alpha, label, labelDy) {
    if (alpha <= 0.01) return;
    ctx.save();
    ctx.globalAlpha *= alpha;
    ctx.beginPath();
    ctx.ellipse(x, y, rx, ry, 0, 0, Math.PI * 2);
    ctx.fillStyle = colour.replace("ALPHA", "0.16");
    ctx.fill();
    ctx.lineWidth = 2;
    ctx.strokeStyle = colour.replace("ALPHA", "0.9");
    ctx.stroke();
    if (label) text(ctx, label, x, y + (labelDy || 0), { size: 13, weight: 700, color: colour.replace("ALPHA", "1") });
    ctx.restore();
  }

  function polyline(ctx, points, colour, width, alpha) {
    if (points.length < 2 || alpha <= 0.01) return;
    ctx.save();
    if (alpha !== undefined) ctx.globalAlpha *= alpha;
    ctx.beginPath();
    ctx.moveTo(points[0][0], points[0][1]);
    for (var i = 1; i < points.length; i++) ctx.lineTo(points[i][0], points[i][1]);
    ctx.strokeStyle = colour;
    ctx.lineWidth = width || 5;
    ctx.lineCap = "round";
    ctx.lineJoin = "round";
    ctx.stroke();
    ctx.restore();
  }

  function arrow(ctx, x1, y1, x2, y2, colour, alpha) {
    if (alpha <= 0.01) return;
    ctx.save();
    ctx.globalAlpha *= alpha === undefined ? 1 : alpha;
    var a = Math.atan2(y2 - y1, x2 - x1);
    ctx.strokeStyle = ctx.fillStyle = colour;
    ctx.lineWidth = 2.5;
    ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(x2, y2);
    ctx.lineTo(x2 - 11 * Math.cos(a - 0.42), y2 - 11 * Math.sin(a - 0.42));
    ctx.lineTo(x2 - 11 * Math.cos(a + 0.42), y2 - 11 * Math.sin(a + 0.42));
    ctx.closePath(); ctx.fill();
    ctx.restore();
  }

  function legend(ctx, items, x, y) {
    items.forEach(function (item) {
      ctx.save();
      ctx.strokeStyle = item[0];
      ctx.lineWidth = 5;
      ctx.lineCap = "round";
      ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x + 26, y); ctx.stroke();
      ctx.restore();
      text(ctx, item[1], x + 36, y, { size: 13, align: "left", color: SOFT, weight: 600 });
      x += 52 + item[1].length * 7.4;
    });
  }

  var HBONDS = { A: 2, T: 2, U: 2, G: 3, C: 3 };

  // ---- replication -------------------------------------------------------------------------------
  function drawReplication(ctx, s, d) {
    var n = d.top.length, dx = 40, x0 = 140, cy = 262, gap = 80, primer = d.primer, i, k;
    var apart = ease(s.apart);
    var swap = clamp(s.swap), ligate = clamp(s.ligate);
    function X(i) { return x0 + i * dx; }
    function open(i) { return ease((s.fork - i) / 3.2); }
    function yTop(i) { return cy - 24 - (gap + 34 * apart) * open(i); }
    function yBot(i) { return cy + 24 + (gap + 34 * apart) * open(i); }
    var busy = 1 - apart; // helper proteins fade when the copies separate

    // template strands (old)
    var topLine = [], botLine = [];
    for (i = 0; i < n; i++) { topLine.push([X(i), yTop(i) - 20]); botLine.push([X(i), yBot(i) + 20]); }
    polyline(ctx, topLine, OLD, 5);
    polyline(ctx, botLine, OLD, 5);
    for (i = 0; i < n; i++) {
      tile(ctx, X(i), yTop(i), d.top[i]);
      tile(ctx, X(i), yBot(i), d.bottom[i]);
      dots(ctx, X(i), cy, HBONDS[d.top[i]], 1 - clamp(open(i) * 5), true);
    }
    text(ctx, "5′", X(0) - 36, yTop(0), { size: 18, weight: 700, color: OLD });
    text(ctx, "3′", X(n - 1) + 36, yTop(n - 1), { size: 18, weight: 700, color: OLD });
    text(ctx, "3′", X(0) - 36, yBot(0), { size: 18, weight: 700, color: OLD });
    text(ctx, "5′", X(n - 1) + 36, yBot(n - 1), { size: 18, weight: 700, color: OLD });

    // a new nucleotide: RNA while it is a primer that has not been replaced yet
    function newTile(i, y, letter, isPrimer, amount, from) {
      var a = clamp(amount);
      if (a <= 0) return;
      var yy = y + from * 34 * (1 - ease(a));
      if (isPrimer && swap < 1) {
        tile(ctx, X(i), yy, letter === "T" ? "U" : letter, { rna: true, alpha: a * (1 - swap) });
        tile(ctx, X(i), yy, letter, { alpha: a * swap });
      } else {
        tile(ctx, X(i), yy, letter, { alpha: a });
      }
    }

    // leading strand: pairs with the bottom template, grows left to right (towards the fork)
    var lead = Math.min(s.lead, n), leadLine = [], primerLine = [];
    for (i = 0; i < n; i++) {
      var amount = lead - i;
      if (amount <= 0) break;
      var y = yBot(i) - 48;
      newTile(i, y, d.top[i], i < primer, amount, -1);
      dots(ctx, X(i), yBot(i) - 24, HBONDS[d.top[i]], clamp(amount), false);
      if (amount >= 0.5) { (i < primer ? primerLine : leadLine).push([X(i), y - 20]); if (i === primer - 1) leadLine.push([X(i) + dx / 2, y - 20]); }
    }
    if (primerLine.length) primerLine.push([primerLine[primerLine.length - 1][0] + (lead > primer ? dx / 2 : 0), primerLine[0][1]]);
    polyline(ctx, leadLine, NEW, 5);
    polyline(ctx, primerLine, swap >= 1 ? NEW : RNA, 5);
    if (lead > 0.5) {
      text(ctx, "5′", X(0) - 30, yBot(0) - 48, { size: 14, weight: 700, color: NEW });
      var tip = Math.min(Math.ceil(lead - 0.5), n) - 1;
      if (lead < n - 0.01 || apart > 0) text(ctx, "3′", X(tip) + 30, yBot(tip) - 48, { size: 14, weight: 700, color: NEW, alpha: lead >= n ? apart : 1 });
    }

    // lagging strand: pairs with the top template, each fragment grows right to left
    var joints = [];
    d.fragments.forEach(function (fragment, index) {
      var progress = s["f" + index] || 0, startI = fragment[0], endI = fragment[1], line = [], pline = [];
      for (var j = endI; j >= startI; j--) {
        var amount = progress - (endI - j);
        if (amount <= 0) break;
        var y = yTop(j) + 48, isPrimer = endI - j < primer;
        newTile(j, y, d.bottom[j], isPrimer, amount, 1);
        dots(ctx, X(j), yTop(j) + 24, HBONDS[d.bottom[j]], clamp(amount), false);
        if (amount >= 0.5) { (isPrimer ? pline : line).push([X(j), y + 20]); if (endI - j === primer - 1) line.push([X(j) - dx / 2, y + 20]); }
      }
      if (pline.length && progress > primer) pline.push([pline[pline.length - 1][0] - dx / 2, pline[0][1]]);
      polyline(ctx, line, NEW, 5);
      polyline(ctx, pline, swap >= 1 ? NEW : RNA, 5);
      if (progress > 0.5) {
        text(ctx, "5′", X(endI) + 11, yTop(endI) + 48 + 34, { size: 12, weight: 700, color: swap >= 1 ? NEW : RNA, alpha: 1 - ligate });
        if (progress >= endI - startI + 1 && index > 0) joints.push([X(startI) - dx / 2, yTop(startI) + 48 + 20, index]);
      }
      if (progress >= endI - startI + 0.9) {
        text(ctx, "Okazaki fragment " + (index + 1), X((startI + endI) / 2), yTop(startI) + 48 + 44, { size: 12.5, color: SOFT, alpha: (1 - ligate) * busy });
      }
      // polymerase or primase working on this fragment
      if (progress > 0 && progress < endI - startI + 1) {
        var at = X(endI - progress + 0.5), yy = yTop(Math.max(startI, Math.round(endI - progress))) + 48;
        if (progress <= primer) blob(ctx, at, yy, 26, 24, "rgba(242,156,203,ALPHA)", busy, "primase", 40);
        else blob(ctx, at, yy, 34, 28, "rgba(90,209,255,ALPHA)", busy, "DNA pol", 42);
      }
    });
    // nicks between fragments: ligase closes them
    joints.forEach(function (joint) {
      polyline(ctx, [[joint[0] - dx / 2, joint[1]], [joint[0] + dx / 2, joint[1]]], NEW, 5, ligate);
      if (s.ligate > 0 && ligate < 1) blob(ctx, joint[0], joint[1] + 2, 20, 16, "rgba(255,230,163,ALPHA)", busy * (swap >= 1 ? 1 : 0.4), "ligase", 30);
    });

    // proteins at the fork
    if (s.fork > 0.01 && s.fork < n + 1.5) {
      var xf = X(s.fork - 0.4), fade = clamp(s.fork) * clamp(n + 1.5 - s.fork) * busy;
      ctx.save();
      ctx.globalAlpha = fade;
      ctx.beginPath();
      for (k = 0; k < 6; k++) {
        var ang = (Math.PI / 3) * k + Math.PI / 6;
        ctx[k ? "lineTo" : "moveTo"](xf + 30 * Math.cos(ang), cy + 30 * Math.sin(ang));
      }
      ctx.closePath();
      ctx.fillStyle = "rgba(242,174,28,0.28)"; ctx.fill();
      ctx.lineWidth = 2.5; ctx.strokeStyle = "#f2ae1c"; ctx.stroke();
      ctx.restore();
      text(ctx, "helicase", xf + 44, cy - 64, { size: 13.5, weight: 700, color: "#f2ae1c", alpha: fade, align: "left" });
      if (s.fork < n) arrow(ctx, xf + 112, cy - 64, xf + 160, cy - 64, "#f2ae1c", fade);
    }
    // single-strand binding proteins on template not yet copied
    for (i = 0; i < n; i++) {
      if (open(i) < 0.75) continue;
      var a = busy * clamp(open(i) * 4 - 3);
      var laggingCopied = d.fragments.some(function (fragment, index) {
        return i >= fragment[0] && i <= fragment[1] && (s["f" + index] || 0) - (fragment[1] - i) > 0.2;
      });
      ctx.save();
      ctx.fillStyle = "rgba(185,163,255,0.9)";
      if (!laggingCopied) { ctx.globalAlpha = a; ctx.beginPath(); ctx.arc(X(i), yTop(i) + 30, 6, 0, 7); ctx.fill(); }
      if (lead - i <= 0.2) { ctx.globalAlpha = a; ctx.beginPath(); ctx.arc(X(i), yBot(i) - 30, 6, 0, 7); ctx.fill(); }
      ctx.restore();
    }
    if (lead > primer && lead < n) blob(ctx, X(lead - 0.5), yBot(Math.max(0, Math.round(lead) - 1)) - 48, 34, 28, "rgba(90,209,255,ALPHA)", busy, "DNA pol", -42);
    else if (lead > 0 && lead <= primer && s.fork < n) blob(ctx, X(lead - 0.5), yBot(0) - 48, 26, 24, "rgba(242,156,203,ALPHA)", busy * clamp((primer + 0.01 - lead) * 50 + (s.f0 < primer ? 1 : 0)), "primase", -40);

    // captions inside the drawing
    var named = clamp(s.fork / 6) * busy;
    if (lead > primer) text(ctx, "leading strand: continuous, towards the fork", X(Math.min(lead, n) / 2), yBot(0) + 52, { size: 14, color: NEW, alpha: named, weight: 700 });
    if ((s.f0 || 0) > primer) text(ctx, "lagging strand: in fragments, away from the fork", X(Math.min(s.fork, n) / 2), yTop(0) - 46, { size: 14, color: NEW, alpha: named, weight: 700 });
    if (apart > 0) {
      text(ctx, "daughter helix 1: old strand + new strand", 600, yTop(0) - 46, { size: 16, weight: 700, alpha: apart });
      text(ctx, "daughter helix 2: old strand + new strand", 600, yBot(0) + 50, { size: 16, weight: 700, alpha: apart });
    }
    legend(ctx, [[OLD, "old (template) strand"], [NEW, "new DNA"], [RNA, "RNA primer"]], 40, 30);
    ctx.save(); ctx.fillStyle = "rgba(185,163,255,0.9)"; ctx.beginPath(); ctx.arc(668, 30, 6, 0, 7); ctx.fill(); ctx.restore();
    text(ctx, "single-strand binding protein", 682, 30, { size: 13, align: "left", color: SOFT });
  }

  // ---- transcription -----------------------------------------------------------------------------
  function drawTranscription(ctx, s, d) {
    var n = d.coding.length, dx = 34, x0 = 600 - ((n - 1) * dx) / 2, cy = 150, flank = d.flank, i, j;
    var total = d.rna.length, pos = Math.min(s.pos, total), hybrid = 5, lift = 120;
    var release = ease(s.release), compare = ease(s.compare);
    function X(i) { return x0 + i * dx; }
    var centre = flank + pos - 2;
    function bump(i) { return s.bubble * ease(1 - Math.abs(i - centre) / 5.2); }
    function yCoding(i) { return cy - 24 - 50 * bump(i); }
    var yTemplate = cy + 24;

    // gene arrow and promoter mark
    arrow(ctx, X(flank) - dx / 2, cy - 104, X(flank + total - 1) + dx / 2, cy - 104, SOFT, 1);
    text(ctx, "gene " + d.gene + " (first " + total + " letters): transcribed region", X(flank + (total - 1) / 2), cy - 122, { size: 13.5, color: SOFT });
    ctx.save();
    ctx.strokeStyle = "#f2ae1c"; ctx.lineWidth = 3;
    ctx.beginPath(); ctx.moveTo(X(0) - dx / 2 - 4, cy - 86); ctx.lineTo(X(0) - dx / 2 - 4, cy - 104); ctx.lineTo(X(flank) - dx / 2 - 8, cy - 104); ctx.stroke();
    ctx.restore();
    text(ctx, "promoter", X(0) - dx / 2 - 6, cy - 122, { size: 13.5, color: "#f2ae1c", weight: 700, align: "left" });

    var codingLine = [], templateLine = [];
    for (i = 0; i < n; i++) { codingLine.push([X(i), yCoding(i) - 20]); templateLine.push([X(i), yTemplate + 20]); }
    polyline(ctx, codingLine, OLD, 5);
    polyline(ctx, templateLine, OLD, 5, 1 - 0.7 * compare);
    for (i = 0; i < n; i++) {
      var transcribed = i >= flank && i < flank + total;
      var ringed = compare > 0 && transcribed && d.coding[i] === "T";
      tile(ctx, X(i), yCoding(i), d.coding[i], { ring: ringed ? "rgba(255,255,255," + compare + ")" : null });
      tile(ctx, X(i), yTemplate, d.template[i], { alpha: 1 - 0.7 * compare });
      dots(ctx, X(i), (yCoding(i) + yTemplate) / 2, HBONDS[d.coding[i]], (1 - clamp(bump(i) * 4)) * (1 - 0.7 * compare), true);
    }
    text(ctx, "5′", X(0) - 34, cy - 24, { size: 17, weight: 700, color: OLD });
    text(ctx, "3′", X(n - 1) + 34, cy - 24, { size: 17, weight: 700, color: OLD });
    text(ctx, "3′", X(0) - 34, yTemplate, { size: 17, weight: 700, color: OLD, alpha: 1 - 0.7 * compare });
    text(ctx, "5′", X(n - 1) + 34, yTemplate, { size: 17, weight: 700, color: OLD, alpha: 1 - 0.7 * compare });
    text(ctx, "coding strand", X(n - 1) + 58, cy - 24, { size: 13.5, align: "left", color: SOFT });
    text(ctx, "template strand", X(n - 1) + 58, yTemplate, { size: 13.5, align: "left", color: SOFT, alpha: 1 - 0.7 * compare });

    // RNA polymerase
    if (s.enzyme > 0.01) {
      var px = X(centre + 0.5);
      ctx.save();
      ctx.globalAlpha = s.enzyme;
      rounded(ctx, px - 205, cy - 92, 410, 216, 96);
      ctx.fillStyle = "rgba(126,160,255,0.13)"; ctx.fill();
      ctx.lineWidth = 2; ctx.strokeStyle = "rgba(160,190,255,0.75)"; ctx.stroke();
      ctx.restore();
      text(ctx, "RNA polymerase", px, cy + 142, { size: 14, weight: 700, color: "#a9c1ff", alpha: s.enzyme });
      if (pos < total) arrow(ctx, px + 214, cy + 60, px + 262, cy + 60, "#a9c1ff", s.enzyme);
    }

    // the RNA: young nucleotides pair with the template, older ones trail away as a single strand
    var rowY = yTemplate + 150, line = [];
    for (j = 0; j < total; j++) {
      var amount = pos - j;
      if (amount <= 0) break;
      var a = clamp(amount), xi = X(flank + j), age = pos - 1 - j - (hybrid - 1);
      var x = xi, y = yTemplate + 48;
      if (age > 0) { x = X(flank + pos - hybrid) - age * dx * 0.86 + dx * 0.0; y = yTemplate + 48 + lift * ease(Math.min(age, 3.5) / 3.5); }
      if (age > 0 && age < 1) x = lerp(xi, x, ease(age));
      y += 30 * (1 - ease(a)); // arriving ribonucleotide
      var xr = lerp(x, xi, release), yr = lerp(y, rowY, release) - 62 * compare;
      var ring = compare > 0 && d.rna[j] === "U" ? "rgba(255,255,255," + compare + ")" : null;
      tile(ctx, xr, yr, d.rna[j], { rna: true, alpha: a, ring: ring });
      if (age <= 0 && release < 0.5) dots(ctx, xi, yTemplate + 24, HBONDS[d.rna[j]], a * (1 - release * 2), false);
      if (a >= 0.5) line.push([xr, yr + 20]);
      if (compare > 0.02) {
        ctx.save(); ctx.globalAlpha = compare * 0.55; ctx.strokeStyle = SOFT; ctx.setLineDash([3, 5]); ctx.lineWidth = 1.5;
        ctx.beginPath(); ctx.moveTo(xi, cy - 6); ctx.lineTo(xi, yr - 18); ctx.stroke(); ctx.restore();
      }
    }
    polyline(ctx, line, RNA, 5);
    if (line.length) {
      var first = line[0], last = line[line.length - 1];
      text(ctx, "5′", first[0] - 30 * (release > 0.5 ? 1 : 0.86), first[1] - 20, { size: 15, weight: 700, color: RNA });
      text(ctx, "3′", last[0] + 30, last[1] - 20, { size: 15, weight: 700, color: RNA, alpha: pos >= total ? release : 1 });
      text(ctx, release > 0.5 ? "messenger RNA (single strand)" : "new RNA", lerp(first[0], last[0], 0.5), Math.max(first[1], last[1]) + 26, { size: 13.5, color: RNA, weight: 700, alpha: clamp(pos / 4) });
    }
    if (compare > 0.02) text(ctx, "same letters as the coding strand, with U where DNA has T", 600, rowY + 34 - 62 * compare + 52, { size: 15, weight: 700, alpha: compare });
    legend(ctx, [[OLD, "DNA"], [RNA, "RNA"]], 40, H - 26);
  }

  // ---- translation -------------------------------------------------------------------------------
  function drawTranslation(ctx, s, d) {
    var codons = d.codons, nAmino = codons.length - 1, cw = 84, cx = 560, my = 392, i, k;
    var cycle = Math.min(s.cycle, nAmino - 1), kP = Math.floor(cycle + 1e-9), frac = cycle - kP;
    if (kP >= nAmino - 1) { kP = nAmino - 1; frac = 0; }
    var p = kP + ease((frac - 0.6) / 0.4);
    var arrive = ease(frac / 0.4), transfer = ease((frac - 0.4) / 0.2);
    var ribo = s.ribo, free = ease(s.free), rf = ease(s.rf);
    function CX(k) { return cx + (k - p) * cw; }

    // ribosome (behind everything)
    var small = clamp(ribo * 2), large = clamp(ribo * 2 - 1);
    if (small > 0.01) {
      ctx.save(); ctx.globalAlpha = small;
      ctx.beginPath(); ctx.ellipse(cx, my + 18, 226, 60, 0, 0, Math.PI * 2);
      ctx.fillStyle = "rgba(126,160,255,0.16)"; ctx.fill();
      ctx.lineWidth = 2; ctx.strokeStyle = "rgba(160,190,255,0.6)"; ctx.stroke(); ctx.restore();
      text(ctx, "small subunit", cx + 150, my + 64, { size: 13, color: "#a9c1ff", align: "left", alpha: small });
    }
    if (large > 0.01) {
      ctx.save(); ctx.globalAlpha = large;
      ctx.beginPath(); ctx.ellipse(cx, my - 150, 232, 138, 0, 0, Math.PI * 2);
      ctx.fillStyle = "rgba(126,160,255,0.11)"; ctx.fill();
      ctx.lineWidth = 2; ctx.strokeStyle = "rgba(160,190,255,0.6)"; ctx.stroke(); ctx.restore();
      text(ctx, "large subunit", cx + 250, my - 150, { size: 13, color: "#a9c1ff", align: "left", alpha: large });
      [["E", -1], ["P", 0], ["A", 1]].forEach(function (site) {
        text(ctx, site[0], cx + site[1] * cw, my + 62, { size: 17, weight: 700, color: "#a9c1ff", alpha: large });
      });
      text(ctx, "tRNA sites", cx - cw - 34, my + 62, { size: 12.5, color: SOFT, align: "right", alpha: large });
    }

    // mRNA tape
    var leftEdge = Math.max(30, CX(0) - cw / 2 - 150), rightEdge = Math.min(W - 30, CX(codons.length - 1) + cw / 2 + 60);
    polyline(ctx, [[leftEdge, my + 20], [rightEdge, my + 20]], RNA, 5);
    if (leftEdge < CX(0) - cw / 2 - 20) {
      text(ctx, "5′", leftEdge - 2, my - 2, { size: 16, weight: 700, color: RNA, align: "left" });
      if (CX(0) - cw / 2 - leftEdge > 120) text(ctx, "untranslated leader", leftEdge + 22, my - 2, { size: 12.5, color: SOFT, align: "left" });
    }
    if (rightEdge < W - 40) text(ctx, "3′", rightEdge + 16, my + 2, { size: 16, weight: 700, color: RNA });
    for (k = 0; k < codons.length; k++) {
      var x = CX(k);
      if (x < -cw || x > W + cw) continue;
      for (i = 0; i < 3; i++) tile(ctx, x + (i - 1) * 26, my, codons[k].codon[i], { rna: true, size: 28, width: 24 });
      var mark = k === 0 ? "start" : k === codons.length - 1 ? "stop" : String(k + 1);
      text(ctx, mark, x, my + 38, { size: 12, color: k === 0 ? "#7fe0b6" : k === codons.length - 1 ? "#ff8f8f" : SOFT, weight: k === 0 || k === codons.length - 1 ? 700 : 600 });
      if (k === 0 || k === codons.length - 1) {
        ctx.save(); ctx.strokeStyle = k === 0 ? "#7fe0b6" : "#ff8f8f"; ctx.lineWidth = 2;
        rounded(ctx, x - 41, my - 18, 82, 36, 9); ctx.stroke(); ctx.restore();
      }
    }

    // where a bead of the chain sits: index 0 is the newest residue, on top of its tRNA
    function chainPoint(anchorX, q) {
      var top = my - 172;
      if (q <= 2) return [anchorX, top - q * 38];
      var col = q - 2, perRow = 13, row = Math.floor((col - 1) / perRow), within = (col - 1) % perRow;
      var xx = row % 2 === 0 ? anchorX - (within + 1) * 38 : anchorX - (perRow - within) * 38;
      return [xx, top - 76 - row * 38];
    }
    function bead(x, y, codon, alpha, big) {
      if (alpha <= 0.01) return;
      ctx.save(); ctx.globalAlpha *= alpha;
      ctx.beginPath(); ctx.arc(x, y, big ? 19 : 17, 0, Math.PI * 2);
      ctx.fillStyle = CLASS_COLOUR[codon.cls] || "#93a4bb"; ctx.fill();
      ctx.lineWidth = 2; ctx.strokeStyle = "rgba(255,255,255,0.85)"; ctx.stroke();
      ctx.fillStyle = "#fff"; ctx.font = "700 12.5px " + FONT; ctx.textAlign = "center"; ctx.textBaseline = "middle";
      ctx.fillText(codon.aa, x, y + 1);
      ctx.restore();
    }
    function tRNA(k, ox, oy, alpha, loaded) {
      if (alpha <= 0.01 || k < 0 || k >= nAmino) return;
      var x = CX(k) + ox, y = my + oy;
      ctx.save(); ctx.globalAlpha *= alpha;
      rounded(ctx, x - 25, y - 150, 50, 114, 20);
      ctx.fillStyle = "rgba(255,255,255,0.13)"; ctx.fill();
      ctx.lineWidth = 2; ctx.strokeStyle = "rgba(220,230,245,0.8)"; ctx.stroke();
      ctx.restore();
      for (var m = 0; m < 3; m++) {
        tile(ctx, x + (m - 1) * 26, y - 34, codons[k].anticodon[m], { rna: true, size: 24, width: 22, alpha: alpha });
        if (Math.abs(ox) < 2 && Math.abs(oy) < 2) dots(ctx, x + (m - 1) * 26, y - 18, HBONDS[codons[k].codon[m]] > 2 ? 3 : 2, alpha, false);
      }
      text(ctx, "tRNA", x, y - 100, { size: 12, color: SOFT, alpha: alpha });
      if (loaded) bead(x, y - 172, codons[k], alpha, true);
    }

    if (ribo > 0.5 || free > 0) {
      var show = clamp(ribo * 2 - 1) * (1 - free);
      // tRNA leaving through the E site
      if (kP >= 1) tRNA(kP - 1, -150 * arrive, -110 * arrive, show * (1 - arrive), false);
      // tRNA in the P site: carries the chain until the peptide bond is made
      tRNA(kP, 0, 0, show, false);
      // tRNA arriving in the A site with the next amino acid
      var hasNext = kP + 1 < nAmino && frac > 0;
      if (hasNext) tRNA(kP + 1, 170 * (1 - arrive), -130 * (1 - arrive), show * clamp(arrive * 1.5), transfer < 0.02);
      // release factor in the A site at the stop codon
      if (rf > 0.01) {
        var rx = CX(nAmino) + 170 * (1 - rf), ry = my - 96 - 130 * (1 - rf);
        ctx.save(); ctx.globalAlpha = rf * (1 - free);
        rounded(ctx, rx - 30, ry - 56, 60, 116, 26);
        ctx.fillStyle = "rgba(232,102,42,0.32)"; ctx.fill(); ctx.lineWidth = 2; ctx.strokeStyle = "#ff9d6b"; ctx.stroke(); ctx.restore();
        text(ctx, "release", rx, ry - 8, { size: 12.5, weight: 700, color: "#ffb48c", alpha: rf * (1 - free) });
        text(ctx, "factor", rx, ry + 8, { size: 12.5, weight: 700, color: "#ffb48c", alpha: rf * (1 - free) });
      }
      // the growing chain
      var newest = kP, anchorBefore = CX(kP), anchorAfter = CX(kP + 1);
      var perLine = 13, gridLeft = 600 - ((perLine - 1) * 42) / 2, nEnd = null;
      for (var m2 = 0; m2 <= newest + (hasNext && transfer >= 0.02 ? 1 : 0); m2++) {
        var before = chainPoint(anchorBefore, newest - m2);
        var after = hasNext ? chainPoint(anchorAfter, newest + 1 - m2) : before;
        if (m2 === newest + 1) before = after;
        var bx = lerp(before[0], after[0], transfer), by = lerp(before[1], after[1], transfer);
        if (free > 0) {
          bx = lerp(bx, gridLeft + (m2 % perLine) * 42, free);
          by = lerp(by, 150 + Math.floor(m2 / perLine) * 46, free);
        }
        if (bx > -30 && bx < W + 30) bead(bx, by, codons[m2], clamp(ribo * 2 - 1) + free, false);
        if (m2 === 0) nEnd = [bx, by];
      }
      if (free > 0.5) {
        text(ctx, "N-terminus", gridLeft - 30, 150, { size: 13, color: SOFT, align: "right", alpha: free });
        var lastCol = (nAmino - 1) % perLine, lastRow = Math.floor((nAmino - 1) / perLine);
        text(ctx, "C-terminus", gridLeft + lastCol * 42 + 30, 150 + lastRow * 46, { size: 13, color: SOFT, align: "left", alpha: free });
        text(ctx, "protein " + d.gene + ": " + d.protein, 600, 150 + (lastRow + 1) * 46 + 16, { size: 15, weight: 700, alpha: free });
      } else if (newest >= 3 && free <= 0) {
        if (nEnd) text(ctx, "N", nEnd[0] - 28, nEnd[1], { size: 12.5, color: SOFT, alpha: show });
      }
      // codon being decoded
      if (show > 0.01) {
        var decode = frac > 0 && kP + 1 < nAmino ? kP + 1 : kP;
        var inA = frac > 0 && kP + 1 < nAmino, isStop = rf > 0.01;
        var c = isStop ? codons[nAmino] : codons[decode];
        text(ctx, (isStop || inA ? "A site" : "P site") + ":  " + c.codon + "  →  " + (isStop ? "stop" : c.aa), W - 40, 40, { size: 17, weight: 700, align: "right", pill: "rgba(255,255,255,0.09)", alpha: show });
        text(ctx, "amino acid " + Math.min(decode + 1, nAmino) + " of " + nAmino, W - 40, 72, { size: 13, color: SOFT, align: "right", alpha: show });
      }
    }
    // amino-acid classes
    var lx = 40;
    [["nonpolar", "nonpolar"], ["polar", "polar"], ["acidic", "acidic"], ["basic", "basic"]].forEach(function (item) {
      ctx.save(); ctx.fillStyle = CLASS_COLOUR[item[0]]; ctx.beginPath(); ctx.arc(lx, H - 24, 7, 0, 7); ctx.fill(); ctx.restore();
      text(ctx, item[1], lx + 13, H - 24, { size: 13, align: "left", color: SOFT });
      lx += 40 + item[1].length * 7.6;
    });
    text(ctx, "amino-acid side chain:", 40, H - 46, { size: 12, align: "left", color: SOFT });
  }

  var SCENES = { replication: drawReplication, transcription: drawTranscription, translation: drawTranslation };

  /* ProcessPlayer(root, scenario, options)
       options.compact: no note line (used inside slides); options.autoplay: start playing */
  function ProcessPlayer(root, scenario, options) {
    var o = options || {}, self = this;
    this.scenario = scenario;
    this.steps = scenario.steps;
    this.starts = [];
    var t = 0;
    this.steps.forEach(function (step) { self.starts.push(t); t += step.dur; });
    this.total = t;
    this.time = 0;
    this.speed = 1;
    this.playing = false;
    this.stopAt = null;
    this.disposed = false;
    root.classList.add("ds-proc");
    this.root = root;
    this.stage = NS.el("div", "ds-proc-stage");
    this.canvas = document.createElement("canvas");
    this.canvas.setAttribute("role", "img");
    this.canvas.setAttribute("aria-label", "Animation: " + scenario.title);
    this.stage.appendChild(this.canvas);
    root.appendChild(this.stage);
    this.ctx = this.canvas.getContext("2d");

    var caption = NS.el("div", "ds-proc-caption");
    this.stepNo = NS.el("div", "ds-step-no", "1");
    this.stepTitle = NS.el("div", "ds-step-title");
    this.stepText = NS.el("div", "ds-step-text");
    caption.appendChild(this.stepNo); caption.appendChild(this.stepTitle); caption.appendChild(this.stepText);
    caption.setAttribute("aria-live", "polite");
    root.appendChild(caption);

    var controls = NS.el("div", "ds-controls");
    function button(label, cls, handler, title) {
      var b = NS.el("button", "ds-btn" + (cls ? " " + cls : ""), label);
      b.type = "button";
      if (title) b.title = title;
      b.addEventListener("click", handler);
      controls.appendChild(b);
      return b;
    }
    this.prevButton = button("◀ Step", "", function () { self.prev(); }, "Back to the start of this step, or to the previous step");
    this.playButton = button("▶ Play", "primary", function () { self.toggle(); }, "Play or pause the whole animation");
    this.nextButton = button("Step ▶", "", function () { self.next(); }, "Play the next step and stop");
    this.dots = NS.el("div", "ds-dots");
    this.steps.forEach(function (step, index) {
      var dot = NS.el("button", "ds-dot");
      dot.type = "button";
      dot.title = (index + 1) + ". " + step.title;
      dot.setAttribute("aria-label", "Step " + (index + 1) + ": " + step.title);
      dot.addEventListener("click", function () { self.goTo(index, true); });
      self.dots.appendChild(dot);
    });
    controls.appendChild(this.dots);
    this.scrub = document.createElement("input");
    this.scrub.type = "range"; this.scrub.className = "ds-scrub";
    this.scrub.min = "0"; this.scrub.max = String(this.total); this.scrub.step = "0.01"; this.scrub.value = "0";
    this.scrub.setAttribute("aria-label", "Position in the animation");
    this.scrub.addEventListener("input", function () { self.pause(); self.time = parseFloat(self.scrub.value); self.draw(); });
    controls.appendChild(this.scrub);
    var speed = NS.el("label", "ds-speed", "Speed ");
    var select = document.createElement("select");
    [["0.5", "0.5×"], ["1", "1×"], ["2", "2×"]].forEach(function (choice) {
      var option = document.createElement("option");
      option.value = choice[0]; option.textContent = choice[1];
      if (choice[0] === "1") option.selected = true;
      select.appendChild(option);
    });
    select.addEventListener("change", function () { self.speed = parseFloat(select.value); });
    speed.appendChild(select);
    controls.appendChild(speed);
    root.appendChild(controls);
    if (!o.compact && scenario.note) root.appendChild(NS.el("div", "ds-legend ds-note", scenario.note));

    this._onKey = function (event) {
      if (event.target && /INPUT|SELECT|TEXTAREA/.test(event.target.tagName)) return;
      if (event.key === " ") { event.preventDefault(); self.toggle(); }
      else if (event.key === "ArrowRight") { event.preventDefault(); self.next(); }
      else if (event.key === "ArrowLeft") { event.preventDefault(); self.prev(); }
    };
    if (!o.noKeys) document.addEventListener("keydown", this._onKey);
    if (global.ResizeObserver) {
      this.observer = new ResizeObserver(function () { self.draw(); });
      this.observer.observe(this.stage);
    }
    this.draw();
    if (o.autoplay) this.play();
  }

  ProcessPlayer.prototype.stepIndex = function (time) {
    var index = 0;
    for (var i = 0; i < this.starts.length; i++) if (time >= this.starts[i] - 1e-6) index = i;
    return index;
  };

  /* stateAt(time): the interpolated state of the scene (exposed for the tests). */
  ProcessPlayer.prototype.stateAt = function (time) {
    time = Math.min(Math.max(time, 0), this.total);
    var index = this.stepIndex(time), step = this.steps[index];
    if (index === 0) return Object.assign({}, step.to);
    var p = clamp((time - this.starts[index]) / step.dur), before = this.steps[index - 1].to, state = {};
    var win = step.win || {};
    Object.keys(step.to).forEach(function (key) {
      var w = win[key], q = w ? clamp((p - w[0]) / (w[1] - w[0])) : p; // a value may change in part of a step
      state[key] = lerp(before[key], step.to[key], q);
    });
    return state;
  };

  ProcessPlayer.prototype.draw = function () {
    if (this.disposed) return;
    var rect = this.stage.getBoundingClientRect(), dpr = Math.min(global.devicePixelRatio || 1, 2);
    var w = Math.max(1, Math.round(rect.width * dpr)), h = Math.max(1, Math.round(rect.height * dpr));
    if (this.canvas.width !== w || this.canvas.height !== h) { this.canvas.width = w; this.canvas.height = h; }
    var ctx = this.ctx, scale = Math.min(w / W, h / H);
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, w, h);
    ctx.setTransform(scale, 0, 0, scale, (w - W * scale) / 2, (h - H * scale) / 2);
    SCENES[this.scenario.kind](ctx, this.stateAt(this.time), this.scenario.data);
    var index = this.stepIndex(this.time), step = this.steps[index];
    if (this.shown !== index) {
      this.shown = index;
      this.stepNo.textContent = String(index + 1);
      this.stepTitle.textContent = step.title;
      this.stepText.innerHTML = step.text; // text written in this project (processes.py), not user input
      Array.prototype.forEach.call(this.dots.children, function (dot, i) {
        dot.setAttribute("aria-current", String(i === index));
        dot.classList.toggle("done", i < index);
      });
    }
    this.scrub.value = String(this.time);
    this.playButton.textContent = this.playing ? "❚❚ Pause" : this.time >= this.total - 1e-6 ? "↺ Replay" : "▶ Play";
    this.prevButton.disabled = this.time <= 0;
    this.nextButton.disabled = this.time >= this.total - 1e-6;
  };

  ProcessPlayer.prototype._tick = function (now) {
    if (!this.playing || this.disposed) return;
    var dt = this.last ? Math.min(0.1, (now - this.last) / 1000) : 0;
    this.last = now;
    this.time += dt * this.speed;
    var limit = this.stopAt === null ? this.total : this.stopAt;
    if (this.time >= limit) { this.time = limit; this.playing = false; this.stopAt = null; }
    this.draw();
    var self = this;
    if (this.playing) global.requestAnimationFrame(function (t) { self._tick(t); });
  };

  ProcessPlayer.prototype.play = function (stopAt) {
    if (this.time >= this.total - 1e-6) this.time = 0;
    this.stopAt = stopAt === undefined ? null : stopAt;
    if (this.playing) return;
    this.playing = true;
    this.last = 0;
    var self = this;
    global.requestAnimationFrame(function (t) { self._tick(t); });
  };
  ProcessPlayer.prototype.pause = function () { this.playing = false; this.stopAt = null; this.draw(); };
  ProcessPlayer.prototype.toggle = function () { if (this.playing) this.pause(); else this.play(); };

  /* next(): play the following step and stop at its end. Returns false at the very end. */
  ProcessPlayer.prototype.next = function () {
    if (this.time >= this.total - 1e-6) return false;
    var index = this.stepIndex(this.time), end = this.starts[index] + this.steps[index].dur;
    if (index === 0 || end - this.time < 0.05) {
      // step 1 is a still picture, and a finished step has nothing left: move on to the next one
      index = Math.min(index + 1, this.steps.length - 1);
      this.time = this.starts[index];
      end = this.starts[index] + this.steps[index].dur;
    }
    this.playing = false;
    this.play(end);
    return true;
  };
  ProcessPlayer.prototype.prev = function () {
    if (this.time <= 0) return false;
    var index = this.stepIndex(this.time);
    if (this.time - this.starts[index] < 0.3 && index > 0) index--;
    this.goTo(index, false);
    return true;
  };
  ProcessPlayer.prototype.goTo = function (index, playIt) {
    this.playing = false;
    this.stopAt = null;
    this.time = this.starts[index];
    if (playIt && index > 0) this.play(this.starts[index] + this.steps[index].dur); else this.draw();
  };
  ProcessPlayer.prototype.dispose = function () {
    this.disposed = true;
    this.playing = false;
    if (this.observer) this.observer.disconnect();
    document.removeEventListener("keydown", this._onKey);
    this.root.textContent = "";
  };

  NS.ProcessPlayer = ProcessPlayer;
  NS.processScenes = SCENES;
})(typeof window !== "undefined" ? window : globalThis);
