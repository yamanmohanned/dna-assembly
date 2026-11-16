/* DNA Studio, lecture deck: shows the slides made by genomeasm.appkit.lecture.
   Keys: right arrow / space / page down = forward (an animation on the slide advances one step
   first), left arrow / page up = back, Home / End, F = full screen, O = overview, P = play/pause.
   3D viewers and animations are created when their slide is shown and freed when it is left. */
(function (global) {
  "use strict";
  var NS = (global.DNAStudio = global.DNAStudio || {});
  var el = NS.el;
  var STAGE_W = 1600, STAGE_H = 900;

  function html(tag, cls, markup) {
    var node = el(tag, cls);
    node.innerHTML = markup; // lecture text is written in this project, never typed by a visitor
    return node;
  }

  function Deck(root, lecture, options) {
    var self = this;
    this.root = root;
    this.lecture = lecture;
    this.options = options || {};
    this.slides = lecture.slides;
    this.index = -1;
    this.live = []; // viewers and players of the current slide
    this.started = 0;
    root.tabIndex = 0;
    this.stage = el("div", "deck-stage");
    root.appendChild(this.stage);
    this.nodes = this.slides.map(function (slide, i) { return self._build(slide, i); });
    this.progress = el("div", "deck-progress");
    this.stage.appendChild(this.progress);
    this._chrome();
    this._overview();
    this.toastNode = el("div", "deck-toast");
    this.stage.appendChild(this.toastNode);

    this._onKey = function (event) { self._key(event); };
    document.addEventListener("keydown", this._onKey);
    global.addEventListener("resize", function () { self.fit(); });
    document.addEventListener("fullscreenchange", function () { global.setTimeout(function () { self.fit(); }, 60); });
    this.fit();
    var start = 0;
    if (!this.options.embedded && /^#\d+$/.test(global.location.hash)) {
      start = Math.min(this.slides.length - 1, Math.max(0, parseInt(global.location.hash.slice(1), 10) - 1));
    }
    this.go(start);
    root.focus();
  }

  Deck.prototype._build = function (slide, i) {
    var self = this, node = el("section", "slide layout-" + slide.layout + (slide.cls ? " " + slide.cls : ""));
    node.setAttribute("aria-label", "Slide " + (i + 1) + ": " + (slide.title || ""));
    if (slide.layout === "title") {
      var grid = el("div", "title-grid"), text = el("div", "title-text");
      text.appendChild(el("div", "title-kicker", slide.kicker || ""));
      text.appendChild(html("h1", "", slide.title));
      text.appendChild(html("div", "title-sub", slide.lead || ""));
      if (slide.meta) text.appendChild(html("div", "title-meta", slide.meta));
      grid.appendChild(text);
      var figure = el("div", "title-figure");
      figure.dataset.figure = "1";
      grid.appendChild(figure);
      node.appendChild(grid);
      return node;
    }
    if (slide.kicker) node.appendChild(el("div", "slide-kicker", slide.kicker));
    node.appendChild(html("h2", "slide-title", slide.title));
    var body = el("div", "slide-body");
    node.appendChild(body);
    function textColumn() {
      var column = el("div", "slide-text");
      if (slide.lead) column.appendChild(html("p", "lead", slide.lead));
      if (slide.points && slide.points.length) {
        var list = el("ul", "points" + (slide.points.length > 4 ? " small" : ""));
        slide.points.forEach(function (point) { list.appendChild(html("li", "", point)); });
        column.appendChild(list);
      }
      if (slide.html) column.appendChild(html("div", "", slide.html));
      if (slide.facts && slide.facts.length) column.appendChild(self._facts(slide.facts));
      if (slide.takeaway && slide.layout !== "figure") column.appendChild(html("div", "takeaway", "<span>" + slide.takeaway + "</span>"));
      return column;
    }
    function figureColumn() {
      var column = el("div", "slide-figure"), box = el("div", "figure-box");
      box.dataset.figure = "1";
      column.appendChild(box);
      if (slide.figure && slide.figure.caption) column.appendChild(html("div", "figure-caption", slide.figure.caption));
      return column;
    }
    if (slide.layout === "split") {
      body.appendChild(textColumn());
      body.appendChild(figureColumn());
    } else if (slide.layout === "figure") {
      body.appendChild(figureColumn());
      if (slide.takeaway) body.appendChild(html("div", "takeaway", "<span>" + slide.takeaway + "</span>"));
    } else if (slide.layout === "outline") {
      var outline = el("div", "outline");
      slide.items.forEach(function (item, n) {
        var card = el("div", "outline-item");
        card.appendChild(el("div", "outline-no", String(n + 1)));
        var words = el("div");
        words.appendChild(el("div", "outline-name", item.name));
        words.appendChild(el("div", "outline-what", item.what));
        card.appendChild(words);
        card.addEventListener("click", function () { self.go(item.slide); });
        outline.appendChild(card);
      });
      body.appendChild(outline);
    } else if (slide.layout === "summary") {
      var summary = el("div", "summary");
      slide.points.forEach(function (point) { summary.appendChild(html("div", "summary-item", point)); });
      body.appendChild(summary);
    } else if (slide.layout === "references") {
      body.appendChild(html("div", "references", slide.html));
    } else {
      body.appendChild(textColumn());
    }
    if (slide.sources && slide.sources.length) {
      node.appendChild(html("div", "slide-sources", "<b>Sources:</b> " + slide.sources.join(" · ")));
    }
    return node;
  };

  Deck.prototype._facts = function (facts) {
    var grid = el("div", "facts");
    facts.forEach(function (fact) {
      var card = el("div", "fact" + (fact.tone ? " " + fact.tone : ""));
      card.appendChild(html("div", "fact-value" + (String(fact.value).length > 12 ? " long" : ""), fact.value));
      card.appendChild(html("div", "fact-label", fact.label));
      if (fact.source) card.appendChild(el("div", "fact-source", fact.source));
      grid.appendChild(card);
    });
    return grid;
  };

  Deck.prototype._chrome = function () {
    var self = this, footer = el("div", "deck-footer");
    function button(label, title, handler) {
      var b = el("button", "", label);
      b.type = "button";
      b.title = title;
      b.setAttribute("aria-label", title);
      b.addEventListener("click", function () { handler(); self.root.focus(); });
      footer.appendChild(b);
      return b;
    }
    this.timer = el("span", "deck-count", "");
    this.timer.title = "Time since the first slide change";
    footer.appendChild(this.timer);
    button("☰", "Overview of all slides (O)", function () { self.toggleOverview(); });
    this.prevButton = button("‹", "Previous (left arrow)", function () { self.prev(); });
    this.count = el("span", "deck-count", "");
    footer.appendChild(this.count);
    this.nextButton = button("›", "Next (right arrow or space)", function () { self.next(); });
    button("⛶", "Full screen (F)", function () { self.fullscreen(); });
    this.stage.appendChild(footer);
    this.footer = footer;
    global.setInterval(function () {
      if (!self.started) return;
      var s = Math.floor((Date.now() - self.started) / 1000);
      self.timer.textContent = Math.floor(s / 60) + ":" + ("0" + (s % 60)).slice(-2);
    }, 1000);
  };

  Deck.prototype._overview = function () {
    var self = this, panel = el("div", "deck-overview");
    panel.appendChild(el("h2", "", this.lecture.title));
    panel.appendChild(el("div", "ov-help", "Click a slide to jump to it. Keys: arrows or space to move, F full screen, P play or pause an animation, O this overview."));
    var grid = el("div", "ov-grid");
    this.slides.forEach(function (slide, i) {
      var item = el("button", "ov-item");
      item.type = "button";
      item.appendChild(el("small", "", (i + 1) + " · " + (slide.section || "")));
      item.appendChild(document.createTextNode(String(slide.title || "").replace(/<[^>]+>/g, "")));
      item.addEventListener("click", function () { self.toggleOverview(false); self.go(i); });
      grid.appendChild(item);
    });
    panel.appendChild(grid);
    this.stage.appendChild(panel);
    this.overview = panel;
    this.overviewItems = grid.children;
  };

  Deck.prototype.toggleOverview = function (open) {
    var show = open === undefined ? !this.overview.classList.contains("open") : open;
    this.overview.classList.toggle("open", show);
    for (var i = 0; i < this.overviewItems.length; i++) this.overviewItems[i].classList.toggle("current", i === this.index);
  };

  Deck.prototype.toast = function (message) {
    var node = this.toastNode;
    node.textContent = message;
    node.classList.add("show");
    global.clearTimeout(this.toastTimer);
    this.toastTimer = global.setTimeout(function () { node.classList.remove("show"); }, 3200);
  };

  Deck.prototype.fit = function () {
    if (this.options.embedded) {
      // inside the application the frame takes the 16:9 shape of the stage (not in full screen)
      try {
        var frame = global.frameElement;
        if (frame && !global.parent.document.fullscreenElement) {
          var wanted = Math.round((frame.clientWidth * STAGE_H) / STAGE_W) + "px";
          if (frame.style.height !== wanted) {
            frame.style.height = wanted;
            var holder = frame.parentElement; // the application wraps the frame in a box of fixed height
            if (holder && holder.children.length === 1) holder.style.height = wanted;
          }
        }
      } catch (e) { /* a frame of another origin keeps its height */ }
    }
    var w = this.root.clientWidth || global.innerWidth, h = this.root.clientHeight || global.innerHeight;
    var scale = Math.min(w / STAGE_W, h / STAGE_H);
    this.stage.style.transform = "scale(" + scale + ") translate(-50%, -50%)";
    this.scale = scale;
    this.live.forEach(function (item) {
      if (item.viewer) { item.options.fontScale = scale; item.viewer.resize(); }
      if (item.draw) item.draw();
    });
  };

  Deck.prototype._mount = function (slide, node) {
    var self = this, box = node.querySelector("[data-figure]"), figure = slide.figure;
    if (!box || !figure) return;
    function molecule(target, sceneName, extra) {
      var options = Object.assign({ toolbar: false, bigLabels: true, selectable: true, fontScale: self.scale }, figure.options || {}, extra || {});
      var view = new NS.MoleculeView(target, self.lecture.scenes[sceneName], options);
      self.live.push(view);
      return view;
    }
    function inner() { // the viewers lay themselves out; give them their own box that fills the figure
      var holder = el("div");
      holder.style.cssText = "position:absolute;inset:0";
      box.appendChild(holder);
      return holder;
    }
    if (figure.type === "viewer") {
      box.classList.add("stage");
      molecule(slide.layout === "title" ? box : inner(), figure.scene);
    } else if (figure.type === "viewers") {
      var pair = el("div", "figure-pair");
      box.appendChild(pair);
      figure.items.forEach(function (item) {
        var cell = el("div");
        cell.appendChild(el("div", "pair-title", item.title));
        var inner = el("div");
        inner.style.cssText = "position:absolute;inset:0";
        cell.appendChild(inner);
        pair.appendChild(cell);
        molecule(inner, item.scene, item.options);
      });
    } else if (figure.type === "process") {
      var player = new NS.ProcessPlayer(inner(), this.lecture.scenarios[figure.scenario], { compact: true, noKeys: true });
      this.live.push(player);
      this.player = player;
    } else if (figure.type === "image") {
      var image = document.createElement("img");
      image.src = figure.src;
      image.alt = figure.alt || "";
      box.appendChild(image);
    } else if (figure.type === "html") {
      box.appendChild(html("div", "html-figure", figure.html));
    }
  };

  Deck.prototype._unmount = function (node) {
    this.live.forEach(function (item) { if (item.dispose) item.dispose(); });
    this.live = [];
    this.player = null;
    var box = node && node.querySelector("[data-figure]");
    if (box) { box.textContent = ""; box.classList.remove("stage"); }
  };

  Deck.prototype.go = function (index) {
    index = Math.min(this.slides.length - 1, Math.max(0, index));
    if (index === this.index) return;
    var old = this.nodes[this.index], next = this.nodes[index], forward = index > this.index;
    if (this.index >= 0 && !this.started) this.started = Date.now();
    if (old) { old.classList.remove("active"); old.classList.toggle("back", forward); this._unmount(old); }
    if (!next.parentNode) this.stage.insertBefore(next, this.progress);
    next.classList.toggle("back", !forward);
    void next.offsetWidth; // start the transition from the shifted position
    next.classList.add("active");
    next.classList.remove("back");
    this.index = index;
    this._mount(this.slides[index], next);
    this.count.textContent = (index + 1) + " / " + this.slides.length;
    this.progress.style.width = (100 * (index + 1)) / this.slides.length + "%";
    this.prevButton.disabled = index === 0;
    this.nextButton.disabled = index === this.slides.length - 1 && !this.player;
    this.root.classList.toggle("on-dark", this.slides[index].layout === "title");
    if (!this.options.embedded) {
      try { global.history.replaceState(null, "", "#" + (index + 1)); } catch (e) { /* file: pages may refuse */ }
    }
  };

  /* next(): one step of the slide's animation if there is one left, otherwise the next slide. */
  Deck.prototype.next = function () {
    if (this.player && this.player.next()) return;
    this.go(this.index + 1);
  };
  Deck.prototype.prev = function () {
    if (this.player && this.player.time > 0) { this.player.prev(); return; }
    this.go(this.index - 1);
  };

  Deck.prototype.fullscreen = function () {
    var self = this, target = document.documentElement, owner = document;
    try {
      if (global.frameElement) { target = global.frameElement; owner = global.parent.document; }
    } catch (e) { /* a frame of another origin: use this document */ }
    if (owner.fullscreenElement) { owner.exitFullscreen(); return; }
    var request = target.requestFullscreen ? target.requestFullscreen() : null;
    if (!request) { this.toast("Press F11 for full screen"); return; }
    if (request.catch) request.catch(function () { self.toast("Full screen was refused by the browser: press F11"); });
  };

  Deck.prototype._key = function (event) {
    if (event.ctrlKey || event.metaKey || event.altKey) return;
    if (event.target && /INPUT|SELECT|TEXTAREA/.test(event.target.tagName)) return;
    var key = event.key, handled = true;
    if (this.overview.classList.contains("open")) {
      if (key === "Escape" || key === "o" || key === "O") this.toggleOverview(false); else handled = false;
    } else if (key === "ArrowRight" || key === "PageDown" || key === " " || key === "Enter") this.next();
    else if (key === "ArrowLeft" || key === "PageUp" || key === "Backspace") this.prev();
    else if (key === "Home") this.go(0);
    else if (key === "End") this.go(this.slides.length - 1);
    else if (key === "f" || key === "F") this.fullscreen();
    else if (key === "o" || key === "O") this.toggleOverview(true);
    else if ((key === "p" || key === "P") && this.player) this.player.toggle();
    else handled = false;
    if (handled) event.preventDefault();
  };

  NS.Deck = Deck;
})(typeof window !== "undefined" ? window : globalThis);
