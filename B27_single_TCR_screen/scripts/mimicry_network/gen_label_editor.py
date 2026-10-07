# -*- coding: utf-8 -*-
"""Generate the label position drag editor (gen_label_editor.py -> label_editor.html).
- Same data source as the formal figures ({tcr}_labelpos_d3.json written by
  plot_mimicry_red.py, including two-line labels)
- Drag a label to move it (cx, cy); double-click a label to flip its
  alignment edge (ha); other interactions match the proofreading tool
- "Export override JSON" only exports labels that were moved ->
  labelpos_override_d3.json; place it in data/mimicry_network/ and rerun
  plot_mimicry_red.py to render the formal figures with your positions
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))           # -> B27_single_TCR_screen
DATA_DIR = os.path.join(REPO, "data", "mimicry_network")
OUT_DIR = os.path.join(REPO, "figures", "mimicry_network")

DATA = {}
for lo in ("tcr59", "tcr63"):
    with open(os.path.join(OUT_DIR, f"{lo}_labelpos_d3.json"), encoding="utf-8") as f:
        DATA[lo.upper()] = json.load(f)

# merge the user's saved override positions and flag them moved=1
# (preserved on re-export)
OV_PATH = os.path.join(DATA_DIR, "labelpos_override_d3.json")
if os.path.exists(OV_PATH):
    ov = json.load(open(OV_PATH, encoding="utf-8"))
    n_merged = 0
    for t, groups in ov.items():
        if t not in DATA:
            continue
        for key, entries in groups.items():
            by_i = {e["i"]: e for e in entries}
            for b in DATA[t].get(key, []):
                if b["i"] in by_i:
                    b.update(cx=by_i[b["i"]]["cx"], cy=by_i[b["i"]]["cy"],
                             ha=by_i[b["i"]]["ha"], moved=1)
                    n_merged += 1
    print(f"override merged: {n_merged} labels from {OV_PATH}")

HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>label position editor</title>
<style>
 body{margin:0;font-family:Arial,sans-serif;background:#f4f5f7}
 #hd{padding:10px 14px;background:#fff;border-bottom:1px solid #ddd;display:flex;
     flex-wrap:wrap;gap:12px;align-items:center}
 #hd b{font-size:16px;color:#333}
 .tbtn{padding:5px 14px;border:1px solid #bbb;border-radius:5px;background:#fff;
       cursor:pointer;font-size:13px}
 .tbtn.on{background:#1a73e8;color:#fff;border-color:#1a73e8}
 .exp{background:#2e7d32;color:#fff;border:none;border-radius:5px;padding:6px 16px;
      font-size:13px;cursor:pointer}
 #stat{font-size:13px;color:#444}
 #wrap{position:relative}
 canvas{display:block;background:#fff;cursor:default}
 #tip{position:absolute;display:none;background:#fffde7;border:1px solid #d9c94a;
      border-radius:4px;padding:8px 10px;font-size:12px;line-height:1.55;color:#222;
      pointer-events:none;box-shadow:2px 2px 6px rgba(0,0,0,.15);max-width:320px;
      z-index:9;white-space:nowrap}
 #tip .tt{font-weight:bold;font-size:13px}
 #tip .rv{color:#c62828;font-weight:bold}
</style>
</head>
<body>
<div id="hd">
 <b>label editor</b>
 <button class="tbtn on" id="b59">TCR59</button>
 <button class="tbtn" id="b63">TCR63</button>
 <button class="tbtn on" id="lvG">With gene (two lines)</button>
 <button class="tbtn" id="lvN">Sequence only</button>
 <button class="exp" id="bExp">Export override JSON</button>
 <span id="stat"></span>
</div>
<div id="wrap"><canvas id="cv"></canvas><div id="tip"></div></div>
<script>
const DATA = __DATA__;
let TCR = "TCR59", LV = "labels_gene";
let view = {s: 1, tx: 0, ty: 0, W: 0, H: 0};
let hover = null, pin = null, dragPan = null, dragLab = null;
const cv = document.getElementById("cv"), ctx = cv.getContext("2d");
const tip = document.getElementById("tip");
const RED = "#d32f2f", LABC = "#b71c1c", LABS = "#1565c0", LINS = "#7d94b5";
const REDLW = {1: 1.6, 2: 1.3, 3: 1.1};
const WEAK = {1: ["#b0b0b0", 0.8], 2: ["#bdbdbd", 0.55],
              3: ["#cccccc", 0.30], 4: ["#d8d8d8", 0.20]};

function fit() {
  const geo = DATA[TCR].geo;
  const W = innerWidth - 4, H = innerHeight - document.getElementById("hd").offsetHeight - 4;
  cv.width = W * devicePixelRatio; cv.height = H * devicePixelRatio;
  cv.style.width = W + "px"; cv.style.height = H + "px";
  const m = 70;
  view.s = Math.min((W - 2 * m) / (2 * geo.mx), (H - 2 * m) / (2 * geo.my));
  view.tx = W / 2; view.ty = H / 2;
  view.W = W; view.H = H;
}
const px = x => view.tx + x * view.s;
const py = y => view.ty - y * view.s;
const kpt = () => view.s / DATA[TCR].geo.S;

function labs() { return DATA[TCR][LV]; }

function labelBoxes() {
  const k = kpt();
  ctx.font = (6 * k) + "px Arial";
  const lh = 6 * 1.12 * k, out = [];
  for (let L = 0; L < labs().length; L++) {
    const b = labs()[L];
    const w = Math.max(...b.lines.map(t => ctx.measureText(t).width));
    const h = b.lines.length * lh;
    const ax = px(b.cx), ay = py(b.cy);
    const x0 = b.ha ? ax : ax - w;
    out.push({L, x0, y0: ay - h / 2, x1: x0 + w, y1: ay + h / 2, b});
  }
  return out;
}

function partnersOf(n) {
  const out = [];
  for (const e of DATA[TCR].edges)
    if (e[3] && (e[0] === n || e[1] === n))
      out.push([e[0] === n ? e[1] : e[0], e[2]]);
  return out;
}

function draw() {
  const d = DATA[TCR], ns = d.nodes, es = d.edges;
  const k = kpt();
  ctx.setTransform(devicePixelRatio, 0, 0, devicePixelRatio, 0, 0);
  ctx.clearRect(0, 0, view.W, view.H);
  const rel = pin ? pin : hover;
  let relNode = -1, relLab = -1;
  if (rel) {
    relNode = rel.type === "label" ? labs()[rel.idx].i : rel.idx;
    relLab = rel.type === "label" ? rel.idx : labs().findIndex(l => l.i === relNode);
  }
  const pset = relNode >= 0 ? new Set(partnersOf(relNode).map(p => p[0])) : null;

  for (const [i, j, dd, red] of es) {
    if (red) continue;
    const w = WEAK[dd] || ["#ddd", 0.15];
    ctx.strokeStyle = w[0]; ctx.lineWidth = Math.max(w[1] * k, 0.3);
    ctx.globalAlpha = (relNode >= 0 && i !== relNode && j !== relNode) ? 0.10 : 1;
    ctx.beginPath(); ctx.moveTo(px(ns[i].x), py(ns[i].y));
    ctx.lineTo(px(ns[j].x), py(ns[j].y)); ctx.stroke();
  }
  for (const [i, j, dd, red] of es) {
    if (!red) continue;
    ctx.strokeStyle = RED; ctx.lineWidth = REDLW[dd] * k;
    ctx.globalAlpha = (relNode >= 0 && i !== relNode && j !== relNode) ? 0.3 : 1;
    ctx.beginPath(); ctx.moveTo(px(ns[i].x), py(ns[i].y));
    ctx.lineTo(px(ns[j].x), py(ns[j].y)); ctx.stroke();
  }
  ctx.globalAlpha = 1;

  for (let i = 0; i < ns.length; i++) {
    const n = ns[i], X = px(n.x), Y = py(n.y), r = Math.max(n.r * k, 1.5);
    ctx.globalAlpha = relNode < 0 ? 1 :
      (i === relNode ? 1 : (pset.has(i) ? 0.75 : 0.18));
    ctx.beginPath();
    if (n.tri) { ctx.moveTo(X, Y - r); ctx.lineTo(X - r * 0.93, Y + r * 0.62);
                 ctx.lineTo(X + r * 0.93, Y + r * 0.62); ctx.closePath(); }
    else ctx.arc(X, Y, r, 0, 6.2832);
    ctx.fillStyle = n.fill; ctx.fill();
    ctx.strokeStyle = n.stroke; ctx.lineWidth = Math.max(n.slw * k, 0.4);
    ctx.stroke();
  }
  ctx.globalAlpha = 1;

  // always-on label leader lines (matching the formal figure: #a85757, dash 3-2)
  const boxes = labelBoxes();
  for (const bx of boxes) {
    const b = bx.b, n = ns[b.i];
    const X = px(n.x), Y = py(n.y);
    const cx = Math.max(bx.x0, Math.min(X, bx.x1));
    const cy = (Y < bx.y0) ? bx.y0 : (Y > bx.y1 ? bx.y1 : py(b.cy));
    ctx.strokeStyle = n.o === "virus" ? "#a85757" : LINS;
    ctx.lineWidth = Math.max(0.7 * k, 0.4);
    ctx.setLineDash([3 * k, 2 * k]);
    ctx.globalAlpha = relLab < 0 ? 0.9 : (bx.L === relLab ? 1 : 0.15);
    ctx.beginPath(); ctx.moveTo(X, Y); ctx.lineTo(cx, cy); ctx.stroke();
  }
  ctx.setLineDash([]);
  ctx.globalAlpha = 1;

  const lh = 6 * 1.12 * k;
  for (const bx of boxes) {
    const b = bx.b;
    ctx.globalAlpha = relLab < 0 ? 1 : (bx.L === relLab ? 1 : 0.18);
    const cy = py(b.cy), n = b.lines.length;
    for (let q = 0; q < n; q++) {
      const ty = cy + (q - (n - 1) / 2) * lh + lh * 0.34;
      ctx.lineWidth = 1.6 * k; ctx.strokeStyle = "rgba(255,255,255,0.95)";
      ctx.strokeText(b.lines[q], bx.x0, ty);
      ctx.fillStyle = ns[b.i].o === "virus" ? LABC : LABS;
      ctx.fillText(b.lines[q], bx.x0, ty);
    }
    if (b.moved) { ctx.strokeStyle = "#2e7d32"; ctx.lineWidth = 1.2;
      ctx.setLineDash([3, 2]);
      ctx.strokeRect(bx.x0 - 3, bx.y0 - 2, bx.x1 - bx.x0 + 6, bx.y1 - bx.y0 + 4);
      ctx.setLineDash([]); }
  }
  ctx.globalAlpha = 1;

  if (relNode >= 0) {
    const n = ns[relNode];
    const X = px(n.x), Y = py(n.y), r = Math.max(n.r * k, 1.5);
    ctx.setLineDash([3, 2.5]); ctx.strokeStyle = RED; ctx.lineWidth = 1.4;
    ctx.beginPath(); ctx.arc(X, Y, r + 4, 0, 6.2832); ctx.stroke();
    ctx.setLineDash([]);
  }
  const moved = {TCR59: 0, TCR63: 0};
  for (const t of ["TCR59", "TCR63"])
    for (const key of ["labels_gene", "labels_nogene"])
      moved[t] += DATA[t][key].filter(l => l.moved).length;
  document.getElementById("stat").innerHTML =
    (dragLab ? "dragging... " : "") + "moved labels: TCR59 " + moved.TCR59 +
    " + TCR63 " + moved.TCR63 +
    " | drag label = move, double-click label = flip alignment, " +
    "drag empty space = pan, wheel = zoom";
}

function pick(mx, my) {
  const boxes = labelBoxes();
  for (const bx of boxes)
    if (mx >= bx.x0 - 4 && mx <= bx.x1 + 4 && my >= bx.y0 - 4 && my <= bx.y1 + 4)
      return {type: "label", idx: bx.L};
  const ns = DATA[TCR].nodes;
  let best = -1, bd = 11;
  for (let i = 0; i < ns.length; i++) {
    const dd = Math.hypot(px(ns[i].x) - mx, py(ns[i].y) - my);
    if (dd < bd) { bd = dd; best = i; }
  }
  return best >= 0 ? {type: "node", idx: best} : null;
}

cv.addEventListener("mousedown", e => {
  const r = cv.getBoundingClientRect();
  const mx = e.clientX - r.left, my = e.clientY - r.top;
  const p = pick(mx, my);
  if (p && p.type === "label") { dragLab = {L: p.idx, ox: mx, oy: my}; pin = p; }
  else { dragPan = {x: mx, y: my}; pin = p; }
  draw();
});
cv.addEventListener("mousemove", e => {
  const r = cv.getBoundingClientRect();
  const mx = e.clientX - r.left, my = e.clientY - r.top;
  if (dragLab) {
    const b = labs()[dragLab.L];
    b.cx = round5((mx - view.tx) / view.s);
    b.cy = round5(-(my - view.ty) / view.s);
    if (!b.moved) b.moved = 1;
    draw(); return;
  }
  if (dragPan) { view.tx += mx - dragPan.x; view.ty += my - dragPan.y;
                 dragPan = {x: mx, y: my}; draw(); return; }
  const p = pick(mx, my);
  if (JSON.stringify(p) !== JSON.stringify(hover)) { hover = p; draw(); }
  if (p) {
    const node = p.type === "label" ? labs()[p.idx].i : p.idx;
    const n = DATA[TCR].nodes[node];
    const ps = partnersOf(node);
    let h = "<div class='tt'>" + n.aa + "</div>gene: " + n.g +
      "<br>origin: " + (n.o === "virus" ? "virus" : "self (human)") +
      "<br>cluster: " + n.cl + (n.mim ? " | mimicry participant" : "");
    if (p.type === "label")
      h += "<br><span class='rv'>&#9654; label points to this node</span>" +
           " (double-click label to flip alignment)";
    if (ps.length) {
      h += "<br><span class='rv'>mimicry partners:</span>";
      for (const [j, dd] of ps) h += "<br>&nbsp;&nbsp;" +
        DATA[TCR].nodes[j].aa + " (" + DATA[TCR].nodes[j].g + ", " +
        DATA[TCR].nodes[j].o + ") d=" + dd;
    }
    tip.innerHTML = h; tip.style.display = "block";
    tip.style.left = Math.min(mx + 16, view.W - 330) + "px";
    tip.style.top = Math.max(my - 10, 4) + "px";
  } else tip.style.display = "none";
});
function round5(v) { return Math.round(v * 1e5) / 1e5; }
addEventListener("mouseup", () => { dragPan = null; dragLab = null; });
cv.addEventListener("dblclick", e => {
  const r = cv.getBoundingClientRect();
  const p = pick(e.clientX - r.left, e.clientY - r.top);
  if (p && p.type === "label") {
    const b = labs()[p.idx];
    b.ha = b.ha ? 0 : 1; b.moved = 1; draw();
  }
});
cv.addEventListener("wheel", e => {
  e.preventDefault();
  const r = cv.getBoundingClientRect();
  const mx = e.clientX - r.left, my = e.clientY - r.top;
  const f = e.deltaY < 0 ? 1.12 : 0.89;
  view.tx = mx - (mx - view.tx) * f; view.ty = my - (my - view.ty) * f;
  view.s *= f; draw();
}, {passive: false});
cv.addEventListener("mouseleave", () => { hover = null; tip.style.display = "none"; draw(); });

document.getElementById("bExp").addEventListener("click", () => {
  const out = {};
  for (const t of ["TCR59", "TCR63"]) {
    out[t] = {};
    for (const key of ["labels_gene", "labels_nogene"])
      out[t][key] = DATA[t][key].filter(l => l.moved)
        .map(l => ({i: l.i, cx: l.cx, cy: l.cy, ha: l.ha}));
  }
  const txt = JSON.stringify(out, null, 1);
  const blob = new Blob([txt], {type: "application/json"});
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = "labelpos_override_d3.json";
  a.click(); URL.revokeObjectURL(a.href);
  let tot = 0;
  for (const t of Object.keys(out))
    for (const key of Object.keys(out[t])) tot += out[t][key].length;
  alert("Exported " + tot + " moved label positions.\n" +
    "Place labelpos_override_d3.json in data/mimicry_network/, then rerun " +
    "plot_mimicry_red.py to render the formal figures with your positions.");
});
function setTCR(t) {
  TCR = t; pin = null; hover = null;
  document.getElementById("b59").classList.toggle("on", t === "TCR59");
  document.getElementById("b63").classList.toggle("on", t === "TCR63");
  fit(); draw();
}
function setLV(v) {
  LV = v; pin = null; hover = null;
  document.getElementById("lvG").classList.toggle("on", v === "labels_gene");
  document.getElementById("lvN").classList.toggle("on", v === "labels_nogene");
  draw();
}
document.getElementById("b59").addEventListener("click", () => setTCR("TCR59"));
document.getElementById("b63").addEventListener("click", () => setTCR("TCR63"));
document.getElementById("lvG").addEventListener("click", () => setLV("labels_gene"));
document.getElementById("lvN").addEventListener("click", () => setLV("labels_nogene"));
addEventListener("resize", () => { fit(); draw(); });
fit(); draw();
</script>
</body>
</html>
"""

html = HTML.replace("__DATA__", json.dumps(DATA, ensure_ascii=False, separators=(",", ":")))
out = os.path.join(OUT_DIR, "label_editor.html")
with open(out, "w", encoding="utf-8") as f:
    f.write(html)
print("written:", out, os.path.getsize(out) // 1024, "KB")
for t, d in DATA.items():
    g2 = sum(1 for b in d["labels_gene"] if len(b["lines"]) == 2)
    vg = sum(1 for b in d["labels_gene"]
             if len(b["lines"]) == 2 and "_" in b["lines"][1])
    print(f"  {t}: gene labels {len(d['labels_gene'])} (two-line {g2}, "
          f"virus gene_species {vg}), nogene {len(d['labels_nogene'])}")
