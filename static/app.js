"use strict";

let STATE = null;

const $ = (id) => document.getElementById(id);
const NS = "http://www.w3.org/2000/svg";

function el(tag, attrs, text) {
  const node = document.createElementNS(NS, tag);
  for (const k in attrs) node.setAttribute(k, attrs[k]);
  if (text !== undefined) node.textContent = text;
  return node;
}

function toast(msg, bad) {
  const t = $("toast");
  t.textContent = msg;
  t.style.borderColor = bad ? "#7a3030" : "#355c31";
  t.style.background = bad ? "#2a1818" : "#1d2a1a";
  t.style.color = bad ? "#ff8585" : "#5dd39e";
  t.style.display = "block";
  clearTimeout(toast._t);
  toast._t = setTimeout(() => (t.style.display = "none"), 3200);
}

async function api(path, opts) {
  const res = await fetch(path, opts);
  if (!res.ok) {
    const txt = await res.text();
    throw new Error(res.status + " " + txt);
  }
  return res.json();
}

function scale(domain, range) {
  const [d0, d1] = domain, [r0, r1] = range;
  const pad = (d1 - d0) * 0.04 || 1;
  const lo = d0 - pad, hi = d1 + pad;
  return {
    lo, hi,
    f: (v) => r0 + ((v - lo) / (hi - lo || 1)) * (r1 - r0),
  };
}

function logScale(domain, range) {
  const l0 = Math.log10(domain[0]), l1 = Math.log10(domain[1]);
  return { lo: domain[0], hi: domain[1],
    f: (v) => range[0] + ((Math.log10(v) - l0) / (l1 - l0 || 1)) * (range[1] - range[0]) };
}

function axes(svg, xs, ys, xLabel, yLabel, logx, logy) {
  const W = svg.clientWidth || 800, H = svg.clientHeight || 180;
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  svg.innerHTML = "";
  const ml = 52, mr = 12, mt = 8, mb = 22;
  const iw = W - ml - mr, ih = H - mt - mb;
  svg.appendChild(el("rect", { x: ml, y: mt, width: iw, height: ih,
    fill: "none", stroke: "#283240" }));
  const xt = logx ? niceLogTicks(xs.lo, xs.hi) : niceTicks(xs.lo, xs.hi);
  xt.forEach((v) => {
    const x = xs.f(v);
    svg.appendChild(el("line", { x1: x, y1: mt, x2: x, y2: mt + ih,
      stroke: "#1c2530" }));
    svg.appendChild(el("text", { x, y: mt + ih + 14, "font-size": 9,
      fill: "#7d8b9c", "text-anchor": "middle" }, fmt(v)));
  });
  const yt = logy ? niceLogTicks(ys.lo, ys.hi) : niceTicks(ys.lo, ys.hi);
  yt.forEach((v) => {
    const y = ys.f(v);
    svg.appendChild(el("line", { x1: ml, y1: y, x2: ml + iw, y2: y,
      stroke: "#1c2530" }));
    svg.appendChild(el("text", { x: ml - 4, y: y + 3, "font-size": 9,
      fill: "#7d8b9c", "text-anchor": "end" }, fmt(v)));
  });
  svg.appendChild(el("text", { x: ml + iw / 2, y: H - 2, "font-size": 10,
    fill: "#93a1b3", "text-anchor": "middle" }, xLabel));
  svg.appendChild(el("text", { x: 12, y: mt + ih / 2, "font-size": 10,
    fill: "#93a1b3", transform: `rotate(-90 12 ${mt + ih / 2})`,
    "text-anchor": "middle" }, yLabel));
  return { ml, mt, iw, ih };
}

function fmt(v) {
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e4 || a < 1e-3) return v.toExponential(1);
  return Math.abs(v - Math.round(v)) < 1e-9 ? String(Math.round(v)) : v.toPrecision(3);
}
function niceTicks(lo, hi) {
  const n = 5, span = hi - lo || 1, step0 = span / n;
  const mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const norm = step0 / mag;
  const step = (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step)
    out.push(+v.toFixed(10));
  return out;
}
function niceLogTicks(lo, hi) {
  const out = [];
  for (let p = Math.floor(Math.log10(lo)); p <= Math.ceil(Math.log10(hi)); p++)
    out.push(Math.pow(10, p));
  return out;
}

function polyline(svg, pts, xs, ys, box, color, width) {
  if (pts.length < 2) return;
  const d = pts.map((p, i) =>
    `${i ? "L" : "M"}${xs.f(p[0]).toFixed(1)},${ys.f(p[1]).toFixed(1)}`).join("");
  svg.appendChild(el("path", { d, fill: "none", stroke: color,
    "stroke-width": width || 1.6 }));
}

const CAND_COLORS = { storage: "#b18cff", radial: "#4da3ff", boundary: "#ff8f6b" };

function drawRateChart() {
  const svg = $("rateChart");
  const steps = STATE.steps;
  const ts = STATE.samples.map((s) => s.time);
  const qs = steps.map((s) => s.rate);
  const xs = scale([Math.min(...ts), Math.max(...ts)], [52, (svg.clientWidth || 800) - 12]);
  const ys = scale([0, Math.max(...qs) * 1.1], [(svg.clientHeight || 200) - 30, 8]);
  const box = axes(svg, xs, ys, "实际时间 t (h)", "流量 (m³/d)");
  // right-continuous step line
  let d = "";
  steps.forEach((s, i) => {
    const x = xs.f(s.time), y = ys.f(s.rate);
    d += (i === 0 ? `M${x},${ys.f(0)}` : `L${x},${ys.f(steps[i - 1].rate)}`) + ` L${x},${y}`;
  });
  d += `L${xs.f(Math.max(...ts))},${ys.f(steps[steps.length - 1].rate)}`;
  svg.appendChild(el("path", { d, fill: "none", stroke: "#5dd39e", "stroke-width": 2 }));
  // markers for steps landing exactly on a sample
  const sampleTimes = new Set(ts.map((t) => +t.toFixed(9)));
  steps.forEach((s) => {
    if (s.time > 0 && sampleTimes.has(+s.time.toFixed(9)))
      svg.appendChild(el("circle", { cx: xs.f(s.time), cy: ys.f(s.rate),
        r: 3.5, fill: "#5dd39e", stroke: "#0f141b" }));
  });
}

function drawPressChart() {
  const svg = $("pressChart");
  const data = STATE.samples;
  const xs = scale([data[0].time, data[data.length - 1].time],
    [52, (svg.clientWidth || 800) - 12]);
  const ps = data.map((s) => s.pressure_corrected);
  const ys = scale([Math.min(...ps), Math.max(...ps)],
    [(svg.clientHeight || 200) - 30, 8]);
  axes(svg, xs, ys, "实际时间 t (h)", "压力 (kPa，换档校正)");
  const seg0 = data.filter((s) => s.segment === 0)
    .map((s) => [s.time, s.pressure_corrected]);
  const seg1 = data.filter((s) => s.segment >= 1)
    .map((s) => [s.time, s.pressure_corrected]);
  polyline(svg, seg0, xs, ys, null, "#4da3ff", 1.6);
  polyline(svg, seg1, xs, ys, null, "#4da3ff", 1.6);
  // raw jump at shifts shown in red
  STATE.shifts.forEach((sh) => {
    svg.appendChild(el("line", { x1: xs.f(sh.time), y1: 8,
      x2: xs.f(sh.time), y2: (svg.clientHeight || 200) - 30,
      stroke: "#ff6b6b", "stroke-dasharray": "4 3", opacity: 0.7 }));
  });
  data.forEach((s) => {
    if (s.duplicate)
      svg.appendChild(el("circle", { cx: xs.f(s.time),
        cy: ys.f(s.pressure_raw), r: 2.6, fill: "#ffb454" }));
  });
}

function drawLogChart() {
  const svg = $("logChart");
  const W = svg.clientWidth || 800, H = svg.clientHeight || 320;
  const data = STATE.samples.filter((s) => s.teq > 0 && s.delta_p > 0);
  const xAll = data.map((s) => s.teq);
  const yDp = data.map((s) => s.delta_p);
  const yDer = data.filter((s) => s.derivative > 0).map((s) => s.derivative);
  const xs = logScale([Math.min(...xAll) * 0.9, Math.max(...xAll) * 1.1], [52, W - 12]);
  const ys = logScale(
    [Math.min(...yDp, ...yDer) * 0.8, Math.max(...yDp, ...yDer) * 1.2],
    [H - 30, 8]);
  axes(svg, xs, ys, "等效时间 t_eq (h，对数)", "Δp / 对数导数 (kPa，对数)", true, true);
  const dpPts = data.map((s) => [s.teq, s.delta_p]);
  // split derivative by segment so a gauge jump is never connected
  STATE.shifts.forEach((sh, si) => {});
  polyline(svg, dpPts, xs, ys, null, "#4da3ff", 1.5);
  const groups = {};
  data.forEach((s) => {
    if (s.derivative > 0)
      (groups[s.segment] ||= []).push([s.teq, s.derivative]);
  });
  Object.values(groups).forEach((pts) =>
    polyline(svg, pts, xs, ys, null, "#ffb454", 2));
  data.forEach((s) => {
    if (s.derivative > 0)
      svg.appendChild(el("circle", { cx: xs.f(s.teq), cy: ys.f(s.derivative),
        r: 2.2, fill: "#ffb454" }));
  });
  // candidate interval vertical spans
  Object.entries(STATE.candidates || {}).forEach(([reg, ab]) => {
    const t0 = teqAt(ab[0]), t1 = teqAt(ab[1]);
    if (t0 && t1)
      svg.appendChild(el("rect", { x: xs.f(t0), y: 8,
        width: xs.f(t1) - xs.f(t0), height: H - 38,
        fill: CAND_COLORS[reg] || "#888", opacity: 0.08 }));
  });
}

function teqAt(t) {
  const hit = STATE.samples.find((s) => Math.abs(s.time - t) < 1e-9 && s.teq > 0);
  return hit ? hit.teq : null;
}

function renderDiag() {
  const d = STATE.diagnostics;
  const chips = [
    [`流量修订版本 r${STATE.rate_version}`, "ok"],
    [`重复时刻 ${d.duplicates}`, d.duplicates ? "warn" : ""],
    [`阶跃对齐样本 ${d.step_aligned}（右连续）`, "ok"],
    [`非正等效时间 ${d.nonpositive_teq}（保留诊断）`, "warn"],
    [`换档分段 ${d.shift_segments}`, "bad"],
    [`换档边界单侧导数 ${d.one_sided}`, "warn"],
  ];
  $("diag").innerHTML = chips
    .map(([t, c]) => `<span class="chip ${c}">${t}</span>`).join("");
  const shiftTxt = STATE.shifts
    .map((s) => `@${s.time}h=${(+s.offset_kpa).toFixed(2)}kPa`).join("，");
  $("meta").textContent =
    `样本 ${STATE.samples.length} 条 · 流量台阶 ${STATE.steps.length} 个 · `
    + `换档 ${shiftTxt || "无"} · 平滑 L=${STATE.smoothing}`;
}

function renderSchemes() {
  const box = $("schemeBox");
  if (!STATE.schemes || !STATE.schemes.length) {
    box.innerHTML = '<p class="muted">尚无方案。新建方案会自动加入存储/径向/边界候选。</p>';
    return;
  }
  box.innerHTML = STATE.schemes.map((s) => {
    const rows = s.intervals.map((iv) => {
      const f = iv.data.fit || {};
      const detail = [
        f.permeability_md != null ? `k=${(+f.permeability_md).toFixed(0)}mD` : "",
        f.storage_si != null ? `C=${(+f.storage_si).toExponential(1)}m³/Pa` : "",
        f.doubling_ratio != null ? `翻倍×${(+f.doubling_ratio).toFixed(2)}` : "",
        f.boundary_distance_m != null ? `L≈${(+f.boundary_distance_m).toFixed(0)}m` : "",
      ].filter(Boolean).join("，");
      const o = iv.data;
      return `<tr><td>${regimeName(o.regime)}</td>
        <td>${o.start_open ? "(" : "["}${o.start_time}–${o.end_time}${o.end_open ? ")" : "]"}</td>
        <td>${f.n_points ?? "-"}</td><td>${detail}</td></tr>`;
    }).join("");
    return `<div style="margin-bottom:8px">
      <b>#${s.id} ${s.name}</b> <span class="muted">L=${s.smoothing}</span>
      <table><tr><th>流态</th><th>开闭端点</th><th>有效样本</th><th>导出参数</th></tr>${rows}</table>
      </div>`;
  }).join("");
}

function regimeName(r) {
  return { storage: "井筒储集", radial: "径向流", boundary: "边界效应" }[r] || r;
}

async function renderRuns() {
  const data = await api("/api/runs");
  $("runBox").innerHTML = data.runs.length
    ? data.runs.map((r) =>
        `<div>运行 #${r.id} <b>${r.label}</b>
         <span class="muted">r${r.rate_version}</span>
         <code>${r.fingerprint}</code></div>`).join("")
    : '<p class="muted">尚无冻结运行。</p>';
}

async function refresh(smoothing) {
  STATE = await api("/api/state?smoothing=" + (smoothing ?? 0));
  drawRateChart();
  drawPressChart();
  drawLogChart();
  renderDiag();
  renderSchemes();
  renderRuns();
}

$("rateBtn").onclick = async () => {
  try {
    const r = await api("/api/rates", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        time: +$("rateTime").value, rate: +$("rateVal").value }),
    });
    toast("流量修订已提交，版本 r" + r.rate_version);
    refresh(+$("schemeL").value || 0);
  } catch (e) { toast(String(e), true); }
};

$("shiftBtn").onclick = async () => {
  try {
    const offRaw = $("shiftOff").value.trim();
    const body = { time: +$("shiftTime").value };
    if (offRaw !== "") body.offset_kpa = +offRaw;
    const r = await api("/api/shifts", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body) });
    toast("换档已标记，分段版本 " + r.shift_version);
    refresh(+$("schemeL").value || 0);
  } catch (e) { toast(String(e), true); }
};

$("schemeBtn").onclick = async () => {
  const L = +$("schemeL").value || 0;
  try {
    const sc = await api("/api/schemes", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ name: $("schemeName").value, smoothing: L }) });
    const cands = STATE.candidates || {};
    const lockRaw = $("lockedPlateau").value.trim();
    for (const [regime, ab] of Object.entries(cands)) {
      const payload = {
        regime, start_time: ab[0], end_time: ab[1], smoothing: L };
      if (regime === "radial" && lockRaw !== "")
        payload.locked_plateau = +lockRaw;
      await api(`/api/schemes/${sc.scheme_id}/intervals?smoothing=${L}`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload) });
    }
    toast("方案已创建并加入 " + Object.keys(cands).length + " 个候选区段");
    refresh(L);
  } catch (e) { toast(String(e), true); }
};

$("runBtn").onclick = async () => {
  try {
    const r = await api("/api/runs?smoothing=" + (+$("schemeL").value || 0), {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ label: $("runLabel").value }) });
    toast("运行已冻结，指纹 " + r.fingerprint + "（含 r" + r.rate_version + "）");
    renderRuns();
  } catch (e) { toast(String(e), true); }
};

$("exportBtn").onclick = async () => {
  const blob = await api("/api/export");
  const url = URL.createObjectURL(new Blob([JSON.stringify(blob, null, 2)],
    { type: "application/json" }));
  const a = document.createElement("a");
  a.href = url; a.download = "drawdown-export.json"; a.click();
  URL.revokeObjectURL(url);
  toast("已导出含指纹的运行记录 JSON");
};

$("resetBtn").onclick = async () => {
  if (!confirm("清空数据库并重新导入固定 fixture？")) return;
  await api("/api/reset", { method: "POST" });
  toast("已清空并重新导入 fixture");
  refresh(0);
};

window.addEventListener("resize", () => {
  if (STATE) { drawRateChart(); drawPressChart(); drawLogChart(); }
});

refresh(0).catch((e) => toast(String(e), true));
