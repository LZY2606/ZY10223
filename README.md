# 压降译井台（Pressure Drawdown Interpretation Bench）

变流量压降试井的本地解释工作台：在**实际、非等间隔时间坐标**上计算多流量等效时间、
压力变化与 Bourdet 对数导数；保留早期井筒储集、径向流、边界效应三类候选区段；
支持修正单个流量阶段、锁定平台斜率、标注仪器换档，多个区段方案可并排且**不覆盖原始读数**。

技术栈：Python + FastAPI + NumPy + SQLite，前端为零依赖的原生 HTML/JS/SVG。

## 安装

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## 演示

```bash
.venv/bin/python -m pytest -q
.venv/bin/uvicorn app:app --host 127.0.0.1 --port 5563
```

浏览器打开 <http://127.0.0.1:5563>，页面标题为 **压降译井台**。
首次进入点“导入固定 fixture”（或调用 `POST /api/import`）。

## 数据口径

- **单位**：时间 h，流量 m³/h，压力 kPa，内部一致的抽象油藏单位；μ=1、B=1、h=10 m。
- **流量历史右连续**：流量在 `start` 时刻立即切换；当压力样本恰好落在阶段边界上时，
  该样本取**阶跃后**流量，且叠加时间存在对数奇异点，故其等效时间记 `teq=0`
  并打标 `rate_boundary_sample` + `teq_nonpositive`，不参与导数。
- **多流量等效时间（标准变流量叠加时间）**：在第 n 阶段内
  ```
  teq(t) = exp( (1/qn) * Σ_{k≤n} (qk − q_{k−1}) · ln(t − t_{k−1}) )
  ```
- **压力变化**：`Δp = p_ref − p_corrected`，压降为正；`p_ref` 取开井基线压力。
- **对数导数**：`dΔp / d(ln teq)`，用非等间隔 Bourdet 三点法在**实际样本的等效时间
  坐标**上直接计算（不做等间隔重采样）。同时保留一份 `dΔp/d(ln t)` 的实际时间导数。
  - 平滑尺度：默认左右各取 **2 个有效邻居**（`smooth_neighbors`），可用 ln 半窗 `smooth_l`
    额外收紧；平滑与导数都**不会跨越仪器换档段**，段边缘退化为单侧并标记 `one_sided`。
  - **重复时刻**与 **teq≤0** 的样本保留在逐样本表中（打标），但不参与导数。
- **仪器换档**：换档点（`index`）及之后样本归入新段；`offset` 为观测值相对地层真值的跳变，
  缺省时用换档点前后局部线性外推估计。`p_corrected = p_raw − 累计 offset`。
  导数严格分段，跳变**不会被平滑器解释成地层响应**；`pressure_raw` 原样保留。
- **候选区段**：自动给出早期井筒储集（der/Δp≈1 的单元斜率）、径向流（导数平台）、
  边界效应（末段导数单调上抬）三类候选；流量阶跃后的瞬态通过“同换档段+同流量块”隔离，
  不会被误判为平台或边界。
- **区段记录**：开闭端点（左/右 open）、有效样本数、平滑尺度、成员样本，以及导出参数：
  - 井筒储集：`C = (dΔp/d ln teq)/q`；
  - 径向流：导数平台、log-log 斜率、`k = 162.6·q·μ·B/(h·平台)`，可锁定平台斜率；
  - 边界：log-log 斜率、单位线截距与数量级边界半径估计。
- **运行指纹**：SHA-256 截断，混入分析器版本、**流量修订版本**、修订后阶段、换档标注、
  平滑尺度、参考压力、区段定义与原始数据摘要；任何一项变化指纹都会变。

## 操作

- 页面上方四图联动：流量历史（右连续阶梯）、压力（原始 vs 校正、换档标记）、
  log-log 诊断图（Δp、导数、候选区段高亮），下方为逐样本诊断表。
- 左侧可编辑任意流量阶段并“保存为新版本”；原始读数永不被覆盖。
- 可按 index/offset 标记仪器换档（offset 留空自动估计）。
- “新建方案”用当前流量版本与换档计算；区段端点、开闭、平台斜率锁定、平滑 L 均可改，
  改动后指纹实时重算。多个方案通过顶部下拉并排切换、互不覆盖。

## 运行记录导出与复核

- `GET /api/runs/{id}/export` 导出包含原始读数、全部流量修订、换档、区段与指纹的 JSON
  （页面“导出运行”）。
- `POST /api/runs/replay` 用导出包**重新计算并核对指纹**，返回 `matches:true/false`
  （页面“重放并核对”）。
- **清空后复核**：`POST /api/reset` 清空全部表，再 `POST /api/import` 从固定 fixture
  重新导入；原始读数来自仓库内 `fixtures/well_test_fixture.json`（确定性，
  可由 `python scripts/make_fixture.py` 重新生成并逐字节比对）。

## 固定 fixture

`fixtures/well_test_fixture.json` 含 44 个样本，刻意植入验收要素：

- 开井前短读（q=0，teq≤0）；
- 两个流量阶跃 `t=10`（0→15）、`t=600`（15→25）**恰落在样本时刻**，且各带一条短读重复；
- 一次仪器换档：`index=37`（t=1360）观测值上跳 18 kPa。

## 目录

```
analysis.py     核心数学：等效时间、换档分段、Bourdet 导数、候选检测、导出参数、指纹
db.py           SQLite 持久层（原始读数 / 流量修订 / 运行方案）
app.py          FastAPI 路由
static/index.html 操作页面
scripts/make_fixture.py 确定性 fixture 生成器
fixtures/       固定验收数据
tests/          pytest 自动化测试（核心数学 + 服务 + 回放）
```

## 主要 API

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/api/import` | 清空后导入固定 fixture |
| POST | `/api/reset` | 清空数据库 |
| GET | `/api/raw` | 原始读数、单位、当前流量、修订列表 |
| POST | `/api/revisions` | 保存一个新的流量修订版本 |
| POST | `/api/runs` | 新建解释方案（含换档/平滑/自动候选） |
| GET | `/api/runs/{id}` | 方案详情（样本、区段、导出参数） |
| PUT | `/api/runs/{id}/intervals` | 更新区段（端点/开闭/锁定/平滑） |
| GET | `/api/runs/{id}/export` | 导出运行记录 |
| POST | `/api/runs/replay` | 重放导出包并核对指纹 |
