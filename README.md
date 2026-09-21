# 压降译井台（Drawdown Interpretation Platform）

变流量条件下的压降试井解释本地服务：在**实际时间**坐标上计算等效时间、
压力变化 Δp 与对数导数，保留早期井筒储集、径向流、边界效应三类区段候选，
支持流量阶段修订、平台斜率锁定与仪器换档标记。多个解释方案并排保存，
**原始压力读数永不被覆盖**。

## 安装

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## 演示 / 验收

```bash
.venv/bin/python -m pytest -q
.venv/bin/uvicorn app:app --host 127.0.0.1 --port 5563
```

浏览器访问 <http://127.0.0.1:5563>，页面标题为 **压降译井台**。

数据库默认落在 `data/drawdown.db`；可用环境变量 `DRAWDOWN_DB` 指向其他
路径（测试套件即用临时库，互不干扰）。

## 数据口径（单位与符号约定）

| 量 | 单位 | 说明 |
|----|------|------|
| 时间 `t` | h | 自测试开始 t=0 的实际经过时间，非等间隔 |
| 流量 `q` | m³/d | 地面产液，生产为正；台阶**右连续** |
| 压力 | kPa | 仪器原始绝对读数；Δp = p_initial − p_corrected |

固定 fixture（`data/fixture.json`，由 `welltest/fixture.py` 确定性生成，
无随机数）的真值参数：渗透率 1000 mD、厚度 20 m、黏度 1 mPa·s、孔隙度
0.2、综合压缩系数 4e-10 1/Pa、井半径 0.1 m、井筒储集无因次 C_D=5、
表皮 0、封闭断层距离 0.8 m（镜像井）。

fixture 刻意包含验收要素：

- **短读重复时刻**两处：`t=0.02 h` 与 `t=150 h`（后者读数偏移 −0.07 kPa）；
- **一次仪器换档跳变 +8.000 kPa**，起点恰好是 `t=48 h` 的压力样本；
- **两个流量阶跃恰落在压力样本上**：`10 h: 0→60 m³/d`、
  `200 h: 60→200 m³/d`；
- 另有一个**不落样本**的阶跃 `2400 h: 200→80 m³/d`。

## 计算口径

### 右连续流量台阶

`rate_at(steps, t)` 规定时刻恰好等于样本时间的阶跃在该样本上**已经生效**。
阶跃瞬间的几何叠加时间是奇异的，该样本标记 `rate_step_aligned`，
`t_eq=0`，不参与导数，但保留在诊断中。

### 变流量等效时间（几何叠加）

```
t_eq = exp( Σ_k (Δq_k / q_N) · ln(t − t_k) ),   t_k < t
```

它满足 `d ln t_eq / dt = S(t)/q_N`，因此径向流平台在不同流量阶段保持恒定，
可在同一双对数图上叠加解释。以下样本**保留诊断但不参与导数拟合**：

- 当前流量非正（`nonpositive_rate`）；
- 阶跃恰好落在样本上（`rate_step_aligned`）；
- 重复时刻（`duplicate_timestamp`，仅最后一条读数有效）；
- 几何叠加为 0 或非正（`non_positive_teq`，骤降流量后可能出现）。

### 实际时间坐标上的对数导数（Bourdet）

- 邻居搜索在**实际时间对数轴** `|ln Δt|` 上进行（非等间隔、非均匀网格），
  割线则取 `x = ln t_eq`；
- 平滑尺度 `L`（默认 0）限制邻居时间比不超过 `e^L`；
- **邻居绝不跨越仪器换档分段**：换档样本及其后第一点为
  `one_sided_segment_boundary`，跳变不可能被三点平滑“抹成”地层响应；
- 当 `ln t_eq` 在实际时间序上非单调（常见于骤降流量后），标记
  `non_monotonic_teq`，不强行桥接触线。

### 区段与导出参数

每个区间记录开闭端点（`[a,b]` / `(a,b)` 等）、有效样本数、平滑尺度与
导出参数：

- **井筒储集**：以累积产量 `Np = ∫q dt` 对 Δp 过原点拟合 `C = Np/Δp`
  （SI: m³/Pa），并给出 log-log 斜率；
- **径向流**：由（可锁定的）平台斜率反算渗透率
  `k = μBq/(4πh·m′)`（mD），中点处给表皮；
- **边界效应**：相对同速率径向平台的翻倍比，及封闭断层近似距离。

`suggest_candidates()` 在“同一换档段 + 同一恒定流量”的连续分组内识别
候选，避免把阶跃重启或换档跳变连成一个区段。

## 页面与接口

页面三图联动：流量历史（右连续台阶）、压力历史（换档分段着色）、
Δp/对数导数双对数诊断（候选区段半透明带）。侧栏可：

- 修正任一流量阶段（生成新的 `rate_version`，原始读数不变）；
- 标记/标定仪器换档（跳变可留空自动估计，或手工给 kPa）；
- 新建并排方案并自动加入三类候选区间（含开闭端点、有效样本、平滑尺度）；
- 冻结运行指纹、导出 JSON、清空并重导 fixture 复核。

主要接口：

| 方法 | 路径 | 作用 |
|------|------|------|
| GET  | `/api/state?smoothing=L` | 全量计算结果与候选 |
| POST | `/api/rates` | 流量阶段修订（版本 +1） |
| POST | `/api/shifts` | 标记/标定仪器换档 |
| POST/GET | `/api/schemes` | 新建/列出并排方案 |
| POST | `/api/schemes/{id}/intervals` | 加区间并返回拟合参数 |
| POST/GET | `/api/runs` | 冻结/列出运行（指纹含流量修订版本） |
| GET  | `/api/runs/{id}` | 取某次运行快照重放 |
| GET  | `/api/export` | 导出含指纹的完整数据包 |
| POST | `/api/import` | 校验指纹并导入（不通过返回 409） |
| POST | `/api/reset` | 清空数据库并重导固定 fixture |

## 运行指纹与清空后复核

指纹对“原始读数 + 折叠后的流量台阶 + 换档偏移 + **rate_version**”做
规范化 SHA-256（取前 16 位）。导出包自带指纹；`/api/import` 会先重算指纹
比对，篡改任何读数或流量都会被拒绝（HTTP 409）。

清空数据库后复核：

```bash
# 页面上点“清空并重导 fixture”，或：
curl -s -X POST http://127.0.0.1:5563/api/reset
# 导出 -> 删除库文件 -> 重启 -> /api/import 回导，指纹应一致
```

## 目录

```
app.py                 # FastAPI 入口与序列化
welltest/engine.py     # 右连续台阶、等效时间、Bourdet 导数、区段拟合
welltest/fixture.py    # 确定性合成 fixture（重复/换档/阶跃/断层）
welltest/store.py      # SQLite：原始读数、修订版本、方案、运行指纹
static/index.html      # 操作页面
static/app.js          # SVG 三图与接口联动
data/fixture.json      # 固化 fixture
tests/                 # 引擎不变量与端到端接口测试（22 项）
```
