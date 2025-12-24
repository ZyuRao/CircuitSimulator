# CircuitSimulator

## 项目概况
CircuitSimulator 是一个简化版电路仿真器，读取类 SPICE 网表并完成电路分析。
当前支持 OP / TRAN / HB，结果输出到 `out/` 目录（CSV，按网表文件名命名）。
.DC / .AC (相关Sweep) 会被解析，但暂未实现，运行时会提示跳过。

## 功能概况

**已支持**

- 网表解析（续行、注释、SPICE 数值后缀）
- 器件：R / C / L / V / I / M(MOS)
    - MOS 为 3 端口（bulk 默认接 source）
- 源：
   - V 源：DC、SIN（可 DC + SIN）
   - I 源：DC
- 分析：
    - `.OP`（工作点）
    - `.TRAN`（瞬态）
    - `.HB`（谐波平衡，Newton）
    - `.PSS`(打靶发求解周期稳态)

- 观测/输出控制：
    - `.PRINT` / `.PROBE`
    - `.PLOTNV` / `.PLOTNC`
- 输出：`.csv` or `.png`（默认 out/）
- 已解析但暂未实现（会跳过）
    - `.DC`
    - `.AC`


## 构建与运行
**依赖**

- C++17 编译器
- Eigen3（Ubuntu/WSL 推荐安装 libeigen3-dev）
- CMake
- Python3 + pandas + matplotlib(绘图脚本)
- OpenMP（可选，仅用于实验的 `CSIM_LU_TRAIL_PAR`）

**构建**
```bash
mkdir -p build
cmake -S . -B build
cmake --build build -j
```

**运行**
```bash
./build/CircuitSimulator ./tests/<filename>.sp
```

## 输出文件说明

- 默认输出目录：./out/

`.csv`文件命名：用 netlist 基名`name`：
- out/name_op.csv
- out/name_dc.csv
- out/name_tran.csv
- out/name_hb.csv
- out/name_shoot.csv

.`PROBE`触发的图片命名：
- out/name_<analysis>_probe.png

## 运行时总耗时输出
程序结束时会输出一次总运行时间：
`[TIME] total_wall_ms=...`

## 网表语法手册（以当前 Parser 实现为准）

### 1) 注释与续行
- **整行注释**：行首（忽略前导空白）为 `*` 或 `;` 的行会被跳过。
- **行内注释**：使用 `$`，`$` 之后内容会被丢弃。
  - 例：`R1 a b 1k  $ this is a comment`
- **续行**：行首（忽略前导空白）为 `+` 表示与上一条逻辑行拼接（中间自动加空格）。
  - 常用于 `.MODEL`、`.PRINT/.PROBE` 写很长的情况。

> 建议：如果要写标题，请写成注释行（例如 `* buffer example`），不要裸写标题文本行。

---

### 2) 数值格式
- 支持科学计数法：`1e-3`
- 支持 SPICE 后缀：如 `1k`, `1u`, `3.3meg` 等（由 `parseSpiceNumber()` 决定）

---

### 3) 器件行（R/C/L/V/I/M）

#### 电阻 / 电容 / 电感
```spice
Rname n+ n- value
Cname n+ n- value
Lname n+ n- value
```

#### 独立电压源 V

支持 DC 与 SIN，允许写法如下（等价形式）：

**DC**：
```spice
Vname np nm <value>
Vname np nm DC <value>
```
**SIN**：
```spice
Vname np nm SIN v0 va freq [phi]
Vname np nm SIN v0 va freq [td phi]
```
- `phi `以 度（**degree**） 输入，内部会转为弧度
- 如果只给 1 个可选参数，按 `phi` 解释（不是 td）

**DC + SIN**（先设直流，再叠加正弦）：
```spice
Vname np nm DC <value> SIN v0 va freq [phi]
Vname np nm <value>    SIN v0 va freq [td phi]
```
#### 独立电流源 I
```spice
Iname np nm <value>
Iname np nm DC <value>
```
#### MOSFET（3 端口，bulk 默认接 source）
支持两种 token 数：

**7-token**
```spice
Mname nd ng ns modelId W L
```

**8-token**（助教/示例网表风格 buffer.sp）
```spice
Mname nd ng ns p|n W L modelId
```

`p|n` 这个 token 目前会被忽略

实际 NMOS/PMOS 由 `.MODEL `里的 `VT` 符号决定：
- $V_T < 0$ 视为 PMOS（内部会转成正值存储）
- $V_T \geq 0$ 视为 NMOS

### 4) 控制卡（.OP / .TRAN / .HB / .PSS / .DC / .AC / .MODEL / .PRINT / .PROBE / .PLOTNV / .PLOTNC）
`.OP`
```spice
.OP
```

`.TRAN`
```spice
.TRAN tstep tstop [tstart]
```
`.HB`
```spice
.HB f0 nHarm
```

`.PSS`
```spice
.PSS f0 tstep
```
`.DC`**（已解析，暂未实现运行）**
```spice
.DC sourceName start stop step
```
`.AC`**（已解析，暂未实现运行）**
```spice
.AC {LIN|DEC|OCT} nPoints fstart fstop
```
`sweepType` 支持：`LIN / OCT / DEC`（未识别时默认 `DEC`）

`.MODEL`
```spice
.MODEL modelId VT <val> MU <val> COX <val> LAMBDA <val> CJ0 <val>
```

参数以 key/value 成对读取（可用续行 `+` 拆开写）

`CJ0/CJO` 都会识别为同一个参数

### 5) 波形表达式（用于 .PRINT / .PROBE）

支持：
- 节点电压：`V(node)`
- 差分电压：`V(n1,n2)`
- 支路电流：`I(elementName)`(例如 `I(V1)`)

表达式里不要插空格（例如不要写 `V( 1 , 0 )`），因为 parser 是按空白切 token 的。

`.PRINT`
```spice
.PRINT <analysis> <expr1> <expr2> ...
```

`analysis`支持：`OP / DC / AC / TRAN / HB`

`.PROBE`
```spice
.PROBE <analysis> <expr1> <expr2> ...
```
语法与 `.PRINT` 一致（用于触发绘图）

`.PLOTNV`**（节点电压）**
```spice
.PLOTNV node1 node2 ...
```

等价于生成 `V(node1) V(node2) ... `的 probe 列表（不带 analysis）

`.PLOTNC`**（器件电流，支持端口信息字符串）**
```spice
.PLOTNC M1(d) V1(+) R1(-) I1(+)
```
`(...)`内内容会被保存为 `elePort` 字段（用于端口/方向描述）

### 完整示例网表：
```spice
* simple netlist example
V1 in 0 DC 5
R1 in out 1k
C1 out 0 1u
M1 out in 0 1 10u 0.35u
.MODEL 1 VT 0.7 MU 0.05 COX 3e-5 LAMBDA 0.02 CJ0 4e-14

.OP
.TRAN 1e-6 1e-3
.PRINT OP V(in) V(out) I(V1)
.PRINT TRAN V(out)
.PLOTNV out
```

## Core
（待补充）
