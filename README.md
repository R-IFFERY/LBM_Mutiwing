# LBM_service

二维 LBM 翼型/多翼型仿真程序，使用 `Fortran + MPI + OpenACC + CUDA` 实现，支持两种运行模式：

- 离线模式：从 `motion_params.dat` 读取运动参数
- 在线模式：通过 ZeroMQ 从 Python 或上位机实时发送控制指令

当前代码默认参数下：

- 翼数量：2
- 物面边界条件：IBM (`bc_num=4`)
- 网格尺寸：`width=40`，`height=20`
- 四层固定嵌套网格，间距从外到内为 `0.08、0.04、0.02、0.01`
- 最外层范围：`[-20,20] × [-10,10]`，最细层范围：`[-8,8] × [-4,4]`
- 流体特征速度：`UU=0.01`
- 最细层仿真时间步：`dt = mesh_min_scale * UU = 0.0001`，运动输入与在线接口使用此时间步

四层范围、层间交换、输出及续算说明见 [NESTED_GRID.md](./NESTED_GRID.md)。

离线运动文件名不再要求包含 `k...st...`。流场在 `t=30、31、32、33、34、35` 自动各保存一次，
时间单位与运动文件第一列相同；运动数据必须覆盖所需保存时刻。当前时间步下对应
`F300000.dat、F310000.dat、F320000.dat、F330000.dat、F340000.dat、F350000.dat`。

## 1. 目录与主要文件

- `main.f90`：主程序入口、全局参数
- `mode.f90`：离线运动参数读取、在线目标构造、输出目录管理
- `driver.f90`：在线 ZeroMQ 驱动
- `zmq_comm.f90`：ZeroMQ 协议解析与回复
- `solution.f90`：主求解流程、力输出、流场输出
- `nested_grid.f90`：嵌套层推进、时空插值、MRT 分布函数转换和细到粗回传
- `nested_io.f90`：去除覆盖区的多层流场输出、四层续算文件
- `airfoil.f90`：边界读取与翼型几何更新
- `boundary_condition.f90`：浸没边界、受力和力矩计算
- `motion_params.dat`：离线模式运动输入
- `boundary.dat`、`boundary1.dat`：第 1/第 2 个翼的边界点

## 2. 编译

项目当前 `Makefile` 使用：

- 编译器：`mpif90`
- OpenACC：`-acc`
- GPU：`-gpu=deepcopy -cuda`
- CUDA 库：`cublas`、`cusolver`
- ZeroMQ：`-lzmq`

编译命令：

```bash
make
```

清理：

```bash
make clean
```

## 3. 运行方式

### 3.1 离线模式

离线模式是默认模式，不需要设置环境变量：

```bash
mpirun -np 1 ./exe
```

离线模式从 `motion_params.dat` 读取运动参数。

注意：

- 现在离线模式下，`motion_params.dat` 中的 `x/y` 位置按绝对位置解释
- 不会再额外叠加初始偏移

### 3.2 在线模式

在线模式通过环境变量启用：

```bash
export LBM_MOTION_SOURCE=ZMQ
export LBM_ZMQ_ENDPOINT=tcp://*:5555
mpirun -np 1 ./exe
```

如果想从在线模式切回离线模式，在当前终端执行：

```bash
unset LBM_MOTION_SOURCE
unset LBM_ZMQ_ENDPOINT
mpirun -np 1 ./exe
```

或仅对单次命令临时取消：

```bash
env -u LBM_MOTION_SOURCE -u LBM_ZMQ_ENDPOINT mpirun -np 1 ./exe
```

## 4. 关键参数

下面几个参数定义在 [main.f90](./main.f90)：

```fortran
parameter (wing_count=2)
parameter (mesh_max_scale=0.08d0)
parameter (wing_initial_x_origin=0.0d0)
parameter (wing_initial_x_spacing=1.5d0)
parameter (wing_initial_y_origin=0.5d0)
parameter (output_stride_x=2,output_stride_y=2)
```

含义：

- `wing_count`：翼数量
- `mesh_max_scale`：最外层网格尺寸
- `mesh_min_scale`：由层数推导的最细层网格尺寸，目前为 `0.01`
- `wing_initial_x_origin`：在线模式下第 1 个翼的基准 `x`
- `wing_initial_x_spacing`：在线模式下相邻翼的基准 `x` 间距
- `wing_initial_y_origin`：在线模式下所有翼共同的基准 `y`
- `output_stride_x/y`：流场输出抽样步长，只影响输出密度，不影响求解

说明：

- 在线模式下，控制输入中的 `sa/sb` 是相对位移，程序会在内部叠加上述基准偏移
- 离线模式下，`motion_params.dat` 里的 `sa/sb` 按绝对位置解释

## 5. 离线输入：motion_params.dat

### 5.1 格式

每一行表示一个时刻的运动参数。

当前 `wing_count=2`，每行共 `19` 列：

```text
time
wing1_angle wing1_omega wing1_omgac wing1_x wing1_xv wing1_xa wing1_y wing1_yv wing1_ya
wing2_angle wing2_omega wing2_omgac wing2_x wing2_xv wing2_xa wing2_y wing2_yv wing2_ya
```

也可以写成单行形式：

```text
t  a1 om1 aa1 x1 xv1 xa1 y1 yv1 ya1  a2 om2 aa2 x2 xv2 xa2 y2 yv2 ya2
```

### 5.2 说明

- 第 1 列：时间
- 每个翼后续 9 列依次为：
  - `angle`
  - `omega`
  - `omgac`
  - `x`
  - `x_velocity`
  - `x_acceleration`
  - `y`
  - `y_velocity`
  - `y_acceleration`

离线模式下：

- `x`、`y` 为绝对位置
- 程序内部会自动减去在线基准偏移，再进入统一几何更新流程

## 6. 在线接口

在线模式通过 ZeroMQ `REQ/REP` 通讯。

- 求解器端：`REP`
- Python/上位机端：`REQ`
- 必须严格一发一收

### 6.1 支持的命令

#### `RESET`

格式：

```text
RESET 0
```

说明：

- 当前实现要求 `RESET` 时间必须为 `0`
- 每个新 episode 都要先发一次 `RESET`

#### `STEP`

当前 `wing_count=2` 时格式：

```text
STEP <target_time> <w1_domega> <w1_dxv> <w1_dyv> <w2_domega> <w2_dxv> <w2_dyv>
```

说明：

- 在线模式传的是速度增量，不是完整状态
- 每个翼 3 个输入量依次是：
  - `domega`
  - `dx_velocity`
  - `dy_velocity`
- 程序会在当前段起点和目标段终点之间构造平滑过渡

时间约束：

- `target_time` 必须是 `0.0001` 的整数倍
- `target_time` 必须不小于当前待推进时刻

例如：

```text
STEP 0.1 0 0 0 0 0 0
```

表示推进到 `t=0.1`，两翼速度增量全为 0。

#### `SAVE`

格式：

```text
SAVE
```

或：

```text
OUTPUT
```

说明：

- `SAVE` 保存当前已经完成计算的最新流场
- 不会中断正在执行的 `STEP`
- `RESET` 后不能立刻保存，必须先完成至少一次有效推进，并且不能处于移网格过渡阶段

#### `STOP`

格式：

```text
STOP
```

### 6.2 回复格式

#### 力回复

`RESET` 和 `STEP` 的回复格式：

```text
FORCE <len> <v1> <v2> ... <vN>
```

当前 `wing_count=2` 时，返回数组顺序为：

```text
1.  time
2.  wing1_fx_avg
3.  wing1_fy_avg
4.  wing1_m_avg
5.  wing2_fx_avg
6.  wing2_fy_avg
7.  wing2_m_avg
8.  total_fx_avg
9.  total_fy_avg
10. total_m_avg
11. wing1_x_last
12. wing1_xv_last
13. wing1_xa_last
14. wing1_y_last
15. wing1_yv_last
16. wing1_ya_last
17. wing1_angle_last
18. wing1_omega_last
19. wing1_omgac_last
20. wing2_x_last
21. wing2_xv_last
22. wing2_xa_last
23. wing2_y_last
24. wing2_yv_last
25. wing2_ya_last
26. wing2_angle_last
27. wing2_omega_last
28. wing2_omgac_last
```

说明：

- 力和力矩列是“最近一次参数更新后的整段平均值”
- 运动参数列是“该段最后一步的值”

#### 状态回复

保存流场成功：

```text
STATUS FLOW_SAVED F1001.dat
```

停止：

```text
STATUS STOPPING
```

错误：

```text
ERROR <message>
```

## 7. 输出文件

所有输出都写到新建的运行目录：

```text
output/run_YYYYMMDD_HHMMSS_mmm/
```

最新一次运行目录会记录在：

```text
output/latest_run.txt
```

### 7.1 离线模式输出

- `force.dat`
- `airfoil_motion.dat`
- 固定时刻触发的流场文件：`F<n>.dat`

### 7.2 在线模式输出

- `force_episode_<episode>.dat`
- `airfoil_motion_episode_<episode>.dat`

如果某个 episode 调用了 `SAVE`，程序会自动创建：

```text
output/run_.../episode_<episode>/
```

然后把该回合的流场文件写入其中，文件名保持原样：

```text
F<n>.dat
```

这样不同 episode 即使保存了相同步号的流场，也不会互相覆盖。

## 8. 输出内容说明

### 8.1 force.dat / force_episode_*.dat

每步写一行，当前列顺序为：

```text
time
[每个翼: Fx Fy M]
[总: Fx Fy M]
[每个翼: x xv xa y yv ya angle omega omgac]
```

这里的位置列写的是实际位置。

### 8.2 airfoil_motion.dat / airfoil_motion_episode_*.dat

这是当前代码保留的辅助输出文件。

当前每个翼只写了 3 列：

```text
time
[每个翼: 0.0  body_extra_v/UU  x_abs]
```

它不是完整的二维运动状态文件，不能替代 `force.dat` 中的运动参数列。

### 8.3 F<n>.dat

Tecplot 风格流场输出，四层组成 13 个矩形流场 zone，粗层被细层覆盖的面积不重复输出，包含：

- `X`
- `Y`
- `U`
- `V`
- `Pressure`
- `Vorticity`
- `body`
- `dvordt`

并在文件末尾追加翼型边界几何。

## 9. Python/上位机接入建议

在线模式推荐交互顺序：

```text
RESET 0
STEP 0.1 ...
如需保存：SAVE
STEP 0.2 ...
如需保存：SAVE
STOP
```

推荐做法：

- 上位机按控制周期发送 `STEP`
- 每次收到 `FORCE` 后再决定是否发送 `SAVE`
- `SAVE` 保存的是该时刻最新完成的流场

## 10. 常见问题

### 10.1 为什么我切不回离线模式？

如果你执行过：

```bash
export LBM_MOTION_SOURCE=ZMQ
export LBM_ZMQ_ENDPOINT=tcp://*:5555
```

那么当前终端会一直保持在线模式，直到你执行：

```bash
unset LBM_MOTION_SOURCE
unset LBM_ZMQ_ENDPOINT
```

### 10.2 为什么 `STEP 0.1` 不是只算一步？

因为内部时间步是：

```text
dt = 0.0001
```

所以 `STEP 0.1` 会连续推进 1000 个小步，然后返回这一整段的平均力结果。

### 10.3 为什么 `SAVE` 有时会失败？

当前代码在移网格后的过渡阶段禁止保存流场。如果收到：

```text
ERROR 当前时刻处于移网格后的过渡阶段，请在下一次 STEP 后再保存流场
```

继续推进一段后再保存即可。

## 11. 当前默认物理/几何设定

当前代码中的默认在线基准位置参数为：

- 第 1 个翼：`x = 0.0`
- 第 2 个翼：`x = 1.5`
- 所有翼：`y = 0.5`

这些参数只作为在线模式的基准位置使用；离线模式以 `motion_params.dat` 中给出的绝对位置为准。
