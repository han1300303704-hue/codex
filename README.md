# 射频受损与量化约束下的近场预测波束跟踪

本仓库实现开题报告中的 MATLAB 仿真框架：面向 300 GHz XL-MIMO 近场移动用户，在公共 CFO、Wiener 相位噪声和有限比特移相器约束下，完成预测波束跟踪、失锁检测和重捕获。

## 快速开始

需要 MATLAB R2019b 或更新版本；工程不依赖任何额外工具箱。

```matlab
run_experiment
```

默认设置为 256 阵元半波长 ULA、300 GHz、10 m 初始距离、20 度初始角度、`v_r=2 m/s`、`v_t=3 m/s`、50 us 跟踪周期、10 dB SNR、2 kHz CFO、100 Hz 相位噪声线宽、3-bit 移相器和 50 次蒙特卡洛试验。运行后在 `results/` 生成：

- `tracking_results.mat`：全部数值结果与配置；
- `beam_gain.fig/.png`：各对照组的归一化波束增益；
- `tracking_error.fig/.png`：位置、距离和角度 RMSE；
- `lock_loss_probability.fig/.png`：固定参数下的失锁概率；
- `robustness_scans.fig/.png`：量化位宽和相位噪声扫描。

结果目录默认不纳入 Git；如需保留某次实验，请将输出复制到其他位置。

## 算法流程

1. 通过球面波 ULA 模型生成逐阵元 LoS 信道和非均匀多普勒。
2. 先进行局部距离-角度波束扫描获得初始焦点，再用初始导频突发完成公共 CFO 粗估计，以阵列投影追踪公共相位。
3. 以 `[r, theta, v_r, v_t, CFO, CPE]` 为状态，使用信息形式 EKF 直接处理复导频的实虚观测。
4. 连续波束或 `B` 位离散相移器波束在预测区域内进行选择；联合方法以不确定性平均增益减去伪焦点守卫区响应为评分。
5. 当增益连续 3 个时隙低于 -3 dB 时，使用宽范围离散码本进行一次测量驱动的重捕获。

固定比较五个场景：理想硬件连续波束、受损未补偿、仅射频补偿、仅量化感知与联合算法。

## 参数配置

所有参数位于 `src/default_config.m`。可在命令行构建覆盖项：

```matlab
addpath('src');
cfg = default_config(struct('n_mc', 10, 'n_slots', 100));
results = run_suite(cfg);
```

量化位宽和相位噪声扫描参数位于 `cfg.scan`；重捕获阈值、局部码本及守卫区惩罚位于 `cfg.beam`。

## 测试

```matlab
addpath('tests');
run_tests
```

测试覆盖目标点聚焦峰值、相移量化级别、理想链路的可用跟踪增益、联合算法相对未补偿算法的固定种子改善，以及全部图和数据产物的生成。

## 目录

```text
run_experiment.m       主入口
src/                   信道、观测、EKF、波束选择与绘图模块
tests/run_tests.m      可执行冒烟测试
docs/model_spec.md     数学模型和性能指标
results/               本地生成的实验结果（已忽略）
```

原始开题报告仅在本地作为需求依据，不会被提交到 GitHub。
