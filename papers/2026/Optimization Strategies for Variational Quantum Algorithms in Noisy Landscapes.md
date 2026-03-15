# Optimization Strategies for Variational Quantum Algorithms in Noisy Landscapes

**読了日:** 2026/03/15
**学会・年:** arXiv 2025年
**引用キー:** novak2025optimization

**keywords:** Variational Quantum Algorithms, VQE, Metaheuristic Optimization, CMA-ES, iL-SHADE, Noisy Landscapes, Barren Plateaus, NISQ

## 論文情報

題名は **Optimization Strategies for Variational Quantum Algorithms in Noisy Landscapes**。著者は **Vojtěch Novák, Ivan Zelinka, Václav Snášel**。arXiv: **2506.01715v2 [quant-ph] 27 Sep 2025**。

## Abstract（要旨）

本論文は、**Variational Quantum Algorithms (VQAs)**、とくに **Variational Quantum Eigensolver (VQE)** における最適化問題を、**noise**, **barren plateaus**, **complex energy landscapes** という観点から扱う。著者らは **50種類超の metaheuristic algorithms** を3段階でベンチマークした。第1段階は **Ising model** 上での初期スクリーニング、第2段階は **9 qubits** までのスケーリング試験、第3段階は **192-parameter Hubbard model** における収束性評価である。ランドスケープ可視化により、**noiseless** では滑らかでほぼ凸な basin が、**finite-shot sampling** 下では歪み、rugged になり、これが gradient-based local methods の失敗理由を説明すると示した。総合的に **CMA-ES** と **iL-SHADE** が最良で、**Simulated Annealing (Cauchy)**、**Harmony Search**、**Symbiotic Organisms Search** も頑健だった。一方で **PSO**, **GA**, 標準的な **DE variants** はノイズ下で大きく悪化した。結論として、**noisy VQE** に対して有望な少数の optimizer 群を特定し、近未来量子デバイス向けの指針を与える。

## I. INTRODUCTION（導入）

量子計算は、古典計算では困難な問題、たとえば量子系シミュレーション、大規模最適化、機械学習に対して大きな可能性を持つ。近未来量子計算で有望なのが **VQAs** であり、量子化学、凝縮系、組合せ最適化などを、量子回路と古典最適化のハイブリッドで扱う。VQAs は **NISQ devices** の代表的パラダイムである。VQE はその代表例で、分子系の基底状態エネルギー近似に使われる。また VQAs は量子相転移、parameterized circuits を使う量子機械学習、**QAOA**、quantum metrology にも応用される。

しかし VQAs の最適化ランドスケープは厳しい。最大の問題は **barren plateau phenomenon** で、qubit 数の増加とともに勾配が指数的に消失し、最適化が事実上困難になることだ。さらに **local minima** と、測定ショット由来の **sampling noise** が訓練を複雑化する。測定ノイズはショット数 **N** に対して **1/√N** でスケールし、誤り訂正があっても精度限界を作る。barren plateau の小さな勾配は統計揺らぎに埋もれるため、勾配法には指数的なショット数が必要になり、現実的でなくなる。

このため著者らは、**Differential Evolution (DE)**, **Particle Swarm Optimization (PSO)**, **Covariance Matrix Adaptation Evolution Strategy (CMA-ES)** のような、局所勾配に依存しにくい **metaheuristic algorithms** に注目する。論文では、1D Ising model ではノイズなしなら滑らかでほぼ凸の basin だが、sampling noise を入れると偽の極小や勾配消失が現れること、Fermi-Hubbard model ではさらに rugged, multimodal, non-convex な surface となり、local descent が破綻しやすいことを示す。IEEE CEC competitions の multimodal・rotated・shifted・compositionally complex functions と VQA の確率的・高次元ランドスケープの類似も指摘し、そこで強い **iL-SHADE** や **CMA-ES** が noisy VQE でも強いのは自然だと位置づける。最後に論文構成を説明する。

## II. VARIATIONAL QUANTUM ALGORITHMS

VQAs は quantum hardware と classical optimization を組み合わせる手法で、qubit 数と回路深さが制約された **NISQ** を前提に設計されている。大きく **problem-driven VQAs**（VQE, QAOA など）と **data-driven QML models** に分けられるが、本論文は前者、とくに parameterized quantum circuits (PQCs) を使う VQE に集中する。PQCs が trial states を生成し、その expectation values を測定し、古典 optimizer がパラメータを更新する。

### A. The Barren Plateau Phenomenon

n-qubit 系の初期状態を **ρ**、parameterized quantum circuit を **U(θ)** とすると、進化後状態は
**ρ(θ) = U(θ)ρU†(θ)**。
測定演算子 **O** に対する loss は
**ℓθ(ρ, O) = Tr[ρ(θ)O]** で定義される。実際には有限ショット **N** で推定されるため、推定量 **ℓ̂θ(ρ, O)** には **1/√N** オーダーの分散がある。最適化の目的は **ℓθ(ρ, O)** を最小化する **θ** を見つけること。

barren plateau は、系サイズ増加に伴って loss またはその勾配が平均の周りに指数的に集中する現象である。勾配成分の分散は
**Varθ[∇θμℓθ(ρ, O)] ∈ O(1/bⁿ)**, **b > 1**。
勾配信号の減衰が統計ノイズ **1/√N** より速いため、下降方向を識別するには指数的ショット数が必要になる。

著者らは barren plateau を2種類に分ける。1つは **Probabilistic concentration with narrow gorges**。この場合、
**Varθ[ℓθ(ρ, O)] ∈ O(1/bⁿ), b > 1**。
Chebyshev により
**Pr(|ℓθ(ρ, O) − Eθ[ℓθ(ρ, O)]| ≥ δ) ∈ O(1/bⁿ)**。
つまり大半のパラメータでは平坦だが、指数的に細い高勾配領域 **narrow gorges** が残る。見つければ最適化経路になるが、見つける確率は指数的に下がる。もう1つは **Deterministic concentration** で、
**|ℓθ(ρ, O) − ℓ₀| ∈ O(1/bⁿ), ∀θ**。
これは landscape 全体が定数付近に押しつぶされる最も深刻な場合で、谷も峡谷も残らず、勾配法はほぼ不可能になる。

原因として、Hilbert 空間の **curse of dimensionality**、過度に表現力の高い circuit によるランダム探索、depolarizing noise による maximally mixed state への駆動が挙げられる。初期状態や observable の選択も影響する。barren plateau は local minima とは違い、真の地形構造ではなく統計的集中に由来する。回路起因の **intrinsic** と、noise 起因の **extrinsic** がある。こうした事情が、正確な勾配を必要としない population-based metaheuristics を使う動機になる。

## III. BENCHMARK MODELS

metaheuristic algorithm の有効性は、multimodality、noise sensitivity、barren plateaus の有無といった landscape complexity に依存する。そこで著者らは、比較可能でありながら難しさを持つ2つの benchmark models を選ぶ。VQAs は local minima, barren plateaus, stochastic measurements に苦しむため、これらは population-based metaheuristics の優位性が期待される条件でもある。

### A. The 1D Ising Model

一次元 transverse-field Ising model（ただし **external magnetic field なし**）を主要ベンチマークに採用する。Hamiltonian は
**H = − Σᵢ₌₁ⁿ⁻¹ σ_z^(i) σ_z^(i+1)**  … (6)
で、**σ_z^(i)** は qubit i 上の Pauli-Z。基底状態は「全スピン上向き」または「全スピン下向き」の **2-fold degenerate ground states** を持ち、symmetric double-well potential 的な構造になる。これは global minima を見つけられるかを試すのに向いている。

qubit 数が増えると Hilbert 空間が指数的に増大し、local minima 問題も悪化する。Ising の multimodal landscape は、evolutionary / swarm-based algorithms と局所法の違いを試すのにちょうどよい。著者らは **3 〜 9 qubits**、すなわち **12 〜 36 parameters** に系サイズを変え、dimensionality に対する性能劣化を調べる。

### B. The Fermi-Hubbard Model

より厳しい benchmark として **Fermi-Hubbard model** を使う。6-site Hubbard Hamiltonian は
**H = −t Σ⟨i,j⟩,s (c†ᵢ,scⱼ,s + c†ⱼ,scᵢ,s) + U Σᵢ nᵢ,↑ nᵢ,↓**  … (7)
で、ここでは **t = U = 1**。**Jordan-Wigner transformation** により fermionic operators を **12-qubit system** に写し、ansatz parameterization 後は **192 parameters** になる。

Hubbard model は strongly-correlated electrons の重要な物理を含む。Ising の対称的 landscape と違い、パラメータ数が増えると energy error が plateau しやすく、高度に frustrated な optimization surface ができ、散在する local minima が早期停止を引き起こす。VQE 自体は 12 sites まで基底エネルギー計算に成功しているが、192-parameter 空間の最適化はかなりきつい。著者らは exact diagonalization を基準に absolute optimizer accuracy を評価する。つまり、Ising が「制御しやすい multimodal test」、Hubbard が「現実的で rugged な correlated test」という二層構成になっている。

## IV. METHODOLOGY

### A. Cost Function Evaluation via Quantum Expectation Values

VQE は parameterized quantum state **|ψ(θ)⟩** に対し、target Hamiltonian **Ĥ** の energy expectation value
**E(θ) = ⟨ψ(θ)|Ĥ|ψ(θ)⟩**  … (8)
を最小化する。variational principle により **E(θ) ≥ E₀** であり、真の基底エネルギーの上界になる。量子実装では Hamiltonian を可測な Pauli operators の線形結合
**Ĥ = Σ_k w_k P̂_k**, **P̂_k = ⊗ⱼ σ_kj**  … (9)
に分解し、期待値は
**⟨Ĥ⟩ = Σ_k w_k⟨P̂_k⟩**  … (10)
と書ける。非ゼロ Pauli term 数が qubit 数に対し高々多項式であることが重要で、Ising と Hubbard は局所相互作用構造ゆえにこの条件を満たす。

### B. Parametrized Quantum Circuits and Ansatz Design

本論文の主 ansatz は **TwoLocal structure**。single-qubit rotation layers と two-qubit entanglement layers を交互に積む：
**|ψ(θ)⟩ = [∏ᵣ₌₁ᴿ U_ent U_rot(θ^(r))] U_rot(θ^(0)) |0⟩^⊗n**  … (11)
ここで **U_rot(θ^(0))** は初期の single-qubit rotations（例：RY）、各 repetition では **controlled-Z gates** からなる **U_ent** の後に **RY, RZ** の parameterized rotation を置く。図1では 3 qubits・1 repetition の回路が示される。

両モデルとも、全 qubit に **RY(π/4)** を入れて equal superposition state から開始する。1D Ising では **Z₂ symmetry** を尊重し、Hubbard では balanced な初期点を与える。entanglement は linear connectivity で、隣接 qubit に controlled-Z をかける。図1の説明では、各 qubit に **R(π/4, π/2)** の初期化回転、その後に **RY**, **RZ** 層と linear topology の controlled-Z 層が続き、barrier が視認性のために入れられている。

ただし著者らは、結論が ansatz 依存である可能性を明示する。**Hamiltonian Variational Ansatz (HVA/VHA)** は Hamiltonian の locality と symmetry を反映して reachable state manifold を制限し、trainability 改善や barren plateau 緩和が報告されている。noisy setting の最近の VHA 研究では、sampling noise が optimizer ranking を変え、noiseless では gradient methods が強い一方、finite-shot noise 下では **CMA-ES** のような population-based optimizers がより頑健で、physics-informed initialization が function-evaluation budget を削る。さらに chemistry-inspired ansatz として **UCCSD** は Trotterization と operator ordering に敏感であり、「Trotterized UCCSD」は一意に定まるものではない。**GUCCSD** を含む **SA-OO-VQE** では、複数の near-degenerate states を扱い、excited states や conical intersections 近傍の quasi-diabatic representation に有利である。したがって TwoLocal の結果は baseline と読むべきで、将来は HVA/VHA や UCCSD/GUCCSD 系でも再評価すべきだと述べる。

### C. Sampling Noise and Statistical Uncertainty Analysis

古典最適化と違い、量子測定は本質的に統計的不確かさを持つ。各 Pauli operator の測定で、**N shots** に対する推定値は
**P̃_k = (1/N) Σ_i=1^N m_i**  … (12)
で、**m_i ∈ {−1, +1}**。各結果の分散は
**Var[m_i] = 1 − ⟨P̂_k⟩²**  … (13)。
中心極限定理より、大きな **N** では
**P̃_k ∼ N(⟨P̂_k⟩, (1 − ⟨P̂_k⟩²)/N)**,
**Var[P̃_k] = (1 − ⟨P̂_k⟩²)/N ≤ 1/N**  … (14)。
全エネルギー推定分散は
**Var[Ẽ] = Σ_k w_k² (1 − ⟨P̂_k⟩²)/N_k**  … (15)。
独立・等ショット数 **N_k ≡ N** を仮定すると
**Var[Ẽ] ≤ (1/N) Σ_k w_k²**  … (16)。

Chebyshev により
**P(|Ẽ − E| ≥ ϵ) ≤ Var[Ẽ]/ϵ² ≤ (Σ_k w_k²)/(Nϵ²)**  … (17)。
さらに Hoeffding / Bernstein により
**P(|Ẽ − E| ≥ ϵ) ≤ 2 exp(−2Nϵ² / Σ_k w_k²)**  … (18)。
したがって、精度 **ϵ** と信頼度 **1 − δ** を達成するのに必要なショット数は
**N ≥ Σ_k w_k² /(δϵ²)**  … (19)。
固定総ショット数 **N_tot** のもとで分散最小の shot allocation は
**N_k ∝ |w_k|√(1 − ⟨P̂_k⟩²)**、
**min Var[Ẽ] = (Σ_k |w_k|√(1 − ⟨P̂_k⟩²))² / N_tot**  … (20)。
このあたり、数式は静かに見えてかなり牙をむいている。

sampling noise は最適化問題の性質を根本的に変える。optimizer は deterministic gradients ではなく stochastic cost landscapes を相手にする。effective noise floor は
**σ_noise ≈ √(Σ_k w_k²) / √N**。
これより小さい energy difference は識別困難で、進捗が隠れる。対象系では noise floor は **10⁻²（low shots）〜 10⁻³（high shots）** 程度。gradient estimation も難しく、parameter shift rule は unbiased だが circuit evaluations を増やしてショットコストを膨らませる。

本研究は **finite-shot estimation が誘起する stochastic objective functions に対する optimizer robustness** の検証として位置づけられる。shot noise は hardware 改善後も残る基本的要素であり、error mitigation により coherent noise が effectively stochastic に圧縮される場合にも対応しやすい。一方で、著者らは限界も認める。**coherent and correlated errors** は landscape をバイアスし、noise-induced barren plateaus や global minimizer の切り替わりすら起こしうる。**randomized compiling**, **ZNE**, data-driven regressors などで partially stochastic regime を回復できても、コストや sample budget は変わりうる。よって本研究の結果は「**statistical fluctuations に対する robustness の下限評価**」と解釈すべきだとする。将来は randomized compiling の on/off や injected coherent over-rotations を入れて rank change を比較すること、grouping や low-rank factorization のような measurement-cost reduction と optimizer stochasticity の相互作用を調べることを提案している。

## V. EXPERIMENTS DESIGN

著者らは多数の metaheuristics を、**bio-based**, **evolutionary-based**, **human-based**, **math-based**, **music-based**, **physics-based**, **swarm-based**, **system-based** に分類して評価した。Table I には具体的なアルゴリズム名が列挙される。まず gradient-based optimizers でも preliminary test を行ったが、**10⁻¹** という比較的緩い tolerance でも global minimum に届かないことが多く、成功率が低かったため、metaheuristics の詳細評価に進んだ。

### Table I に載る optimizer 群

Bio-Based: **BBOA, SMA, BBO, BMO, EOA, IWO, SBO, SOA, SOS, TPO, TSA, VCS, WHO, ABC**。
Evolutionary-Based: **CMA-ES, CRO, EP, DE, GA, FPA, MA, SHADE, HyDE**。
Human-Based: **BRO, IBSO, CA, CHIO, FBIO, GSKA, HBO, HCO, ICA, LCO, QSA, SARO, SPBO, SSDO, TLO, TOA**。
Math-Based: **AOA, CEM, CGO, CircleSA, GBO, HC, PSS, RUN, SCA, SHIO, TS**。
Physics-Based: **ASO, ArchOA, CDO, EO, EVO, FLA, HGSO, MVO, NRO, RIME, TWO, WDO, SA**。
Swarm-Based: **PSO, iSOMA, WOA, ALO, EHO, HHO**。
表そのものが「そんなに試したの？」という気持ちにさせるが、本当に多い。

実験環境は **Intel(R) Core(TM) i5-8400 CPU @ 2.80 GHz**, **16 GB RAM**。ただし比較は runtime ではなく **function evaluations (FEs)** ベースなので、ハード自体は本質ではない。明示がない限り estimator は **5120 shots**。FEs は各 optimizer の **5 runs の平均**。コードは **[https://github.com/VojtechNovak/VQA_metaheuristics](https://github.com/VojtechNovak/VQA_metaheuristics)** で公開され、実装は **Qiskit version 1** ベース。多くの optimizer は **mealpy**、CMA-ES は **cma** module、gradient-based optimizers は **Qiskit + scipy**、iL-SHADE は **pyade** 由来。

### 三段階評価

**Phase 1**: 5-qubit Ising, **20 parameters**、5 independent runs。**10⁻¹** 精度で global minimum に少なくとも1回到達した optimizer だけ次へ。
**Phase 2**: qubit 数 **3〜9** の Ising に対し、global minimum に **10⁻¹** 精度で到達するまでの **FEs** を比較。どれかの run で到達できなければ表では **"—"**。
**Phase 3**: 前段で上位だった optimizer を、**192 parameters** の Hubbard model で比較。**64 shots** と **5120 shots** の2条件で convergence curves を評価。

## VI. OPTIMIZATION LANDSCAPE VISUALIZATION AND ANALYSIS

著者らは optimizer の挙動予測のため、landscape の2次元切り出し可視化を行う。方法は、既に得られた最適パラメータのうち大部分を固定し、2つの selected parameters だけを2D grid 上で動かすもの。Ising では **δ = 10** の global view と **δ = 0.5** の local view を見る。

Ising では ansatz と spin couplings により明確な periodicity が現れる。**noiseless statevector** では最小近傍が smooth で approximately convex なので、gradient-based local optimizers が比較的うまく動く。だが **finite-shot sampling** を入れると、contour は distorted になり、勾配が消え、偽の local minima が出る。ショット数が少ないほどこの歪みは強く、COBYLA や SPSA の効率低下を説明する。一方、**CMA-ES**, **iL-SHADE**, **PSO**, **iSOMA** などの population-based / swarm methods は、滑らかな勾配に頼らず broad sampling で頑健性を保つ。

Hubbard model は別物で、landscape が highly correlated かつ nonconvex、local minima が irregular に散在する。可視化では残り190個を固定して2つだけ動かすため、slice 上のエネルギーは global optimum 付近の **−17** 前後に留まり、見た目の flatness が full space の ruggedness を隠す面はあるが、それでも irregular local structure は見える。local methods は trap されやすく、global metaheuristics の方が shallow minima から脱出しやすい。要するに、**periodicity と convexity は noiseless では局所法を助けるが、noise が入るとその利点は急速に消える**、ということ。

### 図2・図3・図4の読み取り

**Fig. 2** と **Fig. 3** は Ising 6-qubit における **θ₀ vs θ₁**, **θ₁ vs θ₂** の landscape。noiseless・δ=10 では繰り返し構造が見え、δ=0.5 ではほぼ同心円状の convex basin が見える。**64 shots** では local view の輪郭がギザギザに崩れ、局所構造が人工的に増えている。**5120 shots** では再び滑らかさがある程度戻る。
**Fig. 4** は Hubbard 6-site, 192 parameters の **θ₀ vs θ₁**, **θ₁ vs θ₂**。noiseless でも既に不規則だが、64 shots ではさらにざらつきが増し、5120 shots でも Ising のような滑らかな basin にはならない。これは optimizer 選好が model と shot regime で変わる理由の視覚的説明になっている。

## VII. RESULTS

### A. Local Optimizers and Limitations in Noisy Environments

**SPSA** と **COBYLA** は VQE でよく使われる local optimizers だが、noiseless では Sec. VI のとおり smooth convex basin に乗って収束しやすい一方、sampling noise が入ると contour の歪み、広範な gradient vanishing、spurious minima のために reliability を失う。global optimum ではなく shallow local minima や plateau region に落ちやすい。

**Fig. 5** の benchmark によれば、COBYLA の success rate (SR) は qubit 数が大きいと **約20%** に低下し、SPSA でも **約50%** 程度。量子化学で要求される **10⁻⁴** 級の tolerance には遠い。Hubbard の rugged, nonconvex surface は stable gradients をほとんど与えず、noise の有無にかかわらず premature convergence を招く。結論として local optimizers には構造的限界があり、noisy・高次元 VQE では population-based metaheuristics を検討すべきだと述べる。

### B. Metaheuristic algorithms の総論

Sec. VI の可視化が示したように、finite shots は optimization surface を局所法に不利な形へ変形する。そこで著者らは三段階評価を実施した。以下がその結果。

### C. Phase 1: Initial Screening

5-qubit Ising, tolerance **10⁻¹**, sampling noise あり、5 runs の初期スクリーニング。
Bio-inspired では **SOS**, **TPO**, **SOA**, **WOA** が通過し、**SOS** と **TPO** が速かった。
Evolutionary では **CMA-ES**, **DE**, **FPA**, **SHADE**, **HyDE**, **iL-SHADE** が成功し、とくに **CMA-ES** と **iL-SHADE** が最も一貫していた。
Human-inspired では **IBSO**, **ICA**, **LCO**, **FBIO** が進出し、**ICA** が最も効率的。
Math-based では **RUN** と **Harmony Search (HS)** のみ。HS は音楽系起源だが数学的定式化のためここに置かれている。
Physics-based では **FLA**, **NRO**, **SA (Cauchy / Boltzmann / fast schedules)** が通過。
Swarm-based では **PSO** と **iSOMA** のみで、iSOMA が優勢。
System-based では **WCA** のみ通過。
総じて **evolutionary algorithms** が最も成功した。global sampling と parameter adaptation が noise-induced distortions に強いからだと解釈される。

### D. Phase 2: Ising model 上の Function Evaluation 比較

3〜9 qubits の Ising（no external field）について、**10⁻¹** 精度で global minimum に到達するまでの **FEs** を比較した。ノイズ下では stricter tolerance に多くの手法が届かなかったため、この緩い tolerance を採用した。local minima にハマるか、数十万 FEs 改善しない場合は **"—"**。一般傾向として、qubit 数が増えるほど大多数の optimizer は FEs が急増するか、失敗する。これは parameter space の増大と noise による landscape distortion の結果である。

**CMA-ES** は全 qubit 数で明確なトップで、system size が増えても成長は比較的穏やか。**iL-SHADE** も強く、基準コストはやや高いが 9 qubits でも致命的には崩れない。**SA Cauchy** は小〜中規模で効率的だが、9 qubits では FEs が一桁以上増える。**HS** は小規模で良いが、6 qubits 以降で急増。標準 **DE** 変種は mixed で、**DE/best/1/bin** は 6〜7 qubits まではそこそこ、9 qubits では実用的でない FEs、**DE/best/1/exp** や **DE/rand1** はもっと早く崩れる。**iSOMA** や **PSO** など swarm-based は初期は良いが次元増大で急速に悪化する。**ICA**, **LCO**, **VCS**, **ALO**, **HyDE variants** などは大規模で失敗が目立つ。結論として、**CMA-ES**, **iL-SHADE**, 次点で **SA Cauchy** が noisy かつ scalable な少数派として分離される。

### Table II（Ising FE 比較の主要数値）

表の数値をそのまま拾うと、代表的には次の通り。
**CMA-ES**: 3Q **750**, 4Q **1200**, 5Q **1500**, 6Q **1800**, 7Q **2280**, 8Q **2700**, 9Q **3200**。
**SA Cauchy**: **751, 1801, 3500, 4201, 7000, 8001, 9451**。
**iL-SHADE**: **1035, 2039, 3274, 3333, 4368, 4374, 6559**。
**HS**: **356, 1241, 1740, 1726, 5046, 7066, 10586**。
**DE best1bin**: **948, 1520, 4460, 8112, 11676, 14720, 17136**。
**iSOMA**: **1357, 3245, 5552, 15263, 22115, 28144, 33899**。
**SOS**: **1866, 2720, 4880, 13440, 15320, 32000, 37843**。
その下では、**DE best1exp**, **LCO**, **ICA**, **GA**, **DE/rand1**, **CRO**, **FPA**, **SA fast**, **IBSO**, **BBO**, **WOA**, **SOMA T3A**, **HyDE**, **FBIO**, **HyDE-DF**, **PSO**, **VCS**, **SA Boltzmann**, **ALO**, **TPO**, **SOA** が続くが、多くが途中で **—** になっている。表を眺めると、9 qubits 付近ではかなり露骨に勝ち負けが分かれる。まるで就活のESみたいに残酷である。

### E. Phase 3: Hubbard model 上の収束

**Fig. 6** は **64 shots** の Hubbard 6-site, 192-parameter の収束。sampling noise が大きく、estimated energy の揺らぎで true gradient direction が見えにくい。Ising よりはるかに harsh で、correlated fermionic structure, sign constraints, stronger entanglement により rugged, nonconvex, many local traps な landscape になる。最近の研究では shallow traps にいる variational states でも phase information を含みうるが、classical optimizer が escape しにくいことが知られており、本研究もそれに沿う。

この条件で最良だったのは **fine tuned CMA-ES (CMA-ES-ft)**。最少 FEs で exact global minimum に到達した。**iL-SHADE** も高い FEs を要するが exact global minimum に到達。対して他手法は、比較的緩い tolerance に達するのにすら数万 FEs を要した。**SA with multivariate Cauchy** は速いが lower tolerances で global minimum 探索に苦戦。**Harmony Search** と **SOS** がそれに続く。興味深いのは、古典関数最適化では評判のよい **DE** と **iSOMA** が VQE landscape では低調なこと。SciPy 実装の DE は stagnation criterion により predetermined number of FEs まで到達せず早期終了した。つまり「古典で強い」が、そのまま noisy VQE で強いとは限らない。

**Fig. 7** は **5120 shots** の場合。shot 数増加で noise が減り、landscape は smoother になる。それでも傾向自体は変わらず、**fine tuned CMA-ES** が先頭、**iL-SHADE** が有望、**SA Cauchy** と **HS** は似た収束を見せるがやや劣る。SciPy 実装の **DEbest1bin**, **DEbest1exp** は meaningful solution に収束できず、**iSOMA** も高エラーで stagnate。local gradient-based optimizers である **SPSA** と **COBYLA** は noisy conditions で高い局所極小に早期終了しやすく、Hubbard の nonconvexity を裏づける。

## 先行研究との比較（Discussion 相当）

著者らは関連文献と照合して結果を議論する。
a. **DE in VQAs**: noiseless VQE では DE の binomial / exponential crossover が局所法より優れ、hybrid DE も有望とされるが、本研究では noise を入れると標準 DE は早期 stagnation で大きく悪化した。一方 **iL-SHADE** のような advanced adaptive DE は 192-parameter Hubbard でも競争力を維持した。
b. **CEC competitions**: DE 系は state-of-the-art とされるが、本研究では noisy VQE では **CMA-ES** と **SA** の方が広い意味で強く、deterministic benchmarks だけでは VQE の stochastic ruggedness を捉えきれないことが示された。
c. **Swarm intelligence and RL**: **iSOMA** は swarm の中では頑健だが、noise 下では **CMA-ES** と **iL-SHADE** に後れを取る。
d. **Noise-aware classical optimizers**: 先行研究同様、**COBYLA** や **BFGS** のような勾配系は modest noise でも崩れ、本研究の **COBYLA** と **SPSA** の失敗もそれを再現した。
e. **VQE benchmarks across correlated systems**: Hubbard の shallow traps に optimizer が捕まること、しかし trapped solutions も物理情報を持つこと、**CMA-ES** の careful hyperparameter tuning が有効なことなど、近年研究と整合する。また本研究は **SA Cauchy**, **HS**, **SOS** も noise-robust な選択肢であることを広げて示した。総合すると、**noiseless study の ranking は noisy, large-parameter VQE では再評価が必要**だというのが著者らの立場である。

## CONCLUSION（結論）

本研究は **50種類超の metaheuristic algorithms** を対象に、**sampling noise 下での robustness** に重点を置いて VQAs 用 optimizer を系統評価した。手順は、**5-qubit Ising** での初期スクリーニング、**9 qubits** までの Ising スケーリング、**192-parameter Hubbard model** での convergence study の三段階。これと並行して landscape visualization を行い、geometry と noise が optimizer behavior をどう形作るかを明らかにした。

landscape analysis によれば、noiseless Ising では smooth, convex basin が local descent を助けるが、noisy Ising と Hubbard では rugged・noise-distorted・flat regions が gradient-based methods を損なう。これが **COBYLA** と **SPSA** の系統的失敗を説明し、stable gradients に依存しない global, sampling-based metaheuristics の必要性を裏づける。

全 benchmark を通じて、**CMA-ES** が最も信頼できる optimizer だった。低い evaluation count で global minimum に到達し、system size に対して graceful に scaling した。DE 系では **iL-SHADE** だけが noisy conditions で競争力を維持し、adaptive strategies の重要性を示した。**SA Cauchy**, **Harmony Search**, **Symbiotic Organisms Search** も、とくに Hubbard 上で強かった。一方で **PSO**, **GA**, **WOA** など、よく知られた多くの metaheuristics は小規模を越えると stagnation するか、膨大な FEs を要した。Hubbard は最も厳しい試験であり、correlated fermionic structure が rugged, nonconvex surface を生む。ここでも **CMA-ES** と **iL-SHADE** が先導し、shot 数が増えると **SA Cauchy** と **HS** が差を縮めた。

総括すると、**noiseless benchmarks における metaheuristic performance は noisy quantum simulations の結果を予測しない**。deterministic competitions や exact Ising/Hubbard studies で強い手法でも、sampling noise を入れると急激に悪化する。例外的に **CMA-ES**, **iL-SHADE**, そして少数の stochastic heuristics が頑健だった。将来は、**local refinement** と **noise-robust global exploration** を組み合わせる hybrid strategies が、realistic VQE setting における効率と信頼性のギャップを埋める有望な方向だと結論づけている。

## CREDIT AUTHORSHIP CONTRIBUTION STATEMENT

**Vojtěch Novák**: Software, Investigation, Writing – original draft.
**Ivan Zelinka**: Supervision, Validation, Reviewing.
**Václav Snášel**: Conceptualization, Supervision.

## ACKNOWLEDGEMENTS

研究資金として **grant of SGS No. SP2024/008, VSB-Technical University of Ostrava, Czech Republic** が謝辞に記載されている。

## Appendix A: 各アルゴリズムの説明

付録では各 optimizer の由来・基本挙動・主要 hyperparameters を説明している。たとえば、
**CMA-ES** は covariance matrix adaptation により探索分布を適応的に更新する進化戦略。
**iL-SHADE** は success-history based parameter adaptation と linear population reduction を持つ advanced DE variant。
**DE** は mutation / crossover による multi-particle strategy。
**PSO** は swarm communication と personal/global best を使う。
**iSOMA** は leader 方向への migration に基づく stochastic optimization で、本論文では SOMA T3A より高性能な variant を採用。
**ABC** は employed / onlooker / scout bees による foraging 行動を模す。
**SOS** は mutualism / commensalism / parasitism に基づく。
**WOA** は humpback whale の bubble-net hunting。
**SA** は thermodynamic annealing。
**HS** は musician improvisation に着想。
**WCA** は water cycle。
この付録の役目は、「どのアルゴリズムを試したのか」をブラックボックス化せず、探索と exploitation の考え方を説明することにある。

## Appendix B / IX: Algorithm Parameters

この付録では、各 optimizer に対して使った hyperparameters を
**parameter_name=value [typical_range]**
の形式で列挙している。Phase 1 では**原則として default parameter settings** を用いた。default は software package 付属値または literature で一般的な値で、exploration / exploitation のバランスがよい robust starting points とみなされている。

例として、
**CRO**: `pop_size=50, po=0.4 [0.2-0.5], Fb=0.9 [0.6-0.9], Fa=0.1 [0.05-0.3], Fd=0.1 [0.05-0.5], Pd=0.5 [0.1-0.7], GCR=0.1 [0.05-0.2], gamma_min=0.02 [0.01-0.1], gamma_max=0.2 [0.1-0.5], n_trials=5 [2-10]`
**EP**: `pop_size=50, bout_size=0.05 [0.05-0.2]`
**ES**: `pop_size=50, lambda=0.75 [0.5-1.0]`
**FPA**: `pop_size=50, p_s=0.8 [0.5-0.95], levy_multiplier=0.2 [0.0001-1000]`
などが見える。付録全体では各群の多数アルゴリズムについて、こうした設定が並ぶ。

### CMA-ES の fine tuning

著者らは **CMA-ES-ft** だけ追加で fine tuning を行った。理由は、preliminary experiments で一貫して最も速く良い収束を示し、hyperparameter optimization の有力候補だったから。fine tuning には **IRACE** を用い、候補設定を競わせながら generation model を更新した。調整したのは
**population size ∈ [15, 120]**
**initial step-size σ₀ ∈ [0.1, 1.0]**。
**150 evaluations** を hyperparameter tuning に割り当て、surviving configurations を CMA-ES-ft の最適設定とした。一方で **iL-SHADE** は success-history adaptation と linear population reduction により、外部 tuning なしでも largely self-adjusting なので追加調整しなかった。著者らは、他手法も tuning によって改善し得るが、**高次元で swarm-based methods が弱い**という大きな性能階層は変わらないだろうと述べる。

## この論文の本質を一文で言うと

**「VQE の optimizer は、ノイズなしで強いものではなく、ノイズ下の歪んだランドスケープで勾配に頼らず適応的に動けるものが勝つ。そしてその筆頭が CMA-ES と iL-SHADE だった」**という論文です。
