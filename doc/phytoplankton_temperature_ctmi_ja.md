# 植物プランクトン水温応答の改良: CTMI モデル

## 概要

本文書は、ERSEMの一次生産者（植物プランクトン）の水温応答について、現行のQ10方式に加え、**CTMI（Cardinal Temperature Model with Inflection; 基本水温モデル）** を代替手法として導入する提案をまとめたものである。CTMIにより種群ごとに異なる最適水温（Topt）を設定でき、水温ニッチの分化に基づく季節的な植物プランクトン遷移の表現が可能となる。

主な動機は、東京湾のような高度に富栄養化した水域における溶存酸素（DO）再現性の向上である。現行のQ10方式では、冬季の珪藻ブルームや珪藻から渦鞭毛藻への遷移を表現できない。

## 問題の所在

### ERSEMの現行水温応答

`primary_producer.F90`（356行目）の水温関数:

```fortran
et = max(0.0_rk, self%q10**((ETW-10._rk)/10._rk) - self%q10**((ETW-32._rk)/3._rk))
```

数式表現:

```
f(T) = max(0, Q10^((T - 10)/10) - Q10^((T - 32)/3))
```

- 第1項: 10 degCを基準とした指数増加
- 第2項: 約32 degC以上での急激な抑制

標準値Q10 = 2.0（P1-P4共通）における実効的な水温応答:

| 水温 (degC) | f(T)  | 解釈 |
|-------------|-------|------|
| 5           | 0.71  | 冬季最低 |
| 10          | 1.00  | 基準温度 |
| 15          | 1.41  |      |
| 20          | 1.99  |      |
| 25          | 2.78  |      |
| 28          | 3.36  | ピーク付近 |
| 30          | 3.56  | 実効ピーク |
| 32          | 2.83  | 急減 |
| 35          | 0     | 増殖停止 |

**主要な制約:**

1. **約30 degCまで単調増加**: 低水温で最もよく増殖する種を表現する機構がない。10-15 degCに適応した珪藻は、同一曲線を共有する渦鞭毛藻と区別できない。

2. **全機能群が同一の水温応答を共有**: P1-P4はすべてQ10 = 2.0で、基準温度（10 degC）と抑制温度（32 degC）がハードコードされている。各群は最大生産性（`sum`）でのみ区別され、水温ニッチでは区別されない。

3. **季節的遷移を表現不可能**: 富栄養湾で観測される冬季珪藻ブルーム→春夏の渦鞭毛藻への移行→秋季珪藻復帰のパターンは、機能群間の最適水温の違いによって駆動される。単一の単調曲線ではこの交差を捉えられない。

### モデル精度への影響

TB-GOTM最適化研究（ERSEM_GOTM_opt_v3、1200トライアル）では、P1-P4の成長速度を最適化したにもかかわらず、表層DOの相関係数は0.43にとどまった。最適化の結果、P1（珪藻）は下限値に、P4（微小植物プランクトン）は上限値に張り付き、共通の水温関数では冬季と夏季の一次生産を同時に再現できないことが示された。

## 提案手法: CTMI

### CTMI の数式

CTMIの多項式形式を採用している。原型のRosso et al. (1993) 有理関数は
Topt < (Tmin+Tmax)/2 のとき [Tmin, Tmax] 内に特異点を持つため、
同じ4条件 (f(Tmin)=0, f(Tmax)=0, f(Topt)=1, f'(Topt)=0) を満たす
3次多項式に置換した。

```
f(T) = 0                                              (T <= Tmin or T >= Tmax の場合)

f(T) = (T - Tmin)(T - Tmax)(ctmi_a * T + ctmi_b)       (Tmin < T < Tmax の場合)

ここで:
  a = Topt - Tmin
  b = Topt - Tmax
  ctmi_a = -(a + b) / (a * b)^2
  ctmi_b = (a * b + (a + b) * Topt) / (a * b)^2
```

係数 ctmi_a, ctmi_b は初期化時に事前計算される。
実行時に安全クランプ max(0, min(1, f)) を適用する。

**特性:**
- f(Topt) = 1.0（構造上厳密）
- f(Tmin) = f(Tmax) = 0.0
- 単峰型3次多項式、全ての有効なパラメータ (Tmin < Topt < Tmax) で特異点なし
- 全3パラメータが直接的な生物学的意味を持つ:
  - **Tmin**: 増殖可能な最低水温 (degC)
  - **Topt**: 増殖速度が最大となる最適水温 (degC)
  - **Tmax**: 増殖可能な最高水温 (degC)

### 水温応答モデルの比較

| モデル | 数式型 | Topt表現? | 群別設定? | パラメータ |
|--------|--------|----------|----------|-----------|
| ERSM Q10（現行） | 指数+高温抑制 | 不可（単調） | 不可（同一Q10） | Q10 |
| Eppley (1972) | 指数包絡線 | 不可（単調） | 不可 | a, b |
| Eppley-Norberg (2004) | 包絡線×二次ニッチ | 可 | 可 | a, b, z, w |
| **CTMI（提案）** | 基本水温モデル | **可** | **可** | **Tmin, Topt, Tmax** |
| Sharpe-Schoolfield | Boltzmann-Arrhenius | 可 | 可 | Ea, Ed, Topt |

CTMIをEppley-Norbergモデルより推奨する理由:
1. パラメータが増殖実験から直接測定可能
2. Eppley包絡線に依存しない（Bissinger et al. 2008により増殖を過大評価する可能性が指摘）
3. 最もシンプルなパラメータ化（3対4パラメータ）
4. 15種の微細藻類で検証済み（Bernard & Remond 2012）

### `sum` の意味の変化に関する注意

現行Q10方式では、`sum` は「基準温度（10 degC）における最大比生産性」と定義される。CTMIでは、`sum` は「最適水温（Topt）における最大比生産性」となる。両方式ともf = 1.0となる温度で`sum`が適用されるため（Q10ではT=10、CTMIではT=Topt）、`sum`パラメータの単位と概念は一貫性を保つ。ただし、方式切替時に数値の調整が必要な場合がある。

## 実装設計

### 新規パラメータ

`type_ersem_primary_producer` 型定義に追加:

| パラメータ | 型 | YAMLキー | 単位 | デフォルト | 説明 |
|-----------|------|----------|------|---------|------|
| `iswTemp` | 整数 | `iswTemp` | - | 1 | 水温応答スイッチ（1: Q10, 2: CTMI） |
| `Tmin` | 実数 | `Tmin` | degC | 0.0 | 増殖最低水温（CTMI時のみ有効） |
| `Topt` | 実数 | `Topt` | degC | 20.0 | 最適水温（CTMI時のみ有効） |
| `Tmax` | 実数 | `Tmax` | degC | 35.0 | 増殖最高水温（CTMI時のみ有効） |

`iswTemp = 1`（デフォルト）の場合、従来のQ10方式が使用され、Tmin/Topt/Tmaxは無視される。完全な後方互換性が保持される。

### コード変更

修正が必要なファイルは **1つのみ**: `src/primary_producer.F90`

**型定義**（69行目の後）:
```fortran
integer :: iswTemp
real(rk) :: Tmin, Topt, Tmax
```

**パラメータ登録**（106行目の後、`initialize` サブルーチン内）:
```fortran
call self%get_parameter(self%iswTemp, 'iswTemp', '', &
     'temperature response (1: Q10 with high-T suppression, 2: CTMI cardinal temperature model)', &
     default=1, minimum=1, maximum=2)
call self%get_parameter(self%Tmin, 'Tmin', 'degrees_Celsius', &
     'minimum temperature for growth (CTMI only)', default=0.0_rk)
call self%get_parameter(self%Topt, 'Topt', 'degrees_Celsius', &
     'optimal temperature for growth (CTMI only)', default=20.0_rk)
call self%get_parameter(self%Tmax, 'Tmax', 'degrees_Celsius', &
     'maximum temperature for growth (CTMI only)', default=35.0_rk)
```

**水温応答**（`do` サブルーチン内）:
```fortran
! Temperature response
if (self%iswTemp == 1) then
   ! Original Q10 formulation with high-temperature suppression
   et = max(0.0_rk, self%q10**((ETW-10._rk)/10._rk) - self%q10**((ETW-32._rk)/3._rk))
else
   ! Polynomial CTMI (Cardinal Temperature Model with Inflection)
   ! Singularity-free cubic: f(Tmin)=0, f(Tmax)=0, f(Topt)=1, f'(Topt)=0
   ! Coefficients ctmi_a, ctmi_b precomputed in initialize.
   if (ETW <= self%Tmin .or. ETW >= self%Tmax) then
      et = 0.0_rk
   else
      et = (ETW - self%Tmin) * (ETW - self%Tmax) &
         * (self%ctmi_a * ETW + self%ctmi_b)
      et = max(0.0_rk, min(1.0_rk, et))
   end if
end if
```

### 変更不要なファイル

- `ersem_model_library.F90`: 新モデルクラスなし（既存の `primary_producer` を修正）
- `CMakeLists.txt`: 新ソースファイルなし
- ビルドスクリプト: 同一のコンパイル手順

## 東京湾への適用

### 東京湾の植物プランクトン季節的遷移

長期モニタリング（Nakada et al. 2021; Nagai et al. 2022）に基づき、東京湾の植物プランクトンは水温ニッチで分類できる:

| 季節 | 水温 | 優占分類群 | 機能群 |
|------|------|-----------|--------|
| 冬 (12-2月) | 8-15 degC | *S. japonicum*, *S. marinoi-dohrnii* | 低温適応珪藻 |
| 春 (3-5月) | 15-20 degC | *Chaetoceros*, *Thalassiosira* | 春季珪藻 |
| 夏 (6-8月) | 25-30 degC | *Prorocentrum*, *Heterosigma* | 渦鞭毛藻/ラフィド藻 |
| 秋 (9-11月) | 15-25 degC | 珪藻混合群集 | 移行期群集 |

主要な知見:
- *Skeletonema japonicum* は最低水温期（8-10 degC）にブルームし、15 degC以上では消失する（Nagai et al. 2022）
- 夏季の赤潮は *Heterosigma akashiwo* や渦鞭毛藻により、温水（>22 degC）を必要とする
- 移行期（5-6月）に珪藻から渦鞭毛藻への急激な優占種交代が起こる

### 東京湾向けCTMIパラメータの提案

文献値と水温生態学に基づく:

| 群 | ERSEMインスタンス | 代表する分類群 | Tmin (degC) | Topt (degC) | Tmax (degC) |
|----|-----------------|---------------|------|------|------|
| P1 | 珪藻 | 低温珪藻（*S. japonicum*, *Chaetoceros*） | 2 | 15 | 30 |
| P2 | ナノ植物プランクトン | 小型鞭毛藻、クリプト藻 | 5 | 20 | 33 |
| P3 | ピコ植物プランクトン | 小型独立栄養生物、シアノバクテリア型 | 8 | 25 | 35 |
| P4 | 微小植物プランクトン | 渦鞭毛藻（*Prorocentrum*, *Heterocapsa*） | 10 | 25 | 35 |

**P1（珪藻、Topt = 15 degC）の根拠**:
- *S. japonicum* は25 degC以上で不在、8-15 degCでピーク（Nagai et al. 2022）
- チェサピーク湾の珪藻モデルではTopt ~12-18 degCを使用（Cerco & Cole 1993）
- 2月のブルーミングを可能にし、夏季の珪藻を抑制し、観測と一致

**P4（渦鞭毛藻、Topt = 25 degC）の根拠**:
- *P. minimum* ブルームモデルではTopt ~20-25 degCを使用（Li et al. 2021）
- Anderson & Barton (2021) は渦鞭毛藻が全機能型中最も弱いQ10だが最も高いToptを持つことを示した
- 東京湾の夏季赤潮は温水への特殊化を必要とする

### 期待される水温応答曲線

提案CTMIパラメータによる各群の水温応答:

| 水温 | P1（珪藻） | P2（ナノ） | P3（ピコ） | P4（渦鞭毛藻） |
|------|-----------|-----------|-----------|---------------|
| 5 degC | 0.42 | 0 | 0 | 0 |
| 10 degC | 0.86 | 0.53 | 0.11 | 0 |
| 15 degC | **1.00** | 0.88 | 0.48 | 0.44 |
| 20 degC | 0.88 | **1.00** | 0.84 | 0.83 |
| 25 degC | 0.53 | 0.86 | **1.00** | **1.00** |
| 28 degC | 0.23 | 0.64 | 0.93 | 0.92 |
| 30 degC | 0 | 0.42 | 0.78 | 0.78 |

これにより望ましい遷移が実現される: 冬春（5-15 degC）に珪藻が優占、春（20 degC）にナノ植物プランクトンがピーク、夏（25-30 degC）にピコ/渦鞭毛藻が優占。

## 最適化への考慮

### チューニング対象としてのCTMIパラメータ

各群のTopt値はOptuna最適化の有力な候補である:

```python
tuning_targets = [
    # 最大生産性
    ("P1", "sum", "parameters", 0.5, 3.0),
    ("P2", "sum", "parameters", 0.5, 3.0),
    ("P3", "sum", "parameters", 0.5, 3.0),
    ("P4", "sum", "parameters", 0.5, 3.0),
    # 最適水温
    ("P1", "Topt", "parameters", 10.0, 20.0),  # 珪藻: 低温
    ("P4", "Topt", "parameters", 20.0, 30.0),  # 渦鞭毛藻: 高温
]
```

注: TminとTmaxは感度が低く文献値で固定可能だが、Toptは季節ピークのタイミングを直接制御する。

### 他パラメータとの相互作用

CTMIは以下のパラメータと強く相互作用する:
- **a0w**（光減衰）: 有光層深度を制御し、どの種がブルームできるかに影響
- **B1 sum**（細菌分解）: 生産に対するO2消費量を決定
- **R6 rm**（POM沈降速度）: 有機物が水柱で再無機化されるか底泥に到達するかを決定

包括的な最適化（v4）では、CTMIのTopt値をこれらのパラメータと並行してチューニングすべきである。

## ビルドとテスト手順

1. `~/Github/ersem/src/primary_producer.F90` を上記のとおり修正
2. クリーンビルド: `rm -rf ~/build`
3. コンパイル: `~/Github/fabm/src/drivers/gotm/install_ersem_gotm.sh`
4. バイナリをTB-GOTMディレクトリにコピー（シンボリックリンクではなくコピー）
5. 後方互換性テスト: `iswTemp: 1`（デフォルト）で実行し、修正前の出力と比較
6. CTMIテスト: `fabm.yaml` を提案CTMIパラメータに更新し、シミュレーション実行

## 参考文献

### CTMIモデル
- Rosso, L., Lobry, J.R., & Flandrois, J.P. (1993). An unexpected correlation between cardinal temperatures of microbial growth highlighted by a new model. *J. Theor. Biol.*, 162(4), 447-463.
- Bernard, O. & Remond, B. (2012). Validation of a simple model accounting for light and temperature effect on microalgal growth. *Bioresource Technology*, 123, 520-527.

### 水温-増殖関係
- Eppley, R.W. (1972). Temperature and phytoplankton growth in the sea. *Fishery Bulletin*, 70(4), 1063-1085.
- Norberg, J. (2004). Biodiversity and ecosystem functioning: A complex adaptive systems approach. *Limnol. Oceanogr.*, 49(4 part 2), 1269-1277.
- Bissinger, J.E. et al. (2008). Predicting marine phytoplankton maximum growth rates from temperature: Improving on the Eppley curve. *Limnol. Oceanogr.*, 53(2), 487-493.
- Anderson, S.I. & Barton, A.D. (2021). Marine phytoplankton functional types exhibit diverse responses to thermal change. *Nature Communications*, 12, 6413.
- Thomas, M.K. et al. (2012). A global pattern of thermal adaptation in marine phytoplankton. *Science*, 338, 1085-1088.

### 東京湾の植物プランクトン生態
- Nakada, S. et al. (2021). Phytoplankton species abundance in Tokyo Bay (Japan) from 1998 to 2019. *Ecological Research*, 36.
- Nagai, S. et al. (2022). Temporal niche partitioning of *Skeletonema*: seasonal succession in Tokyo Bay. *Aquatic Microbial Ecology*, 89, ame02000.
- Ogawa, Y. & Ichimura, S. (1997). 東京湾における赤潮事象と植物プランクトン群集組成の変遷 1907-1997. *海/La Mer*, 7(3), 159.

### 富栄養湾モデル
- Li, M. et al. (2021). A three-dimensional mechanistic model of *Prorocentrum minimum* blooms in eutrophic Chesapeake Bay. *Sci. Total Environ.*, 769, 144528.
- Cerco, C.F. & Cole, T. (1993). Three-dimensional eutrophication model of Chesapeake Bay. *J. Environ. Eng.*, 119(6), 1006-1025.
- Sohma, A. et al. (2008). A benthic-pelagic coupled ecosystem model for Tokyo Bay. *Ecol. Model.*, 215(1-3).

### ERSEMモデル
- Butenschon, M. et al. (2016). ERSEM 15.06: a generic model for marine biogeochemistry. *Geosci. Model Dev.*, 9, 1293-1339.
- Blackford, J.C. et al. (2004). Ecosystem dynamics at six contrasting sites. *J. Mar. Syst.*, 52, 191-215.

---

*文書作成日: 2026-02-14*
*対象ファイル: `~/Github/ersem/src/primary_producer.F90` 356行目*
