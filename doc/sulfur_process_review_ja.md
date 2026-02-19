# ERSEM 硫黄過程レビュー（文献ベース）

日付: 2026-02-15
対象リポジトリ: `ersem`
対象範囲: 「底質鉛直多層ではない」前提のもとで、現在実装が目的とする酸化還元帯の遷移（好気分解 -> 脱窒帯 -> 硫酸還元帯）をどこまで再現できているかを、コードと文献で評価。

## 指摘事項（重大度順）

1. 高: `testcases/sulfur_cycle_example.yaml` の設定が `pelagic_base` の実装能力と整合していません。
- コード根拠: `testcases/sulfur_cycle_example.yaml:29` は `model: ersem/pelagic_base` に `composition: e` を指定していますが、`src/pelagic_base.F90:145` では `h` のみ追加され、未知成分は `src/pelagic_base.F90:214` で fatal になります。
- 影響: 提示されている硫黄テストケースが、そのままでは再現実行しにくい可能性があります。

2. 高: Layer 3 の硫酸還元に伴うアルカリ度生成が未実装（TODO）で、炭酸系に系統誤差を入れます。
- コード根拠: 硫酸還元源は `src/benthic_sulfur_cycle.F90:315` で計算されていますが、アルカリ度寄与未実装が `src/benthic_sulfur_cycle.F90:461` に明記されています。
- 重要性: 硫酸還元は底質 TA を増加させる主要過程の一つで、未実装だと pH 緩衝能を過小評価し得ます。
- 文献根拠: 硫酸還元に伴うアルカリ度生成は既知（例: Johnston et al., 2014; Norbisrath et al., 2023）。
- 推論注記: 「硫酸還元 1 mol あたり TA +2」は反応式の化学量論からの推論です。

3. 高: 硝酸共役硫化物酸化が「脱窒（N2生成）」に固定され、DNRA 分岐がありません。
- コード根拠: `src/benthic_sulfur_cycle.F90:441` で `H2S_2` 消費、`src/benthic_sulfur_cycle.F90:442` で NO3 消費、`src/benthic_sulfur_cycle.F90:444` で N 生成先は `G4n`（N2）に固定。NH4 生成分岐は未表現です。
- 重要性: 近年の沿岸底質研究では、硫黄-窒素結合過程は条件により DNRA 側へ寄り、固定窒素が NH4 として保持されることが示されています。
- 文献根拠: Bourceau et al. (2023) は硫酸還元と DNRA の強い結びつきを報告。Jones et al. (2017) などでも硫化物が DNRA を促進し得ることが示されています。
- モデル影響: N2 への窒素損失を過大、NH4 保持を過小に見積もるリスクがあります。

4. 中: FeS 生成は「鉄制限なし」の一次反応シンクとして与えられており、反応性鉄の枯渇を表現できません。
- コード根拠: `src/benthic_sulfur_cycle.F90:416`, `src/benthic_sulfur_cycle.F90:420`, `src/benthic_sulfur_cycle.F90:424`, `src/benthic_sulfur_cycle.F90:431`。
- 重要性: FeS/黄鉄鉱埋没は実在シンクですが、容量は鉄鉱物供給と履歴に依存します。
- 文献根拠: 海底質の硫黄循環レビューでは、硫化物の固定割合は鉄循環と密接に結合し可変です（Jorgensen et al., 2019）。

5. 中: 活性硫黄プールは設計上、質量保存されません（硫酸明示なし・直接シンクあり）。
- コード根拠: `src/sulfur_cycle.F90:147`, `src/sulfur_cycle.F90:161`（pelagic S0 sink）、`src/benthic_sulfur_cycle.F90:336`, `src/benthic_sulfur_cycle.F90:454`（benthic S0 burial）。
- 重要性: 現象再現の簡略化としては許容可能ですが、硫黄収支の閉鎖検証ができず、チューニングの相殺誤差を見えにくくします。

6. 中: 酸化還元帯を厳密分離で扱うため、実海域の重なり・高速切替を表現しにくいです。
- コード根拠: Layer 3 の硫酸還元を常時有効とする前提は `src/benthic_sulfur_cycle.F90:310`、ゾーン固定の前提説明は `src/benthic_sulfur_cycle.F90:17`。
- 文献根拠: 典型的ゾーネーション（上層酸化、下層硫酸還元）は妥当な一次近似（Jorgensen, 1982）。一方で沿岸底質では硫酸還元と硝酸還元の同時進行・切替も観測されます（Bourceau et al., 2023）。
- モデル影響: 季節平均の一次近似には有効ですが、短周期の酸化還元振動には限界があります。

7. 低: Layer 2 の `S0_2` は生成項のみで、同層内の明示的消費項がありません。
- コード根拠: `src/benthic_sulfur_cycle.F90:443` は生成のみ。
- 解釈: `benthic_column_dissolved_matter` による輸送で再配分はされますが、同層内の酸化・不均化反応はこのモジュールでは明示されていません。
- 文献根拠: 元素状硫黄は浅い底質で中間プールを形成し、活発に変換されます。

## 目的に対して妥当な点

- 反応配置自体は、目的とする遷移を概ね満たしています。
  - 上層の酸素酸化: `src/benthic_sulfur_cycle.F90:327`
  - 中層の硝酸共役酸化: `src/benthic_sulfur_cycle.F90:339`
  - 深層の硫酸還元源: `src/benthic_sulfur_cycle.F90:310`
- 層間輸送は既存の 3 層拡散フレームワーク（`src/benthic_column_dissolved_matter.F90:282`）を利用しており、非多層制約下では実装として実用的です。

## 総合評価

現在実装は、3 層 ERSEM の制約下で「好気分解 -> 脱窒帯 -> 硫酸還元帯」を一次近似するという目的に対して、方向性は妥当です。特に硫化水素の発生・酸化バリアの表現には有効です。

一方、定量的な底質生物地球化学（とくに窒素分配と炭酸系応答）に使うには、構造的なバイアスが残ります。主要点は、(i) 硫酸還元 TA 項未実装、(ii) DNRA 分岐欠如、(iii) FeS シンクの鉄制約欠落です。

## 参照文献

- Jorgensen, B. B. (1982). Mineralization of organic matter in the sea bed — the role of sulfate reduction. Nature 296, 643-645. https://www.nature.com/articles/296643a0
- Jorgensen, B. B., Findlay, A. J., & Pellerin, A. (2019). The Biogeochemical Sulfur Cycle of Marine Sediments. Frontiers in Microbiology, 10:849. https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2019.00849/full
- Bourceau, O. M., et al. (2023). Simultaneous sulfate and nitrate reduction in coastal sediments. ISME Communications. https://pmc.ncbi.nlm.nih.gov/articles/PMC9992702/
- Jones, Z. L., et al. (2017). Sulfide-Induced Dissimilatory Nitrate Reduction to Ammonium Supports Anammox. Appl. Environ. Microbiol. https://pmc.ncbi.nlm.nih.gov/articles/PMC5514666/
- Finster, K., Liesack, W., & Thamdrup, B. (1998). Elemental sulfur and thiosulfate disproportionation by Desulfocapsa sulfoexigens. Appl. Environ. Microbiol. https://pmc.ncbi.nlm.nih.gov/articles/PMC124681/
- Norbisrath, M., et al. (2023). Alkalinity and nitrate dynamics reveal dominance of anammox in a hyper-turbid estuary. Biogeosciences. https://bg.copernicus.org/articles/20/4307/2023/
- Johnston, S. G., et al. (2014). Alkalinity capture during microbial sulfate reduction and implications for acidification. Geochimica et Cosmochimica Acta. https://doi.org/10.1016/j.gca.2014.01.009
- Sulpis, O., et al. (2022). RADIv1: a non-steady-state early diagenetic model. Geosci. Model Dev. https://gmd.copernicus.org/articles/15/2105/2022/
