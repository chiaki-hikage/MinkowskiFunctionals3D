# MinkowskiFuctionals3D
本コードは3次元空間の密度場(例えば、宇宙に銀河やダークマターの数密度分布など)のミンコフスキー汎関数を計算するコードです。

ミンコフスキー汎関数は構造の形状やトポロジーを定量化したもので、３次元空間では、体積(V0)、表面積(V1)、平均曲率(V2)、オイラー数(V3)の4つがあります。
下のグラフはN体シミュレーションから計算した宇宙のダークマター分布のミンコフスキー汎関数を表しており、
体積(左上パネル)、表面積(右上パネル)、平均曲率(左下パネル)、オイラー数(右下パネル)を
密度のしきい値νの関数として表したものです。

<img width="621" alt="スクリーンショット 2021-10-14 11 28 37" src="https://user-images.githubusercontent.com/86592645/137240401-e46f22c4-402a-49d6-989b-e21a5d8cc627.png">

ミンコフスキー汎関数は密度場のスムージング(高波数成分のフィルタリング)のスケールR_Gによって値が変わります。
スムージングスケールが小さいほど、より小さい非線形スケールの構造を反映し、ミンコフスキー汎関数の形状も非対称になります。

以下の論文はスローン・デジタル・スカイ・サーベイ(SDSS)で観測された銀河データからミンコフスキー汎関数(あるいはジーナス統計)を計算し、
宇宙の標準的な構造形成モデルの予測と一致することや、銀河の明るさや形状などの属性による構造の違いを調べた先駆的な研究です。

"Minkowski Functionals of SDSS galaxies I: Analysis of Excursion Sets"

C. Hikage, J. Schmalzing, T. Buchert, Y. Suto, I. Kayo, A. Taruya,
M. Vogeley, F. Hoyle, J. R. Gott III, J. Brinkmann for the SDSS collaboration

Publ. Astron. Soc. Japan, Vol.55 No.5 (2003), pp.911-931 
