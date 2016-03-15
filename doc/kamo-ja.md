#KAMO
## 概要
*KAMO (Katappashikara Atsumeta data wo Manual yorimoiikanjide Okaeshisuru) system*は，（高分子）単結晶X線回折データの自動処理＆マージのために開発中のプログラムです．
基本的に[XDS package](http://xds.mpimf-heidelberg.mpg.de)のフロントエンドという形態を取りますが，オプションで[DIALS](https://dials.github.io/)や[Aimless](http://www.ccp4.ac.uk/html/aimless.html)も使用可能（になる予定）です．
現在のところ，small wedgeデータ(5-10°/結晶)の自動処理・マージを主眼に開発しています．

SPring-8でのオンラインデータ解析のために設計されていますが，ローカルのデータに対しても使えるようになっています（但し多くのケースが未テストです）．

本マニュアルは2015-12-18現在のものです．

### 依存プログラム・ライブラリ
以下のプログラム・ライブラリを使用しています．

* [CCTBX](http://cctbx.sourceforge.net/) with [CBFlib](http://www.bernstein-plus-sons.com/software/CBF/) (動作上必須)
* [wxPython 2.8](http://www.wxpython.org/), [Matplotlib 1.3](http://matplotlib.org/), [Networkx](https://networkx.github.io/), [Numpy](http://www.numpy.org/) (動作上必須)
* [XDS](http://xds.mpimf-heidelberg.mpg.de)
* [CCP4](http://www.ccp4.ac.uk/) (BLEND, Pointless, Aimless, Ctruncate)
* [R](https://www.r-project.org/) (BLEND, CCクラスタリングに必要) with rjson
* [DIALS](https://dials.github.io/) (完全には未対応)


### 注意
KAMOはまだ発展途上のプログラムです．インタフェース面や機能面などで，多くの欠陥・バグがありますので，何か不具合に気づかれたり，ご要望などありましたら，**決して遠慮せず**開発者までご連絡ください．

このマニュアルもまだ不足してる点が多いですが，とりあえず公開します．

## 使用方法
### GUIの起動
起動した場所の下（サブディレクトリを含む）にあるデータを処理します．

特定のサブディレクトリのみを対象にしたいときは，`include_dir=hoge-1111-\*`という感じでディレクトリ名を指定する(複数指定可)か，あるいは対象ディレクトリ名が書かれたリストファイルを与えます．

* 例1: BL32XUでsmall wedge

	kamo bl=32xu
* 例2: BL41XUで通常(1 or 数結晶でコンプリート)データ測定

	kamo bl=41xu small_wedges=false
* 例3: BL32XUでZoo (自動データ収集)を使ったとき

	kamo bl=32xu mode=zoo
* 例4: ローカルにあるデータを処理（ディレクトリを直接探索．上記の場合はBSSのログからデータセットの場所を探索）

	kamo bl=other

`-h`オプションを付けるとヘルプ＆全パラメータリストが見れます．

過去のデータを処理する際は例えば`date=2015-12-31`とすると，指定日の2日前から指定日までのデータを探索します（デフォルトは`date=today`）．

起動すると自動的にデータセットを探索し，処理を開始します．
処理結果は，`workdir=`で指定した場所(デフォルト:\_kamoproc/)に，同じディレクトリ構造で，*prefix*\_*startimage*-*endimage*/という名前のディレクトリができ，その下がXDS/DIALSの実行場所になっています．たとえば以下のような感じです．

```
~/151126_BL32XU/
│
└ mydata/ <- ここでKAMOを起動したとする
   ├ sample1/
   │  ├ data1_000{001..180}.img
   │  └ data2_000{001..360}.img
   ├ sample2/
   │  ├ data3_000{001..180}.img
   │  └ data4_000{001..360}.img
   │
   └ _kamoproc/ <- これが自動的に作られる (workdir=オプションで変更可)
      ├ sample1/ <- 同じツリー構造が作られる
      │  ├ data1_1-180/ <- これがデータ処理の作業ディレクトリ
      │  └ data2_1-360/
      └ sample2/
         ├ data3_1-180/
         └ data4_1-360/
```

### GUIの説明
テーブルのカラムの説明は以下のとおりです．カラムをクリックすることでソートできます．

カラム名 | 説明
------------ | -------------
Path | データセットの相対パス
Sample ID | Unipuck ID (Well ID)
Wavelength | 波長 (A)
TotalPhi | トータルの振動幅(°)
DeltaPhi | 1枚あたりの振動幅(°)
Cstatus | data **c**ollectionのステータス(never/running/finished)
Pstatus | data **p**rocessingのステータス(never/running/giveup/finished)
Cmpl. | データセットのCompleteness
SG | 空間群
Resn. | 分解能リミットの推定値 (small_wedges=trueの時はピークサーチの結果から，falseのときはCORRECT.LPの分解能100分割のテーブルでI/sigmaが1を切るところ)


### Small wedgeデータのマージ
1. マージ対象のデータにチェックを入れ（失敗してるものを含んでもエラーにはなりません），`Multi merge strategy`のボタンを押し，少し待つ
2. 同じ格子（reindexも考慮して同じ格子になるもの）がグループ化され，データセットの数でソートされる
3. マージするグループと対称性を選ぶ．対称性は，一番頻度の高いものがデフォルトで選択されているが，既知の場合はそれを選ぶ．
5. Proceedボタンを押し，ターミナル画面を見る．指定した対称性と異なる対称で処理されたデータを，指定した対称で処理しなおしている
6. "Do It Yourself!"と表示されたら完了．Reindex operatorが存在する場合は表示されるので，留意する(`kamo.resolve_indexing_ambiguity`を使って解決できます)
7. ターミナルで指示された場所に移動し，スクリプトを修正・実行する．
スクリプト(merge_blend.sh)は以下のようになっている
```
# settings
dmin=2.8 # resolution
anomalous=false # true or false	
# _______/setting
      
kamo.multi_merge \
    workdir=blend_${dmin}A_framecc_b \
    lstin=formerge.lst d_min=${dmin} anomalous=${anomalous} \
    program=xscale xscale.reference=bmin \
    reject_method=framecc+lpstats rejection.lpstats.stats=em.b \
    clustering=blend blend.min_cmpl=90 blend.min_redun=2 blend.max_LCV=None blend.max_aLCV=None
```
8. このスクリプトを実行すると，まずBLENDによる格子定数に基づいた階層的クラスタリングが行われ，見つかったクラスタのうちcompletenessが90%以上・redundancy 2以上になるクラスタすべてについて，マージを試みる．まず単純にxscaleでマージ(run_01/)し，そのマージ結果とのCCを計算することで悪いフレームを見つける．悪いフレームを除いてマージした結果(run_02/)から，error modelの*b*を基準にOutlierを検出し，悪いデータセットを除いてマージした結果がrun_03/に保存される．これが最終結果となる．各結果のディレクトリ/ccp4にはmtzおよびctruncateとphenix.xtirageのログも保存される．
9. 処理完了後，作業ディレクトリ(blend_\*/)以下のreport.htmlをブラウザで表示すれば全最終結果の統計値を一望できる．結果を受けて，場合によっては分解能リミットを変えて再実行する．精密化に使うクラスタの選び方は，だいたいCC1/2が最大になるものを選べば問題ないと思われる（フィードバックお待ちしています）．

### index ambiguityの解消 (kamo.resolve_indexing_ambiguity)
*P*6や*P*4など，あるいは*P*2でも&beta;~90&deg;の場合など，格子の対称性が空間群の対称性よりも高い場合，index ambiguityが生じます．複数の結晶を用いる場合，これを揃えておかなければなりません．

この作業は，上記のマージの前に行う必要があります．ただし，上記の作業中でReindex operatorが表示されなかったときは必要ありません．

```
kamo.resolve_indexing_ambiguity formerge.lst
```
とすると，必要に応じてreindexを行い，formerge_reindexed.lstを出力します．
kamo.multi_mergeの際には`lstin=formerge_reindexed.lst`を指定して下さい．

デフォルトは`method=selective_breeding`で，Kabschによる["Selective Breeding"アルゴリズム](http://dx.doi.org/10.1107/S1399004714013534)を使用します．これはReferenceデータを必要としません．同じくRefererenceデータを必要としない`method=brehm_diederichs` ([Algorithm 2 in Brehm and Diederichs paper](http://dx.doi.org/10.1107/S1399004713025431))も選択可能です．

Referenceデータに合わせたい場合は，`method=reference`を選択し，`reference_file=`にReferenceのMTZファイルを与えてください．



## KAMOは内部で何をやるのか
### データセットの検出
2つのモードがある

* BSSのログファイルからJob情報を検索（オンライン解析用）
* サブディレクトリを掘って探索（ファイルシステムのアクセス速度に依存）

`bl=other`を指定すると後者のモードになる．
前者の場合は，`date=`の指定とあわせて，BSSのログファイルを探索する．よって，ファイルやディレクトリがリネームor移動されているとデータセットが発見できない．

### 各wedgeの処理 (kamo)
XDSを使って処理を行う．基本的にgenerate_XDS.INPと同じ内容のXDS.INPを用いるが，指数付けには全フレームを用いる．
また，指数付け失敗時には，以下のことを検討する．

* IDXREF.LPの差ベクトルクラスタリングの結果から格子が半分になっていると思われる場合は，倍にして試す
* 既知の格子定数＆対称性が与えられている場合は，その利用を試す

対称性の判断は，まずINTEGRATE.HKLに対してpointlessを実行し，次にXDS_ASCII.HKLに対しても実行する．両者で判断が一致しない場合，ISaが大きい方を採用する．

高分解能カットオフはユーザが決めるべきとの立場を取るが，あまりにノイズの多い高角をスケーリングに含めると結果が適切にならないため，細かく切ったシェルでCC1/2が0.5を下回るところで分解能を切り，そのスケール結果をXDS_ASCII.HKL/CORRECT.LPとしている(GUIにはその結果が表示される)．分解能を切ってないXDS_ASCII_fullres.HKL/CORRECT_fullres.LPは別に保存される．

`small_wedges=true`のときは，CORRECTで経験的補正を行わない（つまりスケーリングを行わない；対称反射の数が少なすぎるため）バージョンの出力も作成(XDS_ASCII.HKL_noscale)し，マージの際にはそれが使われる．XDS_ASCII.HKL/XDS_ASCII_fullres.HKLの方は通常どおりスケーリング済みのものになっている．

### マージの準備 (kamoの"Multi-merge strategy"ボタン)
以下の処理が行われます．

1. 選択されたデータを，*P*1の格子でお互いに比較し，等価なものであるかどうかを調べる
2. 等価なもの同士を繋いでいき，グループ分けを行う
3. 各グループごとに，格子対称(格子定数から許される最大の対称性)を調べ，そのサブグループを選択可能な対称性として列挙する
4. 各wedgeごとに推定された対称性の頻度も表示する
5. ユーザがグループ番号と対称性を選択すると，各ファイル(XDS_ASCII.HKL_noscale)をその対称性の格子に変換する(reindexおよび格子定数の変換を行う)．
6. Index ambiguityが存在するかどうかを調べ，存在する場合はユーザに通知する．

### 複数結晶に由来するデータのマージ (kamo.multi_merge)
種々のオプションがありますが，基本的な流れは以下のとおりです．

1. 格子定数またはデータ間相関係数を用いて階層的クラスタリングを行い，見つかった各クラスタに対してマージを試みる．
2. 対象ファイルをXSCALEを使ってスケール＆マージ (run_01)
3. bad frameの検出＆除去 (run_02)
4. bad datasetの検出＆除去 (run_03)

#### データのクラスタリング
BLENDまたはデータセット間CCを用いたクラスタリングが可能です．BLENDについては本家のマニュアル・論文を参照のこと．

`clustering=cc`とすると，CCを用いたクラスタリングを利用できます．さらに，`cc_clustering.b_scale=`でWilson-Bによるスケーリングを行うかどうか，`cc_clustering.use_normalized=`で規格化構造因子を用いるかどうかを選択できます(true/false)．
`cc_clustering.d_min=`でCCの計算に使用する分解能リミットを制限できます．今のところ単純にRのhclust()関数を使っており，他のすべてのデータと共通反射を持つデータしかクラスタリングに用いることができないため，対称性の低い結晶の場合での使用は現実的ではありません．

#### bad frameの検出

全データをマージした結果と，各フレーム上強度とのCCを計算し，CCの悪いフレームを捨てる(`reject_method=framecc`)．
デフォルトでは`rejection.framecc.method=tukey rejection.framecc.iqr_coeff=1.5`になっており，Tukeyの基準(1.5\*IQR)で悪いフレームを検出する．例えば`rejection.framecc.method=abs rejection.framecc.abs_cutoff=0.9`とすることで，絶対的な閾値を設定することも可能．

#### bad datasetの検出
XSCALE.LPの統計値からbad datasetを検出する(`reject_method=lpstats`)．
デフォルトの挙動では，error modelの*b*からTukeyの方法でOutlierを検出する(`rejection.lpstats.stats=em.b rejection.lpstats.iqr_coeff=1.5`)．`rejection.lpstats.stats=`には`em.b`の他，`pairwise_cc`および`bfactor`が選択可能である．`pairwise_cc`はデータセット間のCCを悪くしているデータセットを除くもので，デフォルトでは0.8以下になるデータセットを除くようになっている(`rejection.lpstats.pwcc.method=abs rejection.lpstats.pwcc.abs_cutoff=0.8`)が，Tukeyの方法を選ぶこともできる.

#### scaling referenceの選定
主に，最終的なOverall B-factorに影響を及ぼします．XSCALEでは最初に書いたINPUT_FILE=がリファレンスになる仕様ですが，KAMOのデフォルトでは`xscale.reference=bmin`になっており，*B*が最も小さい，つまり（XSCALEでは）最も分解能に対するfall-offが小さい（高分解能まで強度が出ている）ものをリファレンスにします．


## FAQ
### KAMO
##### チェックボックスを1つ1つチェックしていくのが面倒
全部にチェックを入れたいときは"Check all"ボタンを押すと，全部にチェックが入ります．
また，1つチェックを入れて，Shiftキーを押しながらもう1つにチェックを入れると，間にあるもの全てにチェックが入ります．

#### 格子定数が既知なのでそれを使って欲しい
KAMOを起動するときに，たとえば

```
kamo known.unit_cell=10,20,30,90,90,90 known.space_group=p222
```
という形で格子定数を与えることができます（必ずspace_groupとセットで与えてください）．
格子定数は指数付けの時の事前情報として使われ，また，一致しない場合はデータ処理が行われません．
この方法だと全処理対象に対して同じ格子定数を用いますので，複数種類のデータがあるときは気をつけて下さい．

#### ヘッダが間違っているので正しい値を与えたい
該当イメージが存在するディレクトリに，`kamo_override.config`というファイルを用意すると，処理開始時にそこから情報を読んで使用します．以下の例のように書いて下さい．上書きする必要の無い情報は書かないで下さい．

```
wavelength= 1
distance= 100
orgx= 2500
orgy= 2500
osc_range= 0.1
rotation_axis= -1 0 0
```

### kamo.multi_merge
#### 精密化に使うためのデータはどこ？
run_03/ccp4/xscale.mtzを使って下さい．run_03/が無いときは，run_\*のうち一番数字が大きいディレクトリが最終サイクルです．

## ローカル環境での使用方法
まず，[依存関係](#依存プログラムライブラリ)をすべて導入してください (R, CCP4, XDS)．

次に，[INSTALL.md](../INSTALL.md)に従って環境構築を行いますが，難しければ，以下のPHENIXを使った方法が簡単です．

### PHENIXを利用した環境構築
1. CCP4, R (rjson packageも含め), XDSをインストールする
2. [PHENIX](http://www.phenix-online.org/)-1.10以上をインストールする
3. networkxをphenix.pythonから使えるようにする
   1. https://pypi.python.org/packages/source/n/networkx/networkx-1.11.tar.gz をダウンロード
   2. `tar xvf networkx-1.11.tar.gz`
   3. `cd networkx-1.11; phenix.python setup.py install`
4. 以下のコマンドを実行する(yamtbxをcloneする場所はどこでも良いので，適当に読み替えて下さい)
```
cd $HOME
git clone https://github.com/keitaroyam/yamtbx.git
cd $PHENIX/modules
ln -s ~/yamtbx/yamtbx .
cd ../build
libtbx.configure yamtbx
```

### 起動
基本的には上記と一緒ですが，

```
kamo bl=other log_root=~/kamo-log/ [batch.sge_pe_name=par]
```

のように，`log_root=`に書き込み権限のあるディレクトリを指定して下さい．
また，SGEのparallel environment (qsub -pe の後に書く文字列)を上記のように指定して下さい(デフォルト: par)．SGEの環境が無く，ローカルコンピュータのみで動かすときは，

```
kamo bl=other log_root=~/kamo-log/ batch.engine=sh batch.sh_max_jobs=8
```
として，同時に動かす最大ジョブ数を指定して下さい．

また，reverse-phiでは無いビームライン(SPring-8以外のビームラインは大体該当)では，`reverse_phi=false`を必ず指定して下さい．縦置きゴニオなどには未対応です．
