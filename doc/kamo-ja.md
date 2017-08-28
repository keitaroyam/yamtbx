# KAMO

English version is [here](kamo-en.md).

## 概要
*KAMO (Katappashikara Atsumeta data wo Manual yorimoiikanjide Okaeshisuru) system*は，（高分子）単結晶X線回折データの自動処理＆マージのために開発中のプログラムです．
基本的に[XDS package](http://xds.mpimf-heidelberg.mpg.de)のフロントエンドという形態を取りますが，オプションで[DIALS](https://dials.github.io/)や[Aimless](http://www.ccp4.ac.uk/html/aimless.html)も使用可能です．
現在のところ，small wedgeデータ(5-10°/結晶)の自動処理・マージを主眼に開発しています．

SPring-8でのオンラインデータ解析のために設計されていますが，ローカルのデータに対しても使えるようになっています（但し多くのケースが未テストです）．

本マニュアルは2015-12-18現在のものです．

   * [概要](#概要)
      * [依存プログラム・ライブラリ](#依存プログラムライブラリ)
      * [注意](#注意)
   * [使用方法](#使用方法)
      * [GUIの起動](#guiの起動)
      * [GUIの説明](#guiの説明)
      * [Small wedgeデータのマージ](#small-wedgeデータのマージ)
      * [index ambiguityの解消 (kamo.resolve_indexing_ambiguity)](#index-ambiguityの解消-kamoresolve_indexing_ambiguity)
   * [KAMOは内部で何をやるのか](#kamoは内部で何をやるのか)
      * [データセットの検出](#データセットの検出)
      * [各wedgeの処理 (kamo)](#各wedgeの処理-kamo)
      * [マージの準備 (kamoの"Multi-merge strategy"ボタン)](#マージの準備-kamoのmulti-merge-strategyボタン)
      * [複数結晶に由来するデータのマージ (kamo.multi_merge)](#複数結晶に由来するデータのマージ-kamomulti_merge)
         * [データのクラスタリング](#データのクラスタリング)
         * [bad frameの検出](#bad-frameの検出)
         * [bad datasetの検出](#bad-datasetの検出)
         * [scaling referenceの選定](#scaling-referenceの選定)
   * [FAQ](#faq)
      * [KAMO](#kamo)
         * [チェックボックスを1つ1つチェックしていくのが面倒](#チェックボックスを1つ1つチェックしていくのが面倒)
         * [格子定数が既知なのでそれを使って欲しい](#格子定数が既知なのでそれを使って欲しい)
         * [ヘッダが間違っているので正しい値を与えたい](#ヘッダが間違っているので正しい値を与えたい)
      * [kamo.multi_merge](#kamomulti_merge)
         * [精密化に使うためのデータはどこ？](#精密化に使うためのデータはどこ)
   * [ローカル環境での使用方法](#ローカル環境での使用方法)
      * [DIALSを利用した環境構築](#dialsを利用した環境構築)
      * [KAMOのアップデート](#kamoのアップデート)
      * [起動](#起動)
   * [文献](#文献)
      * [KAMOの引用](#kamoの引用)
      * [KAMOを利用した研究](#kamoを利用した研究)
   * [バージョン履歴](#バージョン履歴)


### 依存プログラム・ライブラリ
以下のプログラム・ライブラリを使用しています．

* [CCTBX](http://cctbx.sourceforge.net/) with [CBFlib](http://www.bernstein-plus-sons.com/software/CBF/) (動作上必須)
* [wxPython 2.8](http://www.wxpython.org/), [Matplotlib 1.3](http://matplotlib.org/), [Networkx 1.x](https://networkx.github.io/), [Numpy](http://www.numpy.org/), [SciPy](https://www.scipy.org/) (動作上必須)
* [XDS](http://xds.mpimf-heidelberg.mpg.de), [xdsstat](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Xdsstat), [H5ToXds](eiger-ja.md#eiger2cbf-h5toxds互換) (EIGERの場合)
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
 - 同じ空間群でも異なる格子定数の候補がある場合があるので，頻度の表示(および既知情報があればそれも参考)に注意して選択して下さい．
 - <s>XDS_ASCII.HKL_noscaleのファイルはreindexの時上書きされるので，複数の候補を並行して用意することはできません．</s> "into the workdir"にチェックを入れておけば各作業ディレクトリにコピーされるので，複数の候補を試すことが可能です．
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
### DIALSを利用した環境構築
注: これまでPHENIXの環境を利用することを推奨していましたが，最新機能の一部がDIALSのモジュールを利用するようになったため，今後はDIALSを利用することを推奨します．基本的にはPHENIX環境でも動かせます(以下のDIALSをPHENIXに読み替えてください)．

1. CCP4, R (rjson packageも含め), XDSをインストールする
   * XDS/XDSSTATのインストールは[XDSwiki/Installation](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Installation)を参照
   * EIGERデータを処理する場合は[H5ToXds](eiger-ja.md#eiger2cbf-h5toxds互換)も必要です
2. [DIALS](https://dials.github.io/installation.html)-1.5以上をインストールする
3. networkxをdials.pythonから使えるようにする
   1. `cd $DIALS/build`
   2. `./bin/libtbx.python -m easy_install networkx==1.11`
4. scipyをdials.pythonから使えるようにする
   1. Macの場合，[gfortran](http://gcc.gnu.org/wiki/GFortranBinaries#MacOS)をインストール. Linuxの場合はパッケージマネージャ(yum等)でblas-develとlapack-develをインストール
   2. `cd $DIALS/build`
   3. `./bin/libtbx.python -m easy_install scipy==0.18.1`
5. 以下のコマンドを実行する(yamtbxをcloneする場所はどこでも良いので，適当に読み替えて下さい)
```
cd $HOME
git clone https://github.com/keitaroyam/yamtbx.git
cd $DIALS/modules
ln -s ~/yamtbx/yamtbx .
cd ../build
./bin/libtbx.configure yamtbx
```

KAMOセットアップ後，
```
kamo.test_installation
```
を実行することで必要なパッケージがインストールされているかどうか確認できます．

##### トラブルシューティング

* scipyの導入時に"as: I don't understand 'm' flag!"というエラーで止まる
   * MacPortsを使用している場合は/opt/local/binを環境変数PATHから外してもう一度試す．[参考URL](https://stackoverflow.com/questions/41542990/while-installing-on-osx-sierra-via-gcc-6-keep-having-fatal-opt-local-bin-l).
* DIALS/PHENIX環境を使っているのに，kamo.test\_installationでwxPythonがNGになる
   * まずDIALS/PHENIXのGUIがちゃんと立ち上がるか確認してください(dials.image\_viewerなど)．Ubuntuの場合，libjpeg62などを導入する必要があるそうです．

### KAMOのアップデート
以下の手順で最新版にアップデートできます
1. yamtbxをcloneした場所へ移動(`cd`)
2. `git pull`
3. `$DIALS/build/bin/libtbx.refresh`
4. `kamo.test_installation`

### 起動
基本的には上記と一緒ですが，イメージファイルをファイルシステムから探すため常に`bl=other`を指定してください．

```
kamo bl=other [batch.sge_pe_name=par]
```

また，SGEのparallel environment (qsub -pe の後に書く文字列)を上記のように指定して下さい(デフォルト: par)．SGEの環境が無く，ローカルコンピュータのみで動かすときは，

```
kamo bl=other batch.engine=sh batch.sh_max_jobs=2
```
として，同時に動かす最大ジョブ数を指定して下さい．

ゴニオメータの回転軸については，ヘッダから取得できない場合，[generate\_XDS.INP](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Generate_XDS.INP)と同様の方法で(つまりヘッダ情報を見て)判断します．
明示的に指示する場合は，`reverse_phi=false` (or `true`)や`rotation_axis=1,0,0`という形で指定してください．
あるいは`use_dxtbx=true`を指定すると，[dxtbx](https://doi.org/10.1107/S1600576714011996)を使って判断します．

## 文献

### KAMOの引用

現在論文準備中につき，当英語版ドキュメントのURL https://github.com/keitaroyam/yamtbx/blob/master/doc/kamo-en.md を引用して頂ますようお願いします．

### KAMOを利用した研究

* Shihoya *et al.* (2017) "X-ray structures of endothelin ET<sub>B</sub> receptor bound to clinical antagonist bosentan and its analog." *Nature Structural & Molecular Biology* doi: [10.1038/nsmb.3450](http://doi.org/10.1038/nsmb.3450) PDB: [5XPR](http://www.rcsb.org/pdb/explore/explore.do?structureId=5XPR) [5X93](http://www.rcsb.org/pdb/explore/explore.do?structureId=5X93)
* Taniguchi *et al.* (2017) "Structural insights into ligand recognition by the lysophosphatidic acid receptor LPA<sub>6</sub>." *Nature* doi: [10.1038/nature23448](http://doi.org/10.1038/nature23448) PDB: [5XSZ](http://www.rcsb.org/pdb/explore/explore.do?structureId=5XSZ)
* Abe *et al.* (2017) "Crystal Engineering of Self-Assembled Porous Protein Materials in Living Cells." *ACS Nano* doi: [10.1021/acsnano.6b06099](http://doi.org/10.1021/acsnano.6b06099) PDB: [5GQM](http://www.rcsb.org/pdb/explore/explore.do?structureId=5GQM) [5GQN](http://www.rcsb.org/pdb/explore/explore.do?structureId=5GQN) [Processig note](https://github.com/keitaroyam/yamtbx/wiki/Processing-Polyhedra-data-(5GQM-&-5GQN))


## バージョン履歴
日付はGitHub公開時

* 2017-07-22
   * MarCCDの拡張子なし形式(hoge.0001など)に対応
* 2017-07-20
   * (new) kamo.multi\_determine\_symmetry: 複数の(small wedge)データから空間群(点群対称のみ)を推定するプログラムを追加
   * KAMO: 既知格子定数の使い方を指定するknown.method=オプションを追加 (デフォルトはnot\_use\_firstで先ず事前情報無しで指数付け．use\_firstを指定すると最初から使用)
   * KAMO: マージ準備時に指定した格子定数にreindexしたHKLファイルをマージ用のディレクトリ以下にコピーするオプションを追加(デフォルトでON)
   * KAMO: ADSC検出器のビームセンターの読み取り方やreverse phiかどうかをヘッダ情報から判断するように変更 (デフォルトでON)
   * KAMO: (experimental) use\_dxtbx=オプションを追加．trueを指定するとdxtbxを使ってビーム・ゴニオメータ・検出器のジオメトリを取得．
   * KAMO: 動作ログを作業ディレクトリにも保存するように変更．log\_root=オプションの指定は任意．
   * GUIを軽量化 (ジョブ状況取得でGUIをブロックしないように変更; 格子定数によるデータのグループ化をバックグラウンドで進行)
   * kamo.resolve\_indexing\_ambiguity: selective breedingで並列計算を実装(nproc>1)
   * 各プログラムで基本的にmax\_delta=5をデフォルトに設定
   * filter\_cell.Rでplotを出力
* 2017-05-24
   * 2017-03-10の修正に含まれていたバグを修正(batchjob実行時の環境変数コピー)
* 2017-05-23
   * KAMO: DIALSのサポート(`engine=dials`)
   * KAMO: multi-merge準備作業の並列化，DIALSでjoint refinementを行うためのファイル出力オプションを追加(DIALS環境へ組み込まれていることが必要)．
   * kamo.multi_merge: 異方性分析のバグ修正
* 2017-04-19
   * KAMO: eiger2cbfで作成されたCBFファイルのサポートを追加．ビームラインでのオンライン処理時にEiger h5ファイルのダウンロードを待つように修正+minisetサポートを追加．
* 2017-03-24
   * kamo.test_installation: XDS.INPのある場所で実行すると走り出してしまう問題を回避．H5ToXdsのテストを追加．
   * KAMO: Multi-merge strategy開始時にP1 cellをCORRECT.LP_noscaleから読むように変更
* 2017-03-10
   * KAMO: 新しいMac (El Capitan以降)でローカル実行できなかった問題を修正.
* 2017-02-16
   * yamtbx.xds_aniso_analysis: 強度の異方性分析のプログラムを追加．kamo.multi_merge内でも実行．
   * KAMO: nprocからDELPHI=を決める際のバグを修正
   * kamo.resolve_indexing_ambiguity: selective-breedingにおけるバグを修正(少数のファイルの時に失敗する場合があった)
   * kamo.multi_merge: add_test_flag=trueのときに同一のテストセットを全結果に対して与えるように変更
* 2017-02-02
   * kamo.auto_multi_mergeを追加（テスト中）．複数サンプルのデータを同時にマージ
   * kamo.multi_merge reference.data= (test flagをコピー), resolution.estimate= (高分解能カットオフを自動決定)のオプションを追加．Pointlessを使用して螺旋を決定
   * kamo: blconfig= を複数指定可に．mode=を両方(zoo+normal)指定可に．
* 2017-01-18
   * yamtbx.beam_direction_plot: 複合格子の場合におかしくなるバグを修正
* 2016-12-26
   * kamo.multi_merge: `space_group=` オプションを追加(マージ時に使用)．指定がない場合Pointlessの結果をmtzに反映
   * kamo.multi_merge: 出力MTZにMULTIPLICITYカラムを追加
   * qsubするスクリプトの先頭でcdするように変更(bashrcなどでcdしている場合に，適切にcdされなかったバグを修正)
* 2016-12-06
   * GUI: `exclude_ice_resolutions=`オプションを追加．プロットがMacで更新されないバグを修正
   * XSCALE実行後の処理を高速化(ファイル名の置換)
   * kamo.resolve_indexing_ambiguity: `reference_label=`が与えられてない場合にクラッシュするバグを修正
   * kamo.test_installation: Adxvのチェックを追加
* 2016-10-05
   * `auto_frame_exclude_spot_based=`オプションを追加．最初と最後でセンタリングが外れているなど反射が写っていないフレームを含むデータに有効（かも）．
* 2016-07-18
   * XDS/XSCALE実行時にRAMディスクまたは一時ディレクトリを使用するオプションを追加(デフォルト)
   * 対称性でグルーピングする際，頻度計算に格子定数も考慮するよう変更（同じ空間群で異なるa,b,cがある場合に区別される）
   * OSX (phenix-1.10.1環境下?)におけるバグ（Multi mergeでProceedしても先に進まない）を修正
   * KAMOのHTML report作成の軽微なバグを修正
   * 非SGE環境でkamo.multi_mergeが実行できないバグを修正
   * LCVとaLCVを実際にマージされた結晶の格子定数から計算するように変更
   * phenix.xtriageのログのパースを修正．Anisotropyをmax(B_cart)-min(B_cart)と定義．
