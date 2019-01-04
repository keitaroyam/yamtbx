# KAMO

Author: [Keitaro Yamashita](https://docs.google.com/forms/d/e/1FAIpQLSdINVTX6HtreMzuyeQk7VLsycKLFAL3SmDdARqQg8zpt46MXw/viewform).
English version is [here](kamo-en.md).

## 概要
*KAMO (Katappashikara Atsumeta data wo Manual yorimoiikanjide Okaeshisuru) system*は，（高分子）単結晶X線回折データの自動処理＆マージのために開発中のプログラムです．
積分・スケーリングには基本的に[XDS package](http://xds.mpimf-heidelberg.mpg.de)を用いることを想定していますが，オプションで[DIALS](https://dials.github.io/)や[Aimless](http://www.ccp4.ac.uk/html/aimless.html)も使用可能です．
現在のところ，small wedgeデータ(5-10°/結晶)の自動処理・マージを主眼に開発しています．

SPring-8でのオンラインデータ解析のために設計されていますが，ローカルのデータに対しても使えるようになっています．

本マニュアルは2018-02-22現在のものです．

   * [KAMO](#kamo)
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
         * [KAMO](#kamo-1)
            * [チェックボックスを1つ1つチェックしていくのが面倒](#チェックボックスを1つ1つチェックしていくのが面倒)
            * [特定のディレクトリだけ処理したい（除きたい）](#特定のディレクトリだけ処理したい除きたい)
            * [格子定数が既知なのでそれを使って欲しい](#格子定数が既知なのでそれを使って欲しい)
            * [ヘッダが間違っているので正しい値を与えたい](#ヘッダが間違っているので正しい値を与えたい)
         * [kamo.multi_merge](#kamomulti_merge)
            * [精密化に使うためのデータはどこ？](#精密化に使うためのデータはどこ)
      * [ローカル環境での使用方法](#ローカル環境での使用方法)
         * [DIALSを利用した環境構築](#dialsを利用した環境構築)
               * [トラブルシューティング](#トラブルシューティング)
         * [KAMOのアップデート](#kamoのアップデート)
         * [起動](#起動)
      * [SPring-8以外の場所でオンライン処理を行う方法](#SPring-8以外の場所でオンライン処理を行う方法)
      * [文献](#文献)
         * [KAMOの引用](#kamoの引用)
         * [KAMOを利用した研究](#kamoを利用した研究)
      * [バージョン履歴](#バージョン履歴)


### 依存プログラム・ライブラリ
以下のプログラム・ライブラリを使用しています．

* [CCTBX](https://github.com/cctbx/cctbx_project/) with [CBFlib](http://www.bernstein-plus-sons.com/software/CBF/)
* [wxPython 2.8](http://www.wxpython.org/), [Matplotlib 1.3](http://matplotlib.org/), [Networkx 1.x](https://networkx.github.io/), [Numpy](http://www.numpy.org/), [SciPy](https://www.scipy.org/)
* [XDS](http://xds.mpimf-heidelberg.mpg.de), [xdsstat](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Xdsstat), [H5ToXds](eiger-ja.md#eiger2cbf-h5toxds互換) (H5ToXdsはEIGERの場合のみ)
* [CCP4](http://www.ccp4.ac.uk/) (BLEND, Pointless, Aimless, Ctruncate)
* [R](https://www.r-project.org/) with rjson (BLENDを使う場合必要)
* [DIALS](https://dials.github.io/) (完全には未対応)


### 注意
KAMOはまだ発展途上のプログラムです．インタフェース面や機能面などで，多くの欠陥・バグがありますので，何か不具合に気づかれたり，ご要望などありましたら，**決して遠慮せず**開発者までご連絡ください．

このマニュアルもまだ不足してる点が多いですが，とりあえず公開します．

## 使用方法
### GUIの起動
起動した場所の下（サブディレクトリを含む）にあるデータを処理します．

特定のサブディレクトリのみを対象にしたいときは，`include_dir=hoge-1111-\*`という感じでディレクトリ名を指定する(複数指定可)か，あるいは対象ディレクトリ名が書かれたリストファイルを与えます．

* 例1: BL32XUで実験中に自動処理を開始

	kamo bl=32xu
* 例2: BL32XUで自動データ収集(ZOO)を使っているときの自動処理

	kamo bl=32xu mode=zoo
* 例3: ローカルにあるデータを処理（ディレクトリを直接探索．上記の場合はBSSのログからデータセットの場所を探索）

	kamo bl=other

`-h`オプションを付けるとヘルプ＆全パラメータリストが見れます．

ビームラインにおいて過去のデータを処理する際は例えば`date=2015-12-31`とすると，指定日の2日前から指定日までのデータを探索します（デフォルトは`date=today`）．`bl=other`指定時は，date指定は関係なく単に指定したディレクトリ以下のデータがすべて処理されます．

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
SG | 推定された空間群
Resn. | 分解能リミットの推定値 (small-wedgeの場合あまりアテにならない！)


### Small wedgeデータのマージ
1. マージ対象のデータにチェックを入れ（失敗してるものを含んでもエラーにはなりません），`Multi merge strategy`のボタンを押し，少し待つ
2. 同じ格子（reindexも考慮して同じ格子になるもの）がグループ化され，データセットの数でソートされる
3. マージするグループと対称性を選ぶ．対称性は，一番頻度の高いものがデフォルトで選択されているが，既知の場合はそれを選ぶ．
5. Proceedボタンを押し，ターミナル画面を見る．指定した対称性と異なる対称で処理されたデータを，指定した対称で処理しなおしている
6. "Do It Yourself!"と表示されたら完了．Reindex operatorが存在する場合は表示されるので，留意する(`kamo.resolve_indexing_ambiguity`を使って解決できます)
7. ターミナルで指示された場所に移動し，スクリプトを修正・実行する．
スクリプト(たとえばmerge\_blend.sh)は以下のようになっている
```
# settings
dmin=2.8 # resolution
anomalous=false # true or false
lstin=formerge.lst # list of XDS_ASCII.HKL files
use_ramdisk=true # set false if there is few memory or few space in /tmp
# _______/setting

kamo.multi_merge \\
        workdir=blend_${dmin}A_framecc_b+B \\
        lstin=${lstin} d_min=${dmin} anomalous=${anomalous} \\
        space_group=None reference.data=None \\
        program=xscale xscale.reference=bmin xscale.degrees_per_batch=None \\
        reject_method=framecc+lpstats rejection.lpstats.stats=em.b+bfactor \\
        clustering=blend blend.min_cmpl=90 blend.min_redun=2 blend.max_LCV=None blend.max_aLCV=None \\
        max_clusters=None xscale.use_tmpdir_if_available=${use_ramdisk} \\
#        batch.engine=sge batch.par_run=merging batch.nproc_each=8 nproc=8 batch.sge_pe_name=par
```
8. このスクリプトを実行すると，まずBLENDによる格子定数に基づいた階層的クラスタリングが行われ，見つかったクラスタのうちcompletenessが90%以上・redundancy 2以上になるクラスタすべてについて，マージを試みる．まず単純にxscaleでマージ(run\_01/)し，そのマージ結果とのCCを計算することで悪いフレームを見つける．悪いフレームを除いてマージした結果(run\_02/)から，スケールの*B*値とerror modelの*b*を基準にOutlierを検出し，悪いデータセットを除いてマージした結果がrun\_03/に保存される．これが最終結果となる．各結果のディレクトリ/ccp4にはmtzおよびctruncateとphenix.xtirageのログも保存される．
9. 処理完了後，作業ディレクトリ(blend_\*/)以下のreport.htmlをブラウザで表示すれば全最終結果の統計値を一望できる．結果を受けて，場合によっては分解能リミットを変えて再実行する．精密化に使うクラスタの選び方は，だいたいCC<sub>1/2</sub>が最大になるものを選べば問題ないと思われる（フィードバックお待ちしています）．

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

対称性の判断は，まずINTEGRATE.HKLに対してpointlessを実行し，次にXDS_ASCII.HKLに対しても実行する．両者で判断が一致しない場合，Probabilityが大きい方を採用する．

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
`cc_clustering.d_min=`でCCの計算に使用する分解能リミットを制限できます．今の方法では，他のすべてのデータと共通反射を持つデータしかクラスタリングに用いることができないため，対称性の低い結晶の場合での使用は現実的ではありません．

#### bad frameの検出

全データをマージした結果と，各フレーム上強度とのCCを計算し，CCの悪いフレームを捨てる(`reject_method=framecc`)．
デフォルトでは`rejection.framecc.method=tukey rejection.framecc.iqr_coeff=1.5`になっており，Tukeyの基準(1.5\*IQR)で悪いフレームを検出する．例えば`rejection.framecc.method=abs rejection.framecc.abs_cutoff=0.9`とすることで，絶対的な閾値を設定することも可能．

#### bad datasetの検出
XSCALE.LPの統計値からbad datasetを検出する(`reject_method=lpstats`)．
デフォルトの挙動では，スケールの*B*とerror modelの*b*からTukeyの方法でOutlierを検出する(`rejection.lpstats.stats=em.b+bfactor rejection.lpstats.iqr_coeff=1.5`)．`rejection.lpstats.stats=`には`pairwise_cc`も選択可能である．`pairwise_cc`はデータセット間のCCを悪くしているデータセットを除くもので，デフォルトでは0.8以下になるデータセットを除くようになっている(`rejection.lpstats.pwcc.method=abs rejection.lpstats.pwcc.abs_cutoff=0.8`)が，Tukeyの方法を選ぶこともできる.

#### scaling referenceの選定
主に，最終的なOverall B-factorに影響を及ぼします．XSCALEでは最初に書いたINPUT_FILE=がリファレンスになる仕様ですが，KAMOのデフォルトでは`xscale.reference=bmin`になっており，*B*が最も小さい，つまり（XSCALEでは）最も分解能に対するfall-offが小さい（高分解能まで強度が出ている）ものをリファレンスにします．


## FAQ
### KAMO
#### チェックボックスを1つ1つチェックしていくのが面倒
全部にチェックを入れたいときは"Check all"ボタンを押すと，全部にチェックが入ります．
また，1つチェックを入れて，Shiftキーを押しながらもう1つにチェックを入れると，間にあるもの全てにチェックが入ります．

#### 特定のディレクトリだけ処理したい（除きたい）
`include_dir=`あるいは`exclude_dir`を使います．複数指定(`include_dir=`を何度も書く)やリストファイル(\*.lst)も指定可能です．

#### 格子定数が既知なのでそれを使って欲しい
KAMOを起動するときに，たとえば

```
kamo known.unit_cell=10,20,30,90,90,90 known.space_group=p222
```
という形で格子定数を与えることができます（必ずspace_groupとセットで与えてください）．
格子定数は指数付けの時の事前情報として使われ，また，一致しない場合はデータ処理が行われません．
この方法だと全処理対象に対して同じ格子定数を用いますので，複数種類のデータがあるときは気をつけて下さい．

また，格子定数を指定する事によって，かえって結果が悪くなるケースも見られますので，慎重に使って下さい．

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

#### report.htmlで樹形図のcompletenessと表の統計値が一致しない
kamo.multi\_mergeでは，クラスタリング処理の直後に仮に全てマージされた場合のcompleteness, multiplicityを計算し，その後にrejection込みのマージ処理を行います．
よって最初に表示されるcompleteness/multiplicityと，最終結果のcompleteness/multiplicityは，rejectionのために一致しません．

## ローカル環境での使用方法
### DIALSを利用した環境構築
注: これまでPHENIXの環境を利用することを推奨していましたが，最新機能の一部がDIALSのモジュールを利用するようになったため，今後はDIALSを利用することを推奨します．基本的にはPHENIX環境でも動かせます(以下のDIALSをPHENIXに読み替えてください)．

DIALS/PHENIXにはCCTBXおよびその依存関係が含まれているため，DIALS/PHENIXを利用することで簡単にKAMOを導入できます(CCTBXを手動で導入する必要がありません)．

1. CCP4, R (rjson packageも含め), XDS, Adxv (任意)をインストールする
   * XDS/XDSSTATのインストールは[XDSwiki/Installation](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Installation)を参照
   * EIGERデータを処理する場合は[H5ToXds](eiger-ja.md#eiger2cbf-h5toxds互換)も必要です
   * rjsonは，Rを導入後，Rをインストールしたユーザ(rootまたはソフトウェア管理用のユーザアカウント)で起動し，`install.packages("rjson")`とタイプすることでインストールできます．サーバを尋ねられた場合は適当に選択します．
2. [DIALS](https://dials.github.io/installation.html)-1.5以上をインストールする
3. networkxをdials.pythonから使えるようにする
   1. `cd $DIALS/build`
   2. `./bin/libtbx.python -m easy_install networkx==1.11`
4. scipyをdials.pythonから使えるようにする (DIALS 1.10以降の環境を使う場合は不要)
   1. Macの場合，App StoreからXcodeを導入後，
   2. `cd $DIALS/build`
   3. Linuxの場合: `./bin/libtbx.python -m easy_install scipy==0.18.1`<br>Macの場合: `./bin/libtbx.python -m pip install scipy==0.18.1`
5. 以下のコマンドを実行する
```
cd $DIALS/modules
git clone https://github.com/keitaroyam/yamtbx.git
cd $DIALS/build
./bin/libtbx.configure yamtbx
```

KAMOセットアップ後，
```
kamo.test_installation
```
を実行することで必要なパッケージがインストールされているかどうか確認できます．

##### トラブルシューティング

* scipyのインストールに失敗する（ビルドが始まった場合）
   * Macの場合，Xcodeのバージョンに合った[Command-line tools](https://developer.apple.com/download/more/)を導入し，[gfortran](http://gcc.gnu.org/wiki/GFortranBinaries#MacOS)をインストールして再度試す．
   * Linuxの場合はパッケージマネージャ(yum等)でblas-develとlapack-develをインストールして再度試す．
* (Mac) scipyのビルド時に`as: I don't understand 'm' flag!`というエラーで止まる
   * 上記の方法でgfortranを導入したか確認する．MacPortsを使用している場合は/opt/local/binを環境変数PATHから外してもう一度試す．[参考URL](https://stackoverflow.com/questions/41542990/while-installing-on-osx-sierra-via-gcc-6-keep-having-fatal-opt-local-bin-l). あるいは/usr/local/gfortran/bin/gfortranが優先利用されるようにPATHの設定を見直す．
* (Mac) scipyのビルド時に`gcc: error: unrecognized command line option ‘-stdlib=libc++’`というエラーで止まる
   * Xcode付属のgccが使われているか確認する(Xcodeおよびcommand-line toolsは導入済みか？他の方法で入れたgccがシステムに無いか)．/usr/bin/g++が使われるようにPATHの設定を見直す．
* DIALS/PHENIX環境を使っているのに，kamo.test\_installationでwxPythonがNGになる
   * まずDIALS/PHENIXのGUIがちゃんと立ち上がるか確認してください(dials.image\_viewerなど)．Ubuntuの場合，libjpeg62などを導入する必要があるそうです．
### KAMOのアップデート
以下の手順で最新版にアップデートできます
1. `cd $DIALS/modules/yamtbx`
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

## SPring-8以外の場所でオンライン処理を行う方法
この項目はビームラインスタッフ向けの記述です．

SPring-8以外の環境でのオンライン処理に対応するため，KAMOの起動オプションに`dataset_paths_txt=`を追加しました．このオプションは`bl=other`と一緒に指定します．
このオプションにはファイル名を指定し，そのテキストファイルは各行にデータセットのテンプレートと開始・終了番号を含みます（カンマ区切り）．
例えば以下のような感じです
```
/hoge/fuga/180730/01/data1_??????.img, 1, 360
/hoge/fuga/180730/02/data2_??????.img, 1, 3600
```
このようなファイルをデータ収集プログラムに（必ずデータ収集が完了してから）出力してもらえれば，KAMOは実験と並行して処理を進めることが可能です．
このテキストファイルは`logwatch_interval=`で指定した間隔でチェックされます(デフォルトは30秒)．

つまり例えば上記テキストファイルの名前がdataset\_paths.txtだとすると，以下のようにしてKAMOを起動します．
```
kamo bl=other dataset_paths_txt=dataset_paths.txt logwatch_interval=10
```

## 文献

### KAMOの引用

以下の論文を引用ください．
* Yamashita, Hirata, and Yamamoto (2018) "KAMO: towards automated data processing for microcrystals." *Acta Cryst. D*__74__ doi: [10.1107/S2059798318004576](https://doi.org/10.1107/S2059798318004576).

論文出版前は，当英語版ドキュメントのURL https://github.com/keitaroyam/yamtbx/blob/master/doc/kamo-en.md を引用して頂いていました．
両方引いて頂いても構いません．

また，XDS, DIALS, POINTLESS, BLENDなど一緒に使ったプログラムの文献も引用して頂くようお願いします．

### KAMOを利用した研究
1. Morimoto *et al.* (2019) "Crystal structure of the endogenous agonist-bound prostanoid receptor EP3." *Nature Chemical Biology* doi: [10.1038/s41589-018-0171-8](https://doi.org/10.1038/s41589-018-0171-8) PDB: [6AK3](http://www.rcsb.org/pdb/explore/explore.do?structureId=6AK3) Raw data: [CXIDB#91](http://www.cxidb.org/id-91.html)
1. Toyoda *et al.* (2019) "Ligand binding to human prostaglandin E receptor EP4 at the lipid-bilayer interface." *Nature Chemical Biology* doi: [10.1038/s41589-018-0131-3](https://doi.org/10.1038/s41589-018-0131-3) PDB: [5YWY](http://www.rcsb.org/pdb/explore/explore.do?structureId=5YWY) [5YHL](http://www.rcsb.org/pdb/explore/explore.do?structureId=5YHL) Raw data: [Zenodo#1173791](https://zenodo.org/record/1173791)
1. Suno *et al.* (2018) "Structural insights into the subtype-selective antagonist binding to the M<sub>2</sub> muscarinic receptor." *Nature Chemical Biology* doi: [10.1038/s41589-018-0152-y](https://doi.org/10.1038/s41589-018-0152-y) PDB: [5ZK8](http://www.rcsb.org/pdb/explore/explore.do?structureId=5ZK8) [5ZKC](http://www.rcsb.org/pdb/explore/explore.do?structureId=5ZKC) [5ZKB](http://www.rcsb.org/pdb/explore/explore.do?structureId=5ZKB) [5ZK3](http://www.rcsb.org/pdb/explore/explore.do?structureId=5ZK3) [5YC8](http://www.rcsb.org/pdb/explore/explore.do?structureId=5YC8) Raw data: [Zenodo#1172266](https://zenodo.org/record/1172266) [Zenodo#1094808](https://zenodo.org/record/1094808)
1. Shihoya *et al.* (2018) "Crystal structures of human ET<sub>B</sub> receptor provide mechanistic insight into receptor activation and partial activation." *Nature Communications* doi: [10.1038/s41467-018-07094-0](https://doi.org/10.1038/s41467-018-07094-0) PDB: [6IGK](http://www.rcsb.org/pdb/explore/explore.do?structureId=6IGK) [6IGL](http://www.rcsb.org/pdb/explore/explore.do?structureId=6IGL) Raw data: [SBDB/611](https://data.sbgrid.org/dataset/611) [SBDB/612](https://data.sbgrid.org/dataset/612) 
1. Oda *et al.* (2018) "Crystal structure of the red light-activated channelrhodopsin Chrimson." *Nature Communications* doi: [10.1038/s41467-018-06421-9](https://doi.org/10.1038/s41467-018-06421-9) PDB: [5ZIH](http://www.rcsb.org/pdb/explore/explore.do?structureId=5ZIH) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-Chrimson-data-(5ZIH))
1. Tanaka *et al.* (2018) "Crystal structure of Escherichia coli YidC revealing all core regions, including flexible C2 loop." *Biochemical and Biophysical Research Communications* doi: [10.1016/j.bbrc.2018.09.043](https://doi.org/10.1016/j.bbrc.2018.09.043) PDB: [6AL2](http://www.rcsb.org/pdb/explore/explore.do?structureId=6AL2)
1. Shimizu *et al.* (2018) "GEF mechanism revealed by the structure of SmgGDS-558 and farnesylated RhoA complex and its implication for a chaperone mechanism." *PNAS* doi: [10.1073/pnas.1804740115](https://doi.org/10.1073/pnas.1804740115)  PDB: [5ZHX](http://www.rcsb.org/pdb/explore/explore.do?structureId=5ZHX) Raw data: [Zenodo#1134209](https://zenodo.org/record/1134209)
1. Kim *et al.* (2018) "Crystal structure of the natural anion-conducting channelrhodopsin GtACR1." *Nature* doi: [10.1038/s41586-018-0511-6](https://doi.org/10.1038/s41586-018-0511-6) PDB: [6CSM](http://www.rcsb.org/pdb/explore/explore.do?structureId=6CSM) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-GtACR1-data-(6CSM))
1. Kato *et al.* (2018) "Structural mechanisms of selectivity and gating in anion channelrhodopsins." *Nature* doi: [10.1038/s41586-018-0504-5](https://doi.org/10.1038/s41586-018-0504-5) PDB: [6CSN](http://www.rcsb.org/pdb/explore/explore.do?structureId=6CSN) [6CSO](http://www.rcsb.org/pdb/explore/explore.do?structureId=6CSO) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-iC++-data-(6CSN-&-6CSO))
1. Ganasen *et al.* (2018) "Structural basis for promotion of duodenal iron absorption by enteric ferric reductase with ascorbate." *Communications Biology* doi: [10.1038/s42003-018-0121-8](https://doi.org/10.1038/s42003-018-0121-8) PDB: [5ZLE](http://www.rcsb.org/pdb/explore/explore.do?structureId=5ZLE) [5ZLG](http://www.rcsb.org/pdb/explore/explore.do?structureId=5ZLG) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-Dcytb-data-(5ZLE-&-5ZLG))
1. Maestre-Reyna  *et al.* (2018) "Twist and turn: a revised structural view on the unpaired bubble of class II CPD photolyase in complex with damaged DNA." *IUCrJ* doi: [10.1107/S205225251800996X](https://doi.org/10.1107/S205225251800996X) PDB: [5ZCW](http://www.rcsb.org/pdb/explore/explore.do?structureId=5ZCW)
1. Franz *et al.* (2018) "Structure of the bifunctional cryptochrome aCRY from *Chlamydomonas reinhardtii*." *Nucleic Acids Research* doi: [10.1093/nar/gky621](https://doi.org/10.1093/nar/gky621) PDB: [5ZM0](http://www.rcsb.org/pdb/explore/explore.do?structureId=5ZM0)
1. Asada *et al.* (2018) "Crystal structure of the human angiotensin II type 2 receptor bound to an angiotensin II analog." *Nature Structural & Molecular Biology*  doi: [10.1038/s41594-018-0079-8](https://doi.org/10.1038/s41594-018-0079-8) PDB: [5XJM](http://www.rcsb.org/pdb/explore/explore.do?structureId=5XJM)
1. Tsuyuguchi *et al.* (2018) "Crystal structures of human CK2α2 in new crystal forms arising from a subtle difference in salt concentration." *Acta Cryst. F*  doi: [10.1107/S2053230X18005204](https://doi.org/10.1107/S2053230X18005204) PDB: [5Y9M](http://www.rcsb.org/pdb/explore/explore.do?structureId=5Y9M)
1. Furukawa *et al.* (2018) "Remote Coupled Drastic β-Barrel to β-Sheet Transition of the Protein Translocation Motor." *Structure*  doi: [10.1016/j.str.2018.01.002](https://doi.org/10.1016/j.str.2018.01.002) PDB: [5YHF](http://www.rcsb.org/pdb/explore/explore.do?structureId=5YHF)
1. Negishi *et al.* (2018) "Supramolecular protein cages constructed from a crystalline protein matrix." *Chemical Communications* doi: [10.1039/C7CC08689J](https://doi.org/10.1039/C7CC08689J) PDB: [5YR1](http://www.rcsb.org/pdb/explore/explore.do?structureId=5YR1) [5YR9](http://www.rcsb.org/pdb/explore/explore.do?structureId=5YR9) [5YRA](http://www.rcsb.org/pdb/explore/explore.do?structureId=5YRA) [5YRB](http://www.rcsb.org/pdb/explore/explore.do?structureId=5YRB) [5YRC](http://www.rcsb.org/pdb/explore/explore.do?structureId=5YRC) [5YRD](http://www.rcsb.org/pdb/explore/explore.do?structureId=5YRD) Raw data: [Zenodo#1470889](https://zenodo.org/record/1470889)
1. Hori *et al.* (2018) "Na<sup>+</sup>-mimicking ligands stabilize the inactive state of leukotriene B<sub>4</sub> receptor BLT1." *Nature Chemical Biology* doi: [10.1038/nchembio.2547](https://doi.org/10.1038/nchembio.2547) PDB: [5X33](http://www.rcsb.org/pdb/explore/explore.do?structureId=5X33) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-BLT1-data-(5X33))
1. Suno *et al.* (2018) "Crystal Structures of Human Orexin 2 Receptor Bound to the Subtype-Selective Antagonist EMPA." *Structure*  doi: [10.1016/j.str.2017.11.005](https://doi.org/10.1016/j.str.2017.11.005) PDB: [5WQC](http://www.rcsb.org/pdb/explore/explore.do?structureId=5WQC) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-OX2R-data-(5WQC))
1. Miyauchi *et al.* (2017) "Structural basis for xenobiotic extrusion by eukaryotic MATE transporter." *Nature Communications* doi: [10.1038/s41467-017-01541-0](https://doi.org/10.1038/s41467-017-01541-0) PDB: [5Y50](http://www.rcsb.org/pdb/explore/explore.do?structureId=5Y50) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-AtDTX14-data-(5Y50))
1. Abe *et al.* (2017) "Structure of in cell protein crystals containing organometallic complexes." *Phys. Chem. Chem. Phys.* doi: [10.1039/C7CP06651A](https://doi.org/10.1039/C7CP06651A) PDB: [5YHA](http://www.rcsb.org/pdb/explore/explore.do?structureId=5YHA) [5YHB](http://www.rcsb.org/pdb/explore/explore.do?structureId=5YHB) Raw data: [Zenodo#1470891](https://zenodo.org/record/1470891)
1. Lee *et al.* (2017) "Structure of the triose-phosphate/phosphate translocator reveals the basis of substrate specificity." *Nature Plants* doi: [10.1038/s41477-017-0022-8](https://doi.org/10.1038/s41477-017-0022-8) PDB: [5Y78](http://www.rcsb.org/pdb/explore/explore.do?structureId=5Y78) [5Y79](http://www.rcsb.org/pdb/explore/explore.do?structureId=5Y79) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-TPT-data-(5Y78-&-5Y79))
1. Tanaka *et al.* (2017) "Crystal Structure of a Plant Multidrug and Toxic Compound Extrusion Family Protein." *Structure* doi: [10.1016/j.str.2017.07.009](https://doi.org/10.1016/j.str.2017.07.009) PDB: [5XJJ](http://www.rcsb.org/pdb/explore/explore.do?structureId=5XJJ)
1. Shihoya *et al.* (2017) "X-ray structures of endothelin ET<sub>B</sub> receptor bound to clinical antagonist bosentan and its analog." *Nature Structural & Molecular Biology* doi: [10.1038/nsmb.3450](https://doi.org/10.1038/nsmb.3450) PDB: [5XPR](http://www.rcsb.org/pdb/explore/explore.do?structureId=5XPR) [5X93](http://www.rcsb.org/pdb/explore/explore.do?structureId=5X93) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-ETBR-bonsentan-data-(5XPR))
1. Taniguchi *et al.* (2017) "Structural insights into ligand recognition by the lysophosphatidic acid receptor LPA<sub>6</sub>." *Nature* doi: [10.1038/nature23448](https://doi.org/10.1038/nature23448) PDB: [5XSZ](http://www.rcsb.org/pdb/explore/explore.do?structureId=5XSZ) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-LPA6-data-(5XSZ))
1. Abe *et al.* (2017) "Crystal Engineering of Self-Assembled Porous Protein Materials in Living Cells." *ACS Nano* doi: [10.1021/acsnano.6b06099](https://doi.org/10.1021/acsnano.6b06099) PDB: [5GQM](http://www.rcsb.org/pdb/explore/explore.do?structureId=5GQM) [5GQN](http://www.rcsb.org/pdb/explore/explore.do?structureId=5GQN) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-Polyhedra-data-(5GQM-&-5GQN))


## バージョン履歴
日付はGitHub公開時

* 2018-07-30
   * KAMO: SP8以外の施設でのオンライン処理に対応するため，dataset\_paths\_txt=オプションを追加．
* 2018-05-22
   * kamo.test\_installation: 2018-04-25の更新によるバグを修正(Numpy 1.11以前でも動くように)
   * マージ準備ルーチンのバグを修正
* 2018-04-25
   * 論文公開にあわせてドキュメントを更新
   * kamo.test\_installation: H5ToXdsの実行テストを追加
* 2018-02-22
   * kamo.multi\_merge: cc\_clustering.min\_common\_refs= オプションを追加(共通反射数の最小値を設定．下回るデータは除外されます)
   * CC<sub>1/2</sub> vs resolutionのカーブフィッティングのパラメータを改善
   * LCV計算におけるバグを修正
* 2018-01-30
   * KAMO: XDS Nov 11, 2017 (BUILT=20171218)のサポートを追加 (このバージョンから精密化に失敗した際にGXPARM.XDSが生成されない)
   * kamo.multi\_prep\_merging: formerge.lst等でファイル名がソートされるようにした
   * kamo.multi\_merge: cc\_clusteringにおけるエラー処理と樹形図サイズの設定を改善
* 2017-12-26
   * KAMO: マージ用スクリプトのデフォルト設定を変更．rejectionをb+Bに．filter\_cell.Rのバグを修正．
   * kamo.multi\_merge: cc\_clusteringでRでは無くSciPyを利用するように変更(SciPy >=0.18.1が必要)．共通反射数が3未満の場合はデータを除外．クラスタリングの方法を選択可能に．
   * kamo.multi\_merge: `resolution.estimate=`をデフォルトでTrueに．また分解能を実際に切るのでは無く，CC<sub>1/2</sub> vs d<sup>-2</sup>のカーブフィッティングの結果からカットオフを推定するように挙動を変更．
* 2017-12-03
   * kamo.multi\_merge cc_clustering時の重大なバグを修正 (共通反射が少ないために除外されるデータが存在する場合に正しく動作していなかった; 特に対称性の低い空間群ではよく起こる)．また，意図していなかったが実質的にはsqrt(1-cc)でWard法を実行していた事を報告します(実際にはsqrt(2-2cc)は距離としての意味を持つので，これは正しい取扱でした)．
   * kamo.multi\_merge: 異方性分析に関するバグ修正．
   * kamo.auto\_multi\_merge: postrefine=パラメータをcell\_method=に名称変更．コード整理．
   * KAMO: batch.engine=sgeが指定されているがqsubが利用できない場合に分かりやすいメッセージを表示．batch.sh\_max\_jobs=Autoをデフォルトに.
   * KAMO: logger登録時のバグを修正 (例外が無視されてしまい，何が悪いのか気づけなかった)．
* 2017-11-02
   * KAMO: 空間群の選択方法を変更 (既知の格子が与えられている場合はスケーリングでそれを必ず使用．未知の場合はPointlessのProbabilityに基づいて選択). GUI起動時に既知格子の値をチェック．HTML reportのバグ修正．
   * kamo.multi\_merge: frame\_ccを元のファイルでは無くxscale.hklを用いて計算するように変更
   * kamo.multi\_merge: Pointlessの空間群がreference settingになっていない場合のバグを修正
   * kamo.auto\_multi\_merge: batch.engine=sh時に動かなかったバグを修正．分解能決定のステップを細かく．
   * kamo.resolve\_indexing\_ambiguity: 異常なファイルが含まれている場合に警告を表示．無視するオプションを追加．
* 2017-10-12
   * kamo: XDSがexpireしている場合に警告を表示
   * kamo.multi\_merge: d\_max=の指定が動くように修正
   * kamo.auto\_multi\_merge: コマンドスクリプトをコミットし忘れてました
* 2017-09-30
   * kamo.multi\_merge: summary, report html等でXSCALE.LPからR値,Completenessを読む際に小数点以下が捨てられていたバグを修正
   * kamo.multi\_merge: Pointlessの実行時間を短縮 (xscale.hklを予めP1でマージして螺旋判定)，全体の実行時間をログに記録
   * KAMO: known.method=オプションにcorrect\_onlyを追加 (CORRECT時にのみ指定した格子定数・対称性を使用．xds.repeat>1時に有用かも)
   * kamo.resolve\_indexing\_ambiguity: reference-based modeで結果の表示法を変更 (selective-breedingと同じに)
* 2017-08-31
   * kamo.multi\_merge: 異方性分析の結果をHTMLレポートに追加．Anisoカラムを廃止して最良／最悪の軸沿いの分解能を表示．
   * kamo.multi\_merge: degrees\_per\_batch=パラメータを追加(frames\_per\_batch=の角度量版)
   * kamo.multi\_merge: cc\_clustering時にCCの平均・最小値を表示
   * kamo.auto\_multi\_merge: root\_dirをCSVの非必須パラメータに(datadir=引数でグローバル設定可能)
   * kamo.auto\_multi\_merge: 最良結果を最終runのみから選択するように変更．基準にOverall CC1/2に加えてOutershell CC1/2も使用．
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
