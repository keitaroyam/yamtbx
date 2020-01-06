# EIGERデータフォーマットとデータ処理

English version is [here](eiger-en.md).

本稿ではSPring-8 BL32XU/BL41XU/BL44XUにおけるEIGER X 9M/16Mに関して記述します．
本内容はEIGERのfirmware更新などによって変更される可能性があります．

  * [EIGER HDF5について](#eiger-hdf5について)
     * [内部データの圧縮(filter)](#内部データの圧縮filter)
     * [SPring-8におけるmaster.h5の改変](#spring-8におけるmasterh5の改変)
     * [SPring-8におけるonlyhits.h5に関して](#spring-8におけるonlyhitsh5に関して)
     * [マルチ(複数データセットのまとまった)HDF5に関して](#マルチ複数データセットのまとまったhdf5に関して)
     * [4M ROI](#4m-roi)
     * [Bad pixelについて](#bad-pixelについて)
  * [ソフトウェアの導入](#ソフトウェアの導入)
     * [HDF5 software](#hdf5-software)
     * [eiger2cbf (H5ToXds互換)](#eiger2cbf-h5toxds互換)
     * [bitshuffle plugin (少し上級編)](#bitshuffle-plugin-少し上級編)
     * [Neggia plugin (XDSでpluginを使った処理を行いたい場合のみ)](#neggia-plugin-xdsでpluginを使った処理を行いたい場合のみ)
        * [コンパイル済みのバイナリを入手する](#コンパイル済みのバイナリを入手する)
        * [自分でビルドする](#自分でビルドする)
  * [イメージの表示](#イメージの表示)
  * [データ処理の方法](#データ処理の方法)
     * [XDS](#xds)
        * [H5ToXdsを使う方法](#h5toxdsを使う方法)
        * [pluginを使う方法](#pluginを使う方法)
        * [autoPROCを使う場合](#autoprocを使う場合)
     * [DIALS](#dials)
     * [iMosflm](#imosflm)
     * [HKL-2000](#hkl-2000)
  * [参考](#参考)

## EIGER HDF5について

hdf5 (.h5)はこれまでのイメージフォーマット(.img, .cbf)と違い，複数のイメージを1つのファイルにまとめて保存することができます．
HDF5はファイルシステムのような階層構造を持っており，それぞれの階層に任意のデータ（イメージデータや実験に関するパラメータなど）を持つことができます．

EIGERの出力するhdf5は，masterとdataの2種類からなります．例えばprefixをsampleとしてデータを収集すると，

* sample\_master.h5 - 実験に関する情報（波長やビームセンターなど）の他，検出器の設定情報
* sample\_data_??????.h5 - 回折イメージ

data h5には複数の回折イメージが含まれており，通常100枚ごとに1ファイル作られます．
つまり100枚以下ならばsample\_data\_000001.h5のみ，200枚以下ならそれに加えてsample\_data\_000002.h5も出力されます．

sample_master.h5にはdata h5へのリンク(hdf5のExternalLinkという機能)も含まれており，master h5を開くことで，透過的に回折データを読めるようになっています．
実際のデータはdata h5の方に入っていますが，それを意識すること無くデータを読めます．master h5を与えるだけでデータ処理が行えるのは，この仕組みからです．

### 内部データの圧縮(filter)

HDF5は内部的にデータを圧縮することができるようになっています('filter'と呼ばれる機能)．EIGERのデータはLZ4アルゴリズムまたはbitshuffle+LZ4 (bslz4)によって圧縮されるようになっており，SPring-8では通常はbslz4です(bslz4はlz4に比べ30-50%ほど小さくなります)．
master.h5は通常圧縮されていませんが，そのままだと巨大なので(flatfieldやpixel maskなどを含む)，SPring-8では圧縮および不要データの削除をするようにしています．

データを読む際には，圧縮されていようがいまいが，同じように読むことができます（HDF5のライブラリが自動的に解凍処理を行います）．ただし非標準のFilter (bslz4も該当)の場合は，プラグインをインストールしなければ読めません．

### SPring-8におけるmaster.h5の改変

SPring-8 MX-BLでは，EIGERから直接出力されるmaster.h5ではなく以下のような改変を行ったものをユーザに提供しています．
（含まれている情報を適切に修正し，かつ，ファイルサイズを削減する目的です）
* 不要データの除去
  * /entry/instrument/detector/detectorSpecific/detectorModule\_*/以下のflatfield,pixel\_mask,trimbitを削除
  * これらのデータはtrimbitを除き，検出器全体での情報が別に保存されている
* サイズの大きな配列データの圧縮(pixel\_maskのみbslz4，他はgzip+shuf)
* 結晶回転軸の情報を付与(/entry/sample/transformations/omega)
* omegaのテーブルを適切に修正(/entry/sample/goniometer/)

### SPring-8におけるonlyhits.h5に関して

回折スキャン(diffraction scan)の結果は巨大であることが多く，通常ほとんどスポットの写っていないフレームであるため，
SPring-8ではスポットの写っているフレームのみを抽出したonlyhits.h5を作成し，それだけを持ち帰れるようにしています．

中身はmaster.h5のコピーに加え，/entry/data以下にフレーム情報を記録しています．
/entry/data/直下にprefix + image numberの名前のgroupが存在し(そのasttributeとしてn\_spotsの情報あり)，その直下にあるdataset `data`がフレームのデータです．

onlyhits.h5は後述の`yamtbx.adxv_eiger`で開くことができます．

### マルチ(複数データセットのまとまった)HDF5に関して

例えば1つのループから多数のsmall-wedgeデータを収集した場合など，すべてのデータが1つのmaster.h5に集約されている場合があります．
詳細は後日更新します．

### 4M ROI

EIGER 9Mおよび16Mには4M ROI (region of interest)の機能があります．
詳細は後日更新します．

### Bad pixelについて

検出器のイメージデータで，ピクセル強度が65535または4294967295になっているピクセルがある場合があります(それぞれ16bitまたは32bitのとき)．
このようなピクセルになる原因として以下の可能性があります (厳密にはEIGERの運用モードにも依存します)
* pixel\_maskに登録されたbad pixelである (dead regionや，既知の故障pixelなど)
* 内部の12bitカウンタ(4M以上は800 Hz)で計数中に一度でも溢れた (overflow)
* countrate\_correction\_count\_cutoffを超えている (count rate補正の範囲外だった)
* その他？

この数字になったpixelがどう処理されるかは処理プログラム依存です．
cbfファイルを介して処理する場合，負値に置き換える挙動が一般的です．
（データ処理時，もともと死んでるピクセルなのか，強度もといcpsが大きすぎたピクセルなのか区別ができないのはちょっと問題な気も…）

## ソフトウェアの導入

以下のソフトウェアの導入を推奨します．

### HDF5 software
generate\_XDS.INPは[HDF5 software](https://www.hdfgroup.org/HDF5/release/obtain5.html)に含まれるh5dumpを利用します．
各種パッケージマネージャでhdf5を導入可能かと思います．MacでMacPortsをお使いの場合は，
```
sudo port install hdf5
```
で導入できます．
[hdfview](https://www.hdfgroup.org/products/java/release/download.html)も導入しておくと便利です(master h5の中身をツリー表示できるGUIです)．

### eiger2cbf (H5ToXds互換)
eiger2cbfはPILATUSでの標準フォーマットであるCBF形式に変換するためのツールです．
[Downloads - Mosflm](http://www.mrc-lmb.cam.ac.uk/harry/imosflm/ver721/downloads.html#Eiger2CBF)の下の方にリンク(Converter for Eiger files)があり，LinuxおよびMac版のバイナリが配布されています．

eiger2cbfは，XDSがh5を処理する際に必要とするH5ToXdsと互換性のあるプログラムです．
H5ToXdsは[Dectrisのサイト](https://www.dectris.com/features/features-eiger-x/h5toxds?path=products/eiger/eiger-x-for-synchrotron/details/eiger-x-4m)の一番下からダウンロードできますが，Linux版のみでありMac版は配布されていません．
またH5ToXdsはCBFにヘッダ情報を付与しません．
これらの理由から，eiger2cbfをH5ToXdsという名前に変えて(またはH5ToXdsという名前のリンクを作って)使うか，以下のようにしてH5ToXdsという名前のシェルスクリプトを作ることを推奨します(eiger2cbfが標準エラー出力に書き出すメッセージが，XDSの出力と重なって見にくいためこれを抑制しています)．
```
cd (任意のPATHが通った場所)
cat <<+ > H5ToXds
#!/bin/sh
eiger2cbf \$@ 2>/dev/null
+
chmod +x H5ToXds
```

(catの下2行の内容をコピーしてH5ToXdsという名前でファイルを保存する方法でも大丈夫ですが，その場合は`$@`の前のバックスラッシュを除いて下さい)

eiger2cbfは，以下のようにmaster.h5を与えて使います．

* `eiger2cbf sample_master.h5` .. 格納されているフレーム数を出力
* `eiger2cbf sample_master.h5 10 sample_10.cbf` .. 10枚目(最初の番号は1)のフレームをsample\_10.cbfとして保存
* `eiger2cbf sample_master.h5 1:100 sample` .. 1-100枚目をsample\_??????.cbfとして保存

BL32XUで2016Aあたりに測定されたデータの一部はflatfield補正が入ってない場合があり，その場合はデータ処理時に補正を行ったほうが良いです．
[修正版eiger2cbf](https://github.com/keitaroyam/eiger2cbf)を使って頂くとflatfield補正がない場合に補正済みcbfを出力できます．

### bitshuffle plugin (少し上級編)
必ずしも必須ではありませんが，[bitshuffle](https://github.com/kiyo-masui/bitshuffle) pluginを導入しておくとh5dumpやhdfview，あるいはadxv経由でbslz4圧縮されたhdf5も開けるようになります．
導入には開発環境が必要です(MacではXcode? [未検証])

```
(sudo) easy_install cython
(sudo) easy_install h5py
cd ~/tmp (任意の場所)
git clone https://github.com/kiyo-masui/bitshuffle.git
cd bitshuffle/bitshuffle
cython ext.pyx
cython h5.pyx
cd ..
(sudo) python setup.py install --h5plugin
```
デフォルトでは/usr/local/hdf5/lib/pluginにインストールされます．変更したい場合は`--h5plugin`の後に`--h5plugin-dir=`を付けてください．
インストールした場所に，環境変数`HDF5_PLUGIN_PATH`を設定する必要があります(bashなら`export HDF5_PLUGIN_PATH=`をcshなら`setenv HDF5_PLUGIN_PATH `を使って設定してください)

現在では[pip経由](https://pypi.org/project/hdf5plugin/)でもインストールできそうです(開発環境もおそらく不要)．

### Neggia plugin (XDSでpluginを使った処理を行いたい場合のみ)
Pluginを使うと，H5ToXdsを使わず（すなわち一時ファイルとしてcbfファイルを出力せずに）直接h5ファイルを処理できるようになります．DECTRISの公式pluginであるNeggiaは以下の方法で入手できます．

#### コンパイル済みのバイナリを入手する
[DECTRISの公式サイト](https://www.dectris.com/company/news/newsroom/news-details/process-eiger-data-with-xds-fast)からMac/Linux用のバイナリを入手できます．ただしユーザ登録が必要です．

#### 自分でビルドする
[dectris/neggia - Github](https://github.com/dectris/neggia)からコードを入手できます．


## イメージの表示
eiger2cbfを用いてcbfに変換すれば，adxvやその他cbfをサポートするビューアで表示できます．

hdf5形式のまま表示するには，以下の方法があります．

* [ALBULA](https://www.dectris.com/products/albula-software)
* [Adxv](http://www.scripps.edu/tainer/arvai/adxv.html) (要[bitshuffle plugin](#bitshuffle-plugin-少し上級編). まずmaster.h5を開いてからdata h5を開くと波長・カメラ長等の情報が反映されます)
* dials.image\_viewer (DIALSプログラムに同梱)
* yamtbx.adxv\_eiger (32XUで標準使用の拙作スクリプト．adxvを使用)

yamtbx.adxv\_eigerは以下の方法で導入できます．[KAMO](kamo-ja.md#ローカル環境での使用方法)またはyamtbxを導入済みの場合は既に使えるようになっています．

1. [PHENIX](http://www.phenix-online.org/)-1.11以上をインストールする
2. 以下のコマンドを実行する(yamtbxをcloneする場所はどこでも良いので，適当に読み替えて下さい)
```
cd $HOME
git clone https://github.com/keitaroyam/yamtbx.git
cd $PHENIX/modules
ln -s ~/yamtbx/yamtbx .
cd ../build
libtbx.configure yamtbx
```


## データ処理の方法

EIGERのデータを処理するには，以下の2通りの方法があります．

1. hdf5のまま処理する (XDS, DIALS, HKL-2000)
2. cbfに変換して処理する (iMosflm)

### XDS

[generate_XDS.INP](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Generate_XDS.INP) (rev. 0.63以降)を用いてXDS.INPを作成できます(要h5dump)．
generata\_XDS.INPには，直接master h5を与えてください．例えば以下のようにします．
```
generate_XDS.INP ../sample_master.h5
```

手動でXDS.INPを用意する場合，`NAME_TEMPLATE_OF_DATA_FRAMES=`にmaster h5の'master'を'??????'に置き換えたものを指定します．上の例だと以下のようになります．
```
NAME_TEMPLATE_OF_DATA_FRAMES= ../sample_??????.h5
```

実際の処理にはH5ToXdsまたはpluginが必要です．

#### H5ToXdsを使う方法
H5ToXdsはXDSの補助プログラムで，cbfを中間ファイルとして使用するものです．上記の通り，eiger2cbfをH5ToXdsとして使用することを推奨します．

#### pluginを使う方法
pluginを使うと，中間ファイルを生成せずに直接h5から処理できるようになるため，特にI/O律速な環境では処理速度の向上(10%程度)が期待できます．
XDS.INPにおいて
```
LIB= /where-you-installed/plugin-name.so
```
という形で，plugin (.soファイル)の場所を指定します．
たとえば上記の[Neggia plugin](#neggia-plugin-xdsでpluginを使った処理を行いたい場合のみ)が使えます．

**注意！** Neggia pluginを使う場合，2016年にBL32XUで収集したデータを処理できません．これはNeggiaが非公式なHDF5再実装を行っていてGZIP+SHUFで内部データが圧縮されたh5ファイルを扱えないためです．
処理したい場合，変換作業が必要です．yamtbxを導入済みの方は，
```
mv yours_master.h5 yours_master.h5.org
yamtbx.eiger_reduce_master yours_master.h5.org h5out=yours_master.h5
```
として，master.h5を変換して下さい．
yamtbxを導入されてない場合は，[こちらのスクリプト](https://github.com/keitaroyam/yamtbx/blob/master/yamtbx/dataproc/command_line/eiger_reduce_master.py)を`phenix.python` (ver. 1.11以降)から立ち上げて下さい．

#### autoPROCを使う場合
`ReversePhi="yes"`オプションを指定する必要があります．
see: [autoPROC wiki](https://www.globalphasing.com/autoproc/wiki/index.cgi?BeamlineSettings#spring8).

**注意！** 2017Aから2018A(の5月頃まで)の間に収集されたデータをautoPROCで処理しようとするとエラーになります．
これは，autoPROCのツールがbslz4で圧縮されたmaster.h5を処理できないためで，2018A(の5月頃)以降はpixel\_mask以外をgzip+shufによる圧縮に変更しました．
この期間中に測定したデータを処理したい場合は，上述の`yamtbx.eiger_reduce_master`を用いて変換してください．

### DIALS
ver. dev-652 (1.dev.193)以降はBL32XUのEIGER hdf5に対応しています．
最新版か，少なくともCCP4 7.0 update 013以降，またはPHENIX-1.11以降に同梱されているバージョンをご利用下さい．
最新版は[DIALS website](https://dials.github.io/installation.html)からダウンロードできます．

データ処理方法は[本家Tutorial](https://dials.github.io/documentation/tutorials/index.html)をご覧ください．ただし，最初のdials.importでは
```
dials.import ../sample_master.h5
```
といった形で，master.h5をimportして処理を始めます．

**注意！** DIALS (ver. 1.5以前)でEIGERのデータを処理すると，I(+),I(-)があべこべ，つまり異常分散差の符号が逆転する問題が起きています(2017年4月リリースのver. 1.5で修正されました)．Bijvoet差を使用しない場合は問題ありませんが，異常分散差フーリエを見たりSAD/MAD/SIRAS/MIRAS法などで位相決定を行う場合，正しい結果が得られません．dials.import後，次の処理に進む前にdatablock.jsonを開き，
```json
        "panels": [
          {
            "origin": [
              119.39999904559726,
              119.92500419409488,
              -119.99999428590824
            ],
            "fast_axis": [
              -1.0,
              0.0,
              0.0
            ],
            "name": "/entry/instrument/detector",
            "raw_image_offset": [
              0,
              0
            ],
            "slow_axis": [
              -0.0,
              -1.0,
              0.0
            ],
```
となっている箇所を見つけて下さい．`fast_axis`を`-1,0,0`から`1,0,0`に変更し，また`origin`の最初の数字の頭にマイナスを付ける必要があります(`origin`の値はビームセンターとカメラ長なのでデータによって異なります)．上の例だと以下のように変わることになります．
```json
        "panels": [
          {
            "origin": [
              -119.39999904559726,
              119.92500419409488,
              -119.99999428590824
            ],
            "fast_axis": [
              1.0,
              0.0,
              0.0
            ],
            "name": "/entry/instrument/detector",
            "raw_image_offset": [
              0,
              0
            ],
            "slow_axis": [
              -0.0,
              -1.0,
              0.0
            ],
```


### iMosflm

cbfに変換することで処理できます．
[eiger2cbf](#eiger2cbf-h5toxds互換)を用いてcbfに変換して下さい．

[HDF5を直接処理できるバージョン](https://www.mrc-lmb.cam.ac.uk/mosflm/mosflm-hdf5/)もできたようです．

### HKL-2000

ver. 714 (Sep 2016)より，hdf5を直接読めるようになりました．それ以前のバージョンをお使いの場合はcbfへの変換が必要です．

## 参考
* [EIGER X series (Dectris公式サイト)](https://www.dectris.com/products/eiger/eiger-x-for-synchrotron)
* [Eiger - XDSwiki](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Eiger)
* [Casanas et al. (2016) EIGER detector: application in macromolecular crystallography. Acta Cryst. D72, 1036-1048](https://doi.org/10.1107/S2059798316012304)
* [HDF5 - The HDF Group](https://www.hdfgroup.org/HDF5/)
* [NeXus](http://www.nexusformat.org/) 将来的にEIGER hdf5もNeXus互換形式になります（現在は限定的）
