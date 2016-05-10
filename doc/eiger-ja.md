# Eigerデータフォーマットとデータ処理

本稿ではBL32XUにおけるEiger X 9Mに関して記述します．
本内容はEigerのfirmware更新などによって変更される可能性があります．

## Eiger HDF5について

hdf5 (.h5)はこれまでのイメージフォーマット(.img, .cbf)と違い，複数のイメージを1つのファイルにまとめて保存することができます．
HDF5はファイルシステムのような階層構造を持っており，それぞれの階層に任意のデータ（イメージデータや実験に関するパラメータなど）を持つことができます．

Eigerの出力するhdf5は，masterとdataの2種類からなります．例えばprefixをsampleとしてデータを収集すると，

* sample_master.h5 - 実験に関する情報（波長やビームセンターなど）の他，検出器の設定情報
* sample_data_??????.h5 - 回折イメージ

data h5には複数の回折イメージが含まれており，通常100枚ごとに1ファイル作られます．
つまり100枚以下ならばsample_data_000001.h5のみ，200枚以下ならそれに加えてsample_data_000002.h5も出力されます．

sample_master.h5にはdata h5へのリンク(hdf5のExternalLinkという機能)も含まれており，master h5を開くことで，透過的に回折データを読めるようになっています．
実際のデータはdata h5の方に入っていますが，それを意識すること無くデータを読めます．master h5を与えるだけでデータ処理が行えるのは，この仕組みからです．

### 内部データの圧縮(filter)

HDF5は内部的にデータを圧縮することができるようになっています('filter'と呼ばれる機能)．EigerのデータはLZ4アルゴリズムまたはbitshuffle+LZ4 (bslz4)によって圧縮されるようになっており，通常はbslz4です(bslz4はlz4に比べ30-50%ほど小さくなります)．
master.h5は通常圧縮されていません(masterにはflatfieldやpixel maskなどそこそこ大きい情報が格納されています)が，BL32XUではbyteshuffle+gzipによって圧縮するようにしています．

データを読む際には，圧縮されていようがいまいが，同じように読むことができます（HDF5のライブラリが自動的に解凍処理を行います）．ただし非標準のFilter (bslz4も該当)の場合は，プラグインをインストールしなければ読めません．


## ソフトウェアの導入

まず以下のソフトウェアを導入します．

### HDF5 software
generate_XDS.INPは[HDF5 software](https://www.hdfgroup.org/HDF5/release/obtain5.html)に含まれるh5dumpを利用します．
各種パッケージマネージャ(Macの場合macportsなど)でhdf5を導入可能かと思います．

[hdfview](https://www.hdfgroup.org/products/java/release/download.html)も導入しておくと便利です(master h5の中身をツリー表示できるGUIです)．

### eiger2cbf (H5ToXds互換)
eiger2cbfはPILATUSでの標準フォーマットであるCBF形式に変換するためのツールです．
[Downloads - Mosflm](http://www.mrc-lmb.cam.ac.uk/harry/imosflm/ver721/downloads.html#Eiger2CBF)の下の方にリンク(Converter for Eiger files)があり，LinuxおよびMac版のバイナリが配布されています．

eiger2cbfは，XDSがh5を処理する際に必要とするH5ToXdsと互換性のあるプログラムです．
H5ToXdsは[Dectrisのサイト](https://www.dectris.com/EIGER_X_Features.html)の一番下からダウンロードできますが，Linux版のみでありMac版は配布されていません．
またH5ToXdsはCBFにヘッダ情報を付与しません．
これらの理由から，eiger2cbfをH5ToXdsという名前に変えて(またはH5ToXdsという名前のリンクを作って)使うことを推奨します．


## イメージの表示
eiger2cbfを用いてcbfに変換すれば，adxvやその他cbfをサポートするビューアで表示できます．

hdf5形式のまま表示するには，以下の方法があります．

* [ALBULA](https://www.dectris.com/Albula_Overview.html)
* [Adxv](http://www.scripps.edu/tainer/arvai/adxv.html) (要bitshuffle plugin. 但し実験情報が読まれないため分解能が正しく表示されない)
* dials.image_viewer (DIALSプログラムに同梱)


## データ処理の方法

Eigerのデータを処理するには，以下の2通りの方法があります．

1. hdf5のまま処理する (XDS, DIALS)
2. cbfに変換して処理する (iMosflm, 他)

### XDS

[generate_XDS.INP](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Generate_XDS.INP) (rev. 0.63以降)を用いてXDS.INPを作成できます(要h5dump)．
generata\_XDS.INPには，直接master h5を与えてください．例えば以下のようにします．
```
generate_XDS.INP ../sample_master.h5
```
実際の処理にはH5ToXdsが必要です．上記の通り，eiger2cbfをH5ToXdsとして使用することを推奨します．

手動でXDS.INPを用意する場合，`NAME_TEMPLATE_OF_DATA_FRAMES=`にmaster h5の'master'を'??????'に置き換えたものを指定します．上の例だと以下のようになります．
```
NAME_TEMPLATE_OF_DATA_FRAMES= ../sample_??????.h5
```

### DIALS

現在のところCCP4に含まれているバージョンは，BL32XUのEiger hdf5に対応していません．
少なくともdev-652以降は対応していますので，[Nightly builds](http://cci.lbl.gov/dials/installers/)からdev-652または最新版をダウンロードしてインストールして下さい．

データ処理方法は[本家Tutorial](http://dials.lbl.gov/documentation/tutorials/index.html)をご覧ください．ただし，最初のdials.importでは
```
dials.import ../sample_master.h5
```
といった形で，master.h5をimportして処理を始めます．

### iMosflm

cbfに変換することで処理できます．
[eiger2cbf](http://www.mrc-lmb.cam.ac.uk/harry/imosflm/ver721/downloads.html#Eiger2CBF)を用いてcbfに変換して下さい．

## 参考
* [EIGER X series (Dectris公式サイト)](https://www.dectris.com/EIGER_X_Detectors.html#main_head_navigation)
* [Eiger - XDSwiki](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Eiger)
* [HDF5 - The HDF Group](https://www.hdfgroup.org/HDF5/)
* [NeXus](http://www.nexusformat.org/) 将来的にEiger hdf5もNeXus互換形式になります（現在は限定的）
