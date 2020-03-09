# KAMO

Author: [Keitaro Yamashita](https://docs.google.com/forms/d/e/1FAIpQLSdINVTX6HtreMzuyeQk7VLsycKLFAL3SmDdARqQg8zpt46MXw/viewform).
Japanese version is [here](kamo-ja.md).

## Overview
KAMO is the program for automated processing and merging of (MX) crystal diffraction data.
Basically KAMO uses [XDS package](http://xds.mpimf-heidelberg.mpg.de), but it supports [DIALS](https://dials.github.io/) and [Aimless](http://www.ccp4.ac.uk/html/aimless.html) optionally.
Currently, the development is focused on processing and merging of small wedge data (5-10°/crystal).

KAMO is designed for online analysis for SPring-8 MX beamlines; however, it can be used offline for local data. Please let me know if you need support for your beamline/detector.

This manual is for 2018-02-22.

*What does KAMO mean?* 'Kamo' itself means a mallard in Japanese. KAMO is an acronym of Japanese words standing for "the system to process thoroughly-collected data in a better way than manually" (Katappashikara Atsumeta data wo Manual yorimoiikanjide Okaeshisuru).

   * [KAMO](#kamo)
      * [Overview](#overview)
         * [Dependencies](#dependencies)
         * [Warning](#warning)
      * [Usage](#usage)
         * [GUI](#gui)
         * [GUI explanations](#gui-explanations)
         * [Merging small wedge data](#merging-small-wedge-data)
         * [Resolution of index ambiguity (kamo.resolve_indexing_ambiguity)](#resolution-of-index-ambiguity-kamoresolve_indexing_ambiguity)
      * [What KAMO does internally](#what-kamo-does-internally)
         * [Dataset detection](#dataset-detection)
         * [wedge data processing (kamo)](#wedge-data-processing-kamo)
         * [Preparing for merging (kamo's "Multi-merge strategy" button)](#preparing-for-merging-kamos-multi-merge-strategy-button)
         * [Merging of multiple wedges (kamo.multi_merge)](#merging-of-multiple-wedges-kamomulti_merge)
            * [Clustering](#clustering)
            * [Bad frame detection](#bad-frame-detection)
            * [Bad dataset detection](#bad-dataset-detection)
            * [Choice of scaling reference](#choice-of-scaling-reference)
      * [FAQ](#faq)
         * [KAMO](#kamo-1)
               * [I'm too lazy to click all checkboxes](#im-too-lazy-to-click-all-checkboxes)
            * [I want to process specified directories only or exclude some directories](#i-want-to-process-specified-directories-only-or-exclude-some-directories)
            * [I want to use prior unit cell information](#i-want-to-use-prior-unit-cell-information)
            * [I want to give the correct experimental parameters as image header is wrong](#i-want-to-give-the-correct-experimental-parameters-as-image-header-is-wrong)
         * [kamo.multi_merge](#kamomulti_merge)
            * [Where is the data for structural refinement?](#where-is-the-data-for-structural-refinement)
      * [How can I use KAMO at home?](#how-can-i-use-kamo-at-home)
         * [Installation](#installation)
               * [Troubleshooting tips](#troubleshooting-tips)
         * [How to update KAMO](#how-to-update-kamo)
         * [Launch](#launch)
      * [For online use at non-SPring-8 site](#for-online-use-at-non-spring-8-site)
      * [Citations](#citations)
         * [How to cite the use of KAMO](#how-to-cite-the-use-of-kamo)
         * [Researches which used KAMO](#researches-which-used-kamo)
      * [Version history](#version-history)

### Dependencies
KAMO uses following programs and libraries.

* [CCTBX](https://github.com/cctbx/cctbx_project/) with [CBFlib](http://www.bernstein-plus-sons.com/software/CBF/) (essential)
* [wxPython 2.8](http://www.wxpython.org/), [Matplotlib 1.3](http://matplotlib.org/), [Networkx 1.x](https://networkx.github.io/), [Numpy](http://www.numpy.org/), [SciPy](https://www.scipy.org/) (essential)
* [XDS](http://xds.mpimf-heidelberg.mpg.de), [xdsstat](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Xdsstat), [H5ToXds](eiger-en.md#eiger2cbf-h5toxds-compatible) (H5ToXds is only required for EIGER data)
* [CCP4](http://www.ccp4.ac.uk/) (BLEND, Pointless, Aimless, Ctruncate)
* [R](https://www.r-project.org/) with rjson (required for BLEND)
* [DIALS](https://dials.github.io/) (not fully-supported yet)


### Warning
KAMO is still under development and may have many many bugs on interface or functions. If you had any troubles, please do not hesitate to contact author.


## Usage
### GUI
The data in the directory where GUI started will be all processed (includes sub-directories).

If you want to process data in specific sub-directory, give the directory names like `include_dir=hoge-1111-\*` or give the list file containing the directory names.

* case 1: online analysis of small wedge at BL32XU

	kamo bl=32xu
* case 2: online analysis of data collected with ZOO (automatic system) at BL32XU

	kamo bl=32xu mode=zoo
* case 3: offline analysis of local data (KAMO will find directories. In above cases, datasets will be found using BSS log files)

	kamo bl=other

Use `-h` option to see help and list of all parameters.

For processing old data on beamline, for example, give `date=2015-12-31` to find data collected from the specified date to two days before (by default `date=today`).
When `bl=other`, the date parameter is ignored and files in specified directory are processed.

KAMO automatically finds datasets and start processing when the program started.
The processing results are saved in the directory specified by `workdir=` (default: \_kamoproc/), which has the same directory as data directory, and in *prefix*\_*startimage*-*endimage*/ XDS/DIALS is run. Example is below:

```
~/151126_BL32XU/
│
└ mydata/ <- Assuming KAMO started here
   ├ sample1/
   │  ├ data1_000{001..180}.img
   │  └ data2_000{001..360}.img
   ├ sample2/
   │  ├ data3_000{001..180}.img
   │  └ data4_000{001..360}.img
   │
   └ _kamoproc/ <- Automatically created (use workdir= option to change)
      ├ sample1/ <- the same structure as data directory
      │  ├ data1_1-180/ <- working directory for data processing program
      │  └ data2_1-360/
      └ sample2/
         ├ data3_1-180/
         └ data4_1-360/
```

### GUI explanations
The meaning of the table column labels. Click the column label to sort.

Column | Explanation
------------ | -------------
Path | Relative path of dataset
Sample ID | Unipuck ID (Well ID)
Wavelength | wavelength (A)
TotalPhi | Total rotation range (°)
DeltaPhi | Rotation width per 1 frame (°)
Cstatus | data **c**ollection status (never/running/finished)
Pstatus | data **p**rocessing status (never/running/giveup/finished)
Cmpl. | Dataset completeness
SG | Space group estimated
Resn. | Estimated resolution limit (unreliable for small-wedge cases)


### Merging small wedge data
1. Select all target data using the checkboxes (no problem if selecting failed ones), click `Multi merge strategy` button, and wait for a while
2. The data with the same lattice (reindexing is taken into account) are grouped, and sorted by the number of datasets
3. Select the group to be merged and its symmetry. The most frequent symmetry is chosen by default. Choose the correct one if known.
5. Click Proceed button, and look at the terminal. The datasets processed with different symmetry are reprocessed with the specified symmetry.
6. Completed when "Do It Yourself!" appeared. Keep in mind the reindex operator(s) if existed. Indexing ambiguity can be resolved using `kamo.resolve_indexing_ambiguity`.
7. Go to the working directory for merging, and modify and start the script. The scrip (for example merge\_blend.sh) is automatically prepared like this:
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
8. When this script is started, first the hierarchical clustering analysis by BLEND is performed. The clusters with high completeness (>90%) are subjected to merging. First simply run xscale to merge (run\_01/), and find the bad frames based on the CC between the merged result and intensities on each frame. From the merging result after excluding the bad frames (run\_02/), bad datasets are detected based on *B* scale and error model's *b* value, and the merging result without bad datasets is saved in run\_03/. This is the final result. In each run\_\*/ccp4 mtz and the log files of ctruncate and phenix.xtirage are saved.
9. When processing is over, report.html in the working directory (blend\_\*/) can be viewed with web browser. All final results can be visually inspected. Based on the results, re-execute the script by changing high resolution limi if needed. For refinement, empirically, the result with the largest CC<sub>1/2</sub> should be appropriate (we welcome your feedback).

### Resolution of index ambiguity (kamo.resolve_indexing_ambiguity)
Like *P*6, *P*4, or even in *P*2 with &beta;~90&deg;, when the lattice symmetry is higher than space group symmetry, the indexing ambiguity problem exists. This must be resolved when merging multiple crystals.

This procedure should be done before the merging command above *Only when* reindex operator is displayed in KAMO GUI.

```
kamo.resolve_indexing_ambiguity formerge.lst
```
performs reindexing and writes formerge_reindexed.lst.
For kamo.multi_merge give `lstin=formerge_reindexed.lst`.

Default is `method=selective_breeding` (Kabsch's ["Selective Breeding" algorithm](http://dx.doi.org/10.1107/S1399004714013534)), which does not need reference data.
You can also use `method=brehm_diederichs` ([Algorithm 2 in Brehm and Diederichs paper](http://dx.doi.org/10.1107/S1399004713025431)), which does not need it either.

When you want to use reference data, choose `method=reference` and give reference MTZ file by `reference_file=`.



## What KAMO does internally
### Dataset detection
Two modes available:

* Find job information from log file of BSS (for online analysis)
* Find datasets from filesystem (may be very slow; up to your filesystem performance)

When `bl=other`, the latter mode works.
In the former mode, `date=` parameter is used to decide the BSS log file names. In this case, any renamed or moved datasets cannot be found.

### wedge data processing (kamo)
Indexing is performed using XDS. Basically XDS.INP is the same as produced by generate_XDS.INP, but all frames are used for indexing.
When indexing failed, the following things are considered:

* Looking at difference vector clustering in IDXREF.LP; double the cell length if it seems to be halved.
* Try to index using the cell and symmetry if previously known

Concerning the decision of symmetry, first INTEGRATE.HKL is analyzed by pointless, and then XDS_ASCII.HKL is analyzed too. If results are not consistent, the one with larger probability is employed.

High resolution cutoff should be determined by user, but inclusion of much noise at higher resolution shell may affect scaling. KAMO cuts resolution at CC1/2 < 0.5 and the scaled result is saved as XDS_ASCII.HKL/CORRECT.LP (displayed on GUI). The result without any resolution cutoff is saved separately as XDS_ASCII_fullres.HKL/CORRECT_fullres.LP.

When `small_wedges=true`, another special version is prepared as XDS_ASCII.HKL_noscale, where the empirical correction in CORRECT job is switched off, which means no scaling because symmetry-related reflections are too few). This result is used in merging process. XDS_ASCII.HKL/XDS_ASCII_fullres.HKL are scaled results as usual.

### Preparing for merging (kamo's "Multi-merge strategy" button)
Internally,

1. For selected datasets, the unit cell in *P*1 is compared to each other, and checked if they are equivalent.
2. Construct a graph where datasets with equivalent cells are connected, and datasets are grouped.
3. For each group, the possible space group symmetries are listed based on its lattice symmetry.
4. The frequency of actually deduced symmetry for wedges is listed as well.
5. When a user chooses group number and symmetry, each result (XDS_ASCII.HKL_noscale) is transformed to the selected symmetry (reindexing and change of unit cell parameters).
 - Sometimes there are multiple candidates sharing the same space group (but different unit cell parameters). Decide using the frequencies and prior knowledge if any.
 - <s>XDS_ASCII.HKL_noscale files will be overwritten when reindexing. You cannot prepare multiple symmetry candidates at the same time (Sorry).</s> You can do that by checking "into the workdir" as files are copied to the working directory.
6. Check and notify if index ambiguity exits.

### Merging of multiple wedges (kamo.multi_merge)
Basically, the program

1. performs hierarchical clustering based on unit cell or pairwise CC, and clusters with high completeness are detected
2. scaling and merging of all datasets using XSCALE (run_01)
3. detects and removes bad frames (run_02)
4. detects and removes bad datasets (run_03)

#### Clustering
BLEND or CC-based clustering is available. See the original document and paper for BLEND.

When `clustering=cc`, CC-based clustering is invoked. Relevant options are `cc_clustering.b_scale=` (if scale data by Wilson-B before CC calculation) and `cc_clustering.use_normalized=` (if normalized structure factor is used).
Use `cc_clustering.d_min=` to limit the resolution. Currently only datasets which have common reflections with all others can be used. Probably not useful for low-symmetry crystals.

#### Bad frame detection

Correlation coefficient is calculated between merged all-data and each frame. Based on the CC, the bad frames are discarded (when `reject_method=framecc`).
By default `rejection.framecc.method=tukey rejection.framecc.iqr_coeff=1.5`, which detects outliers by Tukey's criterion (1.5\*IQR). You can change this; for example `rejection.framecc.method=abs rejection.framecc.abs_cutoff=0.9` to use absolute threshold.

#### Bad dataset detection
Bad datasets are detected using statistics in XSCALE.LP (`reject_method=lpstats`).
By default, *B* scale and error model's *b* value is subjected to Tukey's outlier detection (`rejection.lpstats.stats=em.b+bfactor rejection.lpstats.iqr_coeff=1.5`). You can also give `pairwise_cc` to `rejection.lpstats.stats=`. `pairwise_cc` is to remove datasets which give bad CC (by default CC<0.8; `rejection.lpstats.pwcc.method=abs rejection.lpstats.pwcc.abs_cutoff=0.8`, but optionally Tukey's method can be used).

#### Choice of scaling reference
This affects overall B-factor of merged data. In XSCALE, the first INPUT_FILE= is used as reference. In KAMO, by default `xscale.reference=bmin`, which selects data with smallest *B* (that has smallest intensity fall-off w.r.t. resolution in XDS) as reference.


## FAQ
### KAMO
##### I'm too lazy to click all checkboxes
If you want to check all, click "Check all" button.
Alternatively, click two checkboxes keeping Shift-key pressed to check all in-between items.

#### I want to process specified directories only or exclude some directories
Use `include_dir=` or `exclude_dir`.  Multiple directories can be specified by repeating `include_dir=`. Alternatively a list file (\*.lst) can be specified.

#### I want to use prior unit cell information
You can give it when KAMO starts:

```
kamo known.unit_cell=10,20,30,90,90,90 known.space_group=p222
```
Please make sure to give pace_group as well.
The unit cell parameters are used in indexing as prior information. If not indexed with the cell, the data processing will not proceed.
Note that you cannot use this parameter when you have more than one kind of crystals.

Please carefully use this function since sometimes this leads to worse result.

#### I want to give the correct experimental parameters as image header is wrong
Please prepare the file named `kamo_override.config` in the image file directory, which will be interpreted by the program. Example below (don't write information which is not need to be overridden):

```
wavelength= 1
distance= 100
orgx= 2500
orgy= 2500
osc_range= 0.1
rotation_axis= -1 0 0
```

### kamo.multi_merge
#### Where is the data for structural refinement?
Please use run_03/ccp4/xscale.mtz. If run_03/ is not there, run_\* with the largest number is the final one.

#### In report.html why do completeness values in the dendrogram and table differ?
kamo.multi\_merge first calculates completeness and multiplicity for each cluster just after clustering procedure, and performs merging procedudures including outlier rejections.
The completeness/multiplicity shown first and those of final results are different due to rejections.

## How can I use KAMO at home?

You can easily install KAMO using DIALS/PHENIX environment as DIALS/PHENIX includes CCTBX and its dependencies (No need to install CCTBX by yourself).

### Installation

1. Install CCP4, R (with rjson package), XDS, Adxv (optional)
   * For installation of XDS/XDSSTAT, see [XDSwiki/Installation](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Installation)
   * If you will process EIGER data (h5 files), [H5ToXds](eiger-en.md#eiger2cbf-h5toxds-compatible) is needed
   * rjson can be installed as follows; after installation of R, start R program from user who installed R (root or an account for software installation), and then type `install.packages("rjson")`.
2. Install [DIALS](https://dials.github.io/installation.html)-1.5 or newer (**choose Python2**)
3. Install networkx to dials.python
   1. `cd $DIALS/build`
   2. `./bin/libtbx.python -m easy_install networkx==1.11`
4. Install scipy to dials.python (**Skip this if using DIALS 1.10 or newer**)
   1. If Mac, install Xcode.
   2. `cd $DIALS/build`
   3. If Linux: `./bin/libtbx.python -m easy_install scipy==0.18.1`<br>If Mac: `./bin/libtbx.python -m pip install scipy==0.18.1`
5. Run the following commands
```
cd $DIALS/modules
git clone https://github.com/keitaroyam/yamtbx.git
cd $DIALS/build
./bin/libtbx.configure yamtbx
```

After installation, run
```
kamo.test_installation
```
to check if dependencies are all installed.

##### Troubleshooting tips
* scipy installation fails (in case building started)
   * If Mac, install [Command-line tools](https://developer.apple.com/download/more/) and [gfortran](http://gcc.gnu.org/wiki/GFortranBinaries#MacOS) and try again.
   * If Linux, install blas-devel and lapack-devel using yum or something.
* (Mac) building of scipy stops with error "as: I don't understand 'm' flag!"
   * Make sure you installed gfortran following the way mentioned above. If you are using MacPorts, try excluding /opt/local/bin from environment variable PATH. [Reference URL](https://stackoverflow.com/questions/41542990/while-installing-on-osx-sierra-via-gcc-6-keep-having-fatal-opt-local-bin-l). Alternatively review PATH to use /usr/local/gfortran/bin/gfortran.
* (Mac) building of scipy stops with error `gcc: error: unrecognized command line option ‘-stdlib=libc++’`
   * Make sure gcc installed with Xcode is used (did you install Xcode and command-line tools? Any other gcc installed with different methods?). Review PATH to use /usr/bin/g++.
* Installing with DIALS/PHENIX environment, but kamo.test\_installation claims wxPython is NG
   * Please check if some GUI program of DIALS/PHENIX successfully starts (e.g. dials.image\_viewer)．In case of Ubuntu, you may need to install libjpeg62 package.

### How to update KAMO
1. `cd $DIALS/modules/yamtbx`
2. `git pull`
3. `$DIALS/build/bin/libtbx.refresh`
4. `kamo.test_installation`

### Launch
Basically, the same as above; but specify always `bl=other` to find image files from filesystem.

```
kamo bl=other [batch.sge_pe_name=par]
```

In addition, give SGE's parallel environment (the strings you usually write after qsub -pe; by default par). If no SGE environment and you want to run it only on local computer, give maximum job number instead:

```
kamo bl=other batch.engine=sh batch.sh_max_jobs=2
```

Goniometer rotation axis is, if header does not have that information, recognized in the same way as [generate\_XDS.INP](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Generate_XDS.INP) which uses header information.
If you want to specify the rotation axis, you can give `reverse_phi=false` (or `true`) or `rotation_axis=1,0,0`.
Alternatively, when `use_dxtbx=true`, [dxtbx](https://doi.org/10.1107/S1600576714011996) is used to recognize beamline geometry.

## For online use at non-SPring-8 site
Here an instruction for beamline staff is described.

For the online use at non-SPring-8 sites, an option `dataset_paths_txt=` has been added to KAMO. This works with `bl=other`.
A file path should be given to this option, and the text file should include dataset template, start/end frame numbers in each line (separated by commas).
The example follows:
```
# comment line
/hoge/fuga/180730/01/data1_??????.img, 1, 360
/hoge/fuga/180730/02/data2_??????.img, 1, 3600
```
KAMO can start data processing immediately if a data collection program updates this text file just after the collection.
This text file is regularly checked with an interval specified by `logwatch_interval=` (30 sec by default).

In conclusion, start KAMO like this if the text file name is dataset\_paths.txt:
```
kamo bl=other dataset_paths_txt=dataset_paths.txt logwatch_interval=10
```


## Citations

### How to cite the use of KAMO

Please cite following paper.
* Yamashita, Hirata, and Yamamoto (2018) "KAMO: towards automated data processing for microcrystals." *Acta Cryst. D*__74__ doi: [10.1107/S2059798318004576](https://doi.org/10.1107/S2059798318004576).

You can also cite this documentation's URL: https://github.com/keitaroyam/yamtbx/blob/master/doc/kamo-en.md

Please also cite literatures of internally used programs like XDS, DIALS, POINTLESS, BLEND.

### Researches which used KAMO 

1. Yamaguchi *et al.* (2020) "Crystal structure of Drosophila Piwi" *Nature Communications* doi:[10.1038/s41467-020-14687-1](https://doi.org/10.1038/s41467-020-14687-1) PDB: [6KR6](https://www.rcsb.org/structure/6KR6) Raw data: [Zenodo#3603539](https://doi.org/10.5281/zenodo.3603539)
1. Yoshizawa *et al.* (2020) "Crystal structures of the cell-division protein FtsZ from Klebsiella pneumoniae and Escherichia coli" *Acta Cryst. F* doi:[10.1107/S2053230X2000076X](https://doi.org/10.1107/S2053230X2000076X) PDB: [6LL5](https://www.rcsb.org/structure/6LL5) [6LL6](https://www.rcsb.org/structure/6LL6)
1. Izume *et al.* (2020) "Crystal structure of human endothelin ET<sub>B</sub> receptor in complex with sarafotoxin S6b" *Biochemical and Biophysical Research Communications* doi:[10.1016/j.bbrc.2019.12.091](https://doi.org/10.1016/j.bbrc.2019.12.091) PDB: [6LRY](https://www.rcsb.org/structure/6LRY) Raw data: [Zenodo#3603541](https://doi.org/10.5281/zenodo.3603541)
1. Asada *et al.* (2019) "The Crystal Structure of Angiotensin II Type 2 Receptor with Endogenous Peptide Hormone" *Structure* doi:[10.1016/j.str.2019.12.003](https://doi.org/10.1016/j.str.2019.12.003) PDB: [6JOD](https://www.rcsb.org/structure/6JOD)
1. Sugishima *et al.* (2019) "Crystal structure of phytochromobilin synthase in complex with biliverdin IXα, a key enzyme in the biosynthesis of phytochrome" *Journal of Biological Chemistry* doi:[10.1074/jbc.RA119.011431](https://doi.org/10.1074/jbc.RA119.011431) PDB: [6KME](https://www.rcsb.org/structure/6KME) [6KMD](https://www.rcsb.org/structure/6KMD)
1. Vuckovic *et al.* (2019) "Crystal structure of the M5 muscarinic acetylcholine receptor" *PNAS* doi:[10.1073/pnas.1914446116](https://doi.org/10.1073/pnas.1914446116) PDB: [6OL9](https://www.rcsb.org/structure/6OL9)
1. Shihoya *et al.* (2019) "Crystal structure of heliorhodopsin" *Nature* doi:[10.1038/s41586-019-1604-6](https://doi.org/10.1038/s41586-019-1604-6) PDB: [6IS6](https://www.rcsb.org/structure/6IS6) Raw data: [Zenodo#3333323](https://zenodo.org/record/3333323)
1. Liu *et al.* (2019) "Mechanism of β<sub>2</sub>AR regulation by an intracellular positive allosteric modulator" *Science* doi: [10.1126/science.aaw8981](https://doi.org/10.1126/science.aaw8981) PDB: [6N48](http://www.rcsb.org/pdb/explore/explore.do?structureId=6N48)
1. Nagiri *et al.* (2019) "Crystal structure of human endothelin ET<sub>B</sub> receptor in complex with peptide inverse agonist IRL2500" *Communications Biology* doi: [10.1038/s42003-019-0482-7](https://doi.org/10.1038/s42003-019-0482-7) PDB: [6K1Q](http://www.rcsb.org/pdb/explore/explore.do?structureId=6K1Q) Raw data: [Zenodo#2803553](https://zenodo.org/record/2803553)
1. Liu *et al.* (2019) "Structural Insights into the Process of GPCR-G Protein Complex Formation" *Cell* doi: [10.1016/j.cell.2019.04.021](https://doi.org/10.1016/j.cell.2019.04.021) PDB: [6E67](http://www.rcsb.org/pdb/explore/explore.do?structureId=6E67) [6EG8](http://www.rcsb.org/pdb/explore/explore.do?structureId=6EG8)
1. Nagamura *et al.* (2019) "Structural basis for oligomerization of the prokaryotic peptide transporter PepT<sub>So2</sub>." *Acta Cryst. F*  doi: [10.1107/S2053230X19003546](https://doi.org/10.1107/S2053230X19003546) PDB: [6JKC](http://www.rcsb.org/pdb/explore/explore.do?structureId=6JKC)  [6JKD](http://www.rcsb.org/pdb/explore/explore.do?structureId=6JKD) Raw data: [Zenodo#2533841](https://zenodo.org/record/2533841)
1. Inoue *et al.* (2019) "Structural Basis of Sarco/Endoplasmic Reticulum Ca<sup>2+</sup>-ATPase 2b Regulation via Transmembrane Helix Interplay." *Cell Reports* doi: [10.1016/j.celrep.2019.03.106](https://doi.org/10.1016/j.celrep.2019.03.106) PDB: [5ZTF](http://www.rcsb.org/pdb/explore/explore.do?structureId=5ZTF)
1. Hashimoto *et al.* (2019) "Protein encapsulation in the hollow space of hemocyanin crystals containing a covalently conjugated ligand." *Biochemical and Biophysical Research Communications* doi: [10.1016/j.bbrc.2019.04.062](https://doi.org/10.1016/j.bbrc.2019.04.062)
1. Umeda *et al.* (2019) "Crystallization of the human tetraspanin protein CD9." *Acta Cryst. F*  doi: [10.1107/S2053230X1801840X](https://doi.org/10.1107/S2053230X1801840X)
1. Kato *et al.* (2019) "Crystal structure of plant vacuolar iron transporter VIT1." *Nature Plants* doi:[10.1038/s41477-019-0367-2](https://doi.org/10.1038/s41477-019-0367-2) PDB: [6IU3](https://www.rcsb.org/structure/6IU3) [6IU4](https://www.rcsb.org/structure/6IU4) [6IU5](https://www.rcsb.org/structure/6IU5) [6IU6](https://www.rcsb.org/structure/6IU6) [6IU8](https://www.rcsb.org/structure/6IU8) [6IU9](https://www.rcsb.org/structure/6IU9) Raw data: [Zenodo#2532136](https://zenodo.org/record/2532136) [Zenodo#2532134](https://zenodo.org/record/2532134) [Zenodo#2532138](https://zenodo.org/record/2532138)
1. Terakado-Kimura *et al.* (2019) "Structures of the 5-HT<sub>2A</sub> receptor in complex with the antipsychotics risperidone and zotepine." *Nature Structural & Molecular Biology* doi:[10.1038/s41594-018-0180-z](https://doi.org/10.1038/s41594-018-0180-z) PDB: [6A93](http://www.rcsb.org/pdb/explore/explore.do?structureId=6A93) [6A94](http://www.rcsb.org/pdb/explore/explore.do?structureId=6A94)
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

## Version history
Dates when the code became available on GitHub are shown

* 2020-01-06
   * Fix for XDS BUILT=20191015 where SNRC= parameter was introduced and MINIMUM\_I/SIGMA= was obsoleted in XSCALE.
   * kamo.multi\_prep\_merging: dispatched as a command (the same function as Multi-merge button in KAMO GUI)
   * KAMO: bug fix for split\_data\_by\_deg option when multi-trigger data
* 2019-10-05
   * kamo.auto\_multi\_merge: fix for filtering.choice=cell (not actually excluded unless there is indexing ambiguity)
   * KAMO: tentative support for Eiger2 at SLS (need use\_dxtbx=true and dials 1.14.12 or later)
* 2019-09-05
   * KAMO: change DATA\_RANGE= if first frames have zero or too weak counts (for a rare trouble of EIGER)
   * KAMO: dataset\_paths\_txt mode as default. do not shorten frame numbers automatically. do not wait for all data h5 files if not needed.
   * kamo.multi\_merge: safe guard for figure size when cc\_clustering
   * kamo.auto\_multi\_merge: select group based on reference symmetry
* 2019-06-17
   * kamo.multi\_merge: bug fix for rejection by bfactor (reported by Dr. Shimamura)
   * kamo.auto\_multi\_merge: add cell-based filtering option (filtering.choice=cell)
   * KAMO: allow comment lines (starting with #) in dataset\_paths\_txt=.
   * KAMO: For XDS with EIGER or PILATUS, set SEPMIN=4 CLUSTER\_RADIUS=2 by default.
* 2019-05-20
   * kamo.test\_installation: support new xdsstat
   * KAMO: bug fix in plotting spot numbers, support of scan\_varying (DIALS), support of online use at BL45XU.
* 2019-03-02
   * kamo.auto\_multi\_merge: add option to specify preferred unit cell and symmetry
   * kamo.test\_installation: changed h5 reading test due to [cctbx#282](https://github.com/cctbx/cctbx_project/issues/282)
* 2018-07-30
   * KAMO: dataset\_paths\_txt= option was added for online use at non-SPring-8 sites
* 2018-05-22
   * kamo.test\_installation: fixed a bug introduced in 2018-04-25 (did not work before Numpy 1.11)
   * Bug fixes in file preparation for merging
* 2018-04-25
   * Updated this doc with published KAMO paper
   * kamo.test\_installation: updated H5ToXds test by actually running a program
* 2018-02-22
   * kamo.multi\_merge: added an option cc\_clustering.min\_common\_refs= to set the minimum acceptable number of common reflections.
   * Improved least-square fitting for CC<sub>1/2</sub> vs resolution curve.
   * Bug fix in LCV calculation
* 2018-01-30
   * KAMO: supporting XDS Nov 11, 2017 (BUILT=20171218) that does not create GXPARM.XDS when refinement unsuccessful.
   * kamo.multi\_prep\_merging: sort filenames in formerge.lst etc.
   * kamo.multi\_merge: better error handling and setting minimum figure size in cc\_clustering.
* 2017-12-26
   * KAMO: default rejection parameter in merging script; no rejection is based on b+B. bug fix in filter\_cell.R.
   * kamo.multi\_merge: in cc\_clustering now don't use R, but SciPy (>=0.18.1 is required) for cluster analysis. rejects data if the number of common reflections is less than 3. now can choose cc\_to\_distance functions and linkage methods.
   * kamo.multi\_merge: `resolution.estimate=` is now True by default; changed behaviour of this option to just estimate the cutoff based on curve-fitting of CC<sub>1/2</sub> vs d<sup>-2</sup>, not actually cut.
* 2017-12-03
   * kamo.multi\_merge cc_clustering: severe bug fix when unused data files existed (typical when low space group symmetry). Additionally note that we virtually (unintentionally) did clustering using sqrt(1-cc) with ward.D2 method, which was actually an appropriate way because sqrt(2-2CC)) can be used as a "distance".
   * kamo.multi\_merge: proper error handling and report html fix in anisotropic analysis.
   * kamo.auto\_multi\_merge: postrefine= parameter is now obsolete and cell\_method= should be used. code refactoring.
   * KAMO: user-friendlly message when batch.engine=sge and qsub is not available. batch.sh\_max\_jobs=Auto is now default.
   * KAMO: bug fix in logger registration (exception was ignored before)
* 2017-11-02
   * KAMO: method changed to choose space group (by default use given symmetry for scaling. if not given, choose pointless result with higher probability). sanity check of given symmetry when GUI started. bug fix in html report
   * kamo.multi\_merge: calculate frame\_cc using xscale.hkl instead of original files
   * kamo.multi\_merge: fixed a bug when pointless gave space group in non-reference setting
   * kamo.auto\_multi\_merge: bug fix (did not work when batch.engine=sh), use finer step in resolution determination
   * kamo.resolve\_indexing\_ambiguity: show warning if no reflections left and added option to ignore them
* 2017-10-12
   * kamo: show warning if XDS is expired
   * kamo.multi\_merge: now d\_max= works
   * kamo.auto\_multi\_merge: available now (missed commit)
* 2017-09-30
   * kamo.multi\_merge: Fixed a bug in reading completeness and R-values from table of XSCALE.LP (numbers ignored after decimal point); so those numbers in summary and report html were a bit wrong (sorry)
   * kamo.multi\_merge: Pointless speed-up by giving data merged in P1 instead of xscale.hkl directly.
   * KAMO: add correct\_only choice to known.method= option, which uses given symmetry information in CORRECT only (hopefully useful when xds.repeat>1).
   * kamo.resolve\_indexing\_ambiguity: changed logging format in reference-based mode (now same style as selective-breeding)
* 2017-08-31
   * kamo.multi\_merge: show anisotropic resolution cutoffs in HTML report.
   * kamo.multi\_merge: add degrees\_per\_batch= parameter (degrees version of frames\_per\_batch=).
   * kamo.multi\_merge: show minimum and average value of CC when cc\_clustering.
   * kamo.auto\_multi\_merge: root\_dir now can be omitted in CSV file (specify datadir= for global value)
   * kamo.auto\_multi\_merge: select best result from final runs only. use outershell CC1/2 for decision as well as overall CC1/2.
* 2017-07-22
   * Support MarCCD's no-extension format like foo.0001
* 2017-07-20
   * (new) kamo.multi\_determine\_symmetry: determines space group (point group only) from multiple (small wedge) datasets
   * KAMO: new option known.method= to specify how prior cell knowledge is used (default is not\_use\_first that tries indexing without prior knowledge first; use\_first is to use it first)
   * KAMO: new option on Multi-Merge GUI to copy transformed HKL files into working directory for merging (ON by default)
   * KAMO: now automatically recognizes beam center convention of ADSC detectors and goniometer rotation axis using header information.
   * KAMO: (experimental) new use\_dxtbx= option; when true, dxtbx is used to recognize beam, goniometer, and detector geometries.
   * KAMO: now saves log file in working directory. log\_root= option is not mandatory.
   * KAMO: less stressful GUI (job status acquisition not blocking GUI; grouping datasets using unit cells in the background)
   * kamo.resolve\_indexing\_ambiguity: selective breeding now supports parallel computation (nproc>1)
   * Now basically use max\_delta=5 by default in all programs.
   * filter\_cell.R creates plot
* 2017-05-24
   * Fixed a bug introduced on 2017-03-10 (avoid failure in copying environment variables).
* 2017-05-23
   * KAMO: Add DIALS support (`engine=dials`).
   * KAMO: Multiprocessing of preparation of multi-merge. New option to prepare files for joint refinement using DIALS (requires KAMO built on DIALS environment)
   * kamo.multi_merge: bug fix in anisotropy analysis.
* 2017-04-19
   * KAMO: Added support of eiger2cbf converted CBF files. For online processing at beamline, wait until all Eiger h5 files downloaded, and added support of miniset.
* 2017-03-24
   * kamo.test_installation: Fixed a problem on testing XDS. Add H5ToXds test.
   * KAMO: When multi-merge strategy started, read unit cell in P1 from CORRECT.LP_noscale instead of CORRECT.LP.
* 2017-03-10
   * KAMO: Fixed a Mac (El Capitan or later) specific problem where script couldn't run in local computer.
* 2017-02-16
   * yamtbx.xds_aniso_analysis: New program to perform anisotropy analysis (CC1/2 and I/sigma) for XDS unmerged data (executed in kamo.multi_merge)
   * KAMO: Fixed a silly bug in nproc-based determination of DELPHI= in XDS.
   * kamo.resolve_indexing_ambiguity: fixed a bug in selective-breeding (sometimes failed maybe when small number of files?)
   * kamo.multi_merge: when add_test_flag=true, first generate test set and copy them to all
* 2017-02-02
   * kamo.auto_multi_merge: automatic merging for multiple samples
   * kamo.multi_merge: new options reference.data= (to copy test flags), resolution.estimate= (to automatically decide high resolution cutoff), use pointless to decide screws
   * kamo: now blconfig= can be multiple and mode= can specify both (zoo+normal)
* 2017-01-18
   * yamtbx.beam_direction_plot: fixed a bug in non-primitive space group case
* 2016-12-26
   * kamo.multi_merge: add `space_group=` option (used in merging). use pointless result for mtz if not specified.
   * kamo.multi_merge: add MULTIPLICITY column in mtz
   * bug fix (change directory in qsub script)
* 2016-12-06
   * GUI: add `exclude_ice_resolutions=` option, fixed a bug that plots were not updated on Mac
   * faster string (file name) replacement for XSCALE outputs
   * kamo.resolve_indexing_ambiguity: fixed a bug when no reference_label was given
   * kamo.test_installation: add Adxv test
* 2016-10-05
   * added `auto_frame_exclude_spot_based=` option, which could be useful for processing data including non-spots images
* 2016-07-18
   * use ramdisk/tmpdir for xds/xscale run
   * calculate frequency of crystal symmetry taking unit cell parameters into account
   * bug fix for OSX (with phenix-1.10.1?)
   * bug fix on html report making
   * bug fix for non-sge environment (multi_merge)
   * calculate LCV & aLCV for actual set of parameters
   * bug fix on parsing xtriage. anisotropy is now max(B_cart)-min(B_cart)
