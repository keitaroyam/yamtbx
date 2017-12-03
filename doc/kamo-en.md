# KAMO

Japanese version is [here](kamo-ja.md).

## Overview
*KAMO (in Japanese: Katappashikara Atsumeta data wo Manual yorimoiikanjide Okaeshisuru; to process thoroughly-collected data in better way than manually and give back to user) system* is the program for automated processing and merging of (MX) crystal diffraction data.
Basically KAMO uses [XDS package](http://xds.mpimf-heidelberg.mpg.de), but it supports [DIALS](https://dials.github.io/) and [Aimless](http://www.ccp4.ac.uk/html/aimless.html) optionally.
Currently, the development is focused on processing and merging of small wedge data (5-10°/crystal).

KAMO is designed for online analysis for SPring-8 MX beamlines; however, it can be used offline for local data. Please let me know if you need support for your beamline/detector.

This manual is for 2015-12-18.

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
      * [KAMO](#kamo)
         * [I'm too lazy to click all checkboxes](#im-too-lazy-to-click-all-checkboxes)
         * [I want to use prior unit cell information](#i-want-to-use-prior-unit-cell-information)
         * [I want to give the correct experimental parameters as image header is wrong](#i-want-to-give-the-correct-experimental-parameters-as-image-header-is-wrong)
      * [kamo.multi_merge](#kamomulti_merge)
         * [Where is the data for structural refinement?](#where-is-the-data-for-structural-refinement)
   * [How can I use KAMO at home?](#how-can-i-use-kamo-at-home)
      * [Installation](#installation)
      * [How to update KAMO](#how-to-update-kamo)
      * [Launch](#launch)
   * [Citations](#citations)
      * [How to cite the use of KAMO](#how-to-cite-the-use-of-kamo)
      * [Researches which used KAMO](#researches-which-used-kamo)
   * [Version history](#version-history)


### Dependencies
KAMO uses following programs and libraries.

* [CCTBX](http://cctbx.sourceforge.net/) with [CBFlib](http://www.bernstein-plus-sons.com/software/CBF/) (essential)
* [wxPython 2.8](http://www.wxpython.org/), [Matplotlib 1.3](http://matplotlib.org/), [Networkx 1.x](https://networkx.github.io/), [Numpy](http://www.numpy.org/), [SciPy](https://www.scipy.org/) (essential)
* [XDS](http://xds.mpimf-heidelberg.mpg.de), [xdsstat](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Xdsstat), H5ToXds (for EIGER data)
* [CCP4](http://www.ccp4.ac.uk/) (BLEND, Pointless, Aimless, Ctruncate)
* [R](https://www.r-project.org/) (required for BLEND, CC-based clustering) with rjson
* [DIALS](https://dials.github.io/) (not fully-supported yet)


### Warning
KAMO is still under development and may have many many bugs on interface or functions. If you had any troubles, please do not hesitate to contact author.


## Usage
### GUI
The data in the directory where GUI started will be all processed (includes sub-directories).

If you want to process data in specific sub-directory, give the directory names like `include_dir=hoge-1111-\*` or give the list file containing the directory names.

* case 1: online analysis of small wedge at BL32XU

	kamo bl=32xu
* case 2: online analysis of normal data (1 or a few crystals for complete data) at BL41XU

	kamo bl=41xu small_wedges=false
* case 3: online analysis of data collected with Zoo (automatic system) at BL32XU

	kamo bl=32xu mode=zoo
* case 4: offline analysis of local data (KAMO will find directories. In above cases, datasets will be found using BSS log files)

	kamo bl=other

Use `-h` option to see help and list of all parameters.

For processing old data, for example, give `date=2015-12-31` to find data collected from the specified date to two days before (by default `date=today`).

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
SG | Space group
Resn. | Estimated resolution limit (Based on peak search result when small_wedges=true; based on I/sigma<1 basis using the 100 resolution-bins table in CORRECT.LP)


### Merging small wedge data
1. Select all target data using the checkboxes (no problem if selecting failed ones), click `Multi merge strategy` button, and wait for a while
2. The data with the same lattice (reindexing is taken into account) are grouped, and sorted by the number of datasets
3. Select the group to be merged and its symmetry. The most frequent symmetry is chosen by default. Choose the correct one if known.
5. Click Proceed button, and look at the terminal. The datasets processed with different symmetry are reprocessed with the specified symmetry.
6. Completed when "Do It Yourself!" appeared. Keep in mind the reindex operator(s) if existed. Indexing ambiguity can be resolved using `kamo.resolve_indexing_ambiguity`.
7. Go to the working directory for merging, and modify and start the script. The scrip (merge_blend.sh) is automatically prepared like this:
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
8. When this script is started, first the hierarchical clustering analysis by BLEND is performed. The clusters with high completeness (>90%) are subjected to merging. First simply run xscale to merge (run_01/), and find the bad frames based on the CC between the merged result and intensities on each frame. From the merging result after excluding the bad frames (run_02/), bad datasets are detected based on error model's *b* value, and the merging result without bad datasets is saved in run_03/. This is the final result. In each run_\*/ccp4 mtz and the log files of ctruncate and phenix.xtirage are saved.
9. When processing is over, report.html in the working directory (blend_\*/) can be viewed with web browser. All final results can be visually inspected. Based on the results, re-execute the script by changing high resolution limi if needed. For refinement, empirically, the result with the largest CC1/2 should be appropriate (we welcome your feedback).

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

Concerning the decision of symmetry, first INTEGRATE.HKL is analyzed by pointless, and then XDS_ASCII.HKL is analyzed too. If results are not consistent, the one with larger ISa is employed.

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
Use `cc_clustering.d_min=` to limit the resolution. Currently we are using hclust() function in R and only datasets which have common reflections with all others can be used. Probably not useful for low-symmetry crystals.

#### Bad frame detection

Correlation coefficient is calculated between merged all-data and each frame. Based on the CC, the bad frames are discarded (when `reject_method=framecc`).
By default `rejection.framecc.method=tukey rejection.framecc.iqr_coeff=1.5`, which detects outliers by Tukey's criterion (1.5\*IQR). You can change this; for example `rejection.framecc.method=abs rejection.framecc.abs_cutoff=0.9` to use absolute threshold.

#### Bad dataset detection
Bad datasets are detected using statistics in XSCALE.LP (`reject_method=lpstats`).
By default, error model's *b* value is subjected to Tukey's outlier detection (`rejection.lpstats.stats=em.b rejection.lpstats.iqr_coeff=1.5`). You can give `pairwise_cc` and `bfactor` as well as `em.b` to `rejection.lpstats.stats=`. `pairwise_cc` is to remove datasets which give bad CC (by default CC<0.8; `rejection.lpstats.pwcc.method=abs rejection.lpstats.pwcc.abs_cutoff=0.8`, but optionally Tukey's method can be used).

#### Choice of scaling reference
This affects overall B-factor of merged data. In XSCALE, the first INPUT_FILE= is used as reference. In KAMO, by default `xscale.reference=bmin`, which selects data with smallest *B* (that has smallest intensity fall-off w.r.t. resolution in XDS) as reference.


## FAQ
### KAMO
##### I'm too lazy to click all checkboxes
If you want to check all, click "Check all" button.
Alternatively, click two checkboxes keeping Shift-key pressed to check all in-between items.

#### I want to use prior unit cell information
You can give it when KAMO starts:

```
kamo known.unit_cell=10,20,30,90,90,90 known.space_group=p222
```
Please make sure to give pace_group as well.
The unit cell parameters are used in indexing as prior information. If not indexed with the cell, the data processing will not proceed.
Note that you cannot use this parameter when you have more than one kind of crystals.

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

## How can I use KAMO at home?

You can easily install KAMO using DIALS/PHENIX environment as DIALS/PHENIX includes CCTBX and its dependencies (No need to install CCTBX by yourself).

### Installation

1. Install CCP4, R (with rjson package), XDS
   * For installation of XDS/XDSSTAT, see [XDSwiki/Installation](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Installation)
   * If you will process EIGER data (h5 files), [H5ToXds](eiger-en.md#eiger2cbf-h5toxds-compatible) is needed
   * rjson can be installed as follows; after installation of R, start R program from user who installed R (root or an account for software installation), and then type `install.packages("rjson")`.
2. Install [DIALS](https://dials.github.io/installation.html)-1.5 or newer
3. Install networkx to dials.python
   1. `cd $DIALS/build`
   2. `./bin/libtbx.python -m easy_install networkx==1.11`
4. Install scipy to dials.python
   1. If Mac, install [gfortran](http://gcc.gnu.org/wiki/GFortranBinaries#MacOS). If Linux, install blas-devel and lapack-devel using yum or something.
   2. `cd $DIALS/build`
   3. `./bin/libtbx.python -m easy_install scipy==0.18.1`
5. Run the following commands (yamtbx can be cloned anywhere you like)
```
cd $HOME
git clone https://github.com/keitaroyam/yamtbx.git
cd $DIALS/modules
ln -s ~/yamtbx/yamtbx .
cd ../build
./bin/libtbx.configure yamtbx
```

After installation, run
```
kamo.test_installation
```
to check if dependencies are all installed.

##### Troubleshooting tips
* Installation of scipy stops with error "as: I don't understand 'm' flag!"
   * If you are using MacPorts, try excluding /opt/local/bin from environment variable PATH. [Reference URL](https://stackoverflow.com/questions/41542990/while-installing-on-osx-sierra-via-gcc-6-keep-having-fatal-opt-local-bin-l).
* Installing with DIALS/PHENIX environment, but kamo.test\_installation claims wxPython is NG
   * Please check if some GUI program of DIALS/PHENIX successfully starts (e.g. dials.image\_viewer)．In case of Ubuntu, you may need to install libjpeg62 package.

### How to update KAMO
1. `cd` where-you-cloned-yamtbx
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

## Citations

### How to cite the use of KAMO

As the paper is in preparation, please refer to this documentation URL: https://github.com/keitaroyam/yamtbx/blob/master/doc/kamo-en.md

### Researches which used KAMO 

* Miyauchi *et al.* (2017) "Structural basis for xenobiotic extrusion by eukaryotic MATE transporter." *Nature Communications* doi: [10.1038/s41467-017-01541-0](https://doi.org/10.1038/s41467-017-01541-0) PDB: [5Y50](http://www.rcsb.org/pdb/explore/explore.do?structureId=5Y50) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-AtDTX14-data-(5Y50))
* Lee *et al.* (2017) "Structure of the triose-phosphate/phosphate translocator reveals the basis of substrate specificity." *Nature Plants* doi: [10.1038/s41477-017-0022-8](https://doi.org/10.1038/s41477-017-0022-8) PDB: [5Y78](http://www.rcsb.org/pdb/explore/explore.do?structureId=5Y78) [5Y79](http://www.rcsb.org/pdb/explore/explore.do?structureId=5Y79) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-TPT-data-(5Y78-&-5Y79))
* Tanaka *et al.* (2017) "Crystal Structure of a Plant Multidrug and Toxic Compound Extrusion Family Protein." *Structure* doi: [10.1016/j.str.2017.07.009](https://doi.org/10.1016/j.str.2017.07.009) PDB: [5XJJ](http://www.rcsb.org/pdb/explore/explore.do?structureId=5XJJ)
* Shihoya *et al.* (2017) "X-ray structures of endothelin ET<sub>B</sub> receptor bound to clinical antagonist bosentan and its analog." *Nature Structural & Molecular Biology* doi: [10.1038/nsmb.3450](https://doi.org/10.1038/nsmb.3450) PDB: [5XPR](http://www.rcsb.org/pdb/explore/explore.do?structureId=5XPR) [5X93](http://www.rcsb.org/pdb/explore/explore.do?structureId=5X93) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-ETBR-bonsentan-data-(5XPR))
* Taniguchi *et al.* (2017) "Structural insights into ligand recognition by the lysophosphatidic acid receptor LPA<sub>6</sub>." *Nature* doi: [10.1038/nature23448](https://doi.org/10.1038/nature23448) PDB: [5XSZ](http://www.rcsb.org/pdb/explore/explore.do?structureId=5XSZ) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-LPA6-data-(5XSZ))
* Abe *et al.* (2017) "Crystal Engineering of Self-Assembled Porous Protein Materials in Living Cells." *ACS Nano* doi: [10.1021/acsnano.6b06099](https://doi.org/10.1021/acsnano.6b06099) PDB: [5GQM](http://www.rcsb.org/pdb/explore/explore.do?structureId=5GQM) [5GQN](http://www.rcsb.org/pdb/explore/explore.do?structureId=5GQN) Raw data and processing note: [link](https://github.com/keitaroyam/yamtbx/wiki/Processing-Polyhedra-data-(5GQM-&-5GQN))

## Version history
Dates when the code became available on GitHub are shown

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
