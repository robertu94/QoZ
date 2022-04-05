QoZ: Dynamic Quality Metric Oriented Error Bounded Lossy Compression for Scientific Datasets

=====
(C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory. See COPYRIGHT in top-level directory.

* Major Authors: Sheng Di, Jinyang Liu, Kai Zhao, Xin Liang
* Supervisor: Franck Cappello
* Other Contributors: Robert Underwood, Sihuan Li, Ali M. Gok

## Introduction

This is the source code for SC 22' submisssion: Dynamic Quality Metric Oriented Error Bounded Lossy Compression for Scientific Datasets

## Dependencies

Please Installing the following dependencies before running the artiact evaluation experiments:

* Python >= 3.6
* numpy 
* pandas 
* qcat (from https://github.com/Meso272/qcat, check its readme for installation guides. Make sure the following executables are successfully installed: calculateSSIM and computeErrAutoCorrelation)

## 3rd party libraries/tools

* Zstandard (https://facebook.github.io/zstd/). Zstandard v1.4.5 is included and will be used if libzstd can not be found by
  pkg-config.

## Installation

* mkdir build && cd build
* cmake -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR] ..
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and header files in [INSTALL_DIR]/include

## Single compression/decompression testing Examples

You can use the executable 'qoz' command to do the compression/decompression. Just run "qoz" command to check the instructions for its arguments.
Currently you need to add a configuration file to the argument line (-c) for activating the new features of QoZ. 
The corresponding cofiguration files for each test dataset can be generated by generate_config.py (details shown on following).

## Artifact Evaluation guides

Step 1: Download the dataset from the following links,then unzip them:

* CESM-ATM: https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data/CESM-ATM/SDRBENCH-CESM-ATM-cleared-1800x3600.tar.gz
* Miranda: https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data/S3D/S3D.tar.gz
* Hurricane-ISABEL: https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data/Hurricane-ISABEL/SDRBENCH-Hurricane-ISABEL-100x500x500.tar.gz
* NYX: https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data/EXASKY/NYX/SDRBENCH-EXASKY-NYX-512x512x512.tar.gz
* SCALE-LETKF: https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data/SCALE_LETKF/SDRBENCH-SCALE-98x1200x1200.tar.gz

Step 2: Preprocess the downloaded data with the preprocess_data.py:

* Usage: python preprocess_data -m {MirandaPath} (MirandaPath is the folder of the Miranda Dataset)

Step 3: Run the generate_config.py to generate the configuration files for QoZ:

* python generate_config.py

After that, the configuration files for QoZ test will generated in configs/ folder in the name format of {DatasetName}\_{Target}.config. 

* DatasetName: The name of dataset (cesm, miranda, hurricane, nyx, scale)
* Target: The optimization target, includes cr (compression ratio), psnr (PSNR), ssim (SSIM), ac (Autocorrelation). For each mode of the QoZ, use the configuration file of the corresponding target.

Step 4: Run test_qoz.py to generate the test results.

* Command: test_qoz.py -i {Datapath} -o {OutputPath} -d {DatasetName} -t {Target}
* Datapath: the folder path of the dataset
* OutputPath: the output data file prefix. The output files will be in format of Outputpath_{Metric}.tsv
* DatasetName: See step 3
* Target: See step 3

The output will contain:
* overall_cr: the overall compression ratio under different vr-rel error bounds
* overall_psnr: the overall compression PSNR under different vr-rel error bounds
* cspeed: the compression speed under different vr-rel error bounds (the speed results in the paper are generated with psnr mode)
* dspeed: the decompression speed under different vr-rel error bounds (the speed results in the paper are generated with psnr mode)
* overall_ssim: the overall compression SSIM under different vr-rel error bounds (ssim target only)
* overall_ac: the overall compression Autocorrelation under different vr-rel error bounds (ac target only)

Tips for plotting the rate-distortions:

* Bit rate = 32 / compression ratio



