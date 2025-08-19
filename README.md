# Model Builder for 4MRNA-G4
You can access this website via the following URL. https://s-ando-biophysics.github.io/4MRNA-G4/

The explanation of what 4MRNA-G4 is is provided after the instructions below.

## Announcement
On **August 20, 2025**, the beta version has been released. 

## Instructions

### User Manual
Please refer to the user manual.

- English Version (in preparation)
- Japanese Version (in preparation) 

### How to Use (Summarized Version of the User Manual)

in preparation

### Preparation
Please install and set up the following software in advance.

| Name | Remarks |
| :----- | :----- |
| [Ubuntu](https://apps.microsoft.com/search?query=Ubuntu) | This is only required on Windows. In addition, it is necessary to turn on "**Windows Subsystem for Linux (WSL)**" and **Virtual Machine Platform** in the Windows settings to be able to use shell scripts. |
| [Phaser](https://www.ccp4.ac.uk/download) | The command "phaser" is included in CCP4. If you have not installed CCP4, please install it. [*] |

[*] The following steps are for Windows (Ubuntu). The procedure for macOS is similar.

    # Please change the directory name and CCP4 version as appropriate.
    # Assume that "ccp4-9.0.010-linux64.tar.gz" has been downloaded to "C:\Users\name\Downloads".
    sudo su
    cd /usr/local
    mv /mnt/c/Users/name/Downloads/ccp4-9.0.010-linux64.tar.gz .
    gunzip ccp4-9.0.010-linux64.tar.gz
    tar -xvf ccp4-9.0.010-linux64.tar
    cd ccp4-9
    ./BINARY.setup
    exit
    cd /home/name
    vi .bashrc
    source /usr/local/ccp4-9/bin/ccp4.setup-sh    # Please add to the last line.


### Supported Environment

|  | Operating system | Browser |
| :----- | :----- | :----- |
| **This website** | Windows | Google Chrome, Microsoft Edge, Firefox, Safari |
| **Generated Shell scripts** | Windows (WSL; Ubuntu) | - |

It is expected to work on macOS and Linux (Rocky Linux), but verification is currently in progress.

### Release Notes
- **2025-08-20** The beta version has been released.

## About 4MRNA-G4
**4MRNA** means <ins>**M**</ins>assive <ins>**M**</ins>ulti-type <ins>**M**</ins>odel <ins>**M**</ins>olecular <ins>**R**</ins>eplacement for <ins>**N**</ins>ucleic <ins>**A**</ins>cids. Please refer to [this page](https://github.com/S-Ando-Biophysics/4MRNA/tree/main?tab=readme-ov-file#about-4mrna).

While 4MRNA targets nucleic acid **duplex** structures, 4MRNA-G4 targets **G-quadruplex (G4)** structures. In other words, **4MRNA-G4** is the G4 version of 4MRNA.

== Further explanation is currently being written. ==

## Reference
- in preparation
