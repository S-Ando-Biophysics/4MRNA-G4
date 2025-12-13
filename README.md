# 4MRNA-G4 (beta)
4MRNA-G4 consists of two steps.
The first step is to create many diverse models by varying the way the G-quartet layers stack.
The second step is to perform molecular replacement using the models created.

You can create models on the website and download the code to run the molecular replacement.
Website: https://s-ando-biophysics.github.io/4MRNA-G4/


## Instructions
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
- **2025-11-11** The beta version has been released.

## Reference
- in preparation
