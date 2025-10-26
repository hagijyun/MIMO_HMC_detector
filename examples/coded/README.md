# Prerequisites
## Path
Add the following to your system PATH:
- Rscript.exe
- Rtools (make, gcc, etc.)
- matlab.exe

## Working directory
The following files contain descriptions assuming e:/ as the working directory, so please modify as appropriate:
- _quick_examples.ps1
- MIMO5_para.Rmd
- send_receive2_para.m
- mean_numErrs.m

# Usage
Although execution takes some time, progress can be monitored via the log file parallel_log.txt.

1. Download all files to the working directory.
2. Move 3GPP_LDPC_BG1.Z96 and send_receive2_para.m to MATLAB's default directory when MATLAB starts (on Windows, this is C:\Users\username\Documents\MATLAB).
3. Run WSH with administrator privileges and navigate to the working directory.
4. Execute by running `.\_quick_examples.ps1 1 5`.
5. Run mean_numErrs.m in MATLAB. The initial BER will be saved to numErrs_mean0.csv, and the BER for subsequent iterations will be saved to numErrs_mean1to5.csv.
