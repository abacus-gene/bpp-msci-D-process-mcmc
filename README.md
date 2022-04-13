# bpp-msci-D-process-mcmc

**(0)** 
The included data files are for the Heliconius examples used in the paper.  To run bpp4, using a command 
like the following

bpp --cfile bpp-D-noncod.ctl

This generates an MCMC sample file called h_mel_tim_num_noncod_mcmc_L500_small.txt.

**(1)** 
The program bpp-msci-D-process-mcmc is for processing the MCMC
sample file from bpp4 under the bidirectional introgression (BDI)
model to remove label switching.

To compile the program bpp-msci-D-process-mcmc.c, cd to the src/
folder and try the following

cd src

cl -Ox bpp-msci-D-process-mcmc.c minsub.c

cc -O2 -o bpp-msci-D-process-mcmc bpp-msci-D-process-mcmc.c minsub.c -lm

**(2)** 
To run the program, use the following command.  

 bpp-msci-D-process-mcmc <mcmc-sample-filename> <phi_X-column> <phi_Y-column> <theta_X-column> <theta_Y-column>

 bpp-msci-D-process-mcmc ../data/h_mel_tim_num_noncod_mcmc_L500_small.txt  12 13 6 7

The first command-line argument is the mcmc sample file.  This lists
the sampled parameters in a spreadsheet separated by spaces or tabs.
The four integers specify the columns in the mcmc sample file for the
four parameters involved in the unidentifiability: phi_X, phi_Y,
theta_X, theta_Y.  From the bpp output, column 0 has Gen for MCMC
generations.  When the program runs, it prints out the variable names
for those columns, like the following:

200000 records, 12 variables
phi_x phi_y theta_x theta_y:
phi_X phi_Y theta_4X theta_5Y

Check to make sure that the correct columns are identified and the
printout is correct.  The program creates three files, with
"_processed", "_processed_CoGN", "_processed_CoG0" in the names.
These can be read in Tracer or read by bpp4.5.0 or later to summarize the
posterior, as follows.

bpp --summary --cfile bpp-D-noncod.ctl

The post-processing algorithms have been implemented in bpp since version 4.5.0, so they are automatically run when you 
fit the bidirectional introgression (DBI) model or the double-DBI model.

##References

Yang Z, Flouri T. 2022. Estimation of cross-species introgression
rates using genomic data despite model unidentifiability. Molecular Biology
and Evolution, 10.1093/molbev/msac083
