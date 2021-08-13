          seed =  -1

       seqfile = h_mel_tim_num_noncod.txt
      Imapfile = h_mel_tim_num_imap.txt
       outfile = h_mel_tim_num_noncod_out_L500_small.txt
      mcmcfile = h_mel_tim_num_noncod_mcmc_L500_small.txt

  speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2    * species delimitation rjMCMC algorithm0 and finetune(e)
* speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
*         speciestree = 1 0.1 0.1 0.2    * species tree SPR/SNL

*   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

species&tree = 3 M T N
                 1 1 1
             (((M, Y[&phi=0.15])X, (X[&phi=0.2], T) Y) S, N) R;

         phase = 1 1 1
       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 500 * 31166 autosomal noncoding loci; 2592 on chr1

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = gamma 2 400  # gamma(a, b) for theta (integrated out by default; add E to also sample theta)
      tauprior = gamma 2 400  # gamma(a, b) for root tau
    phiprior = 1 1  # Beta(a, b) for root tau & Dirichlet(a) for other tau's

      finetune =  1: 3 0.002 0.001 0.0001 0.05 0.9 0.054 0.2 0.3 0.4 # finetune for GBtj, GBspr, theta, tau, mix

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 8000
      sampfreq = 5
       nsample = 20000
       threads = 16 1 1
