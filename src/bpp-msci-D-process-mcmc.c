/* bpp-msci-D-process-mcmc.c
* Ziheng, 2021.7.12
* This processes MCMC samples from msci models D or DD, with one or two bidirectional
* introgression events between sister species.
*
* cl -Ox bpp-msci-D-process-mcmc.c minsub.c
* cc -O2 -o bpp-msci-D-process-mcmc bpp-msci-D-process-mcmc.c minsub.c -lm
*
*
* _D0DD1 = 0 (model D,  nmsci=4, nhyper=8):  2 towers so tower is 0 or 1.
*          1 (model DD, nmsci=8, nhyper=16): 4 towers so tower is 0 1 2 3.
* algorithm = 0:beta-gamma (MLE), 1:CoG-N, 2:coG-0.
          files:_processed, _processed_CoG_N, _processed_CoG_0
* sample[n_mcmc_sample*nmsci]: the original sample from the mcmc file, untransformed.
* p[]: the 4 or 8 msci parameters in one mcmc sample.
*
* nmsci  = 4 or 8: phi_x phi_y theta_x theta_y | phi_z phi_w theta_z theta_w
* nhyper = 8 or 16: p q for phi, a b for theta, in the same order as the msci parameters.
*/
#include "minsub.h"

#define msciD_data  0
/* 0: heliconius-D, 1:simulation-D, 2:simulation-DD */

static int _D0DD1;    /* 0 for model D and 1 for model DD */
static int n_mcmc_sample;
static double mp[8], vp[8], phi_sln[4], phi_sln1[4], theta_s[4], theta_sln[4];
extern int noisy;
extern double Small_Diff;

static char* algorithms[3] = { "beta-gamma", "CoG_N" , "CoG_0" };

int data_summaries(double sample[], int towers[], int nmsci);
void move_to_tower(int tower, double p[]);
double lnlike_msci(double x[], int np);
int compare_towers(int towers[], double sample[], double hyper[], int nmsci, int algorithm);

void fit_beta_moments(double pq[2], double m, double v) {
   double z = m * (1 - m) / v - 1;
   if (z <= 0) { printf("v > m(1-m) in fit_beta_moments");  z = 0.01; }
   pq[0] = m * z;  pq[1] = (1 - m) * z;
}
void fit_gamma_moments(double ab[2], double m, double v) {
   ab[0] = m * m / v;  ab[1] = m / v;
}
double ln_beta_ratio(double xnew, double x, double p, double q)
{
   double lnpd = 0;
   if (fabs(xnew - x) > 1e-200)
      lnpd = (p - 1) * log(xnew / x) + (q - 1) * log((1 - xnew) / (1 - x));
   return lnpd;
}
double ln_gamma_ratio(double xnew, double x, double a, double b)
{
   double lnpd = 0;
   if (fabs(xnew - x) > 1e-200)
      lnpd = (a - 1) * log(xnew / x) - b * (xnew - x);
   return lnpd;
}

int data_summaries(double sample[], int towers[], int nmsci)
{
   int i, j;
   double p[8], mold;

   for (i = 0; i < nmsci; i++) mp[i] = vp[i] = 0;
   for (i = 0; i < 4; i++) /* for both D and DD */
      phi_sln[i] = phi_sln1[i] = theta_s[i] = theta_sln[i] = 0;
   for (i = 0; i < n_mcmc_sample; i++) {
      for (j = 0; j < nmsci; j++)
         p[j] = sample[i * nmsci + j];   /* one p from the mcmc sample */
      if (towers[i])  move_to_tower(towers[i], p);
      for (j = 0; j < nmsci; j++) {
         mold = mp[j];
         mp[j] += (p[j] - mold) / (i + 1.0);
         vp[j] += (p[j] - mold) * (p[j] - mp[j]);
      }
      phi_sln[0] += log(p[0]);  phi_sln1[0] += log(1 - p[0]);     /* phi_x */
      phi_sln[1] += log(p[1]);  phi_sln1[1] += log(1 - p[1]);     /* phi_y */
      theta_s[0] += p[2];  theta_sln[0] += log(p[2]);             /* theta_x */
      theta_s[1] += p[3];  theta_sln[1] += log(p[3]);             /* theta_y */
      if (_D0DD1) {
         phi_sln[2] += log(p[4]);  phi_sln1[2] += log(1 - p[4]);  /* phi_z */
         phi_sln[3] += log(p[5]);  phi_sln1[3] += log(1 - p[5]);  /* phi_w */
         theta_s[2] += p[6];  theta_sln[2] += log(p[6]);          /* theta_z */
         theta_s[3] += p[7];  theta_sln[3] += log(p[7]);          /* theta_w */
      }
   }
   abyx(1.0 / (n_mcmc_sample - 1.0), vp, nmsci);
   abyx(1.0 / n_mcmc_sample, phi_sln, 4);  abyx(1.0 / n_mcmc_sample, phi_sln1, 4);
   abyx(1.0 / n_mcmc_sample, theta_s, 4);  abyx(1.0 / n_mcmc_sample, theta_sln, 4);
   return(0);
}

double lnlike_msci(double x[], int np)
{
   /* loglikelihood for fitting beta(p,q) & gamma(a, b) to phi_x, phi_y, theta_x, theta_y, using
      phi_sln[4] = \sum log phi, for phi_x, phi_y, phi_z, phi_w;
      phi_sln1[4] = sum log (1 - phi);
      theta_s[4] = \sum theta, for theta_x, theta_y, theta_z, theta_w;
      theta_sln[4] = \sum log theta.
      x[] are the 8 or 16 hyperparameters.
   */
   int n = n_mcmc_sample;
   double a, b, p, q, lnp = 0, sumx, sumlnx, sumln1x;

   p = x[0]; q = x[1];  sumlnx = phi_sln[0]; sumln1x = phi_sln1[0];  /* phi_x */
   lnp += lgamma(p + q) - lgamma(p) - lgamma(q) + (p - 1) * sumlnx + (q - 1) * sumln1x;
   p = x[2]; q = x[3];  sumlnx = phi_sln[1]; sumln1x = phi_sln1[1];  /* phi_y */
   lnp += lgamma(p + q) - lgamma(p) - lgamma(q) + (p - 1) * sumlnx + (q - 1) * sumln1x;
   a = x[4]; b = x[5];  sumx = theta_s[0]; sumlnx = theta_sln[0];    /* theta_x */
   lnp += a * log(b) - lgamma(a) + (a - 1) * sumlnx - b * sumx;
   a = x[6]; b = x[7];  sumx = theta_s[1]; sumlnx = theta_sln[1];    /* theta_y */
   lnp += a * log(b) - lgamma(a) + (a - 1) * sumlnx - b * sumx;
   if (_D0DD1) {
      p = x[8]; q = x[9];  sumlnx = phi_sln[2]; sumln1x = phi_sln1[2];  /* phi_z */
      lnp += lgamma(p + q) - lgamma(p) - lgamma(q) + (p - 1) * sumlnx + (q - 1) * sumln1x;
      p = x[10]; q = x[11];  sumlnx = phi_sln[3]; sumln1x = phi_sln1[3];  /* phi_w */
      lnp += lgamma(p + q) - lgamma(p) - lgamma(q) + (p - 1) * sumlnx + (q - 1) * sumln1x;
      a = x[12]; b = x[13];  sumx = theta_s[2]; sumlnx = theta_sln[2];    /* theta_z */
      lnp += a * log(b) - lgamma(a) + (a - 1) * sumlnx - b * sumx;
      a = x[14]; b = x[15];  sumx = theta_s[3]; sumlnx = theta_sln[3];    /* theta_w */
      lnp += a * log(b) - lgamma(a) + (a - 1) * sumlnx - b * sumx;
   }
   return (-lnp);
}

void move_to_tower(int tower, double p[])
{
   /* p[] has the 4 or 8 msci parameters (phi, theta) from the mcmc sample file, which
      are reflected if tower is 1 (for mode D) or 1 2 3 (for model DD).
   */
   double t = -1;
   if (tower == 0) error2("should not be here");
   if (_D0DD1 == 0) {  /* phi_x, phi_y, theta_x, theta_y */
      p[0] = 1 - p[0];  p[1] = 1 - p[1];  swap2(p[2], p[3], t);
   }
   else {    /* phi_x, phi_y, theta_x, theta_y, phi_z, phi_w, theta_z, theta_w */
      if (tower == 1) {
         p[4] = 1 - p[4];  p[5] = 1 - p[5];  swap2(p[6], p[7], t);
      }
      else if (tower == 2) {
         p[0] = 1 - p[0];  p[1] = 1 - p[1];  swap2(p[2], p[3], t);
         swap2(p[4], p[5], t);  swap2(p[6], p[7], t);
      }
      else if (tower == 3) {
         p[0] = 1 - p[0];  p[1] = 1 - p[1];  swap2(p[2], p[3], t);
         t = p[4];  p[4] = 1 - p[5];  p[5] = 1 - t;
      }
   }
}

int compare_towers(int towers[], double sample[], double hyper[], int nmsci, int algorithm)
{
   int i, j, tower, tower_best, nchange = 0;
   double lnp, lnp_best, p[8], pnew[8];  /* p is current point, pnew is mirror point */

   for (i = 0; i < n_mcmc_sample; i++) {
      memmove(p, sample + i * nmsci, nmsci * sizeof(double));
      if (towers[i]) move_to_tower(towers[i], p);
      lnp_best = 0; tower_best = towers[i];
      for (tower = 0; tower < (_D0DD1 == 0 ? 2 : 4); tower++) {
         if (tower == towers[i]) continue;
         memmove(pnew, sample + i * nmsci, nmsci * sizeof(double));
         if (tower) move_to_tower(tower, pnew);
         if (algorithm > 0) {
            for (j = 0, lnp = 0; j < nmsci; j++) {
               if (algorithm == 1)   /* 1: CoG-N */
                  lnp += (p[j] - pnew[j]) * (p[j] + pnew[j] - 2 * mp[j]) / (2 * vp[j]);
               else                  /* 2: CoG-0 */
                  lnp += (p[j] - pnew[j]) * (p[j] + pnew[j] - 2 * mp[j]);
            }
         }
         else {                      /* 0: beta-gamma */
            lnp = ln_beta_ratio(pnew[0], p[0], hyper[0], hyper[1]);  /* phi_x */
            lnp += ln_beta_ratio(pnew[1], p[1], hyper[2], hyper[3]);  /* phi_x */
            lnp += ln_gamma_ratio(pnew[2], p[2], hyper[4], hyper[5]);  /* theta_x */
            lnp += ln_gamma_ratio(pnew[3], p[3], hyper[6], hyper[7]);  /* theta_y */
            if (_D0DD1) {
               lnp += ln_beta_ratio(pnew[4], p[4], hyper[8], hyper[9]);  /* phi_z */
               lnp += ln_beta_ratio(pnew[5], p[5], hyper[10], hyper[11]);  /* phi_w */
               lnp += ln_gamma_ratio(pnew[6], p[6], hyper[12], hyper[13]);  /* theta_z */
               lnp += ln_gamma_ratio(pnew[7], p[7], hyper[14], hyper[15]);  /* theta_w */
            }
         }
         if (lnp > lnp_best) {
            if (noisy >= 3) {
               printf("sample %7d: [%d]%9.6f%9.6f%9.6f%9.6f -> [%d]%9.6f%9.6f%9.6f%9.6f %8.3g\n", i + 1,
                  towers[i], p[0], p[1], p[2], p[3], tower, pnew[0], pnew[1], pnew[2], pnew[3], lnp);
               if (_D0DD1) printf(" %18s%9.6f%9.6f%9.6f%9.6f %6s%9.6f%9.6f%9.6f%9.6f\n",
                  "", p[4], p[5], p[6], p[7], "", pnew[4], pnew[5], pnew[6], pnew[7]);
            }
            tower_best = tower;  lnp_best = lnp;
         }
      }
      if (tower_best != towers[i]) {
         towers[i] = tower_best;
         nchange++;
      }
   }
   return nchange;
}

int initialize(double hyper[], int towers[], double sample[], int nmsci, int nhyper)
{
   /* sample[n*nmsci] has the mcmc sample of the msci paramters
   *  phi_x phi_y theta_x theta_y | phi_z phi_w theta_z theta_w
   *  The CoG algorithm needs mp[] vp[] only, and the rest is waste of time.
   */
   double p[8], epsilon = 1e-9;
   int   i, j;

   for (i = 0; i < n_mcmc_sample; i++) {
      for (j = 0; j < nmsci; j++)
         p[j] = sample[i * nmsci + j];

      for (j = 0; j < 2; ++j) {
         if (p[j] < 0 || p[j] > 1) 
            printf("phi_xy in sample %5d = %.6g out of range\n", i + 1, p[j]);
         p[j] = max2(p[j], epsilon);
         p[j] = min2(p[j], 1 - epsilon);         
      }
      if (_D0DD1)
         for (j = 4; j < 4+2; ++j) {
            if (p[j] < 0 || p[j] > 1)
               printf("phi_zw in sample %5d = %.6g out of range\n", i + 1, p[j]);
            p[j] = max2(p[j], epsilon);
            p[j] = min2(p[j], 1 - epsilon);
         }

      if (_D0DD1 == 0)
         towers[i] = (p[0] + p[1] > 1);
      else {
#if(1)
         if (p[0] <= 0.5 && p[4] <= 0.5) towers[i] = 0;
         else if (p[0] <= 0.5 && p[4] > 0.5) towers[i] = 1;
         else if (p[0] > 0.5 && p[5] <= 0.5) towers[i] = 2;
         else if (p[0] > 0.5 && p[5] > 0.5) towers[i] = 3;
#else  /* bad starting points for testing. */
         if (p[0] <= 0.5 && p[4] <= 0.5)
            towers[i] = 0;
         else
            towers[i] = 1 + i % 3;
#endif
      }
   }
   data_summaries(sample, towers, nmsci);
   /* initial values from methods of moments for beta and gamma */
   fit_beta_moments(hyper + 0, mp[0], vp[0]);  /* p q for phi_x */
   fit_beta_moments(hyper + 2, mp[1], vp[1]);  /* p q for phi_y */
   fit_gamma_moments(hyper + 4, mp[2], vp[2]);  /* a b for theta_x */
   fit_gamma_moments(hyper + 6, mp[3], vp[3]);  /* a b for theta_y */
   if (_D0DD1) {
      fit_beta_moments(hyper + 8, mp[4], vp[4]);  /* p q for phi_z */
      fit_beta_moments(hyper + 10, mp[5], vp[5]);  /* p q for phi_w */
      fit_gamma_moments(hyper + 12, mp[6], vp[6]);  /* a b for theta_z */
      fit_gamma_moments(hyper + 14, mp[7], vp[7]);  /* a b for theta_w */
   }
   return(0);
}

void print_towers(int towers[]) {
   for (int i = 0; i < n_mcmc_sample; i++) {
      printf("%d", towers[i]);
      if ((i + 1) % 100 == 0) printf(" [%5d]\n", i + 1);
   }
}

int process_msci_D_mcmc(int D0DD1, int algorithm, char mcmcf[], int columns[], char outf[])
{
   /* This processes mcmc samples from mcmcf to deal with the label-switching problem.
   */
   FILE* fmcmc = gfopen(mcmcf, "r"), * fout;
   FILE* frub = gfopen("rub", "w");
   double* record, * sample, p[8], near0 = 1e-9, near1 = 1 - near0;
   double hyper[16], bounds[16][2], lnL = 0, e = 1e-7, space[10000];
   int   nmsci, nhyper, nvar = -1, i, j, iround, * towers, points_moved;
   static int  lline = 1000000, ifields[MAXNFIELDS], HasHeader = 1;
   char* line;
   static char varstr[MAXNFIELDS][32] = { "" };

   _D0DD1 = D0DD1;
   nmsci = (_D0DD1 == 0 ? 4 : 8);  nhyper = (_D0DD1 == 0 ? 8 : 16);
   if ((line = (char*)malloc(lline * sizeof(char))) == NULL)
      error2("oom ...");
   scanfile(fmcmc, &n_mcmc_sample, &nvar, &HasHeader, line, ifields);
   printf("\n%d records, %d variables\n", n_mcmc_sample, nvar);
   record = (double*)malloc(nvar * sizeof(double));
   sample = (double*)malloc(n_mcmc_sample * nmsci * sizeof(double));
   towers = (int*)malloc(n_mcmc_sample * sizeof(int));
   if (record == NULL || sample == NULL || towers == NULL)  error2("oom");

   if (HasHeader)
      for (i = 0; i < nvar; i++) sscanf(line + ifields[i], "%s", varstr[i]);
   for (i = 0; i < n_mcmc_sample; i++) {
      for (j = 0; j < nvar; j++)   fscanf(fmcmc, "%lf", &record[j]);
      for (j = 0; j < nmsci; j++)  sample[i * nmsci + j] = record[columns[j]];
      if (sample[i * nmsci + 0] == 0) sample[i * nmsci + 0] = near0;      /* phi_X */
      if (sample[i * nmsci + 0] == 1) sample[i * nmsci + 0] = near1;
      if (sample[i * nmsci + 1] == 0) sample[i * nmsci + 1] = near0;      /* phi_Y */
      if (sample[i * nmsci + 1] == 1) sample[i * nmsci + 1] = near1;
      if (nmsci == 8) {
         if (sample[i * nmsci + 4] == 0) sample[i * nmsci + 4] = near0;     /* phi_Z */
         if (sample[i * nmsci + 4] == 1) sample[i * nmsci + 4] = near1;
         if (sample[i * nmsci + 5] == 0) sample[i * nmsci + 5] = near0;     /* phi_W */
         if (sample[i * nmsci + 5] == 1) sample[i * nmsci + 5] = near1;
      }
   }
   if (_D0DD1 == 0)      printf("phi_x phi_y theta_x theta_y:\n");
   else if (_D0DD1 == 1) printf("phi_x phi_y theta_x theta_y phi_z phi_w theta_z theta_w:\n");
   else                  printf("error");
   for (j = 0; j < nmsci; j++) printf("%s ", varstr[columns[j]]);
   printf("\n");

   initialize(hyper, towers, sample, nmsci, nhyper);
   if (n_mcmc_sample <= 10000) print_towers(towers);

   if (algorithm == 0) {
      lnL = lnlike_msci(hyper, nhyper);
      printf("\nsummaries (phi_sln, phi_sln1, theta_s, theta_sln): ");
      matout2(stdout, phi_sln, 1, (_D0DD1 == 0 ? 2 : 4), 9, 5);
      matout2(stdout, phi_sln1, 1, (_D0DD1 == 0 ? 2 : 4), 9, 5);
      matout2(stdout, theta_s, 1, (_D0DD1 == 0 ? 2 : 4), 9, 5);
      matout2(stdout, theta_sln, 1, (_D0DD1 == 0 ? 2 : 4), 9, 5);
      printf("\ninitials: ");
      matout2(stdout, hyper, nhyper / 8, 8, 9, 5);
      printf("\nlnL0 = %12.6f\n", lnL);
   }
   for (i = 0; i < nhyper; i++) { bounds[i][0] = 0.5;  bounds[i][1] = 99999; }
   for (iround = 0; iround < 100; iround++) {
      points_moved = compare_towers(towers, sample, hyper, nmsci, algorithm);
      printf("Round %2d, %2d points moved...\n", iround, points_moved);
      data_summaries(sample, towers, nmsci);
      if (algorithm == 0) {
         lnL = lnlike_msci(hyper, (_D0DD1 == 0 ? 8 : 16));
         printf("lnL = % 12.6f\n", lnL);
      }
      if (points_moved == 0) break;
      if (algorithm == 0) {
         /* get MLEs for phi_x, phi_y, theta_x, theta_y */
         j = ming2(frub, &lnL, lnlike_msci, NULL, hyper, bounds, space, e, nhyper);
         matout(stdout, hyper, 1, nhyper);
      }
      else {
         printf("means: ");
         matout(stdout, mp, 1, nmsci);
      }
   }
   /* print out processed mcmc samples */
   if (n_mcmc_sample <= 10000) print_towers(towers);
   for (i = 0, points_moved = 0; i < n_mcmc_sample; i++)
      points_moved += (towers[i] > 0);
   printf("\nEnd of algorithm: %4d points reflected\n", points_moved);

   if (j) {
      printf("print processed samples into %s\n", outf);
      fout = gfopen(outf, "w");
      rewind(fmcmc);
      fgets(line, lline, fmcmc);
      fprintf(fout, line);
      for (i = 0; i < n_mcmc_sample; i++) {
         for (j = 0; j < nvar; j++) fscanf(fmcmc, "%lf", &record[j]);
         if (towers[i]) {
            for (j = 0; j < nmsci; j++)  p[j] = record[columns[j]];
            move_to_tower(towers[i], p);
            for (j = 0; j < nmsci; j++)  record[columns[j]] = p[j];
         }
         for (j = 0; j < nvar; j++)
            fprintf(fout, (j == 0 ? "%.0f" : (j == nvar - 1 ? "\t%.3f" : "\t%.6f")), record[j]);
         fprintf(fout, "\n");
      }
      fclose(fout);
   }
   fclose(fmcmc);
   free(record);  free(sample);  free(towers);
   return(0);
}

int main(int argc, char* argv[])
{
#if msciD_data == 0    /* heliconius MSci-D */
   int Dcolumns[4] = { 12, 13, 6, 7 };   /* phi_x, phi_y, theta_x, theta_y */
   char mcmcf[4096] = "../data/h_mel_tim_num_noncod.txt";
   //char mcmcf[4096] = "../data/h_mel_tim_num_exonic.txt";
#elif msciD_data == 1  /* simulation-D-L500 & D-L2000 & D-L8000 */
   int Dcolumns[4] = { 9, 10, 4, 5 };   /* phi_x, phi_y, theta_x, theta_y */
   char mcmcf[4096] = "../../A/0ngoing/bppMSci-twintowers/simulation-D/mcmc_L500.txt";
   // char mcmcf[4096] = "../../A/0ngoing/bppMSci-twintowers/simulation-D/mcmc_L2000.txt";
   // char mcmcf[4096] = "../../A/0ngoing/bppMSci-twintowers/simulation-D/mcmc_L8000-r1.txt";
#elif msciD_data == 2  /* simulation-DD-L500 */
   int Dcolumns[8] = { 14, 16, 5, 7, 13, 15, 4, 6 };  /* phi_x, phi_y, theta_x, theta_y, phi_z, phi_w, theta_z, theta_w  */
   char mcmcf[4096] = "../../../A/0ngoing/bppMSci-twintowers/simulation-DD/L500/mcmc.rep2.txt";
#endif
   char outf[4096], * p, * processed[3] = { "_processed", "_processed_CoGN" , "_processed_CoG0" };
   int algorithm;

   if (argc == 6) {
      strcpy(mcmcf, argv[1]);    /* original mcmc file */
      for (int i = 0; i < 4; i++)  sscanf(argv[2 + i], "%d", &Dcolumns[i]);
   }
   else {
      printf("Usage: \n\tbpp-msci-D-process-mcmc <mcmc-sample-filename> <phi_X-column> <phi_Y-column> <theta_X-column> <theta_Y-column>");
   }
   noisy = 1;
   Small_Diff = 1e-6;
   for (algorithm = 0; algorithm < 3; algorithm++) {
      strcpy(outf, mcmcf);
      p = strstr(outf, ".txt");
      if (!p) p = outf + strlen(outf);
      sprintf(p, "%s.txt", processed[algorithm]);
      printf("\n\n*** algorithm %s ***\n: %s --> %s\n", algorithms[algorithm], mcmcf, outf);
      process_msci_D_mcmc((msciD_data == 2), algorithm, mcmcf, Dcolumns, outf);
   }
}
