#include "./ldpc_parm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
//#include "./lib_mat/lib_mat.h"
#include "./lib_rand/lib_rand.h"
#include "./lib_math/fast_math.h"
#include "./bp_dec/bp_dec.h"
#include "./bp_dec/bp_llr.h"
#include "./OSD/OSD.h"
#include <stdbool.h>

#include <math.h>
#include <time.h>
#include <inttypes.h>



/*
        MAIN
                     */
int main(void)
{
  setvbuf(stdout, NULL, _IONBF, 0);  
  puts(">>> BPOSD starting <<<");    

  /* --------------------------------------------------------------------------
                              Core parameters
  * -------------------------------------------------------------------------- */
  int osdw     = OSDW;
  int max_iter = MAX_ITER;

  /* If AFP is your decoding parameter (e.g., attenuation / scaling / damping)*/
  double alpha_a = 0, alpha_b = AFP, alpha= AFP;

  uint8_t *finalResult  = NULL;
  uint8_t *finalResult2 = NULL;

  uint8_t *lastForGF2 = calloc((size_t)N, sizeof(*lastForGF2));
  double  (*probGF2)[2] = calloc((size_t)N, sizeof(*probGF2));

  int32_t iter = 0, b = 0;

  GFQ_t diff[N], zz[M], zz_G[M_G];
  GFQ_t dX[N], dZ[N], zz_Z[N], zz_X[N];

  /* RNG seed state (ensure it is initialized appropriately for your RNG design) */
  uint64_t tx_seed[4] = {0};

  /* Channel / noise parameters */
  double p_nos   = P_NOISE;
  double rnd_val = 0.0;
  double p_bias  = 0.0;
  double afp     = 0.0;   /* if afp is distinct from alpha, keep it; otherwise remove */

  rnd256_init();
  rnd256_init_priv(tx_seed);

  double bias2[2 * N][2];

  /* Load HX/HZ/A/GX/GZ */
  FILE *fpA_HX = fopen(PATH_A_HX, "r");

  a2_matrix A_HX_buf;
  a2_matrix *A_HX = &A_HX_buf;
  load_A2(fpA_HX, A_HX, "HX");
  fclose(fpA_HX); fpA_HX = NULL;

  printf("A_HX loaded, size = %d x %d\n", A_HX->MM, A_HX->NN);

  BP2_Ctl BP_HX_buf;
  BP2_Ctl *bp_HX = &BP_HX_buf;
  alloc_BPC2(A_HX, bp_HX);

  FILE *fpA_HZ = fopen(PATH_A_HZ, "r");
  a2_matrix A_HZ_buf;
  a2_matrix *A_HZ = &A_HZ_buf;
  load_A2(fpA_HZ, A_HZ, "HZ");
  fclose(fpA_HZ); fpA_HZ = NULL;
  printf("A_HZ loaded, size = %d x %d\n", A_HZ->MM, A_HZ->NN);

  BP2_Ctl BP_HZ_buf;
  BP2_Ctl *bp_HZ = &BP_HZ_buf;
  alloc_BPC2(A_HZ, bp_HZ);

  FILE *fpA = fopen(PATH_A, "r");
  a_matrix_GFQ Amtx;
  a_matrix_GFQ *A = &Amtx;
  load_A_GFQ(fpA, A);
  fclose(fpA); fpA = NULL;


  FILE *fpGZ = fopen(PATH_GZ, "r");
  FILE *fpGX = fopen(PATH_GX, "r");
  g_matrix_GFQ GZ_buf, GX_buf;
  g_matrix_GFQ *GZ = &GZ_buf;
  g_matrix_GFQ *GX = &GX_buf;

  load_G2(fpGZ, GZ, "GZ");
  fclose(fpGZ); fpGZ = NULL;
  printf("GZ loaded, size = %u x %u\n", GZ->MM, GZ->NN);

  load_G2(fpGX, GX, "GX");
  fclose(fpGX); fpGX = NULL;
  printf("GX loaded, size = %u x %u\n", GX->MM, GX->NN);

  char st_name[1000] = "";
  char path_BER[1000] ;  
  char path_Last[1000]; 
  sprintf(st_name, "%s%s_it%u_dec%u", st_name, PRE_A, max_iter, BY_DEC);  // care: st_name in both I/O , behavior not guaranteed
  if(USE_GF2_DEC) sprintf(st_name, "%s_%s", st_name, "GF2"); 

  if(osdw >= 0){sprintf(st_name, "%s_OSDW%d", st_name, osdw);}
  else if(osdw == -1){sprintf(st_name, "%s_NoOSD", st_name);}
  else if(osdw == -2){sprintf(st_name, "%s_OSDfull", st_name);}

  sprintf(path_BER,  "Results/BER_%s.txt", st_name);
  FILE *fpBER = fopen ( path_BER , "a+" );

  uint64_t nCW, nErr, nFa, nEb, nIterAcc, nStop=N_ERR_STOP,  nErrQ, nFaQ, nEq;
  uint32_t nErrBit,  nErrQbit;
  double   bler, ber, far,  qbler, qber, qfar,  chk_bler;

  time_t  tStr, tEnd, sec, sec_last=0;  // care that: gcc32 has year 2038 problem (gcc64 no this problem)
  nCW=0; nErr=0; nFa=0; nEb=0; nIterAcc=0; nErrBit=0;  nErrQ=0; nFaQ=0; nEq=0; nErrQbit=0;
  chk_bler=1.0;
  printf("\n\n//==== init p_nos = %g ====// N = %u K = %u   It = %u  BY_DEC %u \n\n", p_nos,N,K,MAX_ITER,BY_DEC);
  tStr = time(NULL);

  OSD osdDecoder, *osdDec = &osdDecoder;
  initOSD(osdDec);
  FILE *fp = fopen( PATH_A , "r" );
  load_A_OSD(osdDec, fp);
  fclose(fp);

  osdDec -> RankH  = N- K;
  osdDec -> osdw = osdw;

  alpha = alpha_a * log10(p_nos) * 10 + alpha_b;
  if(alpha <= 50)    alpha = 0;
  printf("\n\n//=================  Decode ================= //\n\n");
  printf("\nNew alpha: %lf\n\n", alpha);



  /* ================= Decoding ==================== */
  do {
    //- add noise nn_(Nx1), where (Nx1) means a column vector
    GFQ_t nn[N] = {0};
    GFQ_t nn_Xerr[N] = {0};
    GFQ_t nn_Zerr[N] = {0};
    GFQ_t nn_hat[N] = {0};

    for (b = 0; b < N; b++){
      rnd_val = rv_UnifOne();
      if (rnd_val < p_nos / 3){          
        nn[b] = 3; 
        nn_Xerr[b] = 1;
        nn_Zerr[b] = 1;
      }
      else if (rnd_val < 2 * p_nos / 3){ 
        nn[b] = 2; 
        nn_Zerr[b] = 1;
      }
      else if (rnd_val < p_nos){        
        nn[b] = 1; 
        nn_Xerr[b] = 1;
      }
    }  

    const int MAX_OUTER = MAX_CYCLE; 
    const int MAX_INNER = MAX_ITER; 

    double biasZ[N][2];
    double biasX[N][2]; 

    uint8_t target_syndrome_Z[M];
    uint8_t target_syndrome_X[M]; 
    uint8_t target_syndrome_GF4[M*2];

    uint8_t ttZ_prev[N];
    uint8_t ttX_prev[N];
    memset(ttZ_prev, 0xFF, sizeof(ttZ_prev));
    memset(ttX_prev, 0xFF, sizeof(ttX_prev));   

    uint8_t syndrome_ok_Zerror, syndrome_ok_Xerror;
    uint8_t init_syndro_ok_Zerror, init_syndro_ok_Xerror;

    int actual_run_count = 0;

    double sum_qn_Z[N][2]; 
    double sum_qn_X[N][2];
    memset(sum_qn_Z, 0, sizeof(sum_qn_Z));
    memset(sum_qn_X, 0, sizeof(sum_qn_X));

    int iter;

    GenSyndrome_GF2(A_HX, nn_Zerr, target_syndrome_Z);
    GenSyndrome_GF2(A_HZ, nn_Xerr, target_syndrome_X);
    Quan_GenSyndrome(A, nn, target_syndrome_GF4);

    /* ====================initial biasZ =========================== */
    double p_bsc = p_nos * 2.0 / 3; 
    for (int b = 0; b < N; b++)
    {
      biasZ[b][0] = 1.0 - p_bsc;
      biasZ[b][1] = p_bsc;
    }
    #if !useCorrelation
      for (int b = 0; b < N; b++)
      {
        biasX[b][0] = 1.0 - p_bsc;
        biasX[b][1] = p_bsc;
      }
    #endif
    

    /* ====================== Iterative Binary MBP ================================= */
    for (int outer = 0; outer < MAX_OUTER; outer++){
      actual_run_count++; 

      /* ============================================================= */
      /*         Phase 1: Correlated BP Iteration (No OSD)             */
      /* ============================================================= */
      /**********  (I)  Using Hx to decode Z error  *************/
      init_syndro_ok_Zerror = bp2_init22(bp_HX, biasZ, A_HX, target_syndrome_Z);
      iter = 0;
      syndrome_ok_Zerror = init_syndro_ok_Zerror;
      while (!syndrome_ok_Zerror){
        iter++;
        #if BY_DEC == 20
          syndrome_ok_Zerror = bp2_dec20(bp_HX, A_HX);
        #elif BY_DEC == 24
          syndrome_ok_Zerror = bp2_dec24(bp_HX, A_HX);
        #endif

        if ( syndrome_ok_Zerror || iter==MAX_INNER )   break;
      }

        for(int b=0; b<N; b++){
            sum_qn_Z[b][0] += bp_HX->qn[b][0]; 
            sum_qn_Z[b][1] += bp_HX->qn[b][1]; 
        }

      /* -------------- 產生給 Hz 的新 biasX ---------------------- */
      #if useCorrelation
        for (int b = 0; b < N; b++)
        {
          if (bp_HX->tt[b])
          {                    // qubit b 判定有 Z-error 
            biasX[b][0] = 0.5; // erasure(1/2)          
            biasX[b][1] = 0.5;     
          }
          else
          { 
            biasX[b][0] = 1.0 - (p_nos/3.0) / (1.0 - 2.0*p_nos/3.0);
            biasX[b][1] = (p_nos/3.0) / (1.0 - 2.0*p_nos/3.0);
          }
        }
      #endif
      
      /**********  (II) Using HZ to decode X error *************/
      init_syndro_ok_Xerror = bp2_init22(bp_HZ, biasX, A_HZ, target_syndrome_X);
      iter = 0;
      syndrome_ok_Xerror = init_syndro_ok_Xerror;
      while (!syndrome_ok_Xerror){
        iter++;
        #if BY_DEC == 20
          syndrome_ok_Xerror = bp2_dec20(bp_HZ, A_HZ);
        #elif BY_DEC == 24  
          syndrome_ok_Xerror = bp2_dec24(bp_HZ, A_HZ);
        #endif



        if ( syndrome_ok_Xerror || iter==MAX_INNER )   break;
      }

        for(int b=0; b<N; b++){
            sum_qn_X[b][0] += bp_HZ->qn[b][0];
            sum_qn_X[b][1] += bp_HZ->qn[b][1]; 
        }

      /* ********** 兩個 syndrome 都 OK 就停止 ***********/
      if(syndrome_ok_Zerror && syndrome_ok_Xerror) {
          break;
      }
      /**********  與前一輪完全相同就停止 (prev 相等) **********/
      bool sameZ = true, sameX = true;

      /* Z part: 比對 HX->tt vs ttZ_prev */
      for (int b = 0; b < N; ++b) {
          if (bp_HX->tt[b] != ttZ_prev[b]) { sameZ = false; break; }
      }

      /* X part: 比對 HZ->tt vs ttX_prev */
      for (int b = 0; b < N; ++b) {
          if (bp_HZ->tt[b] != ttX_prev[b]) { sameX = false; break; }
      }
      if (sameZ && sameX) {
          break;
      }

      /* 把本輪結果存成 prev，供下一輪比較 */
      memcpy(ttZ_prev, bp_HX->tt, N * sizeof(bp_HX->tt[0]));
      memcpy(ttX_prev, bp_HZ->tt, N * sizeof(bp_HZ->tt[0]));

      /* -------------- 產生給下一輪 Hx 的新 biasZ ----------------*/
      for (int b = 0; b < N; b++)
      {
        if (bp_HZ->tt[b])
        { 
          biasZ[b][0] = 0.5;
          biasZ[b][1] = 0.5;
        }
        else
        {
          biasZ[b][0] = 1.0 - (p_nos/3.0) / (1.0 - 2.0*p_nos/3.0);
          biasZ[b][1] = (p_nos/3.0) / (1.0 - 2.0*p_nos/3.0);
        }
      }
    }
    /* ===============Outer Iteration Finish================= */

    /* ============================================================= */
    /*                        Phase 2: Quaternary OSD                           */
    /* ============================================================= */
    if(osdw != -1 && (!syndrome_ok_Zerror || !syndrome_ok_Xerror)){
        double **probGF4 = calloc(N, sizeof(double*));
        
        for(int b=0; b<N; b++){
          probGF4[b] = calloc(4, sizeof(double));

          double avg_z0 = sum_qn_Z[b][0] / (double)actual_run_count;
          double avg_z1 = sum_qn_Z[b][1] / (double)actual_run_count;   
          double avg_x0 = sum_qn_X[b][0] / (double)actual_run_count;
          double avg_x1 = sum_qn_X[b][1] / (double)actual_run_count;
          double pZ = avg_z1 / (avg_z0 + avg_z1 + 1e-300);
          double pX = avg_x1 / (avg_x0 + avg_x1 + 1e-300);

          pZ = clamp01(pZ);
          pX = clamp01(pX);

          double rho = 1; 
          double sigma_x = sqrt(pX * (1.0 - pX));
          double sigma_z = sqrt(pZ * (1.0 - pZ));
          double y = (pX * pZ) + rho * sigma_x * sigma_z;
          double y_min = (pX + pZ > 1.0) ? (pX + pZ - 1.0) : 0.0;
          double y_max = (pX < pZ) ? pX : pZ;
          if(y < y_min) y = y_min;   
          if(y > y_max) y = y_max;

          double px_only = pX - y;
          double pz_only = pZ - y;
          double pI = 1.0 - px_only - pz_only - y;


          if(pI < 0) pI = 0;
          if(px_only < 0) px_only = 0;
          if(pz_only < 0) pz_only = 0;
          double s = pI + px_only + pz_only + y;
          probGF4[b][0] = pI / s;       // I
          probGF4[b][1] = px_only / s;  // X
          probGF4[b][2] = pz_only / s;  // Z
          probGF4[b][3] = y / s;        // Y
        }
        
        
        finalResult = post_decodeOSDfull(osdDec, target_syndrome_GF4, (const double **)probGF4, lastForGF2, max_iter);
        for(int b = 0; b < N; b++){
          uint8_t val = finalResult[b];
          bp_HX->tt[b] = (val == 2 || val == 3) ? 1 : 0; 
          bp_HZ->tt[b] = (val == 1 || val == 3) ? 1 : 0; 
        }
        for(int i=0; i<N; i++) free(probGF4[i]);
        free(probGF4);
    }
   
    /* ===================== nn_hat ======================== */
    for (int b = 0; b < N; b++)
      nn_hat[b] = bp_HZ->tt[b] | (bp_HX->tt[b] << 1);
   
    finalResult2 = nn_hat;
    nCW++;

    nErrBit = HamDist(nn, finalResult2, N);

    if(nErrBit!=0) {}  // fine only when both syndrome and HamDist pass
    else { nErr++; nEb+=nErrBit;}

    //-- "Quantum" degeneracy check
    nErrQbit = nErrBit;   // start with classical error
    if(nErrQbit!=0) {
      VecDiff(diff, finalResult2, nn, N);    // 或 VecDiff(diff, nn_hat, nn, N) 同義
      for (int b = 0; b < N; ++b) {
          dX[b] = (uint8_t)(diff[b] & 1u);          // X bit
          dZ[b] = (uint8_t)((diff[b] >> 1) & 1u);   // Z bit
      }

      DegSyndrome_GF2(GZ, dX, zz_Z);   
      DegSyndrome_GF2(GX, dZ, zz_X);

      int okZ_degen = 1, okX_degen = 1;
      for (uint32_t m = 0; m < GZ->MM; ++m) if (zz_Z[m]) { okZ_degen = 0; break; }
      for (uint32_t m = 0; m < GX->MM; ++m) if (zz_X[m]) { okX_degen = 0; break; }

      if (okZ_degen && okX_degen) {
          nErrQbit = 0;
      }      
    }
    if(nErrQbit!=0) nErrQ++;
  
    if(nCW % (1000*1) == 0) { 
      tEnd = time(NULL);  sec = tEnd-tStr;
      # if SORTCOMPARE == LASTFORSORT && STORE_LF && OSDW !=-1
      if(nCW % 1000000 == 0){
          fpLF = fopen(path_saveLF, "a");
          for(int i = 1; i<=max_iter+1; i++)    fprintf(fpLF, "%llu ", storeLF[i]);
          fprintf(fpLF, "\n");
          fclose(fpLF);
      }
      #endif
      if(sec >= sec_last+10) {
        tEnd = time(NULL);  sec = tEnd-tStr;  sec_last=sec;
        qbler=(double)nErrQ/nCW;

        printf("    now p %-10g  nCW %-10"PRIu64"  nErrQ %-6"PRIu64"  => qBLER %-10g  @ %10ld sec\n",
              p_nos, nCW, nErrQ, qbler, sec);
        
        if(nErrQ >= N_ERR_STOP_1em7){
          if(1 > qbler && qbler > 1e-1) nStop = N_ERR_ACCURATE;
          if(qbler <= 1e-1)  nStop = N_ERR_STOP_1em1;  // start to accelerate
          if(qbler <= 1e-2)  nStop = N_ERR_STOP_1em2;  // start to accelerate more ...
          if(qbler <= 1e-3)  nStop = N_ERR_STOP_1em3;
          if(qbler <= 1e-4)  nStop = N_ERR_STOP_1em4;
          if(qbler <= 1e-5)  nStop = N_ERR_STOP_1em5;
          if(qbler <= 1e-6)  nStop = N_ERR_STOP_1em6;
          if(qbler <= 1e-7)  nStop = N_ERR_STOP_1em7;
        }
      }
    }

    if (nErrQ >= nStop) {
      tEnd = time(NULL);  sec = tEnd-tStr;
      qbler=(double)nErrQ/nCW;
      chk_bler = qbler;

      printf("p_nos %-10g  qBLER %-10e  for nCW %-10"PRIu64"  nErrQ %-6"PRIu64"    @ %10ld sec \n\n",
              p_nos,qbler, nCW, nErrQ, sec);
      fprintf(fpBER,  "%g\t %10e\n", p_nos,qbler); 
      
      nCW=0; nErr=0; nFa=0; nEb=0; nIterAcc=0; nErrBit=0;  nErrQ=0; nFaQ=0; nEq=0; nErrQbit=0;

      double step;
      (qbler>1e-5)? (step = P_STEP):(step = P_STEP/2);
      p_nos *= pow(10, step/10);
      alpha = alpha_a * log10(p_nos) * 10 + alpha_b;
      if(alpha <= 50)    alpha = 0;
      printf("\nNew alpha: %lf\n\n", alpha);

      if(1 > qbler && qbler > 1e-1) nStop = N_ERR_ACCURATE;
      if(qbler <= 1e-1)  nStop = N_ERR_STOP_1em1;  // start to accelerate
      if(qbler <= 1e-2)  nStop = N_ERR_STOP_1em2;  // start to accelerate more ...
      if(qbler <= 1e-3)  nStop = N_ERR_STOP_1em3;
      if(qbler <= 1e-4)  nStop = N_ERR_STOP_1em4;
      if(qbler <= 1e-5)  nStop = N_ERR_STOP_1em5;
      if(qbler <= 1e-6)  nStop = N_ERR_STOP_1em6;
      if(qbler <= 1e-7)  nStop = N_ERR_STOP_1em7;
    }
  } while (chk_bler>BLER_STOP);

  fprintf(fpBER,  "%s\n", path_BER);
  fclose(fpBER);
	printf("main func done!\n");
  free_OSD(osdDec);
	return 0;
}

//gcc BPOSD.c bp_dec/bp_dec.c bp_dec/bp_llr.c lib_rand/splitmix64.c lib_rand/lib_rand.c lib_rand/xoshiro256starstar.c lib_math/fast_math.c OSD/OSD.c -o BPOSD
//.\BPOSD.exe