#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>

#include "mem.h"
#include "defs.h"
#include "msg.h"
#include "buffer.h"
#include "levels.h"
#include "cache.h"
#include "common.h"
#include "pmodels.h"
#include "context.h"

Parameters *P; // FOR THREAD SHARING

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - R E A D   B A C K W A R D S - - - - - - - - - - - -

char *fgets_backwards(char *str, int size, FILE *F){
  int i;
  if(ftell(F) == 0) return NULL;
  for(i = 0; ftell(F) != 0 && i < size; i++){
    fseek(F, -1, SEEK_CUR);
    str[i] = (char)fgetc(F);
    if(str[i] == '\n' && i != 0) break;
    fseek(F, -1, SEEK_CUR);
    }
  str[i] = '\0';
  return str;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - S I M P L E   R E V E R S E R - - - - - - - - - - - -

void reverse(char *start){
  size_t len = strlen(start);
  for(char *end = &start[len-1] ; start < end ; start++, end--){
    char tmp = start[0];
    start[0] = end[0];
    end[0]   = tmp;
    }
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S O R - - - - - - - - - - - - - -

void CompressTarget(Threads T){

  char        in_name [4096];
  char        out_name[4096];
  sprintf(in_name,  ".dlcr_%u.dna", T.id + 1);
  sprintf(out_name, ".dlcr_%u.inf", T.id + 1);
  FILE        *Reader  = Fopen(in_name,  "r");
  FILE        *Writter = Fopen(out_name, "w");

  uint32_t    n, k, cModel, totModels, idxPos;
  uint64_t    compressed = 0, nSymbols = 0, nBases = 0, i = 0;
  uint8_t     *readBUF, sym, irSym, *pos;
  PModel      **pModel, *MX;
  FloatPModel *PT;
  CMWeight    *WM;
  CBUF        *symbBUF;
  CModel      **cModels;
  CACHE       *C;
  
  // EXTRA MODELS DERIVED FROM EDITS
  totModels = P->nModels;
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].edits != 0)
      totModels += 1;

  C             = CreateCache(P->region, P->nModels); 
  pModel        = (PModel  **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n]   = CreatePModel      (ALPHABET_SIZE);
  MX            = CreatePModel      (ALPHABET_SIZE);
  PT            = CreateFloatPModel (ALPHABET_SIZE);
  WM            = CreateWeightModel (totModels);

  readBUF  = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  symbBUF  = CreateCBuffer(BUFFER_SIZE, BGUARD);

  cModels = (CModel **) Malloc(P->nModels * sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    cModels[n] = CreateCModel(TARGET, P->model[n].ctx, P->model[n].den,
    P->model[n].ir, P->model[n].hashSize, P->model[n].gamma,
    P->model[n].edits, P->model[n].eDen, P->model[n].eGamma);

  // GIVE SPECIFIC GAMMA:
  int pIdx = 0;
  for(n = 0 ; n < P->nModels ; ++n){
    WM->gamma[pIdx++] = cModels[n]->gamma;
    if(P->model[n].edits != 0)
      WM->gamma[pIdx++] = cModels[n]->SUBS.eGamma;
    }

  while((k = fread(readBUF, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){

      sym = readBUF[idxPos];

      // FINAL FILTERING DNA CONTENT
      if(sym != 'A' && sym != 'C' && sym != 'G' && sym != 'T'){
        fprintf(Writter, "2\n"); // FORCE HIGH COMPLEXITY <- UNKNOWN SYMBOL
	++compressed;
        continue;
        }

      symbBUF->buf[symbBUF->idx] = sym = DNASymToNum(sym);
      memset((void *)PT->freqs, 0, ALPHABET_SIZE * sizeof(double));

      n = 0;
      pos = &symbBUF->buf[symbBUF->idx-1];
      for(cModel = 0 ; cModel < P->nModels ; ++cModel){
        CModel *CM = cModels[cModel];
        GetPModelIdx(pos, CM);
        ComputePModel(CM, pModel[n], CM->pModelIdx, CM->alphaDen);
        ComputeWeightedFreqs(WM->weight[n], pModel[n], PT, 4);
        if(CM->edits != 0){
          ++n;
          CM->SUBS.seq->buf[CM->SUBS.seq->idx] = sym;
          CM->SUBS.idx = GetPModelIdxCorr(CM->SUBS.seq->buf+
          CM->SUBS.seq->idx-1, CM, CM->SUBS.idx);
          ComputePModel(CM, pModel[n], CM->SUBS.idx, CM->SUBS.eDen);
          ComputeWeightedFreqs(WM->weight[n], pModel[n], PT, 4);
          }
        ++n;
        }

      ComputeMXProbs(PT, MX, 4);
      fprintf(Writter, "%.3g\n", PModelSymbolNats(MX, sym) / M_LN2);
      CalcDecayment(WM, pModel, sym);

      for(n = 0 ; n < P->nModels ; ++n){
	switch(cModels[n]->ir){
          case 0:
          C->E[C->pos].M[n].idx    = cModels[n]->pModelIdx;
          C->E[C->pos].M[n].s      = sym;
	  break;
	  case 1:
          C->E[C->pos].M[n].idx    = cModels[n]->pModelIdx;
          C->E[C->pos].M[n].s      = sym;
	  C->E[C->pos].M[n].idx_ir = cModels[n]->pModelIdxIR;
          C->E[C->pos].M[n].s_ir   = irSym;
	  break;
	  case 2:
          C->E[C->pos].M[n].idx_ir = cModels[n]->pModelIdxIR;
          C->E[C->pos].M[n].s_ir   = irSym;
	  break;
	  default: 
	  fprintf(stderr, "ERROR: no store action!\n");
	  exit(1);
	  }
        }

      if(compressed >= C->size - 2){
        for(n = 0 ; n < P->nModels ; ++n){
          uint32_t pos = (C->pos == 0) ? C->size - 1 : C->pos - 1;
          switch(cModels[n]->ir){
            case 0:
            UpdateCModelCounter(cModels[n], C->E[pos].M[n].s, C->E[pos].M[n].idx);
            break;
            case 1:
            UpdateCModelCounter(cModels[n], C->E[pos].M[n].s, C->E[pos].M[n].idx);
            UpdateCModelCounter(cModels[n], C->E[pos].M[n].s_ir, C->E[pos].M[n].idx_ir);
            break;
            case 2:
            UpdateCModelCounter(cModels[n], C->E[pos].M[n].s_ir, C->E[pos].M[n].idx_ir);
            break;
            default:
            fprintf(stderr, "ERROR: no update action!\n");
	    exit(1);
            }
          }
        }

      UpdateCache(C);

/*      
      for(n = 0 ; n < P->nModels ; ++n){
        switch(cModels[n]->ir){
          case 0:
          UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);
          break;
          case 1:
          UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);
          irSym = GetPModelIdxIR(symbBUF->buf+symbBUF->idx, cModels[n]);
          UpdateCModelCounter(cModels[n], irSym, cModels[n]->pModelIdxIR);
          break;
          case 2:
          irSym = GetPModelIdxIR(symbBUF->buf+symbBUF->idx, cModels[n]);
          UpdateCModelCounter(cModels[n], irSym, cModels[n]->pModelIdxIR);
          break;
          default:
          UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);
          break;
          }
        }
*/

      RenormalizeWeights(WM);

      n = 0;
      for(cModel = 0 ; cModel < P->nModels ; ++cModel){
        if(cModels[cModel]->edits != 0)
          CorrectCModelSUBS(cModels[cModel], pModel[++n], sym);
        ++n;
        }

      UpdateCBuffer(symbBUF);
      ++compressed;
      }

  fclose(Writter);
  Free(MX);
  
  // for(n = 0 ; n < P->nModels ; ++n)
  //   FreeCModel(cModels[n]);

  RemoveCache(C);
  for(n = 0 ; n < totModels ; ++n){
    Free(pModel[n]->freqs);
    Free(pModel[n]);
    }
  Free(pModel);
  Free(PT);
  Free(readBUF);
  RemoveCBuffer(symbBUF);
  fclose(Reader);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - F   T H R E A D I N G - - - - - - - - - - - - - - -

void *CompressThread(void *Thr){
  Threads *T = (Threads *) Thr;
  CompressTarget(T[0]);
  pthread_exit(NULL);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - C O M P R E S S O R   M A I N - - - - - - - - - - - -

void CompressAction(Threads *T){

  pthread_t t[2];
  uint32_t n;

  pthread_create(&(t[1]), NULL, CompressThread, (void *) &(T[0]));
  pthread_create(&(t[2]), NULL, CompressThread, (void *) &(T[1]));
    
  pthread_join(t[1], NULL);
  pthread_join(t[2], NULL);

  return;
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t main(int argc, char *argv[]){

  char      **p = *&argv, **xargv, *xpl = NULL;
  int32_t   n, xargc = 0;
  uint32_t  k;
  uint64_t  totalBytes, totalSize;
  clock_t   stop = 0, start = clock();
  Threads   *T;

  P = (Parameters *) Malloc(1 * sizeof(Parameters));

  if((P->help = ArgsState(DEFAULT_HELP, p, argc, "-h", "--help")) == 1
  || argc < 2){
    PrintMenuCompression();
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_VERSION, p, argc, "-V", "--version")){
    PrintVersion();
    return EXIT_SUCCESS;
    }

  if(ArgsState(0, p, argc, "-s", "--show-levels")){
    PrintLevels();
    return EXIT_SUCCESS;
    }

  if(ArgsState(0, p, argc, "-p", "--show-parameters")){
    ModelsExplanation();
    return EXIT_SUCCESS;
    }

  P->force     = ArgsState  (DEFAULT_FORCE,   p, argc, "-F", "--force");
  P->verbose   = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v", "--verbose");
  P->threshold = ArgsDouble (1.0,             p, argc, "-t", "--threshold");
  P->weight    = ArgsDouble (0.025,           p, argc, "-w", "--weight");
  P->ignore    = ArgsNum    (10, p, argc, "-i", "--ignore",     0, 999999999);
  P->region    = ArgsNum    (0, p, argc, "-r", "--region-size", 1, 999999999);
  P->level     = ArgsNum    (0, p, argc, "-l", "--level", MIN_LEVEL, MAX_LEVEL);

  P->nModels = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-m") == 0)
      P->nModels += 1;

  if(P->nModels == 0 && P->level == 0)
    P->level = DEFAULT_LEVEL;

  if(P->level != 0){
    xpl = GetLevels(P->level);
    xargc = StrToArgv(xpl, &xargv);
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-m") == 0)
        P->nModels += 1;
    }

  if(P->nModels == 0){
    fprintf(stderr, "Error: at least you need to use a context model!\n");
    return 1;
    }

  P->model = (ModelPar *) Calloc(P->nModels, sizeof(ModelPar));

  k = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-m") == 0)
      P->model[k++] = ArgsUniqModel(argv[n+1], 0);
  if(P->level != 0){
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-m") == 0)
        P->model[k++] = ArgsUniqModel(xargv[n+1], 0);
    }

  P->filename = argv[argc-1];

  if(P->verbose)
    PrintArgs(P);

  // THREADS SETTING
  T = (Threads *) Calloc(2, sizeof(Threads));
  T[0].id = 0;
  T[1].id = 1;

  // COMPRESSING ==============================================================
  //

  if(P->verbose) fprintf(stderr, "Spliting and reversing streams ...\n");

  FILE *IN   = Fopen(P->filename,   "r");
  FILE *OUT1 = Fopen(".dlcr_1.dna", "w");
  FILE *OUT2 = Fopen(".dlcr_2.dna", "w");

  int sym;
  uint64_t nValues = 0;
  while((sym = getc(IN)) != EOF){
    if(sym == '>') // HEADER FOUND!
      while((sym = fgetc(IN)) != EOF && sym != '\n')
        ;
    if(sym == EOF) break;
    if(sym == '\n') continue;
    fprintf(OUT1, "%c", sym);
    ++nValues;
    }
  fclose(IN);
  fclose(OUT1);

  FILE *IN2 = Fopen(".dlcr_1.dna", "r");
  uint64_t nBytes = NBytesInFile(IN2);

  fseek(IN2, -1L, 2);
  while(nBytes--){
    sym = fgetc(IN2); //this moves the logical pointer ahead be 1 char
    fputc(sym, OUT2);
    fseek(IN2, -2L, 1);
    }

  fclose(IN2);
  fclose(OUT2);

  if(P->verbose) fprintf(stderr, "Done!\n");

  // COMPRESSING ==============================================================
  // 
 
  if(P->verbose) fprintf(stderr, "Compressing streams ...\n");
  CompressAction(T);
  if(P->verbose) fprintf(stderr, "Done!\n");

  // GET THE MINIMUM OF LR & RL DIRECTIONS, FILTER & SEGMENT ==================
  //
  
  if(P->verbose) fprintf(stderr, "Filtering and segmenting ...\n");

  FILE *IN_LR = Fopen(".dlcr_1.inf", "r");
  FILE *IN_RL = Fopen(".dlcr_2.inf", "r");

  fseek(IN_RL, 0, SEEK_END );
  char line_LR[1024];
  char line_RL[1024];
  uint64_t idx = 0;
  float smooth;
  int region;
  uint64_t initPos = 0;

  while(fgets_backwards(line_RL, 1024, IN_RL) != NULL){ 
    reverse(line_RL);	  
    if(!fgets(line_LR, 1024, IN_LR)){
      fprintf(stderr, "ERROR: information files have been changed!\n");
      exit(1);
      }

    float RL = atof(line_RL);
    float LR = atof(line_LR);
    float min = RL < LR ? RL : LR;

    if(idx++ == 0){
      smooth = min;
      region = min < P->threshold ? 0 : 1;
      }
    else 
      smooth = smooth - (P->weight * (smooth - min));

    if(smooth >= P->threshold){ // LOW REGION = 0, HIGH REGION = 1
      if(region == 0){              
        region = 1; 
        if(idx - initPos > P->ignore)	
          fprintf(stdout, "%"PRIu64":%"PRIu64"\n", initPos, idx);
        }
      }
    else{ // val < threshold ====> LOW_REGION
      if(region == 1){
        region  = 0;
        initPos = idx;
        }
      }
    }

  if(region == 0){
    if(idx - initPos > P->ignore)	  
    fprintf(stdout, "%"PRIu64":%"PRIu64"\n", initPos, idx);
    }

  fclose(IN_LR);
  fclose(IN_RL);
  
  if(P->verbose) fprintf(stderr, "Done!\n");

  // FINALIZING ===============================================================
  //

  stop = clock();
  fprintf(stderr, "Spent %g sec.\n", ((double) (stop-start)) / CLOCKS_PER_SEC);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
