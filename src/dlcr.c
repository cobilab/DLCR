#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include "mem.h"
#include "defs.h"
#include "msg.h"
#include "buffer.h"
#include "levels.h"
#include "common.h"
#include "pmodels.h"
#include "context.h"

Parameters *P; // FOR THREAD SHARING

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S O R - - - - - - - - - - - - - -

void Compress(uint8_t id){
  FILE        *Reader  = Fopen("forward", "r");
  char        *name    = concatenate("forward", ".co");
  FILE        *Writter = Fopen(name, "w");
  uint32_t    n, k, cModel, totModels, idxPos;
  uint64_t    compressed = 0, nSymbols = 0, nBases = 0;
  uint8_t     *readBUF, sym, irSym, *pos;
  PModel      **pModel, *MX;
  FloatPModel *PT;
  CMWeight    *WM;
  CBUF        *symbBUF;
  CModel      **cModels;
  uint64_t    i = 0;

  if(P->verbose)
    fprintf(stderr, "Analyzing data and creating models ...\n");

  nBases   = NDNASyminFile (Reader);
  nSymbols = NBytesInFile(Reader);

  // EXTRA MODELS DERIVED FROM EDITS
  totModels = P->nModels;
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].edits != 0)
      totModels += 1;

  pModel        = (PModel  **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n]   = CreatePModel(ALPHABET_SIZE);
  MX            = CreatePModel(ALPHABET_SIZE);
  PT            = CreateFloatPModel(ALPHABET_SIZE);
  WM            = CreateWeightModel(totModels);

  readBUF  = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  symbBUF  = CreateCBuffer(BUFFER_SIZE, BGUARD);

  for(n = 0 ; n < P->nModels ; ++n)
    cModels[n] = CreateCModel(TARGET, P->model[n].ctx, P->model[n].den,
    P->model[n].ir, P->model[n].hashSize, P->model[n].gamma,
    P->model[n].edits, P->model[n].eDen, P->model[n].eGamma);

  // GIVE SPECIFIC GAMMA:
  int pIdx = 0;
  for(n = 0 ; n < P->nModels ; ++n){
    WM->gamma[pIdx++] = cModels[n]->gamma;
    if(P->model[n].edits != 0){
      WM->gamma[pIdx++] = cModels[n]->SUBS.eGamma;
      }
    }

  while((k = fread(readBUF, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){

      if(nSymbols > 100) CalcProgress(nSymbols, ++i);

      sym = readBUF[idxPos];

      // FINAL FILTERING DNA CONTENT
      if(sym != 'A' && sym != 'C' && sym != 'G' && sym != 'T'){
        fprintf(Writter, "2\n"); // FORCE HIGH COMPLEXITY <- UNKNOWN SYMBOL
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
  Free(name);

  /*
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == REFERENCE)
      ResetCModelIdx(cModels[n]);
    else
      FreeCModel(cModels[n]);
      */

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
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t main(int argc, char *argv[]){
  char        **p = *&argv, **xargv, *xpl = NULL;
  int32_t     n, xargc = 0;
  uint32_t    k;
  uint64_t    totalBytes, totalSize;
  clock_t     stop = 0, start = clock();

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

  P->force     = ArgsState  (DEFAULT_FORCE,   p, argc, "-F", "--force");
  P->verbose   = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v", "--verbose");
  P->threshold = ArgsDouble (0, p, argc, "-t", "--threshold");
  P->window    = ArgsNum    (0, p, argc, "-w", "--window-size", 1, 999999999);
  P->region    = ArgsNum    (0, p, argc, "-r", "--region-size", 1, 999999999);
  P->level     = ArgsNum    (0, p, argc, "-l", "--level", MIN_LEVEL, MAX_LEVEL);

  P->nModels  = 0;
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

  if(P->verbose)
    PrintArgs(P);

  Compress(1);

  stop = clock();
  fprintf(stdout, "Spent %g sec.\n", ((double)(stop-start))/CLOCKS_PER_SEC);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
