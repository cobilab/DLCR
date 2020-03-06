#include "cache.h"
#include "defs.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CACHE *CreateCache(uint32_t size){

  CACHE *C     = (CACHE *) Calloc(1, sizeof(CACHE));
  C->pos       = 0;
  C->size      = size;
  C->E         = (C_ENTRY *) Calloc(size, sizeof(C_ENTRY));
  C->E->idx    = 0;
  C->E->idx_ir = 0;
  C->E->s      = 'A';
  C->E->s_ir   = 'T';

  return C;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateCache(CACHE *C, uint64_t i, uint64_t i_ir, uint8_t s, uint8_t s_ir){

  C->E[C->pos].idx = i;
  C->E[C->pos].idx_ir = i_ir;
  C->E[C->pos].s = s;
  C->E[C->pos].s_ir = s_ir;

  if(C->pos == C->size - 1) 
    C->pos = 0;
  else
    C->pos++;

  return;	
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void RemoveCache(CACHE *C){

  Free(C->E);
  Free(C);

  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
