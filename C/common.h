#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <sys/time.h>

#include "util_math.h"

/**Property accumulator*/
typedef struct {
  real val, sum, sum2;
} Prop;

#define PropZero(v)  v.sum = v.sum2 = 0.
#define PropAccum(v)  v.sum += v.val, v.sum2 += Sqr (v.val)
#define PropAvg(v, n) \
   v.sum /= n, v.sum2 = sqrt (Max (v.sum2 / n - Sqr (v.sum), 0.))
#define PropEst(v)  v.sum, v.sum2

/**Error handler*/
enum {ERR_NONE, ERR_BOND_SNAPPED, ERR_CHECKPT_READ, ERR_CHECKPT_WRITE,
   ERR_COPY_BUFF_FULL, ERR_EMPTY_EVPOOL, ERR_MSG_BUFF_FULL,
   ERR_OUTSIDE_REGION, ERR_SNAP_READ, ERR_SNAP_WRITE,
   ERR_SUBDIV_UNFIN, ERR_TOO_MANY_CELLS, ERR_TOO_MANY_COPIES,
   ERR_TOO_MANY_LAYERS, ERR_TOO_MANY_LEVELS, ERR_TOO_MANY_MOLS,
   ERR_TOO_MANY_MOVES, ERR_TOO_MANY_NEBRS, ERR_TOO_MANY_REPLICAS};

char *__errorMsg[] = {"", "bond snapped", "read checkpoint data",
   "write checkpoint data", "copy buffer full", "empty event pool",
   "message buffer full", "outside region", "read snap data",
   "write snap data", "subdivision unfinished", "too many cells",
   "too many copied mols", "too many layers", "too many levels",
   "too many mols", "too many moved mols", "too many neighbors",
   "too many replicas"};

#define ErrExit(code) {                                                 \
        printf ("%s:%d:\tError: %s\n",__FILE__,__LINE__,__errorMsg[code]); \
        exit(code);                                                     \
    }

#define MAX_MPEX_ORD  2
#define I(i, j)  ((i) * ((i) + 1) / 2 + (j))
#define c(i, j)  c[I(i, j)]
#define s(i, j)  s[I(i, j)]
typedef struct {
  real c[I(MAX_MPEX_ORD, MAX_MPEX_ORD) + 1], s[I(MAX_MPEX_ORD, MAX_MPEX_ORD) + 1];
} MpTerms;
typedef struct {
  MpTerms le, me;
  int occ;
} MpCell;

typedef struct {
  real time;
  int left, right, up, circAL, circAR, circBL, circBR, idA, idB;
} EvTree;

#define MOL_LIMIT  10000000
#define NHIST  (NDIM + 2)

#define CHAR_MINUS  '-'
#define CHAR_ZERO   '0'

#ifndef NDIM
#define NDIM  3
#endif
