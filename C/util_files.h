#define ReadF(x)       fread  (&x, sizeof (x), 1, fp)
#define WriteF(x)      fwrite (&x, sizeof (x), 1, fp)
#define ReadFN(x, n)   fread  (x, sizeof (x[0]), n, fp)
#define WriteFN(x, n)  fwrite (x, sizeof (x[0]), n, fp)

/**Namelist*/
#define NameVal(x)                                          \
   if (! strncmp (bp, #x, strlen (#x))) {                   \
     bp += strlen (#x);                                     \
     x = strtod (bp, &bp);                                  \
   }

typedef enum {N_I, N_R} VType;

#define NameI(x)  {#x, &x, N_I, sizeof (x) / sizeof (int)}
#define NameR(x)  {#x, &x, N_R, sizeof (x) / sizeof (real)}

typedef struct {
  char *vName;
  void *vPtr;
  VType vType;
  int vLen, vStatus;
} NameList;

#define ValI(x)  {&x, N_I, sizeof (x) / sizeof (int)}
#define ValR(x)  {&x, N_R, sizeof (x) / sizeof (real)}

typedef struct {
  void *vPtr;
  VType vType;
  int vLen;
} ValList;
