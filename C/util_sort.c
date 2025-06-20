#include "var_type.h"

void Sort (real *a, int *seq, int n)
{
  real q;
  int i, ir, ixt, j, k;

  for (j = 0; j < n; j ++) seq[j] = j;
  if (n > 1) {
    k = n / 2;
    ir = n - 1;
    while (1) {
      if (k > 0) {
        -- k;
        ixt = seq[k];
        q = a[ixt];
      } else {
        ixt = seq[ir];
        q = a[ixt];
        seq[ir] = seq[0];
        -- ir;
        if (ir == 0) {
          seq[0] = ixt;
          break;
        }
      }
      i = k;
      j = 2 * k + 1;
      while (j <= ir) {
        if (j < ir && a[seq[j]] < a[seq[j + 1]]) ++ j;
        if (q < a[seq[j]]) {
          seq[i] = seq[j];
          i = j;
          j = 2 * j + 1;
        } else j = ir + 1;
      }
      seq[i] = ixt;
    }
  }
}
