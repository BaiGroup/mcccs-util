#define VWrap(v, t)                                         \
   if (v.t >= 0.5 * region.t)      v.t -= region.t;         \
   else if (v.t < -0.5 * region.t) v.t += region.t

#define VShift(v, t)                                        \
   if (v.t >= 0.5 * region.t)      shift.t -= region.t;     \
   else if (v.t < -0.5 * region.t) shift.t += region.t

#define VShiftWrap(v, t)                                    \
   if (v.t >= 0.5 * region.t) {                             \
     shift.t -= region.t;                                   \
     v.t -= region.t;                                       \
   } else if (v.t < -0.5 * region.t) {                      \
     shift.t += region.t;                                   \
     v.t += region.t;                                       \
   }

#define VCellWrap(t)                                        \
   if (m2v.t >= cells.t) {                                  \
     m2v.t = 0;                                             \
     shift.t = region.t;                                    \
   } else if (m2v.t < 0) {                                  \
     m2v.t = cells.t - 1;                                   \
     shift.t = - region.t;                                  \
   }

#ifndef NDIM
#define NDIM  3
#endif

#if NDIM == 2

#define VWrapAll(v)                                         \
   {VWrap (v, x);                                           \
   VWrap (v, y);}
#define VShiftAll(v)                                        \
   {VShift (v, x);                                          \
   VShift (v, y);}
#define VCellWrapAll()                                      \
   {VCellWrap (x);                                          \
   VCellWrap (y);}
#define OFFSET_VALS                                         \
   {{0,0}, {1,0}, {1,1}, {0,1}, {-1,1}}
#define N_OFFSET  5

#elif NDIM == 3

#define VWrapAll(v)                                         \
   {VWrap (v, x);                                           \
   VWrap (v, y);                                            \
   VWrap (v, z);}
#define VShiftAll(v)                                        \
   {VShift (v, x);                                          \
   VShift (v, y);                                           \
   VShift (v, z);}
#define VCellWrapAll()                                      \
   {VCellWrap (x);                                          \
   VCellWrap (y);                                           \
   VCellWrap (z);}

#define OFFSET_VALS                                           \
   { {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0},            \
     {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}, {-1,1,1}, {-1,0,1},  \
     {-1,-1,1}, {0,-1,1}, {1,-1,1}                            \
   }

#define N_OFFSET  14

#endif
