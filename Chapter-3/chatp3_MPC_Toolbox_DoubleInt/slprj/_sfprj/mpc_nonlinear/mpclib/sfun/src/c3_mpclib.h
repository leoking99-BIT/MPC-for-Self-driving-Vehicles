#ifndef __c3_mpclib_h__
#define __c3_mpclib_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc3_mpclibInstanceStruct
#define typedef_SFc3_mpclibInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c3_sfEvent;
  boolean_T c3_isStable;
  boolean_T c3_doneDoubleBufferReInit;
  uint8_T c3_is_active_c3_mpclib;
  boolean_T c3_isQP;
  real_T c3_nu;
  real_T c3_ny;
  real_T c3_degrees;
  real_T c3_Hinv[49];
  real_T c3_Kx[42];
  real_T c3_Ku1[18];
  real_T c3_Kut[90];
  real_T c3_Kr[60];
  real_T c3_Kv[36];
  real_T c3_Mlim[18];
  real_T c3_Mx[126];
  real_T c3_Mu1[54];
  real_T c3_Mv[108];
  real_T c3_z_degrees[7];
  real_T c3_utarget[15];
  real_T c3_p;
  real_T c3_uoff[3];
  real_T c3_yoff[2];
  real_T c3_maxiter;
  real_T c3_nxQP;
  boolean_T c3_openloopflag;
  real_T c3_lims_inport;
  real_T c3_no_umin;
  real_T c3_no_umax;
  real_T c3_no_ymin;
  real_T c3_no_ymax;
  real_T c3_switch_inport;
  real_T c3_no_switch;
  real_T c3_enable_value;
  real_T c3_return_cost;
  real_T c3_H[49];
  real_T c3_return_sequence;
  real_T c3_blocking_moves[5];
  real_T c3_Linv[49];
  real_T c3_Ac[126];
  real_T c3_no_ywt;
  real_T c3_no_duwt;
  real_T c3_no_rhoeps;
  real_T c3_Wy;
  real_T c3_Wdu;
  real_T c3_Jm;
  real_T c3_SuJm;
  real_T c3_I2JmWuI2Jm;
  real_T c3_Su1;
  real_T c3_I1WuI2Jm;
  real_T c3_Sx;
  real_T c3_Hv;
  real_T c3_Wu;
  real_T c3_I1;
} SFc3_mpclibInstanceStruct;

#endif                                 /*typedef_SFc3_mpclibInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c3_mpclib_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c3_mpclib_get_check_sum(mxArray *plhs[]);
extern void c3_mpclib_method_dispatcher(SimStruct *S, int_T method, void *data);

#endif
