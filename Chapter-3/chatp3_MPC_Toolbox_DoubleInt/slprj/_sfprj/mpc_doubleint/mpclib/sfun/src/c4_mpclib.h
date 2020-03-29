#ifndef __c4_mpclib_h__
#define __c4_mpclib_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc4_mpclibInstanceStruct
#define typedef_SFc4_mpclibInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c4_sfEvent;
  boolean_T c4_isStable;
  boolean_T c4_doneDoubleBufferReInit;
  uint8_T c4_is_active_c4_mpclib;
  real_T c4_nv;
  real_T c4_ny;
  real_T c4_p;
  real_T c4_yoff;
  real_T c4_voff;
  real_T c4_no_md;
  real_T c4_no_ref;
  boolean_T c4_openloopflag;
} SFc4_mpclibInstanceStruct;

#endif                                 /*typedef_SFc4_mpclibInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c4_mpclib_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c4_mpclib_get_check_sum(mxArray *plhs[]);
extern void c4_mpclib_method_dispatcher(SimStruct *S, int_T method, void *data);

#endif
