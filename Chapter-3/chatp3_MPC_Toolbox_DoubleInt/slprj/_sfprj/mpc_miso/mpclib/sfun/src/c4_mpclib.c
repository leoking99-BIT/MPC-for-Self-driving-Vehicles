/* Include files */

#include <stddef.h>
#include "blas.h"
#include "mpclib_sfun.h"
#include "c4_mpclib.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "mpclib_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)
#define c4_b_p                         (10.0)
#define c4_b_ny                        (1.0)
#define c4_b_nv                        (2.0)
#define c4_b_yoff                      (0.0)
#define c4_b_voff                      (0.0)
#define c4_b_no_md                     (0.0)
#define c4_b_no_ref                    (0.0)
#define c4_b_openloopflag              (FALSE)

/* Variable Declarations */

/* Variable Definitions */
static const char * c4_debug_family_names[16] = { "DataType", "nv", "ny", "p",
  "yoff", "voff", "no_md", "no_ref", "openloopflag", "nargin", "nargout", "ref",
  "md", "rseq", "vseq", "v" };

/* Function Declarations */
static void initialize_c4_mpclib(SFc4_mpclibInstanceStruct *chartInstance);
static void initialize_params_c4_mpclib(SFc4_mpclibInstanceStruct *chartInstance);
static void enable_c4_mpclib(SFc4_mpclibInstanceStruct *chartInstance);
static void disable_c4_mpclib(SFc4_mpclibInstanceStruct *chartInstance);
static void c4_update_debugger_state_c4_mpclib(SFc4_mpclibInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c4_mpclib(SFc4_mpclibInstanceStruct
  *chartInstance);
static void set_sim_state_c4_mpclib(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_st);
static void finalize_c4_mpclib(SFc4_mpclibInstanceStruct *chartInstance);
static void sf_c4_mpclib(SFc4_mpclibInstanceStruct *chartInstance);
static void initSimStructsc4_mpclib(SFc4_mpclibInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c4_machineNumber, uint32_T
  c4_chartNumber);
static const mxArray *c4_sf_marshallOut(void *chartInstanceVoid, void *c4_inData);
static void c4_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_b_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static void c4_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_c_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static void c4_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_d_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static real_T c4_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId);
static void c4_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_e_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static const mxArray *c4_f_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static void c4_info_helper(const mxArray **c4_info);
static const mxArray *c4_emlrt_marshallOut(char * c4_u);
static const mxArray *c4_b_emlrt_marshallOut(uint32_T c4_u);
static void c4_b_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_rseq, const char_T *c4_identifier, real_T c4_y[10]);
static void c4_c_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[10]);
static void c4_d_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_vseq, const char_T *c4_identifier, real_T c4_y[22]);
static void c4_e_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[22]);
static void c4_f_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_v, const char_T *c4_identifier, real_T c4_y[2]);
static void c4_g_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[2]);
static const mxArray *c4_g_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static int32_T c4_h_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId);
static void c4_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static boolean_T c4_i_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId);
static void c4_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static uint8_T c4_j_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_b_is_active_c4_mpclib, const char_T *c4_identifier);
static uint8_T c4_k_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId);
static void init_dsm_address_info(SFc4_mpclibInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c4_mpclib(SFc4_mpclibInstanceStruct *chartInstance)
{
  chartInstance->c4_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c4_is_active_c4_mpclib = 0U;
}

static void initialize_params_c4_mpclib(SFc4_mpclibInstanceStruct *chartInstance)
{
  real_T c4_d0;
  real_T c4_d1;
  real_T c4_d2;
  real_T c4_d3;
  real_T c4_d4;
  real_T c4_d5;
  real_T c4_d6;
  real_T c4_d7;
  sf_set_error_prefix_string(
    "Error evaluating data 'nv' in the parent workspace.\n");
  sf_mex_import_named("nv", sf_mex_get_sfun_param(chartInstance->S, 2, 0),
                      &c4_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c4_nv = c4_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'ny' in the parent workspace.\n");
  sf_mex_import_named("ny", sf_mex_get_sfun_param(chartInstance->S, 3, 0),
                      &c4_d1, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c4_ny = c4_d1;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'p' in the parent workspace.\n");
  sf_mex_import_named("p", sf_mex_get_sfun_param(chartInstance->S, 5, 0), &c4_d2,
                      0, 0, 0U, 0, 0U, 0);
  chartInstance->c4_p = c4_d2;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'yoff' in the parent workspace.\n");
  sf_mex_import_named("yoff", sf_mex_get_sfun_param(chartInstance->S, 7, 0),
                      &c4_d3, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c4_yoff = c4_d3;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'voff' in the parent workspace.\n");
  sf_mex_import_named("voff", sf_mex_get_sfun_param(chartInstance->S, 6, 0),
                      &c4_d4, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c4_voff = c4_d4;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'no_md' in the parent workspace.\n");
  sf_mex_import_named("no_md", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c4_d5, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c4_no_md = c4_d5;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'no_ref' in the parent workspace.\n");
  sf_mex_import_named("no_ref", sf_mex_get_sfun_param(chartInstance->S, 1, 0),
                      &c4_d6, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c4_no_ref = c4_d6;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'openloopflag' in the parent workspace.\n");
  sf_mex_import_named("sf_mex_get_sfun_param", sf_mex_get_sfun_param
                      (chartInstance->S, 4, 0), &c4_d7, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c4_openloopflag = (c4_d7 != 0.0);
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c4_mpclib(SFc4_mpclibInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c4_mpclib(SFc4_mpclibInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c4_update_debugger_state_c4_mpclib(SFc4_mpclibInstanceStruct
  *chartInstance)
{
}

static const mxArray *get_sim_state_c4_mpclib(SFc4_mpclibInstanceStruct
  *chartInstance)
{
  const mxArray *c4_st;
  const mxArray *c4_y = NULL;
  int32_T c4_i0;
  real_T c4_u[10];
  const mxArray *c4_b_y = NULL;
  int32_T c4_i1;
  real_T c4_b_u[2];
  const mxArray *c4_c_y = NULL;
  int32_T c4_i2;
  real_T c4_c_u[22];
  const mxArray *c4_d_y = NULL;
  uint8_T c4_hoistedGlobal;
  uint8_T c4_d_u;
  const mxArray *c4_e_y = NULL;
  real_T (*c4_vseq)[22];
  real_T (*c4_v)[2];
  real_T (*c4_rseq)[10];
  c4_v = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c4_vseq = (real_T (*)[22])ssGetOutputPortSignal(chartInstance->S, 2);
  c4_rseq = (real_T (*)[10])ssGetOutputPortSignal(chartInstance->S, 1);
  c4_st = NULL;
  c4_st = NULL;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_createcellarray(4), FALSE);
  for (c4_i0 = 0; c4_i0 < 10; c4_i0++) {
    c4_u[c4_i0] = (*c4_rseq)[c4_i0];
  }

  c4_b_y = NULL;
  sf_mex_assign(&c4_b_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 1, 10), FALSE);
  sf_mex_setcell(c4_y, 0, c4_b_y);
  for (c4_i1 = 0; c4_i1 < 2; c4_i1++) {
    c4_b_u[c4_i1] = (*c4_v)[c4_i1];
  }

  c4_c_y = NULL;
  sf_mex_assign(&c4_c_y, sf_mex_create("y", c4_b_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_setcell(c4_y, 1, c4_c_y);
  for (c4_i2 = 0; c4_i2 < 22; c4_i2++) {
    c4_c_u[c4_i2] = (*c4_vseq)[c4_i2];
  }

  c4_d_y = NULL;
  sf_mex_assign(&c4_d_y, sf_mex_create("y", c4_c_u, 0, 0U, 1U, 0U, 1, 22), FALSE);
  sf_mex_setcell(c4_y, 2, c4_d_y);
  c4_hoistedGlobal = chartInstance->c4_is_active_c4_mpclib;
  c4_d_u = c4_hoistedGlobal;
  c4_e_y = NULL;
  sf_mex_assign(&c4_e_y, sf_mex_create("y", &c4_d_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c4_y, 3, c4_e_y);
  sf_mex_assign(&c4_st, c4_y, FALSE);
  return c4_st;
}

static void set_sim_state_c4_mpclib(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_st)
{
  const mxArray *c4_u;
  real_T c4_dv0[10];
  int32_T c4_i3;
  real_T c4_dv1[2];
  int32_T c4_i4;
  real_T c4_dv2[22];
  int32_T c4_i5;
  real_T (*c4_rseq)[10];
  real_T (*c4_v)[2];
  real_T (*c4_vseq)[22];
  c4_v = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c4_vseq = (real_T (*)[22])ssGetOutputPortSignal(chartInstance->S, 2);
  c4_rseq = (real_T (*)[10])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c4_doneDoubleBufferReInit = TRUE;
  c4_u = sf_mex_dup(c4_st);
  c4_b_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c4_u, 0)),
                        "rseq", c4_dv0);
  for (c4_i3 = 0; c4_i3 < 10; c4_i3++) {
    (*c4_rseq)[c4_i3] = c4_dv0[c4_i3];
  }

  c4_f_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c4_u, 1)), "v",
                        c4_dv1);
  for (c4_i4 = 0; c4_i4 < 2; c4_i4++) {
    (*c4_v)[c4_i4] = c4_dv1[c4_i4];
  }

  c4_d_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c4_u, 2)),
                        "vseq", c4_dv2);
  for (c4_i5 = 0; c4_i5 < 22; c4_i5++) {
    (*c4_vseq)[c4_i5] = c4_dv2[c4_i5];
  }

  chartInstance->c4_is_active_c4_mpclib = c4_j_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c4_u, 3)), "is_active_c4_mpclib");
  sf_mex_destroy(&c4_u);
  c4_update_debugger_state_c4_mpclib(chartInstance);
  sf_mex_destroy(&c4_st);
}

static void finalize_c4_mpclib(SFc4_mpclibInstanceStruct *chartInstance)
{
}

static void sf_c4_mpclib(SFc4_mpclibInstanceStruct *chartInstance)
{
  int32_T c4_i6;
  int32_T c4_i7;
  int32_T c4_i8;
  real_T c4_hoistedGlobal;
  real_T c4_b_hoistedGlobal;
  real_T c4_ref;
  real_T c4_md;
  uint32_T c4_debug_family_var_map[16];
  char_T c4_DataType[6];
  real_T c4_c_nv;
  real_T c4_c_ny;
  real_T c4_c_p;
  real_T c4_c_yoff;
  real_T c4_c_voff;
  real_T c4_c_no_md;
  real_T c4_c_no_ref;
  boolean_T c4_c_openloopflag;
  real_T c4_nargin = 10.0;
  real_T c4_nargout = 3.0;
  real_T c4_rseq[10];
  real_T c4_vseq[22];
  real_T c4_v[2];
  int32_T c4_i9;
  static char_T c4_cv0[6] = { 'd', 'o', 'u', 'b', 'l', 'e' };

  int32_T c4_i10;
  int32_T c4_i11;
  int32_T c4_i12;
  real_T c4_u;
  const mxArray *c4_y = NULL;
  real_T c4_b_u;
  const mxArray *c4_b_y = NULL;
  real_T c4_c_u;
  const mxArray *c4_c_y = NULL;
  real_T c4_d_u;
  const mxArray *c4_d_y = NULL;
  real_T c4_e_u;
  const mxArray *c4_e_y = NULL;
  real_T c4_f_u;
  const mxArray *c4_f_y = NULL;
  real_T c4_g_u;
  const mxArray *c4_g_y = NULL;
  real_T c4_h_u;
  const mxArray *c4_h_y = NULL;
  real_T c4_i_u;
  const mxArray *c4_i_y = NULL;
  boolean_T c4_j_u;
  const mxArray *c4_j_y = NULL;
  const mxArray *c4_b_v = NULL;
  const mxArray *c4_b_vseq = NULL;
  const mxArray *c4_b_rseq = NULL;
  real_T c4_dv3[10];
  int32_T c4_i13;
  real_T c4_dv4[22];
  int32_T c4_i14;
  real_T c4_dv5[2];
  int32_T c4_i15;
  int32_T c4_i16;
  int32_T c4_i17;
  int32_T c4_i18;
  real_T *c4_b_ref;
  real_T *c4_b_md;
  real_T (*c4_c_rseq)[10];
  real_T (*c4_c_vseq)[22];
  real_T (*c4_c_v)[2];
  c4_c_v = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c4_b_md = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c4_b_ref = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  c4_c_vseq = (real_T (*)[22])ssGetOutputPortSignal(chartInstance->S, 2);
  c4_c_rseq = (real_T (*)[10])ssGetOutputPortSignal(chartInstance->S, 1);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 1U, chartInstance->c4_sfEvent);
  for (c4_i6 = 0; c4_i6 < 10; c4_i6++) {
    _SFD_DATA_RANGE_CHECK((*c4_c_rseq)[c4_i6], 0U);
  }

  for (c4_i7 = 0; c4_i7 < 22; c4_i7++) {
    _SFD_DATA_RANGE_CHECK((*c4_c_vseq)[c4_i7], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c4_b_ref, 2U);
  _SFD_DATA_RANGE_CHECK(*c4_b_md, 3U);
  for (c4_i8 = 0; c4_i8 < 2; c4_i8++) {
    _SFD_DATA_RANGE_CHECK((*c4_c_v)[c4_i8], 4U);
  }

  _SFD_DATA_RANGE_CHECK(chartInstance->c4_nv, 5U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c4_ny, 6U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c4_p, 7U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c4_yoff, 8U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c4_voff, 9U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c4_no_md, 10U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c4_no_ref, 11U);
  _SFD_DATA_RANGE_CHECK((real_T)chartInstance->c4_openloopflag, 12U);
  chartInstance->c4_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 1U, chartInstance->c4_sfEvent);
  c4_hoistedGlobal = *c4_b_ref;
  c4_b_hoistedGlobal = *c4_b_md;
  c4_ref = c4_hoistedGlobal;
  c4_md = c4_b_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 16U, 16U, c4_debug_family_names,
    c4_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_DataType, 0U, c4_f_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_c_nv, 1U, c4_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_c_ny, 2U, c4_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_c_p, 3U, c4_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_c_yoff, 4U, c4_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_c_voff, 5U, c4_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_c_no_md, 6U, c4_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_c_no_ref, 7U, c4_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_c_openloopflag, 8U, c4_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_nargin, 9U, c4_d_sf_marshallOut,
    c4_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_nargout, 10U, c4_d_sf_marshallOut,
    c4_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_ref, 11U, c4_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_md, 12U, c4_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_rseq, 13U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_vseq, 14U, c4_b_sf_marshallOut,
    c4_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_v, 15U, c4_sf_marshallOut,
    c4_sf_marshallIn);
  c4_c_openloopflag = c4_b_openloopflag;
  c4_c_no_ref = c4_b_no_ref;
  c4_c_no_md = c4_b_no_md;
  c4_c_voff = c4_b_voff;
  c4_c_yoff = c4_b_yoff;
  c4_c_p = c4_b_p;
  c4_c_ny = c4_b_ny;
  c4_c_nv = c4_b_nv;
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 3);
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 4);
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 5);
  CV_EML_IF(0, 1, 0, TRUE);
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 7);
  for (c4_i9 = 0; c4_i9 < 6; c4_i9++) {
    c4_DataType[c4_i9] = c4_cv0[c4_i9];
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 8);
  for (c4_i10 = 0; c4_i10 < 10; c4_i10++) {
    c4_rseq[c4_i10] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 9);
  for (c4_i11 = 0; c4_i11 < 22; c4_i11++) {
    c4_vseq[c4_i11] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 10);
  for (c4_i12 = 0; c4_i12 < 2; c4_i12++) {
    c4_v[c4_i12] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 11);
  CV_EML_IF(0, 1, 1, TRUE);
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 12);
  c4_u = c4_ref;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_u, 0, 0U, 0U, 0U, 0), FALSE);
  c4_b_u = c4_md;
  c4_b_y = NULL;
  sf_mex_assign(&c4_b_y, sf_mex_create("y", &c4_b_u, 0, 0U, 0U, 0U, 0), FALSE);
  c4_c_u = c4_b_nv;
  c4_c_y = NULL;
  sf_mex_assign(&c4_c_y, sf_mex_create("y", &c4_c_u, 0, 0U, 0U, 0U, 0), FALSE);
  c4_d_u = c4_b_ny;
  c4_d_y = NULL;
  sf_mex_assign(&c4_d_y, sf_mex_create("y", &c4_d_u, 0, 0U, 0U, 0U, 0), FALSE);
  c4_e_u = c4_b_p;
  c4_e_y = NULL;
  sf_mex_assign(&c4_e_y, sf_mex_create("y", &c4_e_u, 0, 0U, 0U, 0U, 0), FALSE);
  c4_f_u = c4_b_yoff;
  c4_f_y = NULL;
  sf_mex_assign(&c4_f_y, sf_mex_create("y", &c4_f_u, 0, 0U, 0U, 0U, 0), FALSE);
  c4_g_u = c4_b_voff;
  c4_g_y = NULL;
  sf_mex_assign(&c4_g_y, sf_mex_create("y", &c4_g_u, 0, 0U, 0U, 0U, 0), FALSE);
  c4_h_u = c4_b_no_md;
  c4_h_y = NULL;
  sf_mex_assign(&c4_h_y, sf_mex_create("y", &c4_h_u, 0, 0U, 0U, 0U, 0), FALSE);
  c4_i_u = c4_b_no_ref;
  c4_i_y = NULL;
  sf_mex_assign(&c4_i_y, sf_mex_create("y", &c4_i_u, 0, 0U, 0U, 0U, 0), FALSE);
  c4_j_u = c4_b_openloopflag;
  c4_j_y = NULL;
  sf_mex_assign(&c4_j_y, sf_mex_create("y", &c4_j_u, 11, 0U, 0U, 0U, 0), FALSE);
  sf_mex_call_debug("mpcblock_refmd_double_mex", 3U, 10U, 14, c4_y, 14, c4_b_y,
                    14, c4_c_y, 14, c4_d_y, 14, c4_e_y, 14, c4_f_y, 14, c4_g_y,
                    14, c4_h_y, 14, c4_i_y, 14, c4_j_y, &c4_b_rseq, &c4_b_vseq,
                    &c4_b_v);
  c4_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_b_rseq), "rseq", c4_dv3);
  for (c4_i13 = 0; c4_i13 < 10; c4_i13++) {
    c4_rseq[c4_i13] = c4_dv3[c4_i13];
  }

  c4_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_b_vseq), "vseq", c4_dv4);
  for (c4_i14 = 0; c4_i14 < 22; c4_i14++) {
    c4_vseq[c4_i14] = c4_dv4[c4_i14];
  }

  c4_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_b_v), "v", c4_dv5);
  for (c4_i15 = 0; c4_i15 < 2; c4_i15++) {
    c4_v[c4_i15] = c4_dv5[c4_i15];
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, -18);
  _SFD_SYMBOL_SCOPE_POP();
  sf_mex_destroy(&c4_b_rseq);
  sf_mex_destroy(&c4_b_vseq);
  sf_mex_destroy(&c4_b_v);
  for (c4_i16 = 0; c4_i16 < 10; c4_i16++) {
    (*c4_c_rseq)[c4_i16] = c4_rseq[c4_i16];
  }

  for (c4_i17 = 0; c4_i17 < 22; c4_i17++) {
    (*c4_c_vseq)[c4_i17] = c4_vseq[c4_i17];
  }

  for (c4_i18 = 0; c4_i18 < 2; c4_i18++) {
    (*c4_c_v)[c4_i18] = c4_v[c4_i18];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, chartInstance->c4_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_mpclibMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void initSimStructsc4_mpclib(SFc4_mpclibInstanceStruct *chartInstance)
{
}

static void init_script_number_translation(uint32_T c4_machineNumber, uint32_T
  c4_chartNumber)
{
}

static const mxArray *c4_sf_marshallOut(void *chartInstanceVoid, void *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i19;
  real_T c4_b_inData[2];
  int32_T c4_i20;
  real_T c4_u[2];
  const mxArray *c4_y = NULL;
  SFc4_mpclibInstanceStruct *chartInstance;
  chartInstance = (SFc4_mpclibInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  for (c4_i19 = 0; c4_i19 < 2; c4_i19++) {
    c4_b_inData[c4_i19] = (*(real_T (*)[2])c4_inData)[c4_i19];
  }

  for (c4_i20 = 0; c4_i20 < 2; c4_i20++) {
    c4_u[c4_i20] = c4_b_inData[c4_i20];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, FALSE);
  return c4_mxArrayOutData;
}

static void c4_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_v;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[2];
  int32_T c4_i21;
  SFc4_mpclibInstanceStruct *chartInstance;
  chartInstance = (SFc4_mpclibInstanceStruct *)chartInstanceVoid;
  c4_v = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_v), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_v);
  for (c4_i21 = 0; c4_i21 < 2; c4_i21++) {
    (*(real_T (*)[2])c4_outData)[c4_i21] = c4_y[c4_i21];
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_b_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i22;
  real_T c4_b_inData[22];
  int32_T c4_i23;
  real_T c4_u[22];
  const mxArray *c4_y = NULL;
  SFc4_mpclibInstanceStruct *chartInstance;
  chartInstance = (SFc4_mpclibInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  for (c4_i22 = 0; c4_i22 < 22; c4_i22++) {
    c4_b_inData[c4_i22] = (*(real_T (*)[22])c4_inData)[c4_i22];
  }

  for (c4_i23 = 0; c4_i23 < 22; c4_i23++) {
    c4_u[c4_i23] = c4_b_inData[c4_i23];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 1, 22), FALSE);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, FALSE);
  return c4_mxArrayOutData;
}

static void c4_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_vseq;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[22];
  int32_T c4_i24;
  SFc4_mpclibInstanceStruct *chartInstance;
  chartInstance = (SFc4_mpclibInstanceStruct *)chartInstanceVoid;
  c4_vseq = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_vseq), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_vseq);
  for (c4_i24 = 0; c4_i24 < 22; c4_i24++) {
    (*(real_T (*)[22])c4_outData)[c4_i24] = c4_y[c4_i24];
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_c_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i25;
  real_T c4_b_inData[10];
  int32_T c4_i26;
  real_T c4_u[10];
  const mxArray *c4_y = NULL;
  SFc4_mpclibInstanceStruct *chartInstance;
  chartInstance = (SFc4_mpclibInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  for (c4_i25 = 0; c4_i25 < 10; c4_i25++) {
    c4_b_inData[c4_i25] = (*(real_T (*)[10])c4_inData)[c4_i25];
  }

  for (c4_i26 = 0; c4_i26 < 10; c4_i26++) {
    c4_u[c4_i26] = c4_b_inData[c4_i26];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 1, 10), FALSE);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, FALSE);
  return c4_mxArrayOutData;
}

static void c4_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_rseq;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[10];
  int32_T c4_i27;
  SFc4_mpclibInstanceStruct *chartInstance;
  chartInstance = (SFc4_mpclibInstanceStruct *)chartInstanceVoid;
  c4_rseq = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_rseq), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_rseq);
  for (c4_i27 = 0; c4_i27 < 10; c4_i27++) {
    (*(real_T (*)[10])c4_outData)[c4_i27] = c4_y[c4_i27];
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_d_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  real_T c4_u;
  const mxArray *c4_y = NULL;
  SFc4_mpclibInstanceStruct *chartInstance;
  chartInstance = (SFc4_mpclibInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_u = *(real_T *)c4_inData;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, FALSE);
  return c4_mxArrayOutData;
}

static real_T c4_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId)
{
  real_T c4_y;
  real_T c4_d8;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), &c4_d8, 1, 0, 0U, 0, 0U, 0);
  c4_y = c4_d8;
  sf_mex_destroy(&c4_u);
  return c4_y;
}

static void c4_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_nargout;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y;
  SFc4_mpclibInstanceStruct *chartInstance;
  chartInstance = (SFc4_mpclibInstanceStruct *)chartInstanceVoid;
  c4_nargout = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_nargout), &c4_thisId);
  sf_mex_destroy(&c4_nargout);
  *(real_T *)c4_outData = c4_y;
  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_e_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  boolean_T c4_u;
  const mxArray *c4_y = NULL;
  SFc4_mpclibInstanceStruct *chartInstance;
  chartInstance = (SFc4_mpclibInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_u = *(boolean_T *)c4_inData;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_u, 11, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, FALSE);
  return c4_mxArrayOutData;
}

static const mxArray *c4_f_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i28;
  char_T c4_b_inData[6];
  int32_T c4_i29;
  char_T c4_u[6];
  const mxArray *c4_y = NULL;
  SFc4_mpclibInstanceStruct *chartInstance;
  chartInstance = (SFc4_mpclibInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  for (c4_i28 = 0; c4_i28 < 6; c4_i28++) {
    c4_b_inData[c4_i28] = (*(char_T (*)[6])c4_inData)[c4_i28];
  }

  for (c4_i29 = 0; c4_i29 < 6; c4_i29++) {
    c4_u[c4_i29] = c4_b_inData[c4_i29];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 10, 0U, 1U, 0U, 2, 1, 6), FALSE);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, FALSE);
  return c4_mxArrayOutData;
}

const mxArray *sf_c4_mpclib_get_eml_resolved_functions_info(void)
{
  const mxArray *c4_nameCaptureInfo = NULL;
  c4_nameCaptureInfo = NULL;
  sf_mex_assign(&c4_nameCaptureInfo, sf_mex_createstruct("structure", 2, 2, 1),
                FALSE);
  c4_info_helper(&c4_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c4_nameCaptureInfo);
  return c4_nameCaptureInfo;
}

static void c4_info_helper(const mxArray **c4_info)
{
  const mxArray *c4_rhs0 = NULL;
  const mxArray *c4_lhs0 = NULL;
  const mxArray *c4_rhs1 = NULL;
  const mxArray *c4_lhs1 = NULL;
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("mtimes"), "name", "name", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363688678U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c4_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c4_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c4_info, c4_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(1363689356U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c4_info, c4_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c4_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c4_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c4_info, sf_mex_duplicatearraysafe(&c4_lhs1), "lhs", "lhs", 1);
  sf_mex_destroy(&c4_rhs0);
  sf_mex_destroy(&c4_lhs0);
  sf_mex_destroy(&c4_rhs1);
  sf_mex_destroy(&c4_lhs1);
}

static const mxArray *c4_emlrt_marshallOut(char * c4_u)
{
  const mxArray *c4_y = NULL;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c4_u)), FALSE);
  return c4_y;
}

static const mxArray *c4_b_emlrt_marshallOut(uint32_T c4_u)
{
  const mxArray *c4_y = NULL;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c4_y;
}

static void c4_b_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_rseq, const char_T *c4_identifier, real_T c4_y[10])
{
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_rseq), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_rseq);
}

static void c4_c_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[10])
{
  real_T c4_dv6[10];
  int32_T c4_i30;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_dv6, 1, 0, 0U, 1, 0U, 1, 10);
  for (c4_i30 = 0; c4_i30 < 10; c4_i30++) {
    c4_y[c4_i30] = c4_dv6[c4_i30];
  }

  sf_mex_destroy(&c4_u);
}

static void c4_d_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_vseq, const char_T *c4_identifier, real_T c4_y[22])
{
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_vseq), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_vseq);
}

static void c4_e_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[22])
{
  real_T c4_dv7[22];
  int32_T c4_i31;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_dv7, 1, 0, 0U, 1, 0U, 1, 22);
  for (c4_i31 = 0; c4_i31 < 22; c4_i31++) {
    c4_y[c4_i31] = c4_dv7[c4_i31];
  }

  sf_mex_destroy(&c4_u);
}

static void c4_f_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_v, const char_T *c4_identifier, real_T c4_y[2])
{
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_v), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_v);
}

static void c4_g_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[2])
{
  real_T c4_dv8[2];
  int32_T c4_i32;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_dv8, 1, 0, 0U, 1, 0U, 1, 2);
  for (c4_i32 = 0; c4_i32 < 2; c4_i32++) {
    c4_y[c4_i32] = c4_dv8[c4_i32];
  }

  sf_mex_destroy(&c4_u);
}

static const mxArray *c4_g_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_u;
  const mxArray *c4_y = NULL;
  SFc4_mpclibInstanceStruct *chartInstance;
  chartInstance = (SFc4_mpclibInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_u = *(int32_T *)c4_inData;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, FALSE);
  return c4_mxArrayOutData;
}

static int32_T c4_h_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId)
{
  int32_T c4_y;
  int32_T c4_i33;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), &c4_i33, 1, 6, 0U, 0, 0U, 0);
  c4_y = c4_i33;
  sf_mex_destroy(&c4_u);
  return c4_y;
}

static void c4_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_b_sfEvent;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  int32_T c4_y;
  SFc4_mpclibInstanceStruct *chartInstance;
  chartInstance = (SFc4_mpclibInstanceStruct *)chartInstanceVoid;
  c4_b_sfEvent = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_b_sfEvent),
    &c4_thisId);
  sf_mex_destroy(&c4_b_sfEvent);
  *(int32_T *)c4_outData = c4_y;
  sf_mex_destroy(&c4_mxArrayInData);
}

static boolean_T c4_i_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId)
{
  boolean_T c4_y;
  boolean_T c4_b0;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), &c4_b0, 1, 11, 0U, 0, 0U, 0);
  c4_y = c4_b0;
  sf_mex_destroy(&c4_u);
  return c4_y;
}

static void c4_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_c_openloopflag;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  boolean_T c4_y;
  SFc4_mpclibInstanceStruct *chartInstance;
  chartInstance = (SFc4_mpclibInstanceStruct *)chartInstanceVoid;
  c4_c_openloopflag = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_c_openloopflag),
    &c4_thisId);
  sf_mex_destroy(&c4_c_openloopflag);
  *(boolean_T *)c4_outData = c4_y;
  sf_mex_destroy(&c4_mxArrayInData);
}

static uint8_T c4_j_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_b_is_active_c4_mpclib, const char_T *c4_identifier)
{
  uint8_T c4_y;
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_k_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c4_b_is_active_c4_mpclib), &c4_thisId);
  sf_mex_destroy(&c4_b_is_active_c4_mpclib);
  return c4_y;
}

static uint8_T c4_k_emlrt_marshallIn(SFc4_mpclibInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId)
{
  uint8_T c4_y;
  uint8_T c4_u0;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), &c4_u0, 1, 3, 0U, 0, 0U, 0);
  c4_y = c4_u0;
  sf_mex_destroy(&c4_u);
  return c4_y;
}

static void init_dsm_address_info(SFc4_mpclibInstanceStruct *chartInstance)
{
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c4_mpclib_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1220731787U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3834945342U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(419509243U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3998385105U);
}

mxArray *sf_c4_mpclib_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("cvpcgBd8lFaroy8z35qPTH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,8,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(1));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,5,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,5,"type",mxType);
    }

    mxSetField(mxData,5,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,6,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,6,"type",mxType);
    }

    mxSetField(mxData,6,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,7,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,7,"type",mxType);
    }

    mxSetField(mxData,7,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(10);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(22);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c4_mpclib_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c4_mpclib_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c4_mpclib(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x4'type','srcId','name','auxInfo'{{M[1],M[21],T\"rseq\",},{M[1],M[19],T\"v\",},{M[1],M[20],T\"vseq\",},{M[8],M[0],T\"is_active_c4_mpclib\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 4, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c4_mpclib_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc4_mpclibInstanceStruct *chartInstance;
    chartInstance = (SFc4_mpclibInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _mpclibMachineNumber_,
           4,
           1,
           1,
           13,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           ssGetPath(S),
           (void *)S);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          init_script_number_translation(_mpclibMachineNumber_,
            chartInstance->chartNumber);
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_mpclibMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _mpclibMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,2,0,1,"rseq");
          _SFD_SET_DATA_PROPS(1,2,0,1,"vseq");
          _SFD_SET_DATA_PROPS(2,1,1,0,"ref");
          _SFD_SET_DATA_PROPS(3,1,1,0,"md");
          _SFD_SET_DATA_PROPS(4,2,0,1,"v");
          _SFD_SET_DATA_PROPS(5,10,0,0,"nv");
          _SFD_SET_DATA_PROPS(6,10,0,0,"ny");
          _SFD_SET_DATA_PROPS(7,10,0,0,"p");
          _SFD_SET_DATA_PROPS(8,10,0,0,"yoff");
          _SFD_SET_DATA_PROPS(9,10,0,0,"voff");
          _SFD_SET_DATA_PROPS(10,10,0,0,"no_md");
          _SFD_SET_DATA_PROPS(11,10,0,0,"no_ref");
          _SFD_SET_DATA_PROPS(12,10,0,0,"openloopflag");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,2,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,856);
        _SFD_CV_INIT_EML_IF(0,1,0,194,225,703,855);
        _SFD_CV_INIT_EML_IF(0,1,1,447,467,577,702);

        {
          unsigned int dimVector[1];
          dimVector[0]= 10;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_c_sf_marshallOut,(MexInFcnForType)
            c4_c_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 22;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_b_sf_marshallOut,(MexInFcnForType)
            c4_b_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_d_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_d_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_sf_marshallOut,(MexInFcnForType)
            c4_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_d_sf_marshallOut,(MexInFcnForType)c4_d_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_d_sf_marshallOut,(MexInFcnForType)c4_d_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_d_sf_marshallOut,(MexInFcnForType)c4_d_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_d_sf_marshallOut,(MexInFcnForType)c4_d_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_d_sf_marshallOut,(MexInFcnForType)c4_d_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(10,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_d_sf_marshallOut,(MexInFcnForType)c4_d_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(11,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_d_sf_marshallOut,(MexInFcnForType)c4_d_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(12,SF_UINT8,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_e_sf_marshallOut,(MexInFcnForType)c4_f_sf_marshallIn);

        {
          real_T *c4_ref;
          real_T *c4_md;
          real_T (*c4_rseq)[10];
          real_T (*c4_vseq)[22];
          real_T (*c4_v)[2];
          c4_v = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
          c4_md = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c4_ref = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          c4_vseq = (real_T (*)[22])ssGetOutputPortSignal(chartInstance->S, 2);
          c4_rseq = (real_T (*)[10])ssGetOutputPortSignal(chartInstance->S, 1);
          _SFD_SET_DATA_VALUE_PTR(0U, *c4_rseq);
          _SFD_SET_DATA_VALUE_PTR(1U, *c4_vseq);
          _SFD_SET_DATA_VALUE_PTR(2U, c4_ref);
          _SFD_SET_DATA_VALUE_PTR(3U, c4_md);
          _SFD_SET_DATA_VALUE_PTR(4U, *c4_v);
          _SFD_SET_DATA_VALUE_PTR(5U, &chartInstance->c4_nv);
          _SFD_SET_DATA_VALUE_PTR(6U, &chartInstance->c4_ny);
          _SFD_SET_DATA_VALUE_PTR(7U, &chartInstance->c4_p);
          _SFD_SET_DATA_VALUE_PTR(8U, &chartInstance->c4_yoff);
          _SFD_SET_DATA_VALUE_PTR(9U, &chartInstance->c4_voff);
          _SFD_SET_DATA_VALUE_PTR(10U, &chartInstance->c4_no_md);
          _SFD_SET_DATA_VALUE_PTR(11U, &chartInstance->c4_no_ref);
          _SFD_SET_DATA_VALUE_PTR(12U, &chartInstance->c4_openloopflag);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _mpclibMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "GD8JDlchWKMURaoIGLPbgC";
}

static void sf_opaque_initialize_c4_mpclib(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc4_mpclibInstanceStruct*) chartInstanceVar)->S,
    0);
  initialize_params_c4_mpclib((SFc4_mpclibInstanceStruct*) chartInstanceVar);
  initialize_c4_mpclib((SFc4_mpclibInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c4_mpclib(void *chartInstanceVar)
{
  enable_c4_mpclib((SFc4_mpclibInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c4_mpclib(void *chartInstanceVar)
{
  disable_c4_mpclib((SFc4_mpclibInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c4_mpclib(void *chartInstanceVar)
{
  sf_c4_mpclib((SFc4_mpclibInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c4_mpclib(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c4_mpclib((SFc4_mpclibInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c4_mpclib();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c4_mpclib(SimStruct* S, const mxArray *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c4_mpclib();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c4_mpclib((SFc4_mpclibInstanceStruct*)chartInfo->chartInstance,
    mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c4_mpclib(SimStruct* S)
{
  return sf_internal_get_sim_state_c4_mpclib(S);
}

static void sf_opaque_set_sim_state_c4_mpclib(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c4_mpclib(S, st);
}

static void sf_opaque_terminate_c4_mpclib(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc4_mpclibInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_mpclib_optimization_info();
    }

    finalize_c4_mpclib((SFc4_mpclibInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc4_mpclib((SFc4_mpclibInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c4_mpclib(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c4_mpclib((SFc4_mpclibInstanceStruct*)(((ChartInfoStruct *)
      ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c4_mpclib(SimStruct *S)
{
  /* Actual parameters from chart:
     no_md no_ref nv ny openloopflag p voff yoff
   */
  const char_T *rtParamNames[] = { "no_md", "no_ref", "nv", "ny", "openloopflag",
    "p", "voff", "yoff" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));

  /* registration for no_md*/
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0], SS_DOUBLE);

  /* registration for no_ref*/
  ssRegDlgParamAsRunTimeParam(S, 1, 1, rtParamNames[1], SS_DOUBLE);

  /* registration for nv*/
  ssRegDlgParamAsRunTimeParam(S, 2, 2, rtParamNames[2], SS_DOUBLE);

  /* registration for ny*/
  ssRegDlgParamAsRunTimeParam(S, 3, 3, rtParamNames[3], SS_DOUBLE);

  /* registration for openloopflag*/
  ssRegDlgParamAsRunTimeParam(S, 4, 4, rtParamNames[4], SS_BOOLEAN);

  /* registration for p*/
  ssRegDlgParamAsRunTimeParam(S, 5, 5, rtParamNames[5], SS_DOUBLE);

  /* registration for voff*/
  ssRegDlgParamAsRunTimeParam(S, 6, 6, rtParamNames[6], SS_DOUBLE);

  /* registration for yoff*/
  ssRegDlgParamAsRunTimeParam(S, 7, 7, rtParamNames[7], SS_DOUBLE);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_mpclib_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      4);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,4,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,4,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,4);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,4,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,4,3);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=3; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 2; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,4);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1190664269U));
  ssSetChecksum1(S,(1166611586U));
  ssSetChecksum2(S,(3158918617U));
  ssSetChecksum3(S,(259070073U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c4_mpclib(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c4_mpclib(SimStruct *S)
{
  SFc4_mpclibInstanceStruct *chartInstance;
  chartInstance = (SFc4_mpclibInstanceStruct *)utMalloc(sizeof
    (SFc4_mpclibInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc4_mpclibInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c4_mpclib;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c4_mpclib;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c4_mpclib;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c4_mpclib;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c4_mpclib;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c4_mpclib;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c4_mpclib;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c4_mpclib;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c4_mpclib;
  chartInstance->chartInfo.mdlStart = mdlStart_c4_mpclib;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c4_mpclib;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->S = S;
  ssSetUserData(S,(void *)(&(chartInstance->chartInfo)));/* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c4_mpclib_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c4_mpclib(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c4_mpclib(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c4_mpclib(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c4_mpclib_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
