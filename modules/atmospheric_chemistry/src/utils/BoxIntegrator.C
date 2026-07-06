//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BoxIntegrator.h"
#include "MCMBoxModel.h"

// -----------------------------------------------------------------------------
// SUNDIALS headers (only pulled in when SUNDIALS is available at build time).
//
// SUNDIALS 7.x requires a SUNContext for every object.  We use a single shared
// context (created on first use) for the lifetime of the simulation.  All
// objects with a SUNContext argument are constructed with this shared context.
// -----------------------------------------------------------------------------
#if defined(HAVE_SUNDIALS)
#include <cstring>                            // for memset
#include <arkode/arkode_arkstep.h>           // ARKStepCreate, ARKStepSetTableNum
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_matrix.h>         // SUNDenseMatrix_*
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#endif // HAVE_SUNDIALS

MooseImplicitIntegrator::MooseImplicitIntegrator(const MCMBoxModel & box_model)
  : _box(box_model)
{
}

Real
MooseImplicitIntegrator::computeResidual(unsigned int species_idx,
                                          const std::vector<Real> & C) const
{
  return _box.getDCdt(species_idx, C);
}

Real
MooseImplicitIntegrator::computeJacobianDiagonal(unsigned int species_idx,
                                                  const std::vector<Real> & C) const
{
  return _box.getJacobianDiagonal(species_idx, C);
}

Real
MooseImplicitIntegrator::computeJacobianOffDiagonal(unsigned int species_idx,
                                                     unsigned int jvar,
                                                     const std::vector<Real> & C) const
{
  return _box.getJacobianOffDiagonal(species_idx, jvar, C);
}

void
MooseImplicitIntegrator::reinit(Real time) const
{
  _box.markDirty();
  _box.setCurrentTime(time);
}

Real
MooseImplicitIntegrator::ppbToMolec() const
{
  return _box.ppbToMolec();
}

// -----------------------------------------------------------------------------
// SUNDIALS CVODE wrapper for box-mode ODE integration
// -----------------------------------------------------------------------------
#if defined(HAVE_SUNDIALS)

/**
 * Analytical Jacobian callback for SUNDIALS.
 *
 * The Jacobian ∂(dC_i/dt)/∂C_j is the cornerstone of stiff ODE convergence.
 * SUNDIALS' built-in finite-difference approximation was catastrophically
 * fragile on the MCM mechanism (600+ species, 30+ orders of magnitude in
 * timescales).  This callback bypasses FD and computes the Jacobian
 * analytically from the reaction-network structure (_iG / _k / StoichMatrix).
 *
 * The chemistry-aware Jacobian formula (see F0AM Jac_eval.m and AtChem2 jfy):
 *
 *   rate_r  = k_r * Π_{m=0}^{2} C[iG[r][m]]    (any reactant == -1 → 1.0)
 *   d rate_r / d C[j] = k_r * (Π_{m: iG[r][m] ≠ j} C[iG[r][m]]) * n_{r,j}
 *
 *   ∂(dC_s/dt) / ∂C_j = Σ_r [ f(r,s) * (∂rate_r / ∂C_j) ]
 *
 * That is exactly what MCMBoxModel::computeJacobianTriplets computes.
 *
 * The dense SUNMatrix J is filled column by column from the triplets.
 * We zero the matrix once, then INSERT_VALUES into (row, col).  Caller
 * avoids ADD-values aggregation to match SUNDIALS dense-linear-solver
 * ownership semantics; PETSc path continues to use ADD_VALUES because it
 * must accommodate repeated calls into the same matrix.
 */
int
SundialsBoxIntegrator::sundialsJacFn(sunrealtype /*t*/,
                                     N_Vector y,
                                     N_Vector /*fy*/,
                                     SUNMatrix J,
                                     void *user_data,
                                     N_Vector /*tmp1*/,
                                     N_Vector /*tmp2*/,
                                     N_Vector /*tmp3*/)
{
  auto * integrator = static_cast<SundialsBoxIntegrator *>(user_data);
  if (!integrator)
    return -1;

  const sunindextype N = N_VGetLength_Serial(y);
  const sunrealtype * yd = N_VGetArrayPointer_Serial(y);
  if (!yd || N == 0)
    return -1;

  // Per-component non-negative clamp (same defensive clamp as RHS).  Applying
  // the clamp here keeps the Jacobian finite when SUNDIALS probes slightly
  // negative y during the Newton update.
  std::vector<Real> C(N);
  for (sunindextype i = 0; i < N; ++i)
    C[i] = std::max(static_cast<Real>(yd[i]), Real(0.0));

  // Get analytical Jacobian triplets from the MCM chemistry kernel.
  std::vector<std::tuple<unsigned int, unsigned int, Real>> Jtriplets;
  integrator->_box.computeJacobianTriplets(C, Jtriplets);

  // Zero the dense SUNMatrix once, then fill with INSERT_VALUES.
  // SUNMatrix_Data(J) returns a contiguous sunrealtype* memory block for
  // a dense SUNMatrix (column-major? no, SUNDense stores row-major).
  sunrealtype * Jdata = SUNDenseMatrix_Data(J);
  if (!Jdata)
    return -1;
  std::memset(Jdata, 0, sizeof(sunrealtype) * N * N);

  // Walk triplets and write into J.  MCM mechanisms are densely
  // interconnected; each species potentially touches every other species.
  for (const auto & [row, col, val] : Jtriplets)
  {
    // Dense SUNDIALS stores J in row-major layout.
    Jdata[row * N + col] = static_cast<sunrealtype>(val);
  }

  return 0;
}

/**
 * SUNDIALS right-hand-side callback (static, bridges C SUNDIALS to C++ Box API).
 *
 * Signature mandated by SUNDIALS: int f(sunrealtype t, N_Vector y,
 *                                        N_Vector ydot, void* user_data).
 * ydot[i] is set to dC[i]/dt via the cached SundialsBoxIntegrator.
 */
int
SundialsBoxIntegrator::sundialsRHSF(sunrealtype /*t*/,
                                     N_Vector y,
                                     N_Vector dy,
                                     void *user_data)
{
  auto * integrator = static_cast<SundialsBoxIntegrator *>(user_data);
  if (!integrator)
    return -1;

  const sunindextype N = N_VGetLength_Serial(y);
  if (N == 0)
    return -1;

  sunrealtype * yd  = N_VGetArrayPointer_Serial(y);
  sunrealtype * dyd = N_VGetArrayPointer_Serial(dy);
  if (!yd || !dyd)
    return -1;

  // IMPORTANT: do NOT call markDirty() here.  The PETSc TS path follows the
  // same convention — evaluateCoefficients() is invoked exactly once per
  // MOOSE timestep in execute() at the step midpoint, and _k is held
  // constant for the entire step.  markDirty() would invalidate the dC/dt
  // cache on every internal CVODE step without triggering a re-evaluation
  // of _k, leading to wrong derivatives and CV_RHSFUNC_FAIL.
  //
  // The callback's job: given y, return dy/dt using the cached parameters
  // (_k, J-values, etc.) already loaded in MCMBoxModel by the wrapper.

  // Per-component non-negative clamp (belt-and-suspenders).
  //
  // CVodeSetConstraints() already enforces C_i ≥ 0 on the solver level, but
  // its enforcement is applied *after* the Newton iteration.  If a predicted
  // y is slightly negative *during* an RHS call, the clamp below prevents
  // negative rate coefficients (sqrt(C_i), log(C_i)) from entering the
  // chemistry kernel one extra time.  This is the same defense KPP, LVODE,
  // and SMVGEAR use.
  std::vector<Real> C(N);
  for (sunindextype i = 0; i < N; ++i)
    C[i] = std::max(static_cast<Real>(yd[i]), Real(0.0));

  for (sunindextype i = 0; i < N; ++i)
    dyd[i] = static_cast<sunrealtype>(integrator->_box.getDCdt(i, C));

  return 0;
}

/**
 * Advance chemistry from t0 to t1 using SUNDIALS ARKODE with an L-stable
 * Diagonally-Implicit Runge-Kutta (DIRK) integrator.
 *
 * Why ARKODE DIRK instead of CVODE BDF:
 *   - CV_BDF (multi-step) 在 MCM 上反复 Newton 不收敛 → 步长被压到 0 →
 *     mxsteps 耗尽 （见 doxygen SUNDIALS 7.x 文档）。
 *   - ARKODE DIRK (single-step, L-stable) 每步只做一次性线性代数
 *     （Newton-like for DIRK stages，但单步式结构让步长控制更鲁棒），
 *     这是 AtChem2 (ROS2/ROS3) 和 KPP 推荐的刚性化学积分路线。
 *   - 这里用的是 ARKODE_ARK324L2SA_DIRK_4_2_3 (4阶, L-stable, 单步,
 *     对刚性化学最鲁棒的 ARKODE 内置表之一)。
 *
 * Uses a dense analytical Jacobian + dense LU linear solver via
 * CVodeSetJacFn/ARKodeSetJacFn.  对于 ≤ ~2000 物种 (dense NJ ² 容许)
 * 表现良好；更大的机制应改用 SUNSparse + KLU 或提供 analytic Jv 乘积。
 *
 * 参见:
 *   KPP drv/cvodemodel.F90 (CVodeSetMaxNumSteps + JacFn 指标),
 *   AtChem2 src/solverFunctions.f90 (jfy),
 *   F0AM Core/Jac_eval.m (Jac = f' * DratesDy).
 *
 * On failure throws mooseError with the SUNDIALS return flag and t0/t1.
 */
void
SundialsBoxIntegrator::solveSundialsCVODE(Real t0, Real t1,
                                           std::vector<Real> & C) const
{
  const sunindextype N = static_cast<sunindextype>(C.size());
  if (N == 0)
    return;

  // One-time SUNContext setup.  SUNDIALS 7.x 要求每个对象一个 SUNContext；
  // 我们每次调用创建一个新的，单线程模式下安全。
  SUNContext ctx = nullptr;
  SUNErrCode ctx_flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (ctx_flag != SUN_SUCCESS)
    mooseError("SundialsBoxIntegrator: SUNContext_Create failed, rc=", ctx_flag);

  // 将浓度向量包装为 SUNDIALS 串行 N_Vector。
  N_Vector y = N_VNew_Serial(N, ctx);
  if (!y)
  {
    SUNContext_Free(&ctx);
    mooseError("SundialsBoxIntegrator: N_VNew_Serial failed (n=", N, ")");
  }
  sunrealtype * yd = N_VGetArrayPointer_Serial(y);
  if (!yd)
  {
    N_VDestroy(y);
    SUNContext_Free(&ctx);
    mooseError("SundialsBoxIntegrator: N_VGetArrayPointer_Serial failed");
  }
  // 初始状态的 floor seed: ARKODE 没有 CVodeSetConstraints API，
  // 所以 IC 中任何 C_i = 0 不做处理会被 solver 当"下降至零"的反常情况
  // (特别是 RHS 用 clamp 时 ARK 内部仍可能对零状态求 J)。
  // 将零 IC 下压到极小的正值 1e-40 molec/cm³ (远低于化学相关性)，
  // 让 ARKODE 能在后续积分中自由地驱动物种归零。
  for (sunindextype i = 0; i < N; ++i)
  {
    Real val = C[i];
    yd[i] = static_cast<sunrealtype>(val > 0.0 ? val : 1.0e-40);
  }

  // 创建 ARKODE 内存块；选择纯隐式 (DIRK)，故 fe = nullptr，fi = 我们的 RHS。
  void * arkode_mem = ARKStepCreate(/* fe */ nullptr,
                                     /* fi = */ &SundialsBoxIntegrator::sundialsRHSF,
                                     static_cast<sunrealtype>(t0),
                                     y, ctx);
  if (!arkode_mem)
  {
    N_VDestroy(y);
    SUNContext_Free(&ctx);
    mooseError("SundialsBoxIntegrator: ARKStepCreate failed");
  }

  // ARKODE 7.x 需要显式注册 user_data（旧版 SUNDIALS 自动传 user_data 的
  // 行为已废弃）。
  int flag = ARKodeSetUserData(arkode_mem,
                                const_cast<void *>(static_cast<const void *>(this)));
  if (flag != ARK_SUCCESS)
  {
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    SUNContext_Free(&ctx);
    mooseError("SundialsBoxIntegrator: ARKodeSetUserData failed, flag=", flag);
  }

  // 设置为纯隐式 (DIRK) 方法 — 因为 fe = nullptr (没有显式部分)。
  // ARKStepSetImplicit 标记只用 fi 回调的 Dirik 表，跳过 ERK 表和 fe 校验。
  // ARKStepSetImplicit 启动纯 DIRK 默认配置 — SUNRISES 使用内置 SDIRK_2_1_2
  // (梯形法则, 2 阶段, 2 阶, L-stable)。不要再用 ARKStepSetTableNum,
  // 因为 SUNRISES 对 ERK+DIRK 表的 stage 数和阶数匹配有严格校验：任何不匹配
  // (包括纯 DIRK 下实际不会使用的 ERK 占位表) 都报 "incompatible Butcher tables"。
  flag = ARKStepSetImplicit(arkode_mem);
  if (flag != ARK_SUCCESS)
  {
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    SUNContext_Free(&ctx);
    mooseError("SundialsBoxIntegrator: ARKStepSetImplicit failed, flag=", flag);
  }
  if (flag != ARK_SUCCESS)
  {
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    SUNContext_Free(&ctx);
    mooseError("SundialsBoxIntegrator: ARKStepSetTableNum failed, flag=", flag);
  }

  // 标量相对和绝对容差 (来自 MOOSE 输入文件 solver_rtol / solver_atol)。
  flag = ARKodeSStolerances(arkode_mem,
                             static_cast<sunrealtype>(_rtol),
                             static_cast<sunrealtype>(_atol));
  if (flag != ARK_SUCCESS)
  {
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    SUNContext_Free(&ctx);
    mooseError("SundialsBoxIntegrator: ARKodeSStolerances failed, flag=", flag);
  }

  // 内步数上限。  KPP 默认 1E6; 我们留 2E6 预算给强刚性 MCM。
  flag = ARKodeSetMaxNumSteps(arkode_mem, 2000000);
  if (flag != ARK_SUCCESS)
  {
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    SUNContext_Free(&ctx);
    mooseError("SundialsBoxIntegrator: ARKodeSetMaxNumSteps failed, flag=", flag);
  }

  // 允许每步内更多次 nonlinear 收敛失败再放弃 — 给步长控制器更多机会缩小 h 重试。
  flag = ARKodeSetMaxConvFails(arkode_mem, 25);
  if (flag != ARK_SUCCESS)
  {
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    SUNContext_Free(&ctx);
    mooseError("SundialsBoxIntegrator: ARKodeSetMaxConvFails failed, flag=", flag);
  }

  // 稠密解析 Jacobian + 稠密 LU 线性求解器。
  SUNMatrix       A  = SUNDenseMatrix(N, N, ctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, ctx);
  if (!A || !LS)
  {
    if (LS) SUNLinSolFree(LS);
    if (A)  SUNMatDestroy(A);
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    SUNContext_Free(&ctx);
    mooseError("SundialsBoxIntegrator: SUNDIALS dense J/LS create failed");
  }

  flag = ARKodeSetLinearSolver(arkode_mem, LS, A);
  if (flag != ARK_SUCCESS)
  {
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    SUNContext_Free(&ctx);
    mooseError("SundialsBoxIntegrator: ARKodeSetLinearSolver failed, flag=", flag);
  }

  // 注册化学-aware 解析 Jacobian 回调。这是 KPP/AtChem2/F0AM 的标准模式 —
  // 避免有限差分在刚性化学中产生的 Newton 不稳定。
  flag = ARKodeSetJacFn(arkode_mem, &SundialsBoxIntegrator::sundialsJacFn);
  if (flag != ARK_SUCCESS)
  {
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    SUNContext_Free(&ctx);
    mooseError("SundialsBoxIntegrator: ARKodeSetJacFn failed, flag=", flag);
  }

  // ---- ARKODE 积分到 t1 (ARK_NORMAL = return at requested time)。----
  sunrealtype t_return = static_cast<sunrealtype>(t0);
  int ark_flag = ARKodeEvolve(arkode_mem, static_cast<sunrealtype>(t1), y,
                               &t_return, ARK_NORMAL);
  if (ark_flag < 0)
  {
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    SUNContext_Free(&ctx);
    mooseError("SundialsBoxIntegrator: ARKodeEvolve failed, flag=", ark_flag,
               " t0=", t0, " t1=", t1);
  }

  // Copy integrated state back into C[]。
  yd = N_VGetArrayPointer_Serial(y);
  for (sunindextype i = 0; i < N; ++i) {
    Real val = static_cast<Real>(yd[i]);
    // 后处理 clamp: ARKODE 没有 CVodeSetConstraints API，但我们希望输出
    // 浓度非负（即便极小的负值也是数值伪迹物理化学不成立）。
    if (val < 0 && val > -1e-20)
      val = 0.0;
    C[i] = val;
  }

  // ---- SUNDIALS step statistics (mirrors PETSc TS diagnostic pattern) ----
  long int nst = 0, nfe = 0, nje = 0, nni = 0, ncfn = 0, netf = 0;
  ARKodeGetNumSteps(arkode_mem, &nst);
  ARKodeGetNumRhsEvals(arkode_mem, /* partition_index */ 0, &nfe);
  ARKodeGetNumJacEvals(arkode_mem, &nje);
  ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  ARKodeGetNumErrTestFails(arkode_mem, &netf);

  // ---- Concentration sanity check (trap negative / NaN / Inf) ----
  Real c_min = std::numeric_limits<Real>::max();
  Real c_max = std::numeric_limits<Real>::lowest();
  std::size_t n_nan = 0, n_neg = 0;
  for (sunindextype i = 0; i < N; ++i)
  {
    if (std::isnan(C[i]) || std::isinf(C[i]))
      ++n_nan;
    else
    {
      if (C[i] < 0.0)
        ++n_neg;
      c_min = std::min(c_min, C[i]);
      c_max = std::max(c_max, C[i]);
    }
  }
  _console << "SUNDIALS ARK step nst=" << nst
           << " nfe=" << nfe << " nje=" << nje << " nni=" << nni
           << " ncfn=" << ncfn << " netf=" << netf
           << " | Cmin=" << c_min << " Cmax=" << c_max
           << " neg=" << n_neg << " nan/inf=" << n_nan << std::endl;

  if (n_nan > 0 || n_neg > 0)
    _console << "SUNDIALS WARNING: " << n_nan << " NaN/Inf + " << n_neg
                   << " negative species detected in mechanism integration"
                   << std::endl;

  // Tear down SUNDIALS objects in reverse order。
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  ARKodeFree(&arkode_mem);
  N_VDestroy(y);
  SUNContext_Free(&ctx);
}

#else  // !HAVE_SUNDIALS

// Stub implementation: this path is taken only when MOOSE is built without
// SUNDIALS.  Calling solveSundialsCVODE() in that case is a user error
// (solver_type=sundials selected but library not linked) — flag it loudly.
void
SundialsBoxIntegrator::solveSundialsCVODE(Real /*t0*/, Real /*t1*/,
                                           std::vector<Real> & /*C*/) const
{
  mooseError("SundialsBoxIntegrator::solveSundialsCVODE called but MOOSE was "
             "built without SUNDIALS.  Recompile with SUNDIALS enabled, or "
             "select a non-sundials solver_type.");
}

#endif // HAVE_SUNDIALS
