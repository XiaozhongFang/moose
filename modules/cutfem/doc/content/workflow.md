# Development Workflow

Standard MOOSE development process for CutFEM module.

## Environment Setup

### 1. Fork MOOSE Repository

```bash
# On GitHub, fork idaholab/moose to your account
# Then clone locally
git clone https://github.com/YOUR_USERNAME/moose.git
cd moose

# Add upstream remote
git remote add upstream https://github.com/idaholab/moose.git
```

### 2. Activate MOOSE Environment

```bash
# Load modules (example for HPC)
module load gcc/11.2.0
module load openmpi/4.1.0
module load cmake/3.22.0

# Or use WSL/Unix shell
export MOOSE_DIR=/path/to/moose
```

## Development Cycle

### Step 1: Create Feature Branch

```bash
git checkout devel
git pull upstream devel
git checkout -b feature/ghost-penalty-stabilization
```

### Step 2: Implement Code

Create files in standard MOOSE locations:

```
modules/cutfem/
├── include/kernels/GhostPenaltyKernel.h
├── src/kernels/GhostPenaltyKernel.C
├── include/CutFEMApp.h
├── src/CutFEMApp.C
├── examples/poisson_with_ghost_penalty.i
└── test/tests/ghost_penalty/test_gp.i
```

### Step 3: Compile

```bash
# Build MOOSE framework first
cd moose
make -j4 METHOD=opt

# Build the module
cd modules/cutfem
make -j4 METHOD=opt
```

Expected output:
```
Compiling C++ (opt)...
Linking module executable...
===========================
CutFEM Module Built Successfully
===========================
```

### Step 4: Run Tests

```bash
cd modules/cutfem/test
../../../moose_test-opt -i tests/ghost_penalty/test_gp.i

# Or via make
cd modules/cutfem
make test
```

### Step 5: Format Code

MOOSE requires specific code formatting:

```bash
cd modules/cutfem
make format
```

Uses `.clang-format` configuration from MOOSE root.

### Step 6: Commit

```bash
git add include/kernels/GhostPenaltyKernel.h
git add src/kernels/GhostPenaltyKernel.C
git add examples/poisson_with_ghost_penalty.i
git add test/tests/ghost_penalty/test_gp.i

git commit -m "Add Ghost Penalty kernel for CutFEM

- Implements face-based gradient jump penalty
- Stabilizes unfitted mesh discretization
- Achieves O(h^-2) condition number bound
- Includes unit tests and example problem"
```

## Pull Request Process

### 1. Push to Your Fork

```bash
git push origin feature/ghost-penalty-stabilization
```

### 2. Create Pull Request

On GitHub:
- Title: "CutFEM Phase 1: Ghost Penalty Stabilization"
- Template (fill completely):
  - **Description**: What does this PR do?
  - **Issue Link**: (if applicable)
  - **Testing**: How was it tested?
  - **Types of Changes**: [✓] New feature, etc.
  - **Checklist**: 
    - [ ] Code compiles without warnings
    - [ ] All tests pass
    - [ ] Code formatted with clang-format
    - [ ] New features documented

### 3. Wait for CIVET Check

MOOSE uses CIVET CI at https://civet.inl.gov

- Automatically triggered on PR
- Runs all tests
- Checks code formatting
- Reports failures

### 4. Address Review Feedback

```bash
# Make requested changes
# (edit files)

# Stage and commit
git add modified_files
git commit -m "Address review feedback: [specific changes]"

# Push to update PR
git push origin feature/ghost-penalty-stabilization
```

### 5. Merge

Once approved:
- GitHub merge button: "Squash and merge"
- This combines all commits into one clean commit on `devel`

## Code Standards

### Header Comments (Doxygen)

```cpp
/**
 * Represents gradient jump penalty stabilization for CutFEM.
 * 
 * This kernel computes:
 * \f$ \int_F c_F h_F^{2(k-1+\gamma)} [D_n^k u] \cdot [D_n^k v] \, dS \f$
 * 
 * where:
 * - \f$ [D_n^k v]_F \f$ is the jump of the k-th normal derivative across face F
 * - \f$ h_F \f$ is the face characteristic length
 * - \f$ \gamma \in [0,1] \f$ is the stabilization parameter
 */
class GhostPenaltyKernel : public InterfaceKernel
{
  // ...
};
```

### Implementation Comments

```cpp
Real
GhostPenaltyKernel::computeQpResidual(Moose::DG_RESIDUAL_TYPE type)
{
  // Penalty coefficient: c_F * h_F^{2(k-1+gamma)}
  Real penalty = _c_F * std::pow(_h_f_squared, _exponent);
  
  // Normal component of gradient
  Real u_jump = computeNormalComponent(_grad_u[_qp]) 
              - computeNormalComponent(_grad_u_neighbor[_qp]);
  
  // Residual contribution
  return penalty * u_jump * _test[_i][_qp];
}
```

### Error Handling

```cpp
if (_gamma < 0.0 || _gamma > 1.0)
  paramError("gamma", "Must be in range [0, 1]");

if (_k < 1)
  paramError("k", "Must be >= 1");
```

## Git Workflow Diagrams

### Branch Strategy

```
upstream/devel  ────────────────────────→ [merged]
                   ↑
                   └─ origin/feature/xxx (your branch)
                          ↑
                      [develop locally]
                          ↓
                   (push) → [PR review] → [merge]
```

### Commit Flow

```
Start Feature
    ↓
Code → Compile → Test → Format → Commit
    ↓
    └─ Fix Issues? → Loop back
    ↓
Push → Review → Address Feedback → Merge
```

## Common Tasks

### Sync with Upstream

```bash
git fetch upstream
git rebase upstream/devel
```

### Undo Last Commit

```bash
git reset --soft HEAD~1
# Files remain staged; re-edit as needed
```

### Check Test Results Locally

```bash
cd modules/cutfem/test
../../../moose_test-opt
```

## Performance Profiling

### Monitor Solver Time

```ini
[Executioner]
  petsc_options = '-log_view'
[]
```

Output shows:
- Matrix assembly time
- Linear solver iterations
- Preconditioner setup time

### PAPI Profiling (Advanced)

```bash
LIBPAPI_EVENTS="PAPI_FP_OPS" moose_opt -i input.i
```

Tracks floating-point operations count.

## Troubleshooting

### Compilation Errors

```bash
# Clean and rebuild
make clobber
make -j4 METHOD=opt 2>&1 | tee build.log

# Check for missing dependencies
grep -i error build.log
```

### Test Failures

```bash
# Run test with verbose output
moose_test-opt -i test.i -v

# Compare with gold standard
diff test_output.csv test_gold.csv
```

### Git Conflicts

```bash
# During rebase
git status  # Shows conflicted files
# Edit files manually
git add file_1.C file_2.h
git rebase --continue
```

## Resources

- **MOOSE Documentation**: https://mooseframework.inl.gov
- **Contributing Guide**: https://mooseframework.inl.gov/contribute/
- **GitHub**: https://github.com/idaholab/moose
- **Discussions**: https://github.com/idaholab/moose/discussions

## Next Steps

1. Implement GhostPenaltyKernel.C (use implementation.md)
2. Create test_gp.i (use example as template)
3. Verify compilation and tests pass
4. Submit PR following template above
