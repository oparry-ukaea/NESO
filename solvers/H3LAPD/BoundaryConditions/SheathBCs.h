#ifndef H3LAPD_SHEATHBCS_H
#define H3LAPD_SHEATHBCS_H

#include "CustomBCs.h"

namespace LR = Nektar::LocalRegions;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;
namespace NESO::Solvers::H3LAPD {

/**
 * @brief Sheath boundary conditions class
 */
class SheathBCs : public CustomBCs {
public:
  friend class MemoryManager<SheathBCs>;

  /// Creates an instance of this class
  static CustomBCsSharedPtr
  create(const LU::SessionReaderSharedPtr &pSession,
         const Array<OneD, MR::ExpListSharedPtr> &pFields,
         const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
         const int field_idx, const int pSpaceDim, const int bcRegion,
         const int cnt, SD::BoundaryConditionShPtr cnd) {
    CustomBCsSharedPtr p = MemoryManager<SheathBCs>::AllocateSharedPtr(
        pSession, pFields, pTraceNormals, field_idx, pSpaceDim, bcRegion, cnt,
        cnd);
    return p;
  }

  /// Name of the class
  static std::string className;

protected:
  virtual void v_Apply(Array<OneD, Array<OneD, NekDouble>> &Fwd,
                       Array<OneD, Array<OneD, NekDouble>> &physarray,
                       const NekDouble &time);

private:
  SheathBCs(const LU::SessionReaderSharedPtr &pSession,
            const Array<OneD, MR::ExpListSharedPtr> &pFields,
            const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
            const int field_idx, const int pSpaceDim, const int bcRegion,
            const int cnt, SD::BoundaryConditionShPtr cnd);

  virtual ~SheathBCs(void){};

  int m_dens_idx;
  int m_mom_idx;
};

} // namespace NESO::Solvers::H3LAPD

#endif
