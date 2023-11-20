#ifndef H3LAPD_OPENBCS_H
#define H3LAPD_OPENBCS_H

#include "CustomBCs.h"

namespace LR = Nektar::LocalRegions;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;
namespace NESO::Solvers::H3LAPD {

/**
 * @brief Open boundary conditions class
 */
class OpenBCs : public CustomBCs {
public:
  friend class MemoryManager<OpenBCs>;

  /// Creates an instance of this class
  static CustomBCsSharedPtr
  create(const LU::SessionReaderSharedPtr &pSession,
         const Array<OneD, MR::ExpListSharedPtr> &pFields,
         const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
         const int field_idx, const int pSpaceDim, const int bcRegion,
         const int cnt, SD::BoundaryConditionShPtr cnd) {
    CustomBCsSharedPtr p = MemoryManager<OpenBCs>::AllocateSharedPtr(
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
  OpenBCs(const LU::SessionReaderSharedPtr &pSession,
          const Array<OneD, MR::ExpListSharedPtr> &pFields,
          const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
          const int field_idx, const int pSpaceDim, const int bcRegion,
          const int cnt, SD::BoundaryConditionShPtr cnd);

  virtual ~OpenBCs(void){};

  int m_dens_idx;
  int m_field_idx;
};

} // namespace NESO::Solvers::H3LAPD

#endif
