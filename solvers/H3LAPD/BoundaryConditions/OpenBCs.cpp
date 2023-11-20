
#include <boost/core/ignore_unused.hpp>

#include "OpenBCs.h"

namespace NESO::Solvers::H3LAPD {

std::string OpenBCs::className = GetCustomBCsFactory().RegisterCreatorFunction(
    "Open", OpenBCs::create,
    "'Open' (extrapolated) BCs for the H3LAPD solver.");

OpenBCs::OpenBCs(const LU::SessionReaderSharedPtr &pSession,
                 const Array<OneD, MR::ExpListSharedPtr> &pFields,
                 const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
                 const int field_idx, const int pSpaceDim, const int bcRegion,
                 const int cnt, SD::BoundaryConditionShPtr cnd)
    : CustomBCs(pSession, pFields, pTraceNormals, field_idx, pSpaceDim,
                bcRegion, cnt, cnd),
      m_field_idx(field_idx) {

  m_dens_idx = m_field_to_index.get_idx("ne");
  ASSERTL0(m_dens_idx >= 0, "Expected a density field called ne - not found")
}

void OpenBCs::v_Apply(Array<OneD, Array<OneD, NekDouble>> &Fwd,
                      Array<OneD, Array<OneD, NekDouble>> &physarray,
                      const NekDouble &time) {
  boost::ignore_unused(time);
}
} // namespace NESO::Solvers::H3LAPD
