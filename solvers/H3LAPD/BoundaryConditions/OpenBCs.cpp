
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
  const Array<OneD, const int> &traceBndMap =
      m_fields[m_dens_idx]->GetTraceBndMap();

  // Loop over all explists in this region
  auto explists = m_fields[m_dens_idx]->GetBndCondExpansions()[m_bcRegion];
  // std::cout << "Setting region " << m_bcRegion << " BCs" << std::endl;
  for (int e = 0; e < explists->GetExpSize(); ++e) {
    // std::cout << "  Setting exp " << e << " vals" << std::endl;
    //  Current explist
    LR::ExpansionSharedPtr explist = explists->GetExp(e);
    // Offset in the field arrays for this explist
    int explist_offset = explists->GetPhys_Offset(e);
    // Offset in the trace map for this explist
    int trace_offset = m_fields[m_dens_idx]->GetTrace()->GetPhys_Offset(
        traceBndMap[m_offset + e]);

    Array<OneD, NekDouble> xc(explist->GetTotPoints());
    explist->GetCoords(xc);

    Array<OneD, NekDouble> target(explist->GetTotPoints());
    evaluate_expression(explist, target);
    // Loop over points in this explist
    for (int idx_in_explist = 0; idx_in_explist < explist->GetTotPoints();
         idx_in_explist++) {
      int pt_idx = trace_offset + idx_in_explist;

      // std::cout << "xcoord = " << xc[idx_in_explist] << std::endl;

      NekDouble current_val = m_fields[m_dens_idx]->GetPhys()[pt_idx];
      NekDouble target_val = target[idx_in_explist];
      (m_fields[m_field_idx]
           ->GetBndCondExpansions()[m_bcRegion]
           ->UpdatePhys())[explist_offset + idx_in_explist] =
          (0.8 * current_val + 0.2 * target_val);
      // std::cout << "Setting " << 0.8 * current_val + 0.2 * target_val
      //           << std::endl;
      // int foo = 1;
    }
  }
}
} // namespace NESO::Solvers::H3LAPD
