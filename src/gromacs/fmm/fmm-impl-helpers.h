/*! \libinternal \file
 *
 * \brief
 * Common utility functions used by various FMM implementations.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \inlibraryapi
 * \ingroup module_fmm
 */
#ifndef GMX_FMM_FMM_HELPERS_H
#define GMX_FMM_FMM_HELPERS_H

#include <cstddef>

#include <vector>

struct gmx_mtop_t;

namespace gmx
{

/*! \brief Create a list of exclusions for every particle in the system
 *
 * This function returns a vector that has as many entries as there are particles
 * in the topology. For each particle the number of the excluded particles are given.
 * The returned vector is excl[i][j] with i = 0..n_atoms and j = 0..n_exclusions(i).
 */
std::vector< std::vector<int> > makeExclusionList(const gmx_mtop_t *mtop);

}      // namespace gmx

#endif // GMX_FMM_FMM_HELPERS_H
