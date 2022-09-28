#include "gmxpre.h"

#include "fmm-impl-helpers.h"

#include "gromacs/topology/topology.h"
#include "gromacs/utility/gmxassert.h"

//! \brief Optionally output debug information
#define DEBUG_PRINT(fmt, ...)
//#define DEBUG_PRINT(fmt, ...) { fprintf(stderr, fmt, __VA_ARGS__); }


namespace gmx
{

std::vector< std::vector<int> > makeExclusionList(const gmx_mtop_t *mtop)
{
    auto exclusions = std::vector< std::vector<int> > (mtop->natoms);

    // Loop over all molecule types
    int globalInd = 0;

    DEBUG_PRINT("Total molecule types: %lu\n", mtop->moltype.size());

    for (auto iMolType = 0u; iMolType < mtop->moltype.size(); iMolType++)
    {
        // Loop over each copy of this molecule
        int nMol = mtop->molblock[iMolType].nmol;
        for (int iCopy = 0; iCopy < nMol; iCopy++)
        {
            const gmx_moltype_t *mt              = &mtop->moltype[iMolType];
            int                  molFirstAtomInd = globalInd;

            DEBUG_PRINT("Molecule number %3d, copy %3d of %d '%s', first atom %d\n", iMolType, iCopy + 1, nMol, *mt->name, molFirstAtomInd);

            // Loop over all atoms iAtom of the molecule
            int *molExcl = mt->excls.index;
            GMX_ASSERT(mt->excls.nr == mt->atoms.nr, "excls.nr != atoms.nr");
            //DEBUG_PRINT("exclusionsnr %4d\n", mt->excls.nr);
            int count = 0;
            for (int iAtom = 0; iAtom < mt->excls.nr; iAtom++)
            {
                int nExcl = molExcl[iAtom+1] - molExcl[iAtom];
                exclusions[globalInd] = std::vector<int> (nExcl);
                //DEBUG_PRINT("particle %4d has %3d exclusions with indices", globalInd, nExcl);

                for (int iExcl = 0; iExcl < nExcl; iExcl++)
                {
                    int excluded = mt->excls.a[count] + molFirstAtomInd;
                    //DEBUG_PRINT(" %3d ", excluded);
                    exclusions[globalInd][iExcl] = excluded;
                    count++;
                }
                //DEBUG_PRINT("%s", "\n");
                globalInd++;
            }
            GMX_ASSERT(count == mt->excls.nra, "Number of exclusions does not match");
        }
    }
    return exclusions;
}


} // namespace gmx
