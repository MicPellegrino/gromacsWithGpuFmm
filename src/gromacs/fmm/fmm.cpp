/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

/*! \internal \file
 *
 * \brief
 * Interface of \Gromacs to the fast multipole method (FMM).
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \ingroup module_fmm
 *
 */

#include "gmxpre.h"

#include "fmm.h"

#include "gromacs/fmm/fmm-impl.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/strconvert.h"

struct ForceProviderInitOptions;

namespace gmx
{

/*! \brief
 *  Class providing basic infrastructure allowing \Gromacs to interface a FMM.
 *
 *  This class handles FMM-related .mdp input file parameters, including forwarding
 *  them to the underlying FMM implementation class, which, as a force provider, needs
 *  to provide a calculateForces() method.
 */
class Fmm final : public IMDModule,
                  public IMdpOptionProvider
{
    public:
        Fmm() { };

        // From IMDModule
        IMdpOptionProvider *mdpOptionProvider() override { return this; }
        IMDOutputProvider  *outputProvider()    override { return nullptr; }

        void initForceProviders(gmx_unused ForceProviders *forceProviders, ForceProviderInitOptions *options) override
        {
            if (EEL_FMM(options->coulombtype_))
            {
#ifdef GMX_WITH_FMM
                this->impl_ = createFmmImpl(&fmmParameters_, options);
                forceProviders->addForceProvider(impl_.get());
#else
                gmx_fatal(FARGS, "FMM electrostatics requested, but GROMACS was compiled without FMM.\n"
                          "Set -DCMAKE_CXX_FLAGS=-DGMX_WITH_FMM when building GROMACS.");
#endif
            }
        }

        // From IMdpOptionProvider
        void initMdpTransform(IKeyValueTreeTransformRules *transform) override;
        void initMdpOptions(IOptionsContainerWithSections *options) override;
        void buildMdpOutput(KeyValueTreeObjectBuilder *builder) const override;

    private:
        FmmInputParameters fmmParameters_;

#ifdef GMX_WITH_FMM
        // Pointer to the underlying FMM implementation (could for example be the ExaFmm or the Juelich FmSolvr library):
        std::unique_ptr<FmmImpl> impl_;
#endif
};


std::unique_ptr<IMDModule> createFastMultipoleModule()
{
    return std::make_unique<Fmm>();
}


void Fmm::initMdpTransform(IKeyValueTreeTransformRules *rules)
{
    rules->addRule().from<std::string>("/fmm-precision"      ).to<real>("/fmm-precision"      ).transformWith(&fromStdString<real>);
    rules->addRule().from<std::string>("/fmm-override-multipole-order").to<int> ("/fmm-override-multipole-order").transformWith(&fromStdString<int>);
    rules->addRule().from<std::string>("/fmm-override-tree-depth"     ).to<int> ("/fmm-override-tree-depth"     ).transformWith(&fromStdString<int>);
    rules->addRule().from<std::string>("/fmm-override-separation"     ).to<int> ("/fmm-override-separation"     ).transformWith(&fromStdString<int>);
}

void Fmm::initMdpOptions(IOptionsContainerWithSections *options)
{
    // Create the FMM-related options
    options->addOption(RealOption   ("fmm-precision"      ).store(&fmmParameters_.precision_));
    options->addOption(IntegerOption("fmm-override-multipole-order").store(&fmmParameters_.override_order_));
    options->addOption(IntegerOption("fmm-override-tree-depth"     ).store(&fmmParameters_.override_depth_));
    options->addOption(IntegerOption("fmm-override-separation"     ).store(&fmmParameters_.override_separation_));
}

void Fmm::buildMdpOutput(KeyValueTreeObjectBuilder *builder) const
{
    const char *const comment[] = {
        "; Parameters for the FMM electrostatics (in case coulombtype = FMM is chosen)\n"
        "; For the fmsolvr FMM, only fmm-precision should be specified,\n"
        "; all other parameters will be automatically optimized for performance.\n"
        "; For fmm-precision <= 0, the provided values for order, depth, and separation will be used."
    };
    builder->addValue<std::string>("comment-fmm", joinStrings(comment, "\n"));
    builder->addValue<real>("fmm-precision", fmmParameters_.precision_);
    builder->addValue<int>("fmm-override-multipole-order", fmmParameters_.override_order_);
    builder->addValue<int>("fmm-override-tree-depth", fmmParameters_.override_depth_);
    builder->addValue<int>("fmm-override-separation", fmmParameters_.override_separation_);
}

} // namespace gmx
