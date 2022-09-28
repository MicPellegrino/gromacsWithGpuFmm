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
   \brief
   Tests that check whether \Gromacs can use the underlying FMM implementation

   \author Carsten Kutzner <ckutzne@gwdg.de>
   \author R. Thomas Ullmann <tullman@gwdg.de>

   \ingroup module_fmm
 */

#include "gmxpre.h"

#include "gromacs/fmm/fmm.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/units.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"
#include "testutils/testoptions.h"

namespace gmx
{

namespace fmm_test
{

//! \cond
//! flag triggering debug output on user-request via commandline switch [no]debug
bool g_debugBool = false;
GMX_TEST_OPTIONS(MyTestOptions, options)
{
    options->addOption(gmx::BooleanOption("debug").store(&g_debugBool)
                           .description("Print additional debug output if available."));
}
//! \endcond

/* \brief Test whether we get the correct Coulomb energy and forces for a simple geometric arrangement
 *
 */
class FmmInterfaceTest : public ::testing::Test
{
    public:
        //! the constructor sets up an FMM object with the same FMM parameters for all tests
        FmmInterfaceTest()
            : fmm_(createFastMultipoleModule())
        {
            setupFmm();
        }

        /*! \brief
         * Tests whether FMM yields correct results for a simple two-charge system.
         *
         * We put some values into the energy and forces that are passed to the FMM to also
         * make sure that the FMM does not simply overwrite these values but adds to them.
         *
         * \param[in]    ePBC            Type of periodic boundary condition
         * \param[in]    multipoleOrder  If > 0, use this multipole order and not the default one
         * \param[in]    treeDepth       If > 0, use this tree depth and not the default one
         * \param[in]    expectedEnergy  Expected value of the Coulomb energy for this charge configuration
         * \param[in]    expectedForce   Expected value of the Coulomb force in x-direction (others are zero)
         */
        void twoChargesANanometerApart(int    ePBC,
                                       int    multipoleOrder,
                                       int    treeDepth,
                                       double expectedEnergy,
                                       double expectedForce);

    private:
        //! the object implementing the fast multipole method
        std::unique_ptr<IMDModule> fmm_;

        //! set up the FMM object with a common parameter set
        void setupFmm()
        {
            // Fill the module as if from .mdp inputs
            KeyValueTreeBuilder     mdpValues;
            constexpr int           multipoleOrder = (std::is_same<real, double>::value ? 18 : 15);
            const std::string       multipoleOrderStr(std::to_string(multipoleOrder));
            mdpValues.rootObject().addValue("fmm-override-multipole-order", multipoleOrderStr);
            mdpValues.rootObject().addValue("fmm-override-tree-depth", std::to_string(1));

            KeyValueTreeTransformer transform;
            transform.rules()->addRule()
                .keyMatchType("/", StringCompareType::CaseAndDashInsensitive);
            fmm_->mdpOptionProvider()->initMdpTransform(transform.rules());
            auto     result = transform.transform(mdpValues.build(), nullptr);
            Options  moduleOptions;
            fmm_->mdpOptionProvider()->initMdpOptions(&moduleOptions);
            assignOptionsFromKeyValueTree(&moduleOptions, result.object(), nullptr);
        }
};


void FmmInterfaceTest::twoChargesANanometerApart(int    ePBC,
                                                 int    multipoleOrder,
                                                 int    treeDepth,
                                                 double expectedEnergy,
                                                 double expectedForce)
{
    // If requested, override one or both default values as given in the setupFmm() initialization routine
    if (multipoleOrder >= 0 || treeDepth >= 0)
    {
        // Fill the module as if from .mdp inputs
        KeyValueTreeBuilder mdpValues;
        if (multipoleOrder >= 0)
        {
            mdpValues.rootObject().addValue("fmm-multipole-order", std::to_string(multipoleOrder));
        }
        if (treeDepth >= 0)
        {
            mdpValues.rootObject().addValue("fmm-tree-depth", std::to_string(treeDepth));
        }
        KeyValueTreeTransformer transform;
        transform.rules()->addRule().keyMatchType("/", StringCompareType::CaseAndDashInsensitive);
        fmm_->mdpOptionProvider()->initMdpTransform(transform.rules());
        auto     result = transform.transform(mdpValues.build(), nullptr);
        Options  moduleOptions;
        fmm_->mdpOptionProvider()->initMdpOptions(&moduleOptions);
        assignOptionsFromKeyValueTree(&moduleOptions, result.object(), nullptr);
    }

#ifndef GMX_WITH_FMM
    // Common message if the tests are run without an underlying FMM implementation
    fprintf(stderr,  "GROMACS compiled without FMM - cannot test FMM energies & forces!\n");
#else
    // Prepare variables needed for the test
    std::unique_ptr<t_mdatoms>  md             = std::make_unique<t_mdatoms>();
    md->homenr                                 = 2;
    md->chargeA                                = new real[md->homenr];
    md->chargeA[0]                             =  1;
    md->chargeA[1]                             = -1;

    matrix                      box            = { { 3.0, 0.0, 0.0 }, { 0.0, 3.0, 0.0 }, { 0.0, 0.0, 3.0 } };
    // Na, Cl from the NaCl.gro test systems in programs/mdrun/tests
    PaddedVector<RVec>          positions      = { { 1.0, 1.5, 1.5 }, { 2.0, 1.5, 1.5 } };

    // offset to check that GROMACS energy and forces are not overwritten by FMM
    PaddedVector<RVec>          preFmmForces   = { { 1, 2, 3 }, { 4, 5, 6 } };
    PaddedVector<RVec>          forces         = preFmmForces;

    CommrecHandle               commrecHandle  = init_commrec();
    t_commrec                  *cr             = commrecHandle.get();
    std::unique_ptr<t_inputrec> ir             = std::make_unique<t_inputrec>();
    ir->coulombtype                            = eelFMM;
    ir->ePBC                                   = ePBC;

    gmx_mtop_t   mtop;
    mtop.natoms  = 2;
    // offset to check that GROMACS energy and forces are not overwritten by FMM
    double                          preFmmEnergy = 123.456;
    // apply offset
    double                          energy       = preFmmEnergy;
    gmx_enerdata_t                  enerd(1, 0);

    // Test the FMM force provider module
    ForceProviders           forceProvider;
    ForceProviderInitOptions initOptions(ir->ePBC, ir->coulombtype, ir->init_t, ir->delta_t, &mtop, g_debugBool);
    ForceWithVirial          forceWithVirial(forces, TRUE);

    fmm_->initForceProviders(&forceProvider, &initOptions);
    ForceProviderInput  forceProviderInput(positions, *md, 0.0, box, nullptr, TRUE, TRUE, *cr);
    ForceProviderOutput forceProviderOutput(&forceWithVirial, &enerd);
    forceProvider.calculateForces(forceProviderInput, &forceProviderOutput);
    energy += enerd.grpp.ener[egCOULSR][0];

    // Subtract pre-Fmm values:
    energy -= preFmmEnergy;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            forceWithVirial.force_[i][j] -= preFmmForces[i][j];
        }
    }

    // Check the Coulomb energy and forces:
    constexpr double energyTolerance = (std::is_same<real, double>::value ? 5e-5 : 3e-4);
    EXPECT_NEAR(expectedEnergy, energy, energyTolerance);

    constexpr double forceTolerance = (std::is_same<real, double>::value ? 0.0005 : 0.005);
    EXPECT_NEAR(-expectedForce, forceWithVirial.force_[0][XX], forceTolerance);
    EXPECT_NEAR(           0.0, forceWithVirial.force_[0][YY], forceTolerance);
    EXPECT_NEAR(           0.0, forceWithVirial.force_[0][ZZ], forceTolerance);

    EXPECT_NEAR(+expectedForce, forceWithVirial.force_[1][XX], forceTolerance);
    EXPECT_NEAR(           0.0, forceWithVirial.force_[1][YY], forceTolerance);
    EXPECT_NEAR(           0.0, forceWithVirial.force_[1][ZZ], forceTolerance);

    // Clean up
    delete [] md->chargeA;
#endif
}


TEST_F(FmmInterfaceTest, OpenBoundariesDefault)
{
    twoChargesANanometerApart(epbcNONE, -1, -1, -ONE_4PI_EPS0, -ONE_4PI_EPS0);
}

// At depth 0, all interactions are in the near field
TEST_F(FmmInterfaceTest, OpenBoundariesDepth0)
{
    twoChargesANanometerApart(epbcNONE, -1, 0, -ONE_4PI_EPS0, -ONE_4PI_EPS0);
}

TEST_F(FmmInterfaceTest, OpenBoundariesDepth3)
{
    twoChargesANanometerApart(epbcNONE, -1, 3, -ONE_4PI_EPS0, -ONE_4PI_EPS0);
}

TEST_F(FmmInterfaceTest, PeriodicBoundariesDefault)
{
    twoChargesANanometerApart(epbcXYZ, -1, -1, -151.5503252, -109.858909);
}

// At depth 0, all interactions are in the near field
TEST_F(FmmInterfaceTest, PeriodicBoundariesDepth0)
{
    twoChargesANanometerApart(epbcXYZ, -1, 0, -151.5503252, -109.858909);
}

TEST_F(FmmInterfaceTest, PeriodicBoundariesDepth3)
{
    twoChargesANanometerApart(epbcXYZ, -1, 3, -151.5503252, -109.858909);
}
} // namespace fmm_test

} // namespace gmx
