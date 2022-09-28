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
 * \brief
 * Tests FMM electrostatics
 *
 * GROMACS needs to be compiled with FMM electrostatics to activate these tests.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \ingroup module_mdrun_integration_tests
 */


#include "gmxpre.h"

#include "gromacs/fmm/fmm.h"  // single spot where GMX_WITH_FMM is defined
#ifdef GMX_WITH_FMM           // Only build these tests if GROMACS was compiled with the FMM

#include <string>
#include <tuple>

#include <gtest/gtest.h>

#include "gromacs/math/units.h"
#include "gromacs/trajectory/energyframe.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"

#include "energyreader.h"
#include "mdruncomparison.h"
#include "moduletest.h"
#include "trajectoryreader.h"
#include "trajectorycomparison.h"


namespace gmx
{
namespace test
{
namespace
{

//! Number of separate FMM ranks to use in these tests
const auto nSeparateFmmRanks = 0;

//! Energy values in .edr files that contain PME or FMM energies
const std::vector<std::string> &energyValuesCoulomb = { "Coulomb (SR)", "Coul. recip." };


//! Helper function to construct an mdp entry for the FMM order
std::string makeFmmOrderMdpStringDependingOnPrecision(int fmmOrderSingle, int fmmOrderDouble)
{
    int order = -1;

    if (std::is_same<real, double>::value)
    {
        // double precision build
        order = fmmOrderDouble;
    }
    else
    {
        // mixed precision build
        order = fmmOrderSingle;
    }

    return "fmm-override-multipole-order = " + std::to_string(order) + '\n';
}

//! Convenience type to provide MD system test parameters as a single argument
typedef std::tuple<const char *, double, double, std::string> SystemAndParameters;

//! Parameterized test fixture for FMM electrostatics
class FmmTest : public MdrunTestFixture,
                public ::testing::WithParamInterface<SystemAndParameters>
{
    public:
        /*! \brief Convenience method to get the first energy entry from the .edr file
         *
         * Any following energy values in the .edr file are ignored!
         */
        real getFirstEnergyValue(const std::string fn, const std::vector<std::string> &namesOfRequiredEnergyFields) const;

        //! \brief Compares all frames of two trajectories
        gmx_bool compareTrajectories(const std::string fn1, const std::string fn2, double absTolerance) const;

        //! \brief Compares energy on a frame-by-frame basis for two energy files
        gmx_bool compareEnergies(const std::string fn1, const std::string fn2, double absTolerance) const;
};


class EnergyDrift : public MdrunTestFixture,
                    public ::testing::WithParamInterface<SystemAndParameters>
{
    public:
        //! \brief Compares drift in total energy for two energy files
        gmx_bool extractEnergyDrift(const std::string fn1, const std::string fn2) const;
};

//! Tests whether FMM and PME yield the same Coulomb energies and forces
//
// Here we start with some simple test systems that should run reasonably fast
INSTANTIATE_TEST_CASE_P(FmmVsPmeFast,
                        FmmTest,
                            ::testing::Values(
                                SystemAndParameters("NaCl",             // Two charges a nanometer apart
#if GMX_DOUBLE
                                                    5e-5, 0.0005,       // Tolerance for energy and forces if compiled in double precision
#else
                                                    2e-4, 0.0015,       // Same for single precision
#endif
                                                    R"(rcoulomb                = 1.4
                                                       rlist                   = 1.4
                                                       rvdw                    = 1.4
                                                       fourierspacing          = 0.025
                                                       ewald-rtol              = 5e-5
                                                       ns_type                 = grid
                                                       nsteps                  = 3
                                                       fmm-override-tree-depth = 0
                                                      )"
                                                    + makeFmmOrderMdpStringDependingOnPrecision(15, 18)
                                                    ),
                                SystemAndParameters("NaCl51",           // A small NaCl system (PME parameters allow forces to be relatively exact)
#if GMX_DOUBLE
                                                    0.01, 0.1,          // double precision energy and forces tolerance
#else
                                                    0.1, 0.15,          // single precision energy and forces tolerance
#endif
                                                    R"(rcoulomb                = 1.2
                                                       rlist                   = 1.2
                                                       rvdw                    = 1.2
                                                       fourierspacing          = 0.08
                                                       pme-order               = 10
                                                       ns_type                 = grid
                                                       nsteps                  = 2
                                                       fmm-override-tree-depth = 1
                                                      )"
                                                    + makeFmmOrderMdpStringDependingOnPrecision(15, 18)
                                                    ),
                                                                // Now the tests with exclusions:
                                SystemAndParameters("spc2",     // Two SPC water molecules
#if GMX_DOUBLE
                                                    1e-5, 1e-3, // double precision energy and forces tolerance
#else
                                                    4e-4, 3e-3, // single precision energy and forces tolerance
#endif
                                                    R"(rcoulomb                = 0.9
                                                       rlist                   = 0.9
                                                       rvdw                    = 0.9
                                                       fourierspacing          = 0.03
                                                       pme-order               = 12
                                                       ns_type                 = grid
                                                       dt                      = 0.005
                                                       nsteps                  = 2
                                                       fmm-override-tree-depth = 0

                                                      )"
                                                    + makeFmmOrderMdpStringDependingOnPrecision(15, 18)
                                                    ),
                                SystemAndParameters("spc5",             // Five SPC water molecules
#if GMX_DOUBLE
                                                    0.5e-3, 0.05,       // double precision energy and forces tolerance
#else
                                                    2.0e-3, 0.05,       // single precision energy and forces tolerance
#endif
                                                    R"(rcoulomb                = 0.9
                                                       rlist                   = 0.9
                                                       rvdw                    = 0.9
                                                       fourierspacing          = 0.03
                                                       pme-order               = 10
                                                       ns_type                 = grid
                                                       dt                      = 0.004
                                                       nsteps                  = 2
                                                       fmm-override-tree-depth = 1
                                                      )"
                                                    + makeFmmOrderMdpStringDependingOnPrecision(15, 18)
                                                    ),
                                SystemAndParameters("spc216",           // 216 SPC water molecules
#if GMX_DOUBLE
                                                    6e-2, 0.25,         // double precision energy and forces tolerance
#else
                                                    0.25, 0.25,         // single precision energy and forces tolerance
#endif
                                                    R"(nstlist                 = 2
                                                       rcoulomb                = 0.9
                                                       rlist                   = 0.9
                                                       rvdw                    = 0.9
                                                       fourierspacing          = 0.07
                                                       pme-order               = 12
                                                       ns_type                 = grid
                                                       dt                      = 0.005
                                                       nsteps                  = 1
                                                       fmm-override-tree-depth = 1

                                                      )"
                                                    + makeFmmOrderMdpStringDependingOnPrecision(15, 18)
                                                    ),
                                SystemAndParameters("spc-and-methanol",
#if GMX_DOUBLE
                                                    2e-5, 0.002,        // double precision energy and forces tolerance
#else
                                                    5e-4, 0.01,         // single precision energy and forces tolerance
#endif
                                                    R"(rcoulomb                = 0.9
                                                       rlist                   = 0.9
                                                       rvdw                    = 0.9
                                                       fourierspacing          = 0.03
                                                       pme-order               = 6
                                                       ns_type                 = grid
                                                       dt                      = 0.004
                                                       constraints             = all-bonds
                                                       nsteps                  = 2
                                                       fmm-override-tree-depth = 2
                                                      )"
                                                    + makeFmmOrderMdpStringDependingOnPrecision(15, 18)
                                                    ),
                                SystemAndParameters("AKE-noWater",
#if GMX_DOUBLE
                                                    0.1, 0.1,           // double precision energy and forces tolerance
#else
                                                    1.0, 1,             // single precision energy and forces tolerance
#endif
                                                    R"(ns_type                 = grid
                                                       rcoulomb                = 1.2
                                                       rlist                   = 1.2
                                                       rvdw                    = 1.2
                                                       fourierspacing          = 0.07
                                                       pme-order               = 12
                                                       constraints             = all-bonds
                                                       dt                      = 0.002
                                                       nsteps                  = 2
                                                       tcoupl                  = Berendsen
                                                       nsttcouple              = 1
                                                       tc-grps                 = System
                                                       tau-t                   = 0.05
                                                       ref-t                   = 300
                                                       fmm-override-tree-depth = 2
                                                      )"
                                                    + makeFmmOrderMdpStringDependingOnPrecision(15, 18)
                                                    )
                                ));


//! Tests whether FMM and PME yield the same Coulomb energies and forces
//
// Here we test more realistic, complex test cases.
INSTANTIATE_TEST_CASE_P(FmmVsPmeComplex,
                        FmmTest,
                            ::testing::Values(
                                SystemAndParameters("AKE",
#if GMX_DOUBLE
                                                    10, 1.0,  // double precision energy and forces tolerance
#else
                                                    400, 3.0, // single precision energy and forces tolerance
#endif
                                                    R"(ns_type                 = grid
                                                       rcoulomb                = 1.2
                                                       rlist                   = 1.2
                                                       rvdw                    = 1.2
                                                       fourierspacing          = 0.07
                                                       pme-order               = 12
                                                       constraints             = all-bonds
                                                       dt                      = 0.002
                                                       nsteps                  = 1
                                                       tcoupl                  = Berendsen
                                                       nsttcouple              = 1
                                                       tc-grps                 = System
                                                       tau-t                   = 0.05
                                                       ref-t                   = 300
                                                       fmm-override-tree-depth = 3
                                                      )"
                                                    + makeFmmOrderMdpStringDependingOnPrecision(10, 14)
                                                    ),
                                SystemAndParameters("TehA",
#if GMX_DOUBLE
                                                    50, 10.0,   // double precision energy and forces tolerance
#else
                                                    4000, 10.0, // single precision energy and forces tolerance
#endif
                                                    R"(fmm-override-multipole-order = 8
                                                       fmm-override-tree-depth      = 4
                                                       ns_type                      = grid
                                                       rcoulomb                     = 1.2
                                                       rlist                        = 1.2
                                                       rvdw                         = 1.2
                                                       fourierspacing               = 0.1
                                                       pme-order                    = 8
                                                       constraints                  = all-bonds
                                                       dt                           = 0.002
                                                       nsteps                       = 1
                                                       tcoupl                       = Berendsen
                                                       nsttcouple                   = 1
                                                       tc-grps                      = System
                                                       tau-t                        = 0.05
                                                       ref-t                        = 300
                                                      )" )
                                )
                        );

TEST_P(FmmTest, CompareEnergyAndForces)
{
    SystemAndParameters sysAndPar = GetParam();

    // Unpack:
    auto mdSystem        = std::get<0> (sysAndPar);
    auto toleranceEnergy = std::get<1> (sysAndPar);
    auto toleranceForces = std::get<2> (sysAndPar);
    auto mdpString       = std::get<3> (sysAndPar);

    auto fn_trr_pme = fileManager_.getTemporaryFilePath("pme.trr");
    auto fn_trr_fmm = fileManager_.getTemporaryFilePath("fmm.trr");
    auto fn_edr_pme = fileManager_.getTemporaryFilePath("pme.edr");
    auto fn_edr_fmm = fileManager_.getTemporaryFilePath("fmm.edr");
    auto fn_log_pme = fileManager_.getTemporaryFilePath("pme.log");
    auto fn_log_fmm = fileManager_.getTemporaryFilePath("fmm.log");
    auto fn_tpr_pme = fileManager_.getTemporaryFilePath("pme.tpr");
    auto fn_tpr_fmm = fileManager_.getTemporaryFilePath("fmm.tpr");

    auto mdpStub =
        R"(coulomb-modifier = None
           pbc              = xyz
           vdw-type         = Cut-off
           nstxout          = 1
           nstvout          = 1
           nstfout          = 1
           nstenergy        = 1
           continuation     = yes
           comm-mode        = Linear
           )";

    runner_.useTopGroAndNdxFromDatabase(mdSystem);

    CommandLine mdrunCaller;
    mdrunCaller.addOption("-npme", nSeparateFmmRanks);

    // PME ---------------------------------------------------------------------
    {
        runner_.fullPrecisionTrajectoryFileName_ = fn_trr_pme;
        runner_.edrFileName_                     = fn_edr_pme;
        runner_.logFileName_                     = fn_log_pme;
        runner_.tprFileName_                     = fn_tpr_pme;
        runner_.useStringAsMdpFile(mdpStub + mdpString + "coulombtype = PME\n");

        ASSERT_EQ(0, runner_.callGrompp());
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
    }

    // FMM ---------------------------------------------------------------------
    {
        runner_.fullPrecisionTrajectoryFileName_ = fn_trr_fmm;
        runner_.edrFileName_                     = fn_edr_fmm;
        runner_.logFileName_                     = fn_log_fmm;
        runner_.tprFileName_                     = fn_tpr_fmm;
        runner_.useStringAsMdpFile(mdpStub + mdpString + "coulombtype = FMM\n");

        ASSERT_EQ(0, runner_.callGrompp());
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
    }

    // Comparison of energies and produced trajectories with forces ------------
    EXPECT_EQ(0, compareEnergies    (fn_edr_fmm, fn_edr_pme, toleranceEnergy));
    EXPECT_EQ(0, compareTrajectories(fn_trr_fmm, fn_trr_pme, toleranceForces));
}


// This test can be used to look at the drift in the total energy with PME
// and FMM over a time span of one nanosecond.
// Run it with
// ./bin/mdrun-test --gtest_filter=FmmVsPmeComplex/EnergyDrift.* -nodelete-temporary-files
// to be able to extract data from the .edr files for plotting.
INSTANTIATE_TEST_CASE_P(FmmVsPmeComplex,
                        EnergyDrift,
                            ::testing::Values(
                                SystemAndParameters("NaClSolution", -1, -1,
                                                    R"(rcoulomb                = 0.9
                                                       rlist                   = 0.9
                                                       rvdw                    = 0.9
                                                       fourierspacing          = 0.113
                                                       pme-order               = 4
                                                       constraints             = all-bonds
                                                       dt                      = 0.002
                                                       nsteps                  = 500000
                                                       lincs-iter              = 2
                                                       DispCorr                = EnerPres
                                                       nstcalcenergy           = 1
                                                       verlet-buffer-tolerance = 1e-5
                                                       pbc                     = xyz
                                                       vdw-type                = Cut-off
                                                       nstenergy               = 10
                                                       comm-mode               = Linear
                                                    )")
                                ));

// Disable this test for now since it runs for quite some time.
TEST_P(EnergyDrift, DISABLED_EnergyDrift)
{
    SystemAndParameters sysAndPar = GetParam();

    // Unpack:
    auto *mdSystem  = std::get<0> (sysAndPar);
    auto  mdpString = std::get<3> (sysAndPar);

    auto  fn_trr_pme = fileManager_.getTemporaryFilePath("pme.trr");
    auto  fn_trr_fmm = fileManager_.getTemporaryFilePath("fmm.trr");
    auto  fn_edr_pme = fileManager_.getTemporaryFilePath("pme.edr");
    auto  fn_edr_fmm = fileManager_.getTemporaryFilePath("fmm.edr");
    auto  fn_log_pme = fileManager_.getTemporaryFilePath("pme.log");
    auto  fn_log_fmm = fileManager_.getTemporaryFilePath("fmm.log");
    auto  fn_tpr_pme = fileManager_.getTemporaryFilePath("pme.tpr");
    auto  fn_tpr_fmm = fileManager_.getTemporaryFilePath("fmm.tpr");
    auto  fn_mdp_pme = fileManager_.getTemporaryFilePath("pme.mdp");
    auto  fn_mdp_fmm = fileManager_.getTemporaryFilePath("fmm.mdp");

    runner_.useTopGroAndNdxFromDatabase(mdSystem);
    // PME ---------------------------------------------------------------------
    {
        runner_.useStringAsMdpFile(mdpString + "coulombtype = PME\n");
        runner_.fullPrecisionTrajectoryFileName_ = fn_trr_pme;
        runner_.edrFileName_                     = fn_edr_pme;
        runner_.logFileName_                     = fn_log_pme;
        runner_.tprFileName_                     = fn_tpr_pme;
        runner_.mdpOutputFileName_               = fn_mdp_pme;


        ASSERT_EQ(0, runner_.callGrompp());
        CommandLine mdrunCaller;
        mdrunCaller.append("mdrun");
        mdrunCaller.append("-v");
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
    }
    // FMM ---------------------------------------------------------------------
    {
        runner_.useStringAsMdpFile(mdpString + "coulombtype = FMM\n");
        runner_.fullPrecisionTrajectoryFileName_ = fn_trr_fmm;
        runner_.edrFileName_                     = fn_edr_fmm;
        runner_.logFileName_                     = fn_log_fmm;
        runner_.tprFileName_                     = fn_tpr_fmm;
        runner_.mdpOutputFileName_               = fn_mdp_fmm;

        ASSERT_EQ(0, runner_.callGrompp());
        CommandLine mdrunCaller;
        mdrunCaller.append("mdrun");
        mdrunCaller.append("-v");
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
    }

    // Comparison of energies and produced trajectories with forces ------------
    ASSERT_EQ(0, extractEnergyDrift(fn_edr_pme, fn_edr_fmm));

}


gmx_bool EnergyDrift::extractEnergyDrift(const std::string fnPme, const std::string fnFmm) const
{
    static gmx_bool bFirst     = TRUE;
    auto            PmeE       = 0.0;
    auto            FmmE       = 0.0;
    auto            PmeE_drift = 0.0;
    auto            FmmE_drift = 0.0;

    fprintf(stderr, "\nComparing Coulomb energies of %s and %s\n", fnPme.c_str(), fnFmm.c_str());
    const std::vector<std::string> &energyFields = { "Total Energy" };

    auto efrPme = openEnergyFileToReadTerms(fnPme, energyFields);
    auto efrFmm = openEnergyFileToReadTerms(fnFmm, energyFields);

    // Loop over both trajectory files and compare the contained data frame-by-frame
    while (efrPme->readNextFrame() && efrFmm->readNextFrame() )
    {
        PmeE = FmmE = 0.0;
        auto frPme = efrPme->frame();
        auto frFmm = efrFmm->frame();
        for (auto &element : energyFields)
        {
            PmeE = frPme.at(element);
        }
        for (auto &element : energyFields)
        {
            FmmE = frFmm.at(element);
        }
        fprintf(stderr, "\r%s total energy (kJ/mol) is     PME: %f   FMM: %f\n", frPme.frameName().c_str(),
                PmeE, FmmE);

        if (bFirst)
        {
            PmeE_drift -= PmeE;
            FmmE_drift -= FmmE;
            bFirst      = FALSE;
        }
    }
    PmeE_drift += PmeE;
    FmmE_drift += FmmE;
    fprintf(stderr, "\nEnergy drift (kJ/mol) is     PME: %f   FMM: %f\n", PmeE_drift, FmmE_drift);

    return 0;
}


real FmmTest::getFirstEnergyValue(const std::string fn, const std::vector<std::string> &namesOfRequiredEnergyFields) const
{
    auto E   = 0.0;
    auto efr = openEnergyFileToReadTerms(fn, namesOfRequiredEnergyFields );
    if (efr->readNextFrame())
    {
        auto fr = efr->frame();
        for (auto &element : namesOfRequiredEnergyFields)
        {
            E += fr.at(element);
        }
    }
    else
    {
        // What TODO if we cannot read the .edr file?
    }
    return E;
}


gmx_bool FmmTest::compareEnergies(const std::string fnPme, const std::string fnFmm, double absTolerance) const
{
    auto efrPme = openEnergyFileToReadTerms(fnPme, energyValuesCoulomb);
    auto efrFmm = openEnergyFileToReadTerms(fnFmm, energyValuesCoulomb);

    // Loop over both trajectory files and compare the contained data frame-by-frame
    while (efrPme->readNextFrame() && efrFmm->readNextFrame() )
    {
        auto PmeE  = 0.0;
        auto FmmE  = 0.0;
        auto frPme = efrPme->frame();
        auto frFmm = efrFmm->frame();
        for (auto &element : energyValuesCoulomb)
        {
            PmeE += frPme.at(element);
            FmmE += frFmm.at(element);
        }
        fprintf(stderr, "\nRead PME frame '%s' with E_PME = %g and FMM frame '%s' with E_FMM = %g (difference %g)\n",
                frPme.frameName().c_str(), PmeE,
                frFmm.frameName().c_str(), FmmE,
                fabs(PmeE-FmmE));

        EXPECT_NEAR(PmeE, FmmE, absTolerance);
    }

    return 0;
}


gmx_bool FmmTest::compareTrajectories(const std::string fn1, const std::string fn2, double absTolerance) const
{

    fprintf(stderr, "\nComparing frames of %s and %s\n", fn1.c_str(), fn2.c_str());

    // Build the manager that will present matching pairs of frames to compare
    FramePairManager<TrajectoryFrameReader>
    trajectoryManager(std::make_unique<TrajectoryFrameReader>(fn1),
                      std::make_unique<TrajectoryFrameReader>(fn2));

    TrajectoryFrameMatchSettings matchSettings;
    matchSettings.mustCompareBox        = true;
    matchSettings.coordinatesComparison = ComparisonConditions::CompareIfReferenceFound;
    matchSettings.velocitiesComparison  = ComparisonConditions::CompareIfReferenceFound;
    matchSettings.forcesComparison      = ComparisonConditions::CompareIfReferenceFound;
    matchSettings.handlePbcIfPossible   = false;
    matchSettings.requirePbcHandling    = false;

    TrajectoryTolerances tolerances(TrajectoryComparison::s_defaultTrajectoryTolerances);
    tolerances.box         = absoluteTolerance(absTolerance);
    tolerances.coordinates = absoluteTolerance(absTolerance);
    tolerances.velocities  = absoluteTolerance(absTolerance);
    tolerances.forces      = absoluteTolerance(absTolerance);

    // Compare the trajectory frames.
    TrajectoryComparison trajectoryComparison(matchSettings, tolerances);
    trajectoryManager.compareAllFramePairs<TrajectoryFrame>(trajectoryComparison);

    return 0;
}


//! Allow fourierspacing to be overwritten in the environment
// to enable some manual experiments without recompilation
inline const char *get_fourierspacing(const char *dflt)
{
    const char *fs = getenv("UNITTEST_FOURIERSPACING");
    return fs && *fs ? fs : dflt;
}

//! Test whether FMM runs produces the expected Coulomb energy for
// a positive and a negative elementary charge at a distance of one nanometer
//
// 2019-04-05 CKu: Disabling this test for now, as the group scheme has
// recently been removed and Verlet not yet supports open boundaries
TEST_F(FmmTest, DISABLED_FmmOpenBoundariesEnergy)
{
    runner_.useStringAsMdpFile(
            R"(cutoff-scheme       = Group
               coulombtype         = FMM
               pbc                 = no
               ns_type             = simple
               rcoulomb            = 1.0e-06
               vdw-type            = Cut-off
               rvdw                = 1.0e-06
               rlist               = 1.0e-06
               continuation        = yes       ; to trigger less NOTEs in grompp output
               comm-mode           = Linear
              )"
            + makeFmmOrderMdpStringDependingOnPrecision(15, 18)
            );
    runner_.useTopGroAndNdxFromDatabase("NaCl");
    ASSERT_EQ(0, runner_.callGrompp());

    // For this test only (as it uses the group scheme) we want to disable FP exceptions.
    // After the test we need to return to the original FP exception setting,
    // to not alter runtime behavior of any of the following tests.
    gmx_fedisableexcept();
    ASSERT_EQ(0, runner_.callMdrun());
    gmx_feenableexcept();

    auto   Ecoul = getFirstEnergyValue(runner_.edrFileName_, energyValuesCoulomb );
#if GMX_DOUBLE
    auto   enerTolerance = 1e-10;
#else
    auto   enerTolerance = 5e-5;
#endif
    EXPECT_REAL_EQ_TOL(-138.935457839084478, Ecoul, absoluteTolerance(enerTolerance));
    EXPECT_REAL_EQ_TOL(       -ONE_4PI_EPS0, Ecoul, absoluteTolerance(enerTolerance));

}


//! Test whether FMM runs with full periodic boundary conditions and yields correct energy
TEST_F(FmmTest, FmmPeriodicBoundariesEnergy)
{
    runner_.useStringAsMdpFile(
            R"(coulombtype                  = FMM
               pbc                          = xyz
               ns_type                      = grid
               rcoulomb                     = 1
               vdw-type                     = Cut-off
               rvdw                         = 1
               rlist                        = 1
               continuation                 = yes
               comm-mode                    = Linear
               fmm-override-multipole-order = 15
               fmm-override-tree-depth      = 2
)");
    runner_.useTopGroAndNdxFromDatabase("NaCl");
    ASSERT_EQ(0, runner_.callGrompp());

    CommandLine mdrunCaller;
    mdrunCaller.addOption("-npme", nSeparateFmmRanks);
    ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));

    auto Ecoul = getFirstEnergyValue(runner_.edrFileName_, energyValuesCoulomb );

    EXPECT_REAL_EQ_TOL(-151.5503252, Ecoul, absoluteTolerance(3e-4));
    // Note: The above result was obtained with FMM, multipole order 80, tree depth 0.
    // Maybe it should be replaced by a more accurate estimate one day.
}


//! Test whether FMM yields a Coulomb energy of exactly zero for a single water molecule
//
// 2019-04-05 CKu: Disabling this test for now, as the group scheme has
// recently been removed and Verlet not yet supports open boundaries
TEST_F(FmmTest, DISABLED_ExclusionsOpen)
{
    auto        EcoulFmm      = -1.0;
    std::string mdpStringBoth =
        R"(cutoff-scheme = Group
           pbc           = no
           ns_type       = simple
           rcoulomb      = 0.9
           vdw-type      = Cut-off
           rvdw          = 0.9
           rlist         = 0.9
           continuation  = yes
           comm-mode     = Linear
)";

    auto       &mdSystem = "spc1a";

    runner_.useTopGroAndNdxFromDatabase(mdSystem);
    runner_.useStringAsMdpFile(mdpStringBoth + "coulombtype = FMM\n" + makeFmmOrderMdpStringDependingOnPrecision(15, 18));
    ASSERT_EQ(0, runner_.callGrompp());
    ASSERT_EQ(0, runner_.callMdrun());
    EcoulFmm = getFirstEnergyValue(runner_.edrFileName_, energyValuesCoulomb );

    // For a single water molecule the energy should be exactly zero, but due to
    // how exactly exclusions are treated it is most likely not
    EXPECT_NEAR(0.0, EcoulFmm, 1e-10);

    fprintf(stderr, "Ecoulomb for %s via FMM is %g (should be exactly zero)\n", mdSystem, EcoulFmm);
}



//! Test whether FMM and PME yield the same Coulomb energies and forces for a single water
// molecule at different x-positions in the box. The water molecule from "spc1c"
// additionally crosses a box boundary.
TEST_F(FmmTest, ExclusionsPeriodic)
{
#if GMX_DOUBLE
    double enerTolerance = 5e-5;
    double trajTolerance = 5e-4;
#else
    double enerTolerance = 1e-3;
    double trajTolerance = 2e-2;
#endif

    auto        EcoulPme      =  0.0;
    auto        EcoulFmm      = -1.0;
    std::string mdpStringBoth =
        R"(pbc           = xyz
           ns_type       = grid
           rcoulomb      = 0.9
           vdw-type      = Cut-off
           rvdw          = 0.9
           rlist         = 0.9
           nstxout       = 0       ; don't compare x in the output!
           nstvout       = 1
           nstfout       = 1
           nstenergy     = 1
           continuation  = yes
           nsteps        = 2
           comm-mode     = Linear
)";

    auto fn_trr_pme = fileManager_.getTemporaryFilePath("pme.trr");
    auto fn_trr_fmm = fileManager_.getTemporaryFilePath("fmm.trr");
    auto fn_edr_pme = fileManager_.getTemporaryFilePath("pme.edr");
    auto fn_edr_fmm = fileManager_.getTemporaryFilePath("fmm.edr");

    // Use PME as the reference here ------------------------------------------------
    runner_.useTopGroAndNdxFromDatabase("spc1a");
    runner_.useStringAsMdpFile(mdpStringBoth +
                               R"(coulomb-modifier = None
                                  coulombtype      = PME
                                  pme-order        = 12
                                  fourierspacing   = 0.015)");
    runner_.fullPrecisionTrajectoryFileName_ = fn_trr_pme;
    runner_.edrFileName_                     = fn_edr_pme;

    ASSERT_EQ(0, runner_.callGrompp());
    ASSERT_EQ(0, runner_.callMdrun());
    EcoulPme = getFirstEnergyValue(runner_.edrFileName_, energyValuesCoulomb ); // should be -0.0706237
    EXPECT_NEAR(EcoulPme, -0.0706237, enerTolerance);

    // Test FMM on several MD systems which differ in the x-position of the water molecule:
    for (auto &mdSystem : { "spc1a", "spc1b", "spc1c" })
    {
        fprintf(stderr, "\n--------------------------------------------------------\n");
        fprintf(stderr, "Testing: %s\n", mdSystem);

        runner_.useTopGroAndNdxFromDatabase(mdSystem);
        runner_.useStringAsMdpFile(mdpStringBoth +
                                   R"(coulombtype            = FMM
                                     fmm-override-tree-depth = 1
                                     )" +
                                   makeFmmOrderMdpStringDependingOnPrecision(16, 18));
        runner_.fullPrecisionTrajectoryFileName_ = fn_trr_fmm;
        runner_.edrFileName_                     = fn_edr_fmm;
        ASSERT_EQ(0, runner_.callGrompp());
        ASSERT_EQ(0, runner_.callMdrun());
        EcoulFmm = getFirstEnergyValue(runner_.edrFileName_, energyValuesCoulomb );

        // Comparison --------------------------------------------------------------
        EXPECT_NEAR(EcoulPme, EcoulFmm, enerTolerance);

        // Comparison of energies and produced trajectories with forces ------------
        EXPECT_EQ(0, compareEnergies    (fn_edr_fmm, fn_edr_pme, enerTolerance));
        EXPECT_EQ(0, compareTrajectories(fn_trr_fmm, fn_trr_pme, trajTolerance));

        fprintf(stderr, "\nEcoulomb for %s via PME is %g\n", mdSystem, EcoulPme);
        fprintf(stderr, "Ecoulomb for %s via FMM is %g\n", mdSystem, EcoulFmm);
    }
}


//! Test whether FMM and PME yield the same Coulomb energies without exclusions
//
// System TehA-noexcl, with deactivated exclusions:
// Ecoulomb for TehA-noexcl via PME is -8.36623e+07
// Ecoulomb for TehA-noexcl via FMM is -8.36622e+07
// FMM with EXCLUSIONS=0:
// Ecoulomb for TehA-noexcl via FMM is -8.36622e+07 (ok, should make no difference)
//
// Comparison to 'normal' system TehA with activated exclusions:
// Ecoulomb for TehA via PME is -3.40738e+06
// Ecoulomb for TehA via FMM is -3.40727e+06
// FMM with EXCLUSIONS=0:
// Ecoulomb for TehA via FMM is -8.36622e+07
//
TEST_F(FmmTest, DISABLED_NoExclusions)
{
    auto        fn_trr_pme = fileManager_.getTemporaryFilePath("pme.trr");
    auto        fn_trr_fmm = fileManager_.getTemporaryFilePath("fmm.trr");
    auto        EcoulPme   =  0.0;
    auto        EcoulFmm   = -1.0;

    std::string mdpStringBoth =
        R"(pbc                          = xyz
           ns_type                      = grid
           rcoulomb                     = 1.2
           vdw-type                     = Cut-off
           rvdw                         = 1.2
           rlist                        = 1.2
           continuation                 = yes
           comm-mode                    = Linear
           nstxout                      = 1
           nstfout                      = 1
           constraints                  = none         ; switch off LINCS
           constraint-algorithm         = Shake
           dt                           = 0.0005       ; time step does not matter for now, we just compute energies + forces once
           define                       = -DFLEXIBLE   ; workaround to switch off settle
           fmm-override-multipole-order = 10
           fmm-override-tree-depth      = 4
)";

    auto &mdSystem = "TehA-noexcl";

    fprintf(stderr, "\n--------------------------------------------------------\n");

    fprintf(stderr, "Testing: %s\n", mdSystem);

    runner_.useTopGroAndNdxFromDatabase(mdSystem);

    // PME ---------------------------------------------------------------------
    {
        runner_.fullPrecisionTrajectoryFileName_ = fn_trr_pme;;
        runner_.useStringAsMdpFile(mdpStringBoth +
                                   "coulomb-modifier = None\n"
                                   "coulombtype      = PME\n"
                                   "fourierspacing   = " + get_fourierspacing("0.12") + "\n");
        ASSERT_EQ(0, runner_.callGrompp());
        ASSERT_EQ(0, runner_.callMdrun());
        EcoulPme = getFirstEnergyValue(runner_.edrFileName_, energyValuesCoulomb );
    }

    // FMM ---------------------------------------------------------------------
    {
        runner_.fullPrecisionTrajectoryFileName_ = fn_trr_fmm;
        runner_.useStringAsMdpFile(mdpStringBoth + "coulombtype = FMM\n");
        ASSERT_EQ(0, runner_.callGrompp());
        ASSERT_EQ(0, runner_.callMdrun());
        EcoulFmm = getFirstEnergyValue(runner_.edrFileName_, energyValuesCoulomb );
    }

    // Comparison of energies and produced trajectories with forces ------------
#if GMX_DOUBLE
    auto tolerance = 2.5;
#else
    auto tolerance = 100.0;
#endif
    EXPECT_EQ(0, compareTrajectories(fn_trr_fmm, fn_trr_pme, tolerance));
    EXPECT_REAL_EQ_TOL(EcoulPme, -83662314, relativeToleranceAsFloatingPoint(1, 1e-5));
    EXPECT_REAL_EQ_TOL(EcoulPme, EcoulFmm, relativeToleranceAsFloatingPoint(1, 1e-5));
    fprintf(stderr, "\nEcoulomb for %s via PME is %g\n", mdSystem, EcoulPme);
    fprintf(stderr, "Ecoulomb for %s via FMM is %g\n", mdSystem, EcoulFmm);
}

// Test energies for A and B state
TEST_F(FmmTest, DISABLED_EnergiesForStatesAB) // Disabling for now, since the respective functionality is not yet implemented
{
    auto       &mdSystem = "cytochrome-AB";

    std::string mdpStringBoth =
        R"(dt                    = 0.004
           nsteps                = 0
           comm_mode             = Angular
           nstcomm               = 1
           comm_grps             = Protein
           nstenergy             = 1
           nstcalcenergy         = 1
           cutoff-scheme         = Verlet
           nstlist               = 5
           pbc                   = xyz
           rlist                 = 1
           rcoulomb              = 1
           vdwtype               = Cut-off
           rvdw                  = 1
           pme_order             = 4
           Tcoupl                = nose-hoover
           nh-chain-length       = 1
           tc-grps               = Protein Solution
           tau_t                 = 0.5     0.5
           ref_t                 = 300     300
           Pcoupl                = Parrinello-Rahman
           tau_p                 = 5
           compressibility       = 4.5e-5
           ref_p                 = 1.0
           refcoord-scaling      = All
           gen_vel               = yes
           gen_temp              = 300.0
           gen_seed              = 1972
           constraints           = all-bonds
           constraint-algorithm  = Lincs
           lincs-order           = 6
)";

    std::string mdpStringFreeEnergy =
        R"(free_energy           = yes
           couple-lambda0        = vdw-q
           couple-lambda1        = vdw-q
           couple-intramol       = no
           init-lambda-state     = -1
           delta_lambda          = 0
           nstdhdl               = 50
)";

    runner_.useTopGroAndNdxFromDatabase(mdSystem);

    CommandLine gromppCaller;
    gromppCaller.addOption("-maxwarn", 7);

    auto referenceECoulStateA =  -433225.3; // computed with PME
    auto referenceECoulStateB =  -426014.9; // computed with PME
    auto toleranceInEnergy    = relativeToleranceAsFloatingPoint(1, 1e-5);

    // PME ---------------------------------------------------------------------
    {
        // PME, A-State --------------------------------------------------------
        runner_.useStringAsMdpFile(mdpStringBoth + mdpStringFreeEnergy +
                                   "init_lambda              = 0.0              \n"   // A-State
                                   "coulombtype              = PME              \n"
                                   "fourierspacing = " + get_fourierspacing("0.09") + "\n");
        ASSERT_EQ(0, runner_.callGrompp(gromppCaller));
        ASSERT_EQ(0, runner_.callMdrun());
        auto EcoulPmeA = getFirstEnergyValue(runner_.edrFileName_, energyValuesCoulomb );
        EXPECT_REAL_EQ_TOL(EcoulPmeA, referenceECoulStateA, toleranceInEnergy);

        // PME, B-State ---------------------------------------------------------
        runner_.useStringAsMdpFile(mdpStringBoth + mdpStringFreeEnergy +
                                   "init_lambda              = 1.0              \n"   // B-State
                                   "coulombtype              = PME              \n"
                                   "fourierspacing = " + get_fourierspacing("0.09") + "\n");
        ASSERT_EQ(0, runner_.callGrompp(gromppCaller));
        ASSERT_EQ(0, runner_.callMdrun());
        auto EcoulPmeB = getFirstEnergyValue(runner_.edrFileName_, energyValuesCoulomb);
        EXPECT_REAL_EQ_TOL(EcoulPmeB, referenceECoulStateB, toleranceInEnergy);
    }

    // FMM ---------------------------------------------------------------------
    {
        // FMM, A-State --------------------------------------------------------
        runner_.useStringAsMdpFile(mdpStringBoth +
                                   "init_lambda              = 0.0              \n"   // A-State
                                   "coulombtype              = FMM              \n");
        ASSERT_EQ(0, runner_.callGrompp(gromppCaller));
        ASSERT_EQ(0, runner_.callMdrun());
        auto EcoulFmmA = getFirstEnergyValue(runner_.edrFileName_, energyValuesCoulomb);
        EXPECT_REAL_EQ_TOL(EcoulFmmA, referenceECoulStateA, toleranceInEnergy);

        // FMM, B-State --------------------------------------------------------
        runner_.useStringAsMdpFile(mdpStringBoth +
                                   "init_lambda              = 1.0              \n"   // B-State
                                   "coulombtype              = FMM              \n");
        ASSERT_EQ(0, runner_.callGrompp(gromppCaller));
        ASSERT_EQ(0, runner_.callMdrun());
        auto EcoulFmmB = getFirstEnergyValue(runner_.edrFileName_, energyValuesCoulomb );
        EXPECT_REAL_EQ_TOL(EcoulFmmB, referenceECoulStateB, toleranceInEnergy);
    }

}


}      // namespace
}      // namespace test
}      // namespace gmx

#endif // GMX_WITH_FMM
