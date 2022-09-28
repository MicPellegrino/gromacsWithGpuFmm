/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares gmx::IMDModule.
 *
 * See \ref page_mdmodules for an overview of this and associated interfaces.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_IMDMODULE_H
#define GMX_MDTYPES_IMDMODULE_H
#include <string>

#include "gromacs/utility/basedefinitions.h"

struct ForceProviders;
struct gmx_mtop_t;
struct t_inputrec;

class t_state;
struct ObservablesHistory;

namespace gmx
{

class IMDOutputProvider;
class IMdpOptionProvider;
struct MdrunOptions;

/*! \libinternal\brief
    Temporary structure with data that some force providers need during init
 */
struct ForceProviderInitOptions
{
    //! type of periodic boundary conditions
    const int                     ePBC_                  = 0;
    //! type of electrostatics treatment
    const int                     coulombtype_           = 0;
    //! initial time
    const double                  initTime_              = 0.0;              
    //! the integration timestep
    const double                  timeStep_              = 0;                
    //! global, complete system topology struct
    const gmx_mtop_t * const      mtop_                  = nullptr;          
    //! mdrun options, stating e.g, whether we restarted from a checkpoint file
    bool                          beVerbose_             = false;            

    /*! \brief constructor

        \param[in]   ePBC                    \copybrief ePBC_
        \param[in]   coulombtype             \copybrief coulombtype_
        \param[in]   it                      \copybrief initTime_
        \param[in]   dt                      \copybrief timeStep_
        \param[in]   mtop                    \copybrief mtop_
        \param[in]   beVerbose               \copybrief beVerbose_
     */
    ForceProviderInitOptions(int ePBC, int coulombtype, double it, double dt,
                             const gmx_mtop_t *mtop,
                             bool beVerbose)
        : ePBC_(ePBC), coulombtype_(coulombtype),
          initTime_(it), timeStep_(dt),
          mtop_(mtop), beVerbose_(beVerbose)
        { }

    //! default constructor has no initializer list and thus keeps default values
    ForceProviderInitOptions() noexcept {}
};


/*! \libinternal \brief
 * Extension module for \Gromacs simulations.
 *
 * The methods that return other interfaces can in the future return null for
 * those interfaces that the module does not need to implement, but currently
 * the callers are not prepared to generically handle various cases.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
class IMDModule
{
    public:
        virtual ~IMDModule() {}

        //! Returns an interface for handling mdp input (and tpr I/O).
        virtual IMdpOptionProvider *mdpOptionProvider() = 0;
        //! Returns an interface for handling output files during simulation.
        virtual IMDOutputProvider *outputProvider()     = 0;
        //! Initializes force providers from this module.
        virtual void initForceProviders(ForceProviders *forceProviders, ForceProviderInitOptions *options) = 0;
        //! finishes force providers from this module, default implementation does nothing
        virtual void finishForceProviders(gmx_unused ForceProviders *forceProviders, gmx_unused ForceProviderInitOptions *options) { }
};

} // namespace gmx

#endif
