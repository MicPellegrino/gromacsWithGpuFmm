/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "pdb2gmx.h"

#include <cctype>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxlib/conformation_utilities.h"
#include "gromacs/gmxpreprocess/fflibutil.h"
#include "gromacs/gmxpreprocess/genhydro.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/h_db.h"
#include "gromacs/gmxpreprocess/hizzie.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/gmxpreprocess/pgutil.h"
#include "gromacs/gmxpreprocess/specbond.h"
#include "gromacs/gmxpreprocess/ter_db.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/gmxpreprocess/xlate.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/dir_separator.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/utility/stringutil.h"

#include "hackblock.h"
#include "resall.h"

struct RtpRename{
    RtpRename(const char *newGmx, const char *newMain,
              const char *newNter, const char *newCter,
              const char *newBter) :
        gmx(newGmx), main(newMain), nter(newNter), cter(newCter), bter(newBter)
    {}
    std::string gmx;
    std::string main;
    std::string nter;
    std::string cter;
    std::string bter;
};

static const char *res2bb_notermini(const std::string             &name,
                                    gmx::ArrayRef<const RtpRename> rr)
{
    /* NOTE: This function returns the main building block name,
     *       it does not take terminal renaming into account.
     */
    auto found = std::find_if(rr.begin(), rr.end(), [&name](const auto &rename)
                              { return gmx::equalCaseInsensitive(name, rename.gmx); });
    return found != rr.end() ? found->main.c_str() : name.c_str();
}

static const char *select_res(int nr, int resnr,
                              const char *name[], const char *expl[],
                              const char *title,
                              gmx::ArrayRef<const RtpRename> rr)
{
    printf("Which %s type do you want for residue %d\n", title, resnr+1);
    for (int sel = 0; (sel < nr); sel++)
    {
        printf("%d. %s (%s)\n",
               sel, expl[sel], res2bb_notermini(name[sel], rr));
    }
    printf("\nType a number:"); fflush(stdout);

    int userSelection;
    if (scanf("%d", &userSelection) != 1)
    {
        gmx_fatal(FARGS, "Answer me for res %s %d!", title, resnr+1);
    }

    return name[userSelection];
}

static const char *get_asptp(int resnr, gmx::ArrayRef<const RtpRename> rr)
{
    enum {
        easp, easpH, easpNR
    };
    const char *lh[easpNR]   = { "ASP", "ASPH" };
    const char *expl[easpNR] = {
        "Not protonated (charge -1)",
        "Protonated (charge 0)"
    };

    return select_res(easpNR, resnr, lh, expl, "ASPARTIC ACID", rr);
}

static const char *get_glutp(int resnr, gmx::ArrayRef<const RtpRename> rr)
{
    enum {
        eglu, egluH, egluNR
    };
    const char *lh[egluNR]   = { "GLU", "GLUH" };
    const char *expl[egluNR] = {
        "Not protonated (charge -1)",
        "Protonated (charge 0)"
    };

    return select_res(egluNR, resnr, lh, expl, "GLUTAMIC ACID", rr);
}

static const char *get_glntp(int resnr, gmx::ArrayRef<const RtpRename> rr)
{
    enum {
        egln, eglnH, eglnNR
    };
    const char *lh[eglnNR]   = { "GLN", "QLN" };
    const char *expl[eglnNR] = {
        "Not protonated (charge 0)",
        "Protonated (charge +1)"
    };

    return select_res(eglnNR, resnr, lh, expl, "GLUTAMINE", rr);
}

static const char *get_lystp(int resnr, gmx::ArrayRef<const RtpRename> rr)
{
    enum {
        elys, elysH, elysNR
    };
    const  char *lh[elysNR]   = { "LYSN", "LYS" };
    const char  *expl[elysNR] = {
        "Not protonated (charge 0)",
        "Protonated (charge +1)"
    };

    return select_res(elysNR, resnr, lh, expl, "LYSINE", rr);
}

static const char *get_argtp(int resnr, gmx::ArrayRef<const RtpRename> rr)
{
    enum {
        earg, eargH, eargNR
    };
    const  char *lh[eargNR]   = { "ARGN", "ARG" };
    const char  *expl[eargNR] = {
        "Not protonated (charge 0)",
        "Protonated (charge +1)"
    };

    return select_res(eargNR, resnr, lh, expl, "ARGININE", rr);
}

static const char *get_histp(int resnr, gmx::ArrayRef<const RtpRename> rr)
{
    const char *expl[ehisNR] = {
        "H on ND1 only",
        "H on NE2 only",
        "H on ND1 and NE2",
        "Coupled to Heme"
    };

    return select_res(ehisNR, resnr, hh, expl, "HISTIDINE", rr);
}

static void read_rtprename(const char *fname, FILE *fp,
                           std::vector<RtpRename> *rtprename)
{
    char         line[STRLEN], buf[STRLEN];

    int          ncol = 0;
    while (get_a_line(fp, line, STRLEN))
    {
        /* line is NULL-terminated and length<STRLEN, so final arg cannot overflow.
         * For other args, we read up to 6 chars (so we can detect if the length is > 5).
         * Note that the buffer length has been increased to 7 to allow this,
         * so we just need to make sure the strings have zero-length initially.
         */
        char gmx[STRLEN];
        char main[STRLEN];
        char nter[STRLEN];
        char cter[STRLEN];
        char bter[STRLEN];
        gmx[0] = '\0'; main[0] = '\0'; nter[0] = '\0'; cter[0] = '\0'; bter[0] = '\0';
        int  nc            = sscanf(line, "%6s %6s %6s %6s %6s %s",
                                    gmx, main, nter,
                                    cter, bter, buf);
        RtpRename newEntry(gmx, main, nter, cter, bter);
        if (ncol == 0)
        {
            if (nc != 2 && nc != 5)
            {
                gmx_fatal(FARGS, "Residue renaming database '%s' has %d columns instead of %d or %d", fname, ncol, 2, 5);
            }
            ncol = nc;
        }
        else if (nc != ncol)
        {
            gmx_fatal(FARGS, "A line in residue renaming database '%s' has %d columns, while previous lines have %d columns", fname, nc, ncol);
        }

        if (nc == 2)
        {
            /* This file does not have special termini names, copy them from main */
            newEntry.nter = newEntry.main;
            newEntry.cter = newEntry.main;
            newEntry.bter = newEntry.main;
        }
        rtprename->push_back(newEntry);
    }
}

static std::string search_resrename(gmx::ArrayRef<const RtpRename> rr,
                                    const char *name,
                                    bool bStart, bool bEnd,
                                    bool bCompareFFRTPname)
{
    auto found = std::find_if(rr.begin(), rr.end(), [&name, &bCompareFFRTPname](const auto &rename)
                              { return ((!bCompareFFRTPname && (name == rename.gmx)) ||
                                        (bCompareFFRTPname && (name == rename.main))); });

    std::string newName;
    /* If found in the database, rename this residue's rtp building block,
     * otherwise keep the old name.
     */
    if (found != rr.end())
    {
        if (bStart && bEnd)
        {
            newName  = found->bter;
        }
        else if (bStart)
        {
            newName = found->nter;
        }
        else if (bEnd)
        {
            newName = found->cter;
        }
        else
        {
            newName = found->main;
        }

        if (newName[0] == '-')
        {
            gmx_fatal(FARGS, "In the chosen force field there is no residue type for '%s'%s", name, bStart ? ( bEnd ? " as a standalone (starting & ending) residue" : " as a starting terminus") : (bEnd ? " as an ending terminus" : ""));
        }
    }

    return newName;
}

static void rename_resrtp(t_atoms                       *pdba,
                          int                            nterpairs,
                          gmx::ArrayRef<const int>       r_start,
                          gmx::ArrayRef<const int>       r_end,
                          gmx::ArrayRef<const RtpRename> rr,
                          t_symtab                      *symtab,
                          bool                           bVerbose)
{
    bool bFFRTPTERRNM = (getenv("GMX_NO_FFRTP_TER_RENAME") == nullptr);

    for (int r = 0; r < pdba->nres; r++)
    {
        bool bStart = false;
        bool bEnd   = false;
        for (int j = 0; j < nterpairs; j++)
        {
            if (r == r_start[j])
            {
                bStart = true;
            }
        }
        for (int j = 0; j < nterpairs; j++)
        {
            if (r == r_end[j])
            {
                bEnd = true;
            }
        }

        std::string newName = search_resrename(rr, *pdba->resinfo[r].rtp, bStart, bEnd, false);

        if (bFFRTPTERRNM && newName.empty() && (bStart || bEnd))
        {
            /* This is a terminal residue, but the residue name,
             * currently stored in .rtp, is not a standard residue name,
             * but probably a force field specific rtp name.
             * Check if we need to rename it because it is terminal.
             */
            newName = search_resrename(rr,
                                       *pdba->resinfo[r].rtp, bStart, bEnd, true);
        }

        if (!newName.empty() && newName != *pdba->resinfo[r].rtp)
        {
            if (bVerbose)
            {
                printf("Changing rtp entry of residue %d %s to '%s'\n",
                       pdba->resinfo[r].nr, *pdba->resinfo[r].name, newName.c_str());
            }
            pdba->resinfo[r].rtp = put_symtab(symtab, newName.c_str());
        }
    }
}

static void pdbres_to_gmxrtp(t_atoms *pdba)
{
    int i;

    for (i = 0; (i < pdba->nres); i++)
    {
        if (pdba->resinfo[i].rtp == nullptr)
        {
            pdba->resinfo[i].rtp = pdba->resinfo[i].name;
        }
    }
}

static void rename_pdbres(t_atoms *pdba, const char *oldnm, const char *newnm,
                          bool bFullCompare, t_symtab *symtab)
{
    char *resnm;
    int   i;

    for (i = 0; (i < pdba->nres); i++)
    {
        resnm = *pdba->resinfo[i].name;
        if ((bFullCompare && (gmx::equalCaseInsensitive(resnm, oldnm))) ||
            (!bFullCompare && strstr(resnm, oldnm) != nullptr))
        {
            /* Rename the residue name (not the rtp name) */
            pdba->resinfo[i].name = put_symtab(symtab, newnm);
        }
    }
}

static void rename_bb(t_atoms *pdba, const char *oldnm, const char *newnm,
                      bool bFullCompare, t_symtab *symtab)
{
    char *bbnm;
    int   i;

    for (i = 0; (i < pdba->nres); i++)
    {
        /* We have not set the rtp name yes, use the residue name */
        bbnm = *pdba->resinfo[i].name;
        if ((bFullCompare && (gmx::equalCaseInsensitive(bbnm, oldnm))) ||
            (!bFullCompare && strstr(bbnm, oldnm) != nullptr))
        {
            /* Change the rtp builing block name */
            pdba->resinfo[i].rtp = put_symtab(symtab, newnm);
        }
    }
}

static void rename_bbint(t_atoms *pdba, const char *oldnm,
                         const char *gettp(int, gmx::ArrayRef<const RtpRename>),
                         bool bFullCompare,
                         t_symtab *symtab,
                         gmx::ArrayRef<const RtpRename> rr)
{
    int         i;
    const char *ptr;
    char       *bbnm;

    for (i = 0; i < pdba->nres; i++)
    {
        /* We have not set the rtp name yet, use the residue name */
        bbnm = *pdba->resinfo[i].name;
        if ((bFullCompare && (strcmp(bbnm, oldnm) == 0)) ||
            (!bFullCompare && strstr(bbnm, oldnm) != nullptr))
        {
            ptr                  = gettp(i, rr);
            pdba->resinfo[i].rtp = put_symtab(symtab, ptr);
        }
    }
}

static void check_occupancy(t_atoms *atoms, const char *filename, bool bVerbose)
{
    int i, ftp;
    int nzero   = 0;
    int nnotone = 0;

    ftp = fn2ftp(filename);
    if (!atoms->pdbinfo || ((ftp != efPDB) && (ftp != efBRK) && (ftp != efENT)))
    {
        fprintf(stderr, "No occupancies in %s\n", filename);
    }
    else
    {
        for (i = 0; (i < atoms->nr); i++)
        {
            if (atoms->pdbinfo[i].occup != 1)
            {
                if (bVerbose)
                {
                    fprintf(stderr, "Occupancy for atom %s%d-%s is %f rather than 1\n",
                            *atoms->resinfo[atoms->atom[i].resind].name,
                            atoms->resinfo[ atoms->atom[i].resind].nr,
                            *atoms->atomname[i],
                            atoms->pdbinfo[i].occup);
                }
                if (atoms->pdbinfo[i].occup == 0)
                {
                    nzero++;
                }
                else
                {
                    nnotone++;
                }
            }
        }
        if (nzero == atoms->nr)
        {
            fprintf(stderr, "All occupancy fields zero. This is probably not an X-Ray structure\n");
        }
        else if ((nzero > 0) || (nnotone > 0))
        {
            fprintf(stderr,
                    "\n"
                    "WARNING: there were %d atoms with zero occupancy and %d atoms with\n"
                    "         occupancy unequal to one (out of %d atoms). Check your pdb file.\n"
                    "\n",
                    nzero, nnotone, atoms->nr);
        }
        else
        {
            fprintf(stderr, "All occupancies are one\n");
        }
    }
}

static void write_posres(const char *fn, t_atoms *pdba, real fc)
{
    FILE *fp;
    int   i;

    fp = gmx_fio_fopen(fn, "w");
    fprintf(fp,
            "; In this topology include file, you will find position restraint\n"
            "; entries for all the heavy atoms in your original pdb file.\n"
            "; This means that all the protons which were added by pdb2gmx are\n"
            "; not restrained.\n"
            "\n"
            "[ position_restraints ]\n"
            "; %4s%6s%8s%8s%8s\n", "atom", "type", "fx", "fy", "fz"
            );
    for (i = 0; (i < pdba->nr); i++)
    {
        if (!is_hydrogen(*pdba->atomname[i]) && !is_dummymass(*pdba->atomname[i]))
        {
            fprintf(fp, "%6d%6d  %g  %g  %g\n", i+1, 1, fc, fc, fc);
        }
    }
    gmx_fio_fclose(fp);
}

static int read_pdball(const char *inf, bool bOutput, const char *outf, char **title,
                       t_atoms *atoms, rvec **x,
                       int *ePBC, matrix box, bool bRemoveH,
                       t_symtab *symtab, ResidueType *rt, const char *watres,
                       AtomProperties *aps, bool bVerbose)
/* Read a pdb file. (containing proteins) */
{
    int natom, new_natom, i;

    /* READ IT */
    printf("Reading %s...\n", inf);
    readConfAndAtoms(inf, symtab, title, atoms, ePBC, x, nullptr, box);
    natom           = atoms->nr;
    if (atoms->pdbinfo == nullptr)
    {
        snew(atoms->pdbinfo, atoms->nr);
    }
    if (fn2ftp(inf) == efPDB)
    {
        get_pdb_atomnumber(atoms, aps);
    }
    if (bRemoveH)
    {
        new_natom = 0;
        for (i = 0; i < atoms->nr; i++)
        {
            if (!is_hydrogen(*atoms->atomname[i]))
            {
                atoms->atom[new_natom]     = atoms->atom[i];
                atoms->atomname[new_natom] = atoms->atomname[i];
                atoms->pdbinfo[new_natom]  = atoms->pdbinfo[i];
                copy_rvec((*x)[i], (*x)[new_natom]);
                new_natom++;
            }
        }
        atoms->nr = new_natom;
        natom     = new_natom;
    }

    printf("Read");
    if (title[0])
    {
        printf(" '%s',", *title);
    }
    printf(" %d atoms\n", natom);

    /* Rename residues */
    rename_pdbres(atoms, "HOH", watres, false, symtab);
    rename_pdbres(atoms, "SOL", watres, false, symtab);
    rename_pdbres(atoms, "WAT", watres, false, symtab);

    rename_atoms("xlateat.dat", nullptr,
                 atoms, symtab, {}, true,
                 rt, true, bVerbose);

    if (natom == 0)
    {
        return 0;
    }
    if (bOutput)
    {
        write_sto_conf(outf, *title, atoms, *x, nullptr, *ePBC, box);
    }

    return natom;
}

static void process_chain(t_atoms *pdba, gmx::ArrayRef<gmx::RVec> x,
                          bool bTrpU, bool bPheU, bool bTyrU,
                          bool bLysMan, bool bAspMan, bool bGluMan,
                          bool bHisMan, bool bArgMan, bool bGlnMan,
                          real angle, real distance, t_symtab *symtab,
                          gmx::ArrayRef<const RtpRename> rr)
{
    /* Rename aromatics, lys, asp and histidine */
    if (bTyrU)
    {
        rename_bb(pdba, "TYR", "TYRU", false, symtab);
    }
    if (bTrpU)
    {
        rename_bb(pdba, "TRP", "TRPU", false, symtab);
    }
    if (bPheU)
    {
        rename_bb(pdba, "PHE", "PHEU", false, symtab);
    }
    if (bLysMan)
    {
        rename_bbint(pdba, "LYS", get_lystp, false, symtab, rr);
    }
    if (bArgMan)
    {
        rename_bbint(pdba, "ARG", get_argtp, false, symtab, rr);
    }
    if (bGlnMan)
    {
        rename_bbint(pdba, "GLN", get_glntp, false, symtab, rr);
    }
    if (bAspMan)
    {
        rename_bbint(pdba, "ASP", get_asptp, false, symtab, rr);
    }
    else
    {
        rename_bb(pdba, "ASPH", "ASP", false, symtab);
    }
    if (bGluMan)
    {
        rename_bbint(pdba, "GLU", get_glutp, false, symtab, rr);
    }
    else
    {
        rename_bb(pdba, "GLUH", "GLU", false, symtab);
    }

    if (!bHisMan)
    {
        set_histp(pdba, gmx::as_rvec_array(x.data()), symtab, angle, distance);
    }
    else
    {
        rename_bbint(pdba, "HIS", get_histp, true, symtab, rr);
    }

    /* Initialize the rtp builing block names with the residue names
     * for the residues that have not been processed above.
     */
    pdbres_to_gmxrtp(pdba);

    /* Now we have all rtp names set.
     * The rtp names will conform to Gromacs naming,
     * unless the input pdb file contained one or more force field specific
     * rtp names as residue names.
     */
}

/* struct for sorting the atoms from the pdb file */
typedef struct {
    int  resnr;  /* residue number               */
    int  j;      /* database order index         */
    int  index;  /* original atom number         */
    char anm1;   /* second letter of atom name   */
    char altloc; /* alternate location indicator */
} t_pdbindex;

static bool pdbicomp(const t_pdbindex &a, const t_pdbindex &b)
{
    int d = (a.resnr - b.resnr);
    if (d == 0)
    {
        d = (a.j - b.j);
        if (d == 0)
        {
            d = (a.anm1 - b.anm1);
            if (d == 0)
            {
                d = (a.altloc - b.altloc);
            }
        }
    }
    return d < 0;
}

static void sort_pdbatoms(gmx::ArrayRef<const PreprocessResidue> restp_chain,
                          int natoms,
                          t_atoms **pdbaptr,
                          t_atoms **newPdbAtoms,
                          std::vector<gmx::RVec> *x,
                          t_blocka *block, char ***gnames)
{
    t_atoms               *pdba = *pdbaptr;
    std::vector<gmx::RVec> xnew;
    t_pdbindex            *pdbi;
    char                  *atomnm;

    natoms = pdba->nr;
    snew(pdbi, natoms);

    for (int i = 0; i < natoms; i++)
    {
        atomnm = *pdba->atomname[i];
        const PreprocessResidue *localPpResidue   = &restp_chain[pdba->atom[i].resind];
        auto                     found            = std::find_if(localPpResidue->atomname.begin(), localPpResidue->atomname.end(),
                                                                 [&atomnm](char** it){return gmx::equalCaseInsensitive(atomnm, *it); });
        if (found == localPpResidue->atomname.end())
        {
            char buf[STRLEN];

            sprintf(buf,
                    "Atom %s in residue %s %d was not found in rtp entry %s with %d atoms\n"
                    "while sorting atoms.\n%s", atomnm,
                    *pdba->resinfo[pdba->atom[i].resind].name,
                    pdba->resinfo[pdba->atom[i].resind].nr,
                    localPpResidue->resname.c_str(),
                    localPpResidue->natom(),
                    is_hydrogen(atomnm) ?
                    "\nFor a hydrogen, this can be a different protonation state, or it\n"
                    "might have had a different number in the PDB file and was rebuilt\n"
                    "(it might for instance have been H3, and we only expected H1 & H2).\n"
                    "Note that hydrogens might have been added to the entry for the N-terminus.\n"
                    "Remove this hydrogen or choose a different protonation state to solve it.\n"
                    "Option -ignh will ignore all hydrogens in the input." : ".");
            gmx_fatal(FARGS, "%s", buf);
        }
        /* make shadow array to be sorted into indexgroup */
        pdbi[i].resnr  = pdba->atom[i].resind;
        pdbi[i].j      = std::distance(localPpResidue->atomname.begin(), found);
        pdbi[i].index  = i;
        pdbi[i].anm1   = atomnm[1];
        pdbi[i].altloc = pdba->pdbinfo[i].altloc;
    }
    std::sort(pdbi, pdbi+natoms, pdbicomp);

    /* pdba is sorted in pdbnew using the pdbi index */
    std::vector<int> a(natoms);
    srenew(*newPdbAtoms, 1);
    init_t_atoms((*newPdbAtoms), natoms, true);
    (*newPdbAtoms)->nr   = pdba->nr;
    (*newPdbAtoms)->nres = pdba->nres;
    srenew((*newPdbAtoms)->resinfo, pdba->nres);
    std::copy(pdba->resinfo, pdba->resinfo + pdba->nres, (*newPdbAtoms)->resinfo);
    for (int i = 0; i < natoms; i++)
    {
        (*newPdbAtoms)->atom[i]     = pdba->atom[pdbi[i].index];
        (*newPdbAtoms)->atomname[i] = pdba->atomname[pdbi[i].index];
        (*newPdbAtoms)->pdbinfo[i]  = pdba->pdbinfo[pdbi[i].index];
        xnew.push_back(x->at(pdbi[i].index));
        /* make indexgroup in block */
        a[i] = pdbi[i].index;
    }
    /* clean up */
    done_atom(pdba);
    sfree(pdba);
    /* copy the sorted pdbnew back to pdba */
    *pdbaptr = *newPdbAtoms;
    *x       = xnew;
    add_grp(block, gnames, a, "prot_sort");
    sfree(pdbi);
}

static int remove_duplicate_atoms(t_atoms *pdba, gmx::ArrayRef<gmx::RVec> x, bool bVerbose)
{
    int        i, j, oldnatoms, ndel;
    t_resinfo *ri;

    printf("Checking for duplicate atoms....\n");
    oldnatoms    = pdba->nr;
    ndel         = 0;
    /* NOTE: pdba->nr is modified inside the loop */
    for (i = 1; (i < pdba->nr); i++)
    {
        /* compare 'i' and 'i-1', throw away 'i' if they are identical
           this is a 'while' because multiple alternate locations can be present */
        while ( (i < pdba->nr) &&
                (pdba->atom[i-1].resind == pdba->atom[i].resind) &&
                (strcmp(*pdba->atomname[i-1], *pdba->atomname[i]) == 0) )
        {
            ndel++;
            if (bVerbose)
            {
                ri = &pdba->resinfo[pdba->atom[i].resind];
                printf("deleting duplicate atom %4s  %s%4d%c",
                       *pdba->atomname[i], *ri->name, ri->nr, ri->ic);
                if (ri->chainid && (ri->chainid != ' '))
                {
                    printf(" ch %c", ri->chainid);
                }
                if (pdba->pdbinfo)
                {
                    if (pdba->pdbinfo[i].atomnr)
                    {
                        printf("  pdb nr %4d", pdba->pdbinfo[i].atomnr);
                    }
                    if (pdba->pdbinfo[i].altloc && (pdba->pdbinfo[i].altloc != ' '))
                    {
                        printf("  altloc %c", pdba->pdbinfo[i].altloc);
                    }
                }
                printf("\n");
            }
            pdba->nr--;
            /* We can not free, since it might be in the symtab */
            /* sfree(pdba->atomname[i]); */
            for (j = i; j < pdba->nr; j++)
            {
                pdba->atom[j]     = pdba->atom[j+1];
                pdba->atomname[j] = pdba->atomname[j+1];
                if (pdba->pdbinfo)
                {
                    pdba->pdbinfo[j]  = pdba->pdbinfo[j+1];
                }
                copy_rvec(x[j+1], x[j]);
            }
            srenew(pdba->atom,     pdba->nr);
            /* srenew(pdba->atomname, pdba->nr); */
            srenew(pdba->pdbinfo,  pdba->nr);
        }
    }
    if (pdba->nr != oldnatoms)
    {
        printf("Now there are %d atoms. Deleted %d duplicates.\n", pdba->nr, ndel);
    }

    return pdba->nr;
}

static void
checkResidueTypeSanity(t_atoms     *pdba,
                       int          r0,
                       int          r1,
                       ResidueType *rt)
{
    std::string startResidueString = gmx::formatString("%s%d", *pdba->resinfo[r0].name, pdba->resinfo[r0].nr);
    std::string endResidueString   = gmx::formatString("%s%d", *pdba->resinfo[r1-1].name, pdba->resinfo[r1-1].nr);

    // Check whether all residues in chain have the same chain ID.
    bool         allResiduesHaveSameChainID = true;
    char         chainID0                   = pdba->resinfo[r0].chainid;
    char         chainID;
    std::string  residueString;

    for (int i = r0 + 1; i < r1; i++)
    {
        chainID = pdba->resinfo[i].chainid;
        if (chainID != chainID0)
        {
            allResiduesHaveSameChainID  = false;
            residueString               = gmx::formatString("%s%d", *pdba->resinfo[i].name, pdba->resinfo[i].nr);
            break;
        }
    }

    if (!allResiduesHaveSameChainID)
    {
        gmx_fatal(FARGS,
                  "The chain covering the range %s--%s does not have a consistent chain ID. "
                  "The first residue has ID '%c', while residue %s has ID '%c'.",
                  startResidueString.c_str(), endResidueString.c_str(),
                  chainID0, residueString.c_str(), chainID);
    }

    // At this point all residues have the same ID. If they are also non-blank
    // we can be a bit more aggressive and require the types match too.
    if (chainID0 != ' ')
    {
        bool        allResiduesHaveSameType = true;
        std::string restype;
        std::string restype0 = rt->typeOfNamedDatabaseResidue(*pdba->resinfo[r0].name);

        for (int i = r0 + 1; i < r1; i++)
        {
            restype = rt->typeOfNamedDatabaseResidue(*pdba->resinfo[i].name);
            if (!gmx::equalCaseInsensitive(restype, restype0))
            {
                allResiduesHaveSameType = false;
                residueString           = gmx::formatString("%s%d", *pdba->resinfo[i].name, pdba->resinfo[i].nr);
                break;
            }
        }

        if (!allResiduesHaveSameType)
        {
            gmx_fatal(FARGS,
                      "The residues in the chain %s--%s do not have a consistent type. "
                      "The first residue has type '%s', while residue %s is of type '%s'. "
                      "Either there is a mistake in your chain, or it includes nonstandard "
                      "residue names that have not yet been added to the residuetypes.dat "
                      "file in the GROMACS library directory. If there are other molecules "
                      "such as ligands, they should not have the same chain ID as the "
                      "adjacent protein chain since it's a separate molecule.",
                      startResidueString.c_str(), endResidueString.c_str(),
                      restype0.c_str(), residueString.c_str(), restype.c_str());
        }
    }
}

static void find_nc_ter(t_atoms *pdba, int r0, int r1, int *r_start, int *r_end,
                        ResidueType *rt)
{
    int         i;
    gmx::compat::optional<std::string> startrestype;

    *r_start = -1;
    *r_end   = -1;

    int startWarnings = 0;
    int endWarnings   = 0;
    int ionNotes      = 0;

    // Check that all residues have the same chain identifier, and if it is
    // non-blank we also require the residue types to match.
    checkResidueTypeSanity(pdba, r0, r1, rt);

    // If we return correctly from checkResidueTypeSanity(), the only
    // remaining cases where we can have non-matching residue types is if
    // the chain ID was blank, which could be the case e.g. for a structure
    // read from a GRO file or other file types without chain information.
    // In that case we need to be a bit more liberal and detect chains based
    // on the residue type.

    // If we get here, the chain ID must be identical for all residues
    char chainID = pdba->resinfo[r0].chainid;

    /* Find the starting terminus (typially N or 5') */
    for (i = r0; i < r1 && *r_start == -1; i++)
    {
        startrestype = rt->optionalTypeOfNamedDatabaseResidue(*pdba->resinfo[i].name);
        if (!startrestype)
        {
            continue;
        }
        if (gmx::equalCaseInsensitive(*startrestype, "Protein") ||
            gmx::equalCaseInsensitive(*startrestype, "DNA") ||
            gmx::equalCaseInsensitive(*startrestype, "RNA") )
        {
            printf("Identified residue %s%d as a starting terminus.\n", *pdba->resinfo[i].name, pdba->resinfo[i].nr);
            *r_start = i;
        }
        else if (gmx::equalCaseInsensitive(*startrestype, "Ion"))
        {
            if (ionNotes < 5)
            {
                printf("Residue %s%d has type 'Ion', assuming it is not linked into a chain.\n", *pdba->resinfo[i].name, pdba->resinfo[i].nr);
            }
            if (ionNotes == 4)
            {
                printf("Disabling further notes about ions.\n");
            }
            ionNotes++;
        }
        else
        {
            // Either no known residue type, or one not needing special handling
            if (startWarnings < 5)
            {
                if (chainID == ' ')
                {
                    printf("\nWarning: Starting residue %s%d in chain not identified as Protein/RNA/DNA.\n"
                           "This chain lacks identifiers, which makes it impossible to do strict\n"
                           "classification of the start/end residues. Here we need to guess this residue\n"
                           "should not be part of the chain and instead introduce a break, but that will\n"
                           "be catastrophic if they should in fact be linked. Please check your structure,\n"
                           "and add %s to residuetypes.dat if this was not correct.\n\n",
                           *pdba->resinfo[i].name, pdba->resinfo[i].nr, *pdba->resinfo[i].name);
                }
                else
                {
                    printf("\nWarning: No residues in chain starting at %s%d identified as Protein/RNA/DNA.\n"
                           "This makes it impossible to link them into a molecule, which could either be\n"
                           "correct or a catastrophic error. Please check your structure, and add all\n"
                           "necessary residue names to residuetypes.dat if this was not correct.\n\n",
                           *pdba->resinfo[i].name, pdba->resinfo[i].nr);
                }
            }
            if (startWarnings == 4)
            {
                printf("Disabling further warnings about unidentified residues at start of chain.\n");
            }
            startWarnings++;
        }
    }

    if (*r_start >= 0)
    {
        /* Go through the rest of the residues, check that they are the same class, and identify the ending terminus. */
        for (int i = *r_start; i < r1; i++)
        {
            gmx::compat::optional<std::string> restype =
                rt->optionalTypeOfNamedDatabaseResidue(*pdba->resinfo[i].name);
            if (!restype)
            {
                continue;
            }
            if (gmx::equalCaseInsensitive(*restype, *startrestype) && endWarnings == 0)
            {
                *r_end = i;
            }
            else if (gmx::equalCaseInsensitive(*startrestype, "Ion"))
            {
                if (ionNotes < 5)
                {
                    printf("Residue %s%d has type 'Ion', assuming it is not linked into a chain.\n", *pdba->resinfo[i].name, pdba->resinfo[i].nr);
                }
                if (ionNotes == 4)
                {
                    printf("Disabling further notes about ions.\n");
                }
                ionNotes++;
            }
            else
            {
                // Either no known residue type, or one not needing special handling.
                GMX_RELEASE_ASSERT(chainID == ' ', "Chain ID must be blank");
                // Otherwise the call to checkResidueTypeSanity() will
                // have caught the problem.
                if (endWarnings < 5)
                {
                    printf("\nWarning: Residue %s%d in chain has different type ('%s') from\n"
                           "residue %s%d ('%s'). This chain lacks identifiers, which makes\n"
                           "it impossible to do strict classification of the start/end residues. Here we\n"
                           "need to guess this residue should not be part of the chain and instead\n"
                           "introduce a break, but that will be catastrophic if they should in fact be\n"
                           "linked. Please check your structure, and add %s to residuetypes.dat\n"
                           "if this was not correct.\n\n",
                           *pdba->resinfo[i].name, pdba->resinfo[i].nr, restype->c_str(),
                           *pdba->resinfo[*r_start].name, pdba->resinfo[*r_start].nr, startrestype->c_str(), *pdba->resinfo[i].name);
                }
                if (endWarnings == 4)
                {
                    printf("Disabling further warnings about unidentified residues at end of chain.\n");
                }
                endWarnings++;
            }
        }
    }

    if (*r_end >= 0)
    {
        printf("Identified residue %s%d as a ending terminus.\n", *pdba->resinfo[*r_end].name, pdba->resinfo[*r_end].nr);
    }
}

/* enum for chain separation */
enum ChainSepType {
    enChainSep_id_or_ter, enChainSep_id_and_ter, enChainSep_ter,
    enChainSep_id, enChainSep_interactive
};
static const char *ChainSepEnum[] = {"id_or_ter", "id_and_ter", "ter", "id", "interactive"};
static const char *ChainSepInfoString[] =
{
    "Splitting chemical chains based on TER records or chain id changing.\n",
    "Splitting chemical chains based on TER records and chain id changing.\n",
    "Splitting chemical chains based on TER records only (ignoring chain id).\n",
    "Splitting chemical chains based on changing chain id only (ignoring TER records).\n",
    "Splitting chemical chains interactively.\n"
};

static void
modify_chain_numbers(t_atoms *       pdba,
                     ChainSepType    enumChainSep)
{
    int           i;
    char          old_prev_chainid;
    char          old_this_chainid;
    int           old_prev_chainnum;
    int           old_this_chainnum;
    t_resinfo    *ri;
    char          select[STRLEN];
    int           new_chainnum;
    int           this_atomnum;
    int           prev_atomnum;
    const char *  prev_atomname;
    const char *  this_atomname;
    const char *  prev_resname;
    const char *  this_resname;
    int           prev_resnum;
    int           this_resnum;
    char          prev_chainid;
    char          this_chainid;

    /* The default chain enumeration is based on TER records only */
    printf("%s", ChainSepInfoString[enumChainSep]);

    old_prev_chainid  = '?';
    old_prev_chainnum = -1;
    new_chainnum      = -1;

    this_atomname       = nullptr;
    this_atomnum        = -1;
    this_resname        = nullptr;
    this_resnum         = -1;
    this_chainid        = '?';

    for (i = 0; i < pdba->nres; i++)
    {
        ri                 = &pdba->resinfo[i];
        old_this_chainid   = ri->chainid;
        old_this_chainnum  = ri->chainnum;

        prev_atomname      = this_atomname;
        prev_atomnum       = this_atomnum;
        prev_resname       = this_resname;
        prev_resnum        = this_resnum;
        prev_chainid       = this_chainid;

        this_atomname      = *(pdba->atomname[i]);
        this_atomnum       = (pdba->pdbinfo != nullptr) ? pdba->pdbinfo[i].atomnr : i+1;
        this_resname       = *ri->name;
        this_resnum        = ri->nr;
        this_chainid       = ri->chainid;

        switch (enumChainSep)
        {
            case enChainSep_id_or_ter:
                if (old_this_chainid != old_prev_chainid || old_this_chainnum != old_prev_chainnum)
                {
                    new_chainnum++;
                }
                break;

            case enChainSep_id_and_ter:
                if (old_this_chainid != old_prev_chainid && old_this_chainnum != old_prev_chainnum)
                {
                    new_chainnum++;
                }
                break;

            case enChainSep_id:
                if (old_this_chainid != old_prev_chainid)
                {
                    new_chainnum++;
                }
                break;

            case enChainSep_ter:
                if (old_this_chainnum != old_prev_chainnum)
                {
                    new_chainnum++;
                }
                break;
            case enChainSep_interactive:
                if (old_this_chainid != old_prev_chainid || old_this_chainnum != old_prev_chainnum)
                {
                    if (i > 0)
                    {
                        printf("Split the chain (and introduce termini) between residue %s%d (chain id '%c', atom %d %s)\
\n"
                               "and residue %s%d (chain id '%c', atom %d %s) ? [n/y]\n",
                               prev_resname, prev_resnum, prev_chainid, prev_atomnum, prev_atomname,
                               this_resname, this_resnum, this_chainid, this_atomnum, this_atomname);

                        if (nullptr == fgets(select, STRLEN-1, stdin))
                        {
                            gmx_fatal(FARGS, "Error reading from stdin");
                        }
                    }
                    if (i == 0 || select[0] == 'y')
                    {
                        new_chainnum++;
                    }
                }
                break;
            default:
                gmx_fatal(FARGS, "Internal inconsistency - this shouldn't happen...");
        }
        old_prev_chainid  = old_this_chainid;
        old_prev_chainnum = old_this_chainnum;

        ri->chainnum = new_chainnum;
    }
}

struct t_pdbchain {
    char             chainid   = ' ';
    char             chainnum  = ' ';
    int              start     = -1;
    int              natom     = -1;
    bool             bAllWat   = false;
    int              nterpairs = -1;
    std::vector<int> chainstart;
};

struct t_chain {
    char                                 chainid   = ' ';
    int                                  chainnum  = ' ';
    bool                                 bAllWat   = false;
    int                                  nterpairs = -1;
    std::vector<int>                     chainstart;
    std::vector<MoleculePatchDatabase *> ntdb;
    std::vector<MoleculePatchDatabase *> ctdb;
    std::vector<int>                     r_start;
    std::vector<int>                     r_end;
    t_atoms                             *pdba;
    std::vector<gmx::RVec>               x;
};

// TODO make all enums into scoped enums
/* enum for vsites */
enum VSitesType {
    enVSites_none, enVSites_hydrogens, enVSites_aromatics
};
static const char *VSitesEnum[] = {"none", "hydrogens", "aromatics"};

/* enum for water model */
enum WaterType {
    enWater_select, enWater_none, enWater_spc, enWater_spce,
    enWater_tip3p, enWater_tip4p, enWater_tip5p, enWater_tips3p
};
static const char *WaterEnum[] = {
    "select", "none", "spc", "spce",
    "tip3p", "tip4p", "tip5p", "tips3p"
};

/* enum for merge */
enum MergeType {
    enMerge_no, enMerge_all, enMerge_interactive
};
static const char *MergeEnum[] = {"no", "all", "interactive"};

namespace gmx
{

namespace
{

class pdb2gmx : public ICommandLineOptionsModule
{
    public:
        pdb2gmx() :
            bVsites_(FALSE), bPrevWat_(FALSE), bVsiteAromatics_(FALSE),
            enumChainSep_(enChainSep_id_or_ter),
            enumVSites_(enVSites_none),
            enumWater_(enWater_select),
            enumMerge_(enMerge_no),
            itp_file_(nullptr),
            mHmult_(0)
        {
        }

        // From ICommandLineOptionsModule
        void init(CommandLineModuleSettings * /*settings*/) override
        {
        }

        void initOptions(IOptionsContainer                 *options,
                         ICommandLineOptionsModuleSettings *settings) override;

        void optionsFinished() override;

        int run() override;

    private:
        bool                     bNewRTP_;
        bool                     bInter_;
        bool                     bCysMan_;
        bool                     bLysMan_;
        bool                     bAspMan_;
        bool                     bGluMan_;
        bool                     bHisMan_;
        bool                     bGlnMan_;
        bool                     bArgMan_;
        bool                     bTerMan_;
        bool                     bUnA_;
        bool                     bHeavyH_;
        bool                     bSort_;
        bool                     bAllowMissing_;
        bool                     bRemoveH_;
        bool                     bDeuterate_;
        bool                     bVerbose_;
        bool                     bChargeGroups_;
        bool                     bCmap_;
        bool                     bRenumRes_;
        bool                     bRTPresname_;
        bool                     bIndexSet_;
        bool                     bOutputSet_;
        bool                     bVsites_;
        bool                     bWat_;
        bool                     bPrevWat_;
        bool                     bITP_;
        bool                     bVsiteAromatics_;
        real                     angle_;
        real                     distance_;
        real                     posre_fc_;
        real                     long_bond_dist_;
        real                     short_bond_dist_;

        std::string              indexOutputFile_;
        std::string              outputFile_;
        std::string              topologyFile_;
        std::string              includeTopologyFile_;
        std::string              outputConfFile_;
        std::string              inputConfFile_;
        std::string              outFile_;
        std::string              ff_;

        ChainSepType             enumChainSep_;
        VSitesType               enumVSites_;
        WaterType                enumWater_;
        MergeType                enumMerge_;

        FILE                    *itp_file_;
        char                     forcefield_[STRLEN];
        char                     ffdir_[STRLEN];
        char                    *ffname_;
        char                    *watermodel_;
        std::vector<std::string> incls_;
        std::vector<t_mols>      mols_;
        real                     mHmult_;
};

void pdb2gmx::initOptions(IOptionsContainer                 *options,
                          ICommandLineOptionsModuleSettings *settings)
{
    const char *desc[] = {
        "[THISMODULE] reads a [REF].pdb[ref] (or [REF].gro[ref]) file, reads",
        "some database files, adds hydrogens to the molecules and generates",
        "coordinates in GROMACS (GROMOS), or optionally [REF].pdb[ref], format",
        "and a topology in GROMACS format.",
        "These files can subsequently be processed to generate a run input file.",
        "[PAR]",
        "[THISMODULE] will search for force fields by looking for",
        "a [TT]forcefield.itp[tt] file in subdirectories [TT]<forcefield>.ff[tt]",
        "of the current working directory and of the GROMACS library directory",
        "as inferred from the path of the binary or the [TT]GMXLIB[tt] environment",
        "variable.",
        "By default the forcefield selection is interactive,",
        "but you can use the [TT]-ff[tt] option to specify one of the short names",
        "in the list on the command line instead. In that case [THISMODULE] just looks",
        "for the corresponding [TT]<forcefield>.ff[tt] directory.",
        "[PAR]",
        "After choosing a force field, all files will be read only from",
        "the corresponding force field directory.",
        "If you want to modify or add a residue types, you can copy the force",
        "field directory from the GROMACS library directory to your current",
        "working directory. If you want to add new protein residue types,",
        "you will need to modify [TT]residuetypes.dat[tt] in the library directory",
        "or copy the whole library directory to a local directory and set",
        "the environment variable [TT]GMXLIB[tt] to the name of that directory.",
        "Check Chapter 5 of the manual for more information about file formats.",
        "[PAR]",

        "Note that a [REF].pdb[ref] file is nothing more than a file format, and it",
        "need not necessarily contain a protein structure. Every kind of",
        "molecule for which there is support in the database can be converted.",
        "If there is no support in the database, you can add it yourself.[PAR]",

        "The program has limited intelligence, it reads a number of database",
        "files, that allow it to make special bonds (Cys-Cys, Heme-His, etc.),",
        "if necessary this can be done manually. The program can prompt the",
        "user to select which kind of LYS, ASP, GLU, CYS or HIS residue is",
        "desired. For Lys the choice is between neutral (two protons on NZ) or",
        "protonated (three protons, default), for Asp and Glu unprotonated",
        "(default) or protonated, for His the proton can be either on ND1,",
        "on NE2 or on both. By default these selections are done automatically.",
        "For His, this is based on an optimal hydrogen bonding",
        "conformation. Hydrogen bonds are defined based on a simple geometric",
        "criterion, specified by the maximum hydrogen-donor-acceptor angle",
        "and donor-acceptor distance, which are set by [TT]-angle[tt] and",
        "[TT]-dist[tt] respectively.[PAR]",

        "The protonation state of N- and C-termini can be chosen interactively",
        "with the [TT]-ter[tt] flag.  Default termini are ionized (NH3+ and COO-),",
        "respectively.  Some force fields support zwitterionic forms for chains of",
        "one residue, but for polypeptides these options should NOT be selected.",
        "The AMBER force fields have unique forms for the terminal residues,",
        "and these are incompatible with the [TT]-ter[tt] mechanism. You need",
        "to prefix your N- or C-terminal residue names with \"N\" or \"C\"",
        "respectively to use these forms, making sure you preserve the format",
        "of the coordinate file. Alternatively, use named terminating residues",
        "(e.g. ACE, NME).[PAR]",

        "The separation of chains is not entirely trivial since the markup",
        "in user-generated PDB files frequently varies and sometimes it",
        "is desirable to merge entries across a TER record, for instance",
        "if you want a disulfide bridge or distance restraints between",
        "two protein chains or if you have a HEME group bound to a protein.",
        "In such cases multiple chains should be contained in a single",
        "[TT]moleculetype[tt] definition.",
        "To handle this, [THISMODULE] uses two separate options.",
        "First, [TT]-chainsep[tt] allows you to choose when a new chemical chain should",
        "start, and termini added when applicable. This can be done based on the",
        "existence of TER records, when the chain id changes, or combinations of either",
        "or both of these. You can also do the selection fully interactively.",
        "In addition, there is a [TT]-merge[tt] option that controls how multiple chains",
        "are merged into one moleculetype, after adding all the chemical termini (or not).",
        "This can be turned off (no merging), all non-water chains can be merged into a",
        "single molecule, or the selection can be done interactively.[PAR]",

        "[THISMODULE] will also check the occupancy field of the [REF].pdb[ref] file.",
        "If any of the occupancies are not one, indicating that the atom is",
        "not resolved well in the structure, a warning message is issued.",
        "When a [REF].pdb[ref] file does not originate from an X-ray structure determination",
        "all occupancy fields may be zero. Either way, it is up to the user",
        "to verify the correctness of the input data (read the article!).[PAR]",

        "During processing the atoms will be reordered according to GROMACS",
        "conventions. With [TT]-n[tt] an index file can be generated that",
        "contains one group reordered in the same way. This allows you to",
        "convert a GROMOS trajectory and coordinate file to GROMOS. There is",
        "one limitation: reordering is done after the hydrogens are stripped",
        "from the input and before new hydrogens are added. This means that",
        "you should not use [TT]-ignh[tt].[PAR]",

        "The [REF].gro[ref] and [TT].g96[tt] file formats do not support chain",
        "identifiers. Therefore it is useful to enter a [REF].pdb[ref] file name at",
        "the [TT]-o[tt] option when you want to convert a multi-chain [REF].pdb[ref] file.",
        "[PAR]",

        "The option [TT]-vsite[tt] removes hydrogen and fast improper dihedral",
        "motions. Angular and out-of-plane motions can be removed by changing",
        "hydrogens into virtual sites and fixing angles, which fixes their",
        "position relative to neighboring atoms. Additionally, all atoms in the",
        "aromatic rings of the standard amino acids (i.e. PHE, TRP, TYR and HIS)",
        "can be converted into virtual sites, eliminating the fast improper dihedral",
        "fluctuations in these rings (but this feature is deprecated).",
        "[BB]Note[bb] that in this case all other hydrogen",
        "atoms are also converted to virtual sites. The mass of all atoms that are",
        "converted into virtual sites, is added to the heavy atoms.[PAR]",
        "Also slowing down of dihedral motion can be done with [TT]-heavyh[tt]",
        "done by increasing the hydrogen-mass by a factor of 4. This is also",
        "done for water hydrogens to slow down the rotational motion of water.",
        "The increase in mass of the hydrogens is subtracted from the bonded",
        "(heavy) atom so that the total mass of the system remains the same."
    };

    settings->setHelpText(desc);

    options->addOption(BooleanOption("newrtp")
                           .store(&bNewRTP_).defaultValue(false).hidden()
                           .description("Write the residue database in new format to [TT]new.rtp[tt]"));
    options->addOption(RealOption("lb")
                           .store(&long_bond_dist_).defaultValue(0.25).hidden()
                           .description("Long bond warning distance"));
    options->addOption(RealOption("sb")
                           .store(&short_bond_dist_).defaultValue(0.05).hidden()
                           .description("Short bond warning distance"));
    options->addOption(EnumOption<ChainSepType>("chainsep").enumValue(ChainSepEnum)
                           .store(&enumChainSep_)
                           .description("Condition in PDB files when a new chain should be started (adding termini)"));
    options->addOption(EnumOption<MergeType>("merge").enumValue(MergeEnum)
                           .store(&enumMerge_)
                           .description("Merge multiple chains into a single [moleculetype]"));
    options->addOption(StringOption("ff")
                           .store(&ff_).defaultValue("select")
                           .description("Force field, interactive by default. Use [TT]-h[tt] for information."));
    options->addOption(EnumOption<WaterType>("water")
                           .store(&enumWater_).enumValue(WaterEnum)
                           .description("Water model to use"));
    options->addOption(BooleanOption("inter")
                           .store(&bInter_).defaultValue(false)
                           .description("Set the next 8 options to interactive"));
    options->addOption(BooleanOption("ss")
                           .store(&bCysMan_).defaultValue(false)
                           .description("Interactive SS bridge selection"));
    options->addOption(BooleanOption("ter")
                           .store(&bTerMan_).defaultValue(false)
                           .description("Interactive termini selection, instead of charged (default)"));
    options->addOption(BooleanOption("lys")
                           .store(&bLysMan_).defaultValue(false)
                           .description("Interactive lysine selection, instead of charged"));
    options->addOption(BooleanOption("arg")
                           .store(&bArgMan_).defaultValue(false)
                           .description("Interactive arginine selection, instead of charged"));
    options->addOption(BooleanOption("asp")
                           .store(&bAspMan_).defaultValue(false)
                           .description("Interactive aspartic acid selection, instead of charged"));
    options->addOption(BooleanOption("glu")
                           .store(&bGluMan_).defaultValue(false)
                           .description("Interactive glutamic acid selection, instead of charged"));
    options->addOption(BooleanOption("gln")
                           .store(&bGlnMan_).defaultValue(false)
                           .description("Interactive glutamine selection, instead of charged"));
    options->addOption(BooleanOption("his")
                           .store(&bHisMan_).defaultValue(false)
                           .description("Interactive histidine selection, instead of checking H-bonds"));
    options->addOption(RealOption("angle")
                           .store(&angle_).defaultValue(135.0)
                           .description("Minimum hydrogen-donor-acceptor angle for a H-bond (degrees)"));
    options->addOption(RealOption("dist")
                           .store(&distance_).defaultValue(0.3)
                           .description("Maximum donor-acceptor distance for a H-bond (nm)"));
    options->addOption(BooleanOption("una")
                           .store(&bUnA_).defaultValue(false)
                           .description("Select aromatic rings with united CH atoms on phenylalanine, tryptophane and tyrosine"));
    options->addOption(BooleanOption("sort")
                           .store(&bSort_).defaultValue(true).hidden()
                           .description("Sort the residues according to database, turning this off is dangerous as charge groups might be broken in parts"));
    options->addOption(BooleanOption("ignh")
                           .store(&bRemoveH_).defaultValue(false)
                           .description("Ignore hydrogen atoms that are in the coordinate file"));
    options->addOption(BooleanOption("missing")
                           .store(&bAllowMissing_).defaultValue(false)
                           .description("Continue when atoms are missing and bonds cannot be made, dangerous"));
    options->addOption(BooleanOption("v")
                           .store(&bVerbose_).defaultValue(false)
                           .description("Be slightly more verbose in messages"));
    options->addOption(RealOption("posrefc")
                           .store(&posre_fc_).defaultValue(1000)
                           .description("Force constant for position restraints"));
    options->addOption(EnumOption<VSitesType>("vsite")
                           .store(&enumVSites_).enumValue(VSitesEnum)
                           .description("Convert atoms to virtual sites"));
    options->addOption(BooleanOption("heavyh")
                           .store(&bHeavyH_).defaultValue(false)
                           .description("Make hydrogen atoms heavy"));
    options->addOption(BooleanOption("deuterate")
                           .store(&bDeuterate_).defaultValue(false)
                           .description("Change the mass of hydrogens to 2 amu"));
    options->addOption(BooleanOption("chargegrp")
                           .store(&bChargeGroups_).defaultValue(true)
                           .description("Use charge groups in the [REF].rtp[ref] file"));
    options->addOption(BooleanOption("cmap")
                           .store(&bCmap_).defaultValue(true)
                           .description("Use cmap torsions (if enabled in the [REF].rtp[ref] file)"));
    options->addOption(BooleanOption("renum")
                           .store(&bRenumRes_).defaultValue(false)
                           .description("Renumber the residues consecutively in the output"));
    options->addOption(BooleanOption("rtpres")
                           .store(&bRTPresname_).defaultValue(false)
                           .description("Use [REF].rtp[ref] entry names as residue names"));
    options->addOption(FileNameOption("f")
                           .legacyType(efSTX).inputFile()
                           .store(&inputConfFile_).required()
                           .defaultBasename("protein").defaultType(efPDB)
                           .description("Structure file"));
    options->addOption(FileNameOption("o")
                           .legacyType(efSTO).outputFile()
                           .store(&outputConfFile_).required()
                           .defaultBasename("conf")
                           .description("Structure file"));
    options->addOption(FileNameOption("p")
                           .legacyType(efTOP).outputFile()
                           .store(&topologyFile_).required()
                           .defaultBasename("topol")
                           .description("Topology file"));
    options->addOption(FileNameOption("i")
                           .legacyType(efITP).outputFile()
                           .store(&includeTopologyFile_).required()
                           .defaultBasename("posre")
                           .description("Include file for topology"));
    options->addOption(FileNameOption("n")
                           .legacyType(efNDX).outputFile()
                           .store(&indexOutputFile_).storeIsSet(&bIndexSet_)
                           .defaultBasename("index")
                           .description("Index file"));
    options->addOption(FileNameOption("q")
                           .legacyType(efSTO).outputFile()
                           .store(&outFile_).storeIsSet(&bOutputSet_)
                           .defaultBasename("clean").defaultType(efPDB)
                           .description("Structure file"));
}

void pdb2gmx::optionsFinished()
{
    if (inputConfFile_.empty())
    {
        GMX_THROW(InconsistentInputError("You must supply an input file"));
    }
    if (bInter_)
    {
        /* if anything changes here, also change description of -inter */
        bCysMan_ = true;
        bTerMan_ = true;
        bLysMan_ = true;
        bArgMan_ = true;
        bAspMan_ = true;
        bGluMan_ = true;
        bGlnMan_ = true;
        bHisMan_ = true;
    }

    if (bHeavyH_)
    {
        mHmult_ = 4.0;
    }
    else if (bDeuterate_)
    {
        mHmult_ = 2.0;
    }
    else
    {
        mHmult_ = 1.0;
    }

    /* Force field selection, interactive or direct */
    choose_ff(strcmp(ff_.c_str(), "select") == 0 ? nullptr : ff_.c_str(),
              forcefield_, sizeof(forcefield_),
              ffdir_, sizeof(ffdir_));

    if (strlen(forcefield_) > 0)
    {
        ffname_    = forcefield_;
        ffname_[0] = std::toupper(ffname_[0]);
    }
    else
    {
        gmx_fatal(FARGS, "Empty forcefield string");
    }
}

int pdb2gmx::run()
{
    char                       select[STRLEN];
    std::vector<DisulfideBond> ssbonds;

    int                        this_atomnum;
    int                        prev_atomnum;
    const char                *prev_atomname;
    const char                *this_atomname;
    const char                *prev_resname;
    const char                *this_resname;
    int                        prev_resnum;
    int                        this_resnum;
    char                       prev_chainid;
    char                       this_chainid;
    int                        prev_chainnumber;
    int                        this_chainnumber;
    int                        this_chainstart;
    int                        prev_chainstart;

    printf("\nUsing the %s force field in directory %s\n\n",
           ffname_, ffdir_);

    choose_watermodel(WaterEnum[enumWater_], ffdir_, &watermodel_);

    switch (enumVSites_)
    {
        case enVSites_none:
            bVsites_         = false;
            bVsiteAromatics_ = false;
            break;
        case enVSites_hydrogens:
            bVsites_         = true;
            bVsiteAromatics_ = false;
            break;
        case enVSites_aromatics:
            bVsites_         = true;
            bVsiteAromatics_ = true;
            break;
        default:
            gmx_fatal(FARGS, "Internal inconsistency: VSitesEnum='%s'", VSitesEnum[enumVSites_]);
    }       /* end switch */

    /* Open the symbol table */
    t_symtab symtab;
    open_symtab(&symtab);

    /* Residue type database */
    ResidueType rt;

    /* Read residue renaming database(s), if present */
    std::vector<std::string> rrn = fflib_search_file_end(ffdir_, ".r2b", FALSE);

    std::vector<RtpRename>   rtprename;
    for (const auto &filename : rrn)
    {
        printf("going to rename %s\n", filename.c_str());
        FILE *fp = fflib_open(filename);
        read_rtprename(filename.c_str(), fp, &rtprename);
        gmx_ffclose(fp);
    }

    /* Add all alternative names from the residue renaming database to the list
       of recognized amino/nucleic acids. */
    for (const auto &rename : rtprename)
    {
        /* Only add names if the 'standard' gromacs/iupac base name was found */
        if (auto restype = rt.optionalTypeOfNamedDatabaseResidue(rename.gmx))
        {
            rt.addResidue(rename.main, *restype);
            rt.addResidue(rename.nter, *restype);
            rt.addResidue(rename.cter, *restype);
            rt.addResidue(rename.bter, *restype);
        }
    }

    matrix      box;
    const char *watres;
    clear_mat(box);
    if (watermodel_ != nullptr && (strstr(watermodel_, "4p") ||
                                   strstr(watermodel_, "4P")))
    {
        watres = "HO4";
    }
    else if (watermodel_ != nullptr && (strstr(watermodel_, "5p") ||
                                        strstr(watermodel_, "5P")))
    {
        watres = "HO5";
    }
    else
    {
        watres = "HOH";
    }

    AtomProperties aps;
    char          *title = nullptr;
    int            ePBC;
    t_atoms        pdba_all;
    rvec          *pdbx;
    int            natom = read_pdball(inputConfFile_.c_str(), bOutputSet_, outFile_.c_str(),
                                       &title, &pdba_all, &pdbx, &ePBC, box, bRemoveH_,
                                       &symtab, &rt, watres, &aps, bVerbose_);

    if (natom == 0)
    {
        std::string message = formatString("No atoms found in pdb file %s\n", inputConfFile_.c_str());
        GMX_THROW(InconsistentInputError(message));
    }

    printf("Analyzing pdb file\n");
    int nwaterchain = 0;

    modify_chain_numbers(&pdba_all, enumChainSep_);

    int nchainmerges        = 0;

    this_atomname       = nullptr;
    this_atomnum        = -1;
    this_resname        = nullptr;
    this_resnum         = -1;
    this_chainid        = '?';
    this_chainnumber    = -1;
    this_chainstart     = 0;
    /* Keep the compiler happy */
    prev_chainstart     = 0;

    int                     numChains = 0;
    std::vector<t_pdbchain> pdb_ch;

    t_resinfo              *ri;
    bool                    bMerged = false;
    for (int i = 0; (i < natom); i++)
    {
        ri = &pdba_all.resinfo[pdba_all.atom[i].resind];

        /* TODO this should live in a helper object, and consolidate
           that with code in modify_chain_numbers */
        prev_atomname      = this_atomname;
        prev_atomnum       = this_atomnum;
        prev_resname       = this_resname;
        prev_resnum        = this_resnum;
        prev_chainid       = this_chainid;
        prev_chainnumber   = this_chainnumber;
        if (!bMerged)
        {
            prev_chainstart    = this_chainstart;
        }

        this_atomname      = *pdba_all.atomname[i];
        this_atomnum       = (pdba_all.pdbinfo != nullptr) ? pdba_all.pdbinfo[i].atomnr : i+1;
        this_resname       = *ri->name;
        this_resnum        = ri->nr;
        this_chainid       = ri->chainid;
        this_chainnumber   = ri->chainnum;

        bWat_ = gmx::equalCaseInsensitive(*ri->name, watres);

        if ((i == 0) || (this_chainnumber != prev_chainnumber) || (bWat_ != bPrevWat_))
        {
            GMX_RELEASE_ASSERT(pdba_all.pdbinfo, "Must have pdbinfo from reading a PDB file if chain number is changing");
            this_chainstart = pdba_all.atom[i].resind;
            bMerged         = false;
            if (i > 0 && !bWat_)
            {
                if (!strncmp(MergeEnum[enumMerge_], "int", 3))
                {
                    printf("Merge chain ending with residue %s%d (chain id '%c', atom %d %s) and chain starting with\n"
                           "residue %s%d (chain id '%c', atom %d %s) into a single moleculetype (keeping termini)? [n/y]\n",
                           prev_resname, prev_resnum, prev_chainid, prev_atomnum, prev_atomname,
                           this_resname, this_resnum, this_chainid, this_atomnum, this_atomname);

                    if (nullptr == fgets(select, STRLEN-1, stdin))
                    {
                        gmx_fatal(FARGS, "Error reading from stdin");
                    }
                    bMerged = (select[0] == 'y');
                }
                else if (!strncmp(MergeEnum[enumMerge_], "all", 3))
                {
                    bMerged = true;
                }
            }

            if (bMerged)
            {
                pdb_ch[numChains-1].chainstart[pdb_ch[numChains-1].nterpairs] =
                    pdba_all.atom[i].resind - prev_chainstart;
                pdb_ch[numChains-1].nterpairs++;
                pdb_ch[numChains-1].chainstart.resize(pdb_ch[numChains-1].nterpairs+1);
                nchainmerges++;
            }
            else
            {
                /* set natom for previous chain */
                if (numChains > 0)
                {
                    pdb_ch[numChains-1].natom = i-pdb_ch[numChains-1].start;
                }
                if (bWat_)
                {
                    nwaterchain++;
                    ri->chainid = ' ';
                }
                /* check if chain identifier was used before */
                for (int j = 0; (j < numChains); j++)
                {
                    if (pdb_ch[j].chainid != ' ' && pdb_ch[j].chainid == ri->chainid)
                    {
                        printf("WARNING: Chain identifier '%c' is used in two non-sequential blocks.\n"
                               "They will be treated as separate chains unless you reorder your file.\n",
                               ri->chainid);
                    }
                }
                t_pdbchain newChain;
                newChain.chainid  = ri->chainid;
                newChain.chainnum = ri->chainnum;
                newChain.start    = i;
                newChain.bAllWat  = bWat_;
                if (bWat_)
                {
                    newChain.nterpairs = 0;
                }
                else
                {
                    newChain.nterpairs = 1;
                }
                newChain.chainstart.resize(newChain.nterpairs+1);
                /* modified [numChains] to [0] below */
                newChain.chainstart[0] = 0;
                pdb_ch.push_back(newChain);
                numChains++;
            }
        }
        bPrevWat_ = bWat_;
    }
    pdb_ch.back().natom = natom-pdb_ch.back().start;

    /* set all the water blocks at the end of the chain */
    std::vector<int> swap_index(numChains);
    int              j = 0;
    for (int i = 0; i < numChains; i++)
    {
        if (!pdb_ch[i].bAllWat)
        {
            swap_index[j] = i;
            j++;
        }
    }
    for (int i = 0; i < numChains; i++)
    {
        if (pdb_ch[i].bAllWat)
        {
            swap_index[j] = i;
            j++;
        }
    }
    if (nwaterchain > 1)
    {
        printf("Moved all the water blocks to the end\n");
    }

    t_atoms             *pdba;
    std::vector<t_chain> chains(numChains);
    /* copy pdb data and x for all chains */
    for (int i = 0; (i < numChains); i++)
    {
        int si                   = swap_index[i];
        chains[i].chainid    = pdb_ch[si].chainid;
        chains[i].chainnum   = pdb_ch[si].chainnum;
        chains[i].bAllWat    = pdb_ch[si].bAllWat;
        chains[i].nterpairs  = pdb_ch[si].nterpairs;
        chains[i].chainstart = pdb_ch[si].chainstart;
        chains[i].ntdb.clear();
        chains[i].ctdb.clear();
        chains[i].r_start.resize(pdb_ch[si].nterpairs);
        chains[i].r_end.resize(pdb_ch[si].nterpairs);

        snew(chains[i].pdba, 1);
        init_t_atoms(chains[i].pdba, pdb_ch[si].natom, true);
        for (j = 0; j < chains[i].pdba->nr; j++)
        {
            chains[i].pdba->atom[j]     = pdba_all.atom[pdb_ch[si].start+j];
            chains[i].pdba->atomname[j] =
                put_symtab(&symtab, *pdba_all.atomname[pdb_ch[si].start+j]);
            chains[i].pdba->pdbinfo[j] = pdba_all.pdbinfo[pdb_ch[si].start+j];
            chains[i].x.emplace_back(pdbx[pdb_ch[si].start+j]);
        }
        /* Re-index the residues assuming that the indices are continuous */
        int k                = chains[i].pdba->atom[0].resind;
        int nres             = chains[i].pdba->atom[chains[i].pdba->nr-1].resind - k + 1;
        chains[i].pdba->nres = nres;
        for (int j = 0; j < chains[i].pdba->nr; j++)
        {
            chains[i].pdba->atom[j].resind -= k;
        }
        srenew(chains[i].pdba->resinfo, nres);
        for (int j = 0; j < nres; j++)
        {
            chains[i].pdba->resinfo[j]      = pdba_all.resinfo[k+j];
            chains[i].pdba->resinfo[j].name = put_symtab(&symtab, *pdba_all.resinfo[k+j].name);
            /* make all chain identifiers equal to that of the chain */
            chains[i].pdba->resinfo[j].chainid = pdb_ch[si].chainid;
        }
    }

    if (nchainmerges > 0)
    {
        printf("\nMerged chains into joint molecule definitions at %d places.\n\n",
               nchainmerges);
    }

    printf("There are %d chains and %d blocks of water and "
           "%d residues with %d atoms\n",
           numChains-nwaterchain, nwaterchain,
           pdba_all.nres, natom);

    printf("\n  %5s  %4s %6s\n", "chain", "#res", "#atoms");
    for (int i = 0; (i < numChains); i++)
    {
        printf("  %d '%c' %5d %6d  %s\n",
               i+1, chains[i].chainid ? chains[i].chainid : '-',
               chains[i].pdba->nres, chains[i].pdba->nr,
               chains[i].bAllWat ? "(only water)" : "");
    }
    printf("\n");

    check_occupancy(&pdba_all, inputConfFile_.c_str(), bVerbose_);

    /* Read atomtypes... */
    PreprocessingAtomTypes atype = read_atype(ffdir_, &symtab);

    /* read residue database */
    printf("Reading residue database... (%s)\n", forcefield_);
    std::vector<std::string>       rtpf  = fflib_search_file_end(ffdir_, ".rtp", true);
    std::vector<PreprocessResidue> rtpFFDB;
    for (const auto &filename : rtpf)
    {
        readResidueDatabase(filename, &rtpFFDB, &atype, &symtab, false);
    }
    if (bNewRTP_)
    {
        /* Not correct with multiple rtp input files with different bonded types */
        FILE *fp = gmx_fio_fopen("new.rtp", "w");
        print_resall(fp, rtpFFDB, atype);
        gmx_fio_fclose(fp);
    }

    /* read hydrogen database */
    std::vector<MoleculePatchDatabase> ah;
    read_h_db(ffdir_, &ah);

    /* Read Termini database... */
    std::vector<MoleculePatchDatabase>                  ntdb;
    std::vector<MoleculePatchDatabase>                  ctdb;
    std::vector<MoleculePatchDatabase *>                tdblist;
    int                                                 nNtdb = read_ter_db(ffdir_, 'n', &ntdb, &atype);
    int                                                 nCtdb = read_ter_db(ffdir_, 'c', &ctdb, &atype);

    FILE                                               *top_file = gmx_fio_fopen(topologyFile_.c_str(), "w");

    print_top_header(top_file, topologyFile_.c_str(), FALSE, ffdir_, mHmult_);

    t_chain                *cc;
    std::vector<gmx::RVec>  x;
    /* new pdb datastructure for sorting. */
    t_atoms               **sortAtoms  = nullptr;
    t_atoms               **localAtoms = nullptr;
    snew(sortAtoms, numChains);
    snew(localAtoms, numChains);
    for (int chain = 0; (chain < numChains); chain++)
    {
        cc = &(chains[chain]);

        /* set pdba, natom and nres to the current chain */
        pdba     = cc->pdba;
        x        = cc->x;
        natom    = cc->pdba->nr;
        int nres = cc->pdba->nres;

        if (cc->chainid && ( cc->chainid != ' ' ) )
        {
            printf("Processing chain %d '%c' (%d atoms, %d residues)\n",
                   chain+1, cc->chainid, natom, nres);
        }
        else
        {
            printf("Processing chain %d (%d atoms, %d residues)\n",
                   chain+1, natom, nres);
        }

        process_chain(pdba, x, bUnA_, bUnA_, bUnA_, bLysMan_, bAspMan_, bGluMan_,
                      bHisMan_, bArgMan_, bGlnMan_, angle_, distance_, &symtab,
                      rtprename);

        cc->chainstart[cc->nterpairs] = pdba->nres;
        j = 0;
        for (int i = 0; i < cc->nterpairs; i++)
        {
            find_nc_ter(pdba, cc->chainstart[i], cc->chainstart[i+1],
                        &(cc->r_start[j]), &(cc->r_end[j]), &rt);

            if (cc->r_start[j] >= 0 && cc->r_end[j] >= 0)
            {
                j++;
            }
        }
        cc->nterpairs = j;
        if (cc->nterpairs == 0)
        {
            printf("Problem with chain definition, or missing terminal residues.\n"
                   "This chain does not appear to contain a recognized chain molecule.\n"
                   "If this is incorrect, you can edit residuetypes.dat to modify the behavior.\n");
        }

        /* Check for disulfides and other special bonds */
        ssbonds = makeDisulfideBonds(pdba, gmx::as_rvec_array(x.data()), bCysMan_, bVerbose_);

        if (!rtprename.empty())
        {
            rename_resrtp(pdba, cc->nterpairs, cc->r_start, cc->r_end, rtprename,
                          &symtab, bVerbose_);
        }

        for (int i = 0; i < cc->nterpairs; i++)
        {
            /* Set termini.
             * We first apply a filter so we only have the
             * termini that can be applied to the residue in question
             * (or a generic terminus if no-residue specific is available).
             */
            /* First the N terminus */
            if (nNtdb > 0)
            {
                tdblist = filter_ter(ntdb,
                                     *pdba->resinfo[cc->r_start[i]].name);
                if (tdblist.empty())
                {
                    printf("No suitable end (N or 5') terminus found in database - assuming this residue\n"
                           "is already in a terminus-specific form and skipping terminus selection.\n");
                    cc->ntdb.push_back(nullptr);
                }
                else
                {
                    if (bTerMan_ && !tdblist.empty())
                    {
                        sprintf(select, "Select start terminus type for %s-%d",
                                *pdba->resinfo[cc->r_start[i]].name,
                                pdba->resinfo[cc->r_start[i]].nr);
                        cc->ntdb.push_back(choose_ter(tdblist, select));
                    }
                    else
                    {
                        cc->ntdb.push_back(tdblist[0]);
                    }

                    printf("Start terminus %s-%d: %s\n",
                           *pdba->resinfo[cc->r_start[i]].name,
                           pdba->resinfo[cc->r_start[i]].nr,
                           (cc->ntdb[i])->name.c_str());
                    tdblist.clear();
                }
            }
            else
            {
                cc->ntdb.push_back(nullptr);
            }

            /* And the C terminus */
            if (nCtdb > 0)
            {
                tdblist = filter_ter(ctdb,
                                     *pdba->resinfo[cc->r_end[i]].name);
                if (tdblist.empty())
                {
                    printf("No suitable end (C or 3') terminus found in database - assuming this residue\n"
                           "is already in a terminus-specific form and skipping terminus selection.\n");
                    cc->ctdb.push_back(nullptr);
                }
                else
                {
                    if (bTerMan_ && !tdblist.empty())
                    {
                        sprintf(select, "Select end terminus type for %s-%d",
                                *pdba->resinfo[cc->r_end[i]].name,
                                pdba->resinfo[cc->r_end[i]].nr);
                        cc->ctdb.push_back(choose_ter(tdblist, select));
                    }
                    else
                    {
                        cc->ctdb.push_back(tdblist[0]);
                    }
                    printf("End terminus %s-%d: %s\n",
                           *pdba->resinfo[cc->r_end[i]].name,
                           pdba->resinfo[cc->r_end[i]].nr,
                           (cc->ctdb[i])->name.c_str());
                    tdblist.clear();
                }
            }
            else
            {
                cc->ctdb.push_back(nullptr);
            }
        }
        std::vector<MoleculePatchDatabase> hb_chain;
        /* lookup hackblocks and rtp for all residues */
        std::vector<PreprocessResidue>     restp_chain;
        get_hackblocks_rtp(&hb_chain, &restp_chain,
                           rtpFFDB, pdba->nres, pdba->resinfo,
                           cc->nterpairs, &symtab, cc->ntdb, cc->ctdb, cc->r_start, cc->r_end,
                           bAllowMissing_);
        /* ideally, now we would not need the rtp itself anymore, but do
           everything using the hb and restp arrays. Unfortunately, that
           requires some re-thinking of code in gen_vsite.c, which I won't
           do now :( AF 26-7-99 */

        rename_atoms(nullptr, ffdir_,
                     pdba, &symtab, restp_chain, false, &rt, false, bVerbose_);

        match_atomnames_with_rtp(restp_chain, hb_chain, pdba, &symtab, x, bVerbose_);

        if (bSort_)
        {
            char     **gnames;
            t_blocka  *block = new_blocka();
            snew(gnames, 1);
            sort_pdbatoms(restp_chain, natom, &pdba, &sortAtoms[chain], &x, block, &gnames);
            remove_duplicate_atoms(pdba, x, bVerbose_);
            if (bIndexSet_)
            {
                if (bRemoveH_)
                {
                    fprintf(stderr, "WARNING: with the -remh option the generated "
                            "index file (%s) might be useless\n"
                            "(the index file is generated before hydrogens are added)",
                            indexOutputFile_.c_str());
                }
                write_index(indexOutputFile_.c_str(), block, gnames, false, 0);
            }
            for (int i = 0; i < block->nr; i++)
            {
                sfree(gnames[i]);
            }
            sfree(gnames);
            done_blocka(block);
            sfree(block);
        }
        else
        {
            fprintf(stderr, "WARNING: "
                    "without sorting no check for duplicate atoms can be done\n");
        }

        /* Generate Hydrogen atoms (and termini) in the sequence */
        printf("Generating any missing hydrogen atoms and/or adding termini.\n");
        add_h(&pdba, &localAtoms[chain], &x, ah, &symtab,
              cc->nterpairs, cc->ntdb, cc->ctdb, cc->r_start, cc->r_end, bAllowMissing_);
        printf("Now there are %d residues with %d atoms\n",
               pdba->nres, pdba->nr);

        /* make up molecule name(s) */

        int         k = (cc->nterpairs > 0 && cc->r_start[0] >= 0) ? cc->r_start[0] : 0;

        std::string restype = rt.typeOfNamedDatabaseResidue(*pdba->resinfo[k].name);

        std::string molname;
        std::string suffix;
        if (cc->bAllWat)
        {
            molname = "Water";
        }
        else
        {
            this_chainid = cc->chainid;

            /* Add the chain id if we have one */
            if (this_chainid != ' ')
            {
                suffix.append(formatString("_chain_%c", this_chainid));
            }

            /* Check if there have been previous chains with the same id */
            int nid_used = 0;
            for (int k = 0; k < chain; k++)
            {
                if (cc->chainid == chains[k].chainid)
                {
                    nid_used++;
                }
            }
            /* Add the number for this chain identifier if there are multiple copies */
            if (nid_used > 0)
            {
                suffix.append(formatString("%d", nid_used+1));
            }

            if (suffix.length() > 0)
            {
                molname.append(restype);
                molname.append(suffix);
            }
            else
            {
                molname = restype;
            }
        }
        std::string itp_fn   = topologyFile_;;
        std::string posre_fn = includeTopologyFile_;
        if ((numChains-nwaterchain > 1) && !cc->bAllWat)
        {
            bITP_  = true;
            printf("Chain time...\n");
            //construct the itp file name
            itp_fn = stripSuffixIfPresent(itp_fn, ".top");
            itp_fn.append("_");
            itp_fn.append(molname);
            itp_fn.append(".itp");
            //now do the same for posre
            posre_fn = stripSuffixIfPresent(posre_fn, ".itp");
            posre_fn.append("_");
            posre_fn.append(molname);
            posre_fn.append(".itp");
            if (posre_fn == itp_fn)
            {
                posre_fn = Path::concatenateBeforeExtension(posre_fn, "_pr");
            }
            incls_.emplace_back();
            incls_.back()    = itp_fn;
            itp_file_        = gmx_fio_fopen(itp_fn.c_str(), "w");
        }
        else
        {
            bITP_ = false;
        }

        mols_.emplace_back();
        if (cc->bAllWat)
        {
            mols_.back().name = "SOL";
            mols_.back().nr   = pdba->nres;
        }
        else
        {
            mols_.back().name = molname;
            mols_.back().nr   = 1;
        }

        if (bITP_)
        {
            print_top_comment(itp_file_, itp_fn.c_str(), ffdir_, true);
        }

        FILE *top_file2;
        if (cc->bAllWat)
        {
            top_file2 = nullptr;
        }
        else if (bITP_)
        {
            top_file2 = itp_file_;
        }
        else
        {
            top_file2 = top_file;
        }

        pdb2top(top_file2, posre_fn.c_str(), molname.c_str(), pdba, &x, &atype, &symtab,
                rtpFFDB,
                restp_chain, hb_chain,
                bAllowMissing_,
                bVsites_, bVsiteAromatics_, ffdir_,
                mHmult_, ssbonds,
                long_bond_dist_, short_bond_dist_, bDeuterate_, bChargeGroups_, bCmap_,
                bRenumRes_, bRTPresname_);

        if (!cc->bAllWat)
        {
            write_posres(posre_fn.c_str(), pdba, posre_fc_);
        }

        if (bITP_)
        {
            gmx_fio_fclose(itp_file_);
        }

        /* pdba and natom have been reassigned somewhere so: */
        cc->pdba = pdba;
        cc->x    = x;
    }

    if (watermodel_ == nullptr)
    {
        for (int chain = 0; chain < numChains; chain++)
        {
            if (chains[chain].bAllWat)
            {
                auto message = formatString("You have chosen not to include a water model, "
                                            "but there is water in the input file. Select a "
                                            "water model or remove the water from your input file.");
                GMX_THROW(InconsistentInputError(message));
            }
        }
    }
    else
    {
        std::string waterFile = formatString("%s%c%s.itp", ffdir_, DIR_SEPARATOR, watermodel_);
        if (!fflib_fexist(waterFile))
        {
            auto message = formatString("The topology file '%s' for the selected water "
                                        "model '%s' can not be found in the force field "
                                        "directory. Select a different water model.",
                                        waterFile.c_str(), watermodel_);
            GMX_THROW(InconsistentInputError(message));
        }
    }

    print_top_mols(top_file, title, ffdir_, watermodel_, incls_, mols_);
    gmx_fio_fclose(top_file);

    /* now merge all chains back together */
    natom     = 0;
    int nres  = 0;
    for (int i = 0; (i < numChains); i++)
    {
        natom += chains[i].pdba->nr;
        nres  += chains[i].pdba->nres;
    }
    t_atoms *atoms;
    snew(atoms, 1);
    init_t_atoms(atoms, natom, false);
    for (int i = 0; i < atoms->nres; i++)
    {
        sfree(atoms->resinfo[i].name);
    }
    atoms->nres = nres;
    srenew(atoms->resinfo, nres);
    x.clear();
    int k = 0;
    int l = 0;
    for (int i = 0; (i < numChains); i++)
    {
        if (numChains > 1)
        {
            printf("Including chain %d in system: %d atoms %d residues\n",
                   i+1, chains[i].pdba->nr, chains[i].pdba->nres);
        }
        for (int j = 0; (j < chains[i].pdba->nr); j++)
        {
            atoms->atom[k]         = chains[i].pdba->atom[j];
            atoms->atom[k].resind += l; /* l is processed nr of residues */
            atoms->atomname[k]     = chains[i].pdba->atomname[j];
            atoms->resinfo[atoms->atom[k].resind].chainid = chains[i].chainid;
            x.push_back(chains[i].x[j]);
            k++;
        }
        for (int j = 0; (j < chains[i].pdba->nres); j++)
        {
            atoms->resinfo[l] = chains[i].pdba->resinfo[j];
            if (bRTPresname_)
            {
                atoms->resinfo[l].name = atoms->resinfo[l].rtp;
            }
            l++;
        }
        done_atom(chains[i].pdba);
    }

    if (numChains > 1)
    {
        fprintf(stderr, "Now there are %d atoms and %d residues\n", k, l);
        print_sums(atoms, true);
    }

    rvec box_space;
    fprintf(stderr, "\nWriting coordinate file...\n");
    clear_rvec(box_space);
    if (box[0][0] == 0)
    {
        make_new_box(atoms->nr, gmx::as_rvec_array(x.data()), box, box_space, false);
    }
    write_sto_conf(outputConfFile_.c_str(), title, atoms, gmx::as_rvec_array(x.data()), nullptr, ePBC, box);

    done_symtab(&symtab);
    done_atom(&pdba_all);
    done_atom(atoms);
    for (int chain = 0; chain < numChains; chain++)
    {
        sfree(sortAtoms[chain]);
        sfree(localAtoms[chain]);
    }
    sfree(sortAtoms);
    sfree(localAtoms);
    sfree(atoms);
    sfree(title);
    sfree(pdbx);
    printf("\t\t--------- PLEASE NOTE ------------\n");
    printf("You have successfully generated a topology from: %s.\n",
           inputConfFile_.c_str());
    if (watermodel_ != nullptr)
    {
        printf("The %s force field and the %s water model are used.\n",
               ffname_, watermodel_);
        sfree(watermodel_);
    }
    else
    {
        printf("The %s force field is used.\n",
               ffname_);
    }
    printf("\t\t--------- ETON ESAELP ------------\n");

    return 0;
}

}   // namespace

const char pdb2gmxInfo::name[]             = "pdb2gmx";
const char pdb2gmxInfo::shortDescription[] =
    "Convert coordinate files to topology and FF-compliant coordinate files";
ICommandLineOptionsModulePointer pdb2gmxInfo::create()
{
    return std::make_unique<pdb2gmx>();
}

} // namespace gmx
