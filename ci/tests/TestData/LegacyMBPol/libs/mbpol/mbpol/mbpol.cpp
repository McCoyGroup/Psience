#ifdef HAVE_CONFIG_H
#   include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>

#define not_VERBOSE 1

#ifdef VERBOSE
#   include <iostream>
#endif /* VERBOSE */

#include "ps.h"
#include "mbpol.h"

#include "x2b-v9x.h"
#include "x3b-v2x.h"
#include "x2b-dispersion.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//Changed after publication. Old value 100kcal/mol (too high). 
const double E1max =  60.0; // kcal/mol

}

////////////////////////////////////////////////////////////////////////////////

namespace x2o {

//----------------------------------------------------------------------------//

double mbpol::operator()(size_t nw, const double* pos) const
{
    assert(nw > 0 && pos);

    double E1(0);

    bool good = true;

    for (size_t i = 0; i < nw; ++i) {
        const size_t i9 = 9*i;

        const double this_E1 = ps::pot_nasa(pos + i9, 0);
        E1 += this_E1;

        if (this_E1 > E1max)
            good = false;
    }

    double Eelec(0), Eind(0), Edisp(0), E2poly(0), E3poly(0);

    if (good) {
        m_ttm4(nw, pos, Eelec, 0, Eind, 0);

        for (size_t i = 0; i < nw; ++i) {
            const size_t i9 = 9*i;

            for (size_t j = i + 1; j < nw; ++j) {
                const size_t j9 = 9*j;
                Edisp += x2b_dispersion::eval(pos + i9, pos + j9);
                E2poly += x2b_v9x::eval(pos + i9, pos + j9);

                for (size_t k = j + 1; k < nw; ++k) {
                    const size_t k9 = 9*k;
                    E3poly += x3b_v2x::eval(pos + i9, pos + j9, pos + k9);
                }
            }
        }
    } // good

#   ifdef VERBOSE
    std::cout << "\n    E1 = " << E1
              << "\n Eelec = " << Eelec
              << "\n  Eind = " << Eind
              << "\n Edisp = " << Edisp
              << "\nE2poly = " << E2poly
              << "\nE3poly = " << E3poly
              << std::endl;
#   endif /* VERBOSE */

    return E1 + Eelec + Eind + Edisp + E2poly + E3poly;
}

//----------------------------------------------------------------------------//

double mbpol::operator()(size_t nw, const double* pos, double* grd) const
{
    assert(nw > 0 && pos && grd);

    double Etot(0);
    bool good(true);

    for (size_t i = 0; i < nw; ++i) {
        const size_t i9 = 9*i;

        const double E1 = ps::pot_nasa(pos + i9, grd + i9);
        if (E1 > E1max)
            good = false;

        Etot += E1;
    }

    if (good) {
        double Eelec(0), Eind(0), gEelec[9*nw], gEind[9*nw];

        m_ttm4(nw, pos, Eelec, gEelec, Eind, gEind);
        for (size_t i = 0; i < 9*nw; ++i)
            grd[i] += gEelec[i] + gEind[i];

        Etot += Eelec + Eind;

        for (size_t i = 0; i < nw; ++i) {
            const size_t i9 = 9*i;

            for (size_t j = i + 1; j < nw; ++j) {
                const size_t j9 = 9*j;
                Etot += x2b_v9x::eval(pos + i9, pos + j9, grd + i9, grd + j9)
                      + x2b_dispersion::eval(pos + i9, pos + j9,
                                             grd + i9, grd + j9);

                for (size_t k = j + 1; k < nw; ++k) {
                    const size_t k9 = 9*k;
                    Etot += x3b_v2x::eval(pos + i9, pos + j9, pos + k9,
                                          grd + i9, grd + j9, grd + k9);
                }
            }
        }
    } // good

    return Etot;
}

//----------------------------------------------------------------------------//

} // namespace x2o

////////////////////////////////////////////////////////////////////////////////
