#ifdef HAVE_CONFIG_H
#   include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>
#include <algorithm>

#include "macros.h"

#include "ps.h"
#include "ttm4-es.h"
#include "ttm4-smear.h"

#include "constants.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

const double CHARGECON = constants::CHARGECON;

//----------------------------------------------------------------------------//

// polarization parameters

const double polarO = 1.310;
const double polarH = 0.294;
const double polarM = 0.0;

const double polar_sqrt[4] = {
    std::sqrt(polarO),
    std::sqrt(polarH),
    std::sqrt(polarH),
    std::sqrt(polarM)
};

const double polfacO = polarO;
const double polfacH = polarH;
const double polfacM = polarO;

const double AA[4] = {
    std::pow(polfacO, 1.0/6.0),
    std::pow(polfacH, 1.0/6.0),
    std::pow(polfacH, 1.0/6.0),
    std::pow(polfacM, 1.0/6.0)
};

const double aCC = 0.4;
const double aCD = 0.4;

const double aDD_inter    = 0.055;
const double aDD_intra_HH = 0.055;
const double aDD_intra_OH = 0.626;

//----------------------------------------------------------------------------//

// M-site positioning (TTM2.1-F)
const double gammaM = 0.426706882;

const double gamma1 = 1.0 - gammaM;
const double gamma2 = gammaM/2;

//----------------------------------------------------------------------------//

inline void compute_M_site_crd
    (const double* RESTRICT O,
     const double* RESTRICT H1,
     const double* RESTRICT H2,
     double* RESTRICT M)
{
    for (size_t i = 0; i < 3; ++i)
        M[i] = gamma1*O[i] + gamma2*(H1[i] + H2[i]);
}

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace ttm {

//
// for CC-pol monomer geometry
//

const double ttm4_es::molecular_polarizability = 1.43016; // A^3
const double ttm4_es::molecular_dipole = 1.867895838174; // D

//----------------------------------------------------------------------------//

ttm4_es::ttm4_es()
: m_nw(0),
  m_mem(0)
{
}

//----------------------------------------------------------------------------//

ttm4_es::~ttm4_es()
{
    if (m_nw > 0)
        delete[] m_mem;
}

//----------------------------------------------------------------------------//

void ttm4_es::allocate(size_t nw)
{
    if (m_nw > nw)
        return;

    delete[] m_mem;

    // see operator() below
    m_mem = new double[95*nw + 144*nw*nw];

    m_nw = nw;
}

//----------------------------------------------------------------------------//

void ttm4_es::operator()
    (size_t nw, const double* RESTRICT crd,
     double& E_elec, double* RESTRICT gE_elec,
     double& E_ind, double* RESTRICT gE_ind)
{
    assert(crd || nw == 0);

    allocate(nw);

    // setup pointers

    const size_t natom = 4*nw; // O H H M
    const size_t natom3 = 3*natom;

    double* xyz    = m_mem;           // 12*nw
    double* charge = xyz + natom3;    //  4*nw
    double* dipole = charge + natom;  // 12*nw
    double* diag   = dipole + natom3; // 12*nw
    double* grdq   = diag + natom3;   // 12*nw
    double* grd4   = grdq + 27*nw;    // 27*nw
    double* Efq    = grd4 + natom3;   // 12*nw
    double* phi    = Efq + natom3;    //  4*nw
    double* ddt    = phi + natom;     // 144*nw*nw

    // zero out Efq/phi

    std::fill(Efq, Efq + natom3 + natom + natom3*natom3, 0.0);

    // compute M-sites and charges

    for (size_t n = 0; n < nw; ++n) {
        const size_t n4 = 4*n;
        const size_t n9 = 9*n;
        const size_t n12 = 12*n;

        std::copy(crd + n9, crd + n9 + 9, xyz + n12);
        compute_M_site_crd(crd + n9, crd + n9 + 3, crd + n9 + 6,
                           xyz + n12 + 9);

        double q3[3], dq3[27];
        if (gE_elec != 0 || gE_ind != 0)
            ps::dms_nasa(0.0, 0.0, 0.0,
                         crd + n9, q3, dq3, false);
        else
            ps::dms_nasa(0.0, 0.0, 0.0,
                         crd + n9, q3, 0, false);

        // TTM2.1-F assignment
        const double tmp = gamma2/gamma1;

        charge[n4 + 0] = 0.0;                                     // O
        charge[n4 + 1] = CHARGECON*(q3[1] + tmp*(q3[1] + q3[2])); // H1
        charge[n4 + 2] = CHARGECON*(q3[2] + tmp*(q3[1] + q3[2])); // H2
        charge[n4 + 3] = CHARGECON*q3[0]/gamma1;                  // M

#define DQ3(i,j,k) dq3[k + 3*(j + 3*i)]
#define GRDQ(i,j,k) grdq[k + 3*(j + 3*(i + 3*n))]

        if (gE_elec == 0 && gE_ind == 0)
            continue;

        for (size_t k = 0; k < 3; ++k) {
            GRDQ(0, 0, k) = DQ3(0, 0, k) + tmp*(DQ3(0, 0, k) + DQ3(0, 1, k));
            GRDQ(1, 0, k) = DQ3(1, 0, k) + tmp*(DQ3(1, 0, k) + DQ3(1, 1, k));
            GRDQ(2, 0, k) = DQ3(2, 0, k) + tmp*(DQ3(2, 0, k) + DQ3(2, 1, k));

            GRDQ(0, 1, k) = DQ3(0, 1, k) + tmp*(DQ3(0, 1, k) + DQ3(0, 0, k));
            GRDQ(1, 1, k) = DQ3(1, 1, k) + tmp*(DQ3(1, 1, k) + DQ3(1, 0, k));
            GRDQ(2, 1, k) = DQ3(2, 1, k) + tmp*(DQ3(2, 1, k) + DQ3(2, 0, k));

            GRDQ(0, 2, k) = DQ3(0, 2, k) - 2*tmp*(DQ3(0, 0, k) + DQ3(0, 1, k));
            GRDQ(1, 2, k) = DQ3(1, 2, k) - 2*tmp*(DQ3(1, 0, k) + DQ3(1, 1, k));
            GRDQ(2, 2, k) = DQ3(2, 2, k) - 2*tmp*(DQ3(2, 0, k) + DQ3(2, 1, k));
        }
    }

    if (gE_elec != 0 || gE_ind != 0) {
        for (size_t n = 0; n < 27*nw; ++n)
            grdq[n] *= CHARGECON;
    }

    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;

        for (size_t j = i + 1; j < natom; ++j) {
            const size_t j3 = 3*j;

            double Rij[3], Rsq(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = xyz[i3 + k] - xyz[j3 + k];
                Rsq += Rij[k]*Rij[k];
            }

            // charge-charge

            const bool ij_from_same_water = (i/4 == j/4);

            if (!ij_from_same_water) {
                double ts0, ts1;
                ttm4::smear01(std::sqrt(Rsq), AA[i%4]*AA[j%4], aCC, ts0, ts1);

                phi[i] += ts0*charge[j];
                phi[j] += ts0*charge[i];

                for (size_t k = 0; k < 3; ++k) {
                    Efq[i3 + k] += charge[j]*ts1*Rij[k];
                    Efq[j3 + k] -= charge[i]*ts1*Rij[k];
                }
            }

            // dipole-dipole tensor

            const double aDD = ij_from_same_water ?
                (i%4 == 0 ? aDD_intra_OH : aDD_intra_HH) : aDD_inter;

            double ts1, ts2;
            ttm4::smear2(std::sqrt(Rsq), AA[i%4]*AA[j%4], aDD, ts1, ts2);

            double dd3[3][3];

            dd3[0][0] = 3.0*ts2*Rij[0]*Rij[0] - ts1;
            dd3[1][1] = 3.0*ts2*Rij[1]*Rij[1] - ts1;
            dd3[2][2] = 3.0*ts2*Rij[2]*Rij[2] - ts1;
            dd3[0][1] = 3.0*ts2*Rij[0]*Rij[1];
            dd3[0][2] = 3.0*ts2*Rij[0]*Rij[2];
            dd3[1][2] = 3.0*ts2*Rij[1]*Rij[2];
            dd3[1][0] = dd3[0][1];
            dd3[2][0] = dd3[0][2];
            dd3[2][1] = dd3[1][2];

            const double aiaj = polar_sqrt[i%4]*polar_sqrt[j%4];
            for (size_t k = 0; k < 3; ++k)
                for (size_t l = 0; l < 3; ++l)
                    ddt[natom3*(i3 + k) + j3 + l] = -aiaj*dd3[k][l];
        }
    }

    // populate the diagonal of the dipole-dipole tensor

    for (size_t i = 0; i < natom3; ++i)
        ddt[natom3*i + i] = 1.0;

    // perform Cholesky decomposition ddt = L*L^T
    // storing L in the lower triangle of ddt and
    // its diagonal in diag

    for (size_t i = 0; i < natom3; ++i) {
        for (size_t j = i; j < natom3; ++j) {

            double sum = ddt[i*natom3 + j];
            for (size_t k = 0; k < i; ++k)
                sum -= ddt[i*natom3 + k]*ddt[j*natom3 + k];

            if (i == j) {
                assert(sum > 0.0);
                diag[i] = std::sqrt(sum);
            } else {
                ddt[j*natom3 + i] = sum/diag[i];
            }
        }
    }

    // solve L*y = sqrt(a)*Efq storing y in dipole[]

    for (size_t i = 0; i < natom3; ++i) {
        double sum = polar_sqrt[(i/3)%4]*Efq[i];
        for (size_t k = 0; k < i; ++k)
            sum -= ddt[i*natom3 + k]*dipole[k];

        dipole[i] = sum/diag[i];
    }

    // solve L^T*x = y

    for (size_t i = natom3; i > 0; --i) {
        const size_t i1 = i - 1;
        double sum = dipole[i1];
        for (size_t k = i; k < natom3; ++k)
            sum -= ddt[k*natom3 + i1]*dipole[k];

        dipole[i1] = sum/diag[i1];
    }

    for (size_t i = 0; i < natom3; ++i)
        dipole[i] *= polar_sqrt[(i/3)%4];

    // compute the energies

    E_elec = 0.0;
    for (size_t i = 0; i < natom; ++i)
        E_elec += phi[i]*charge[i];

    E_elec /= 2;

    E_ind = 0.0;
    for (size_t i = 0; i < natom3; ++i)
        E_ind -= dipole[i]*Efq[i];

    E_ind /= 2;

    // *****************
    // *** gradients ***
    // *****************

    // charge-charge interactions

    if (gE_elec) {
        for (size_t i = 0; i < natom3; ++i)
            grd4[i] = -charge[i/3]*Efq[i];

        for (size_t n = 0; n < nw; ++n) {

            const size_t io  = 9*n + 0;
            const size_t ih1 = 9*n + 3;
            const size_t ih2 = 9*n + 6;

            // O H H
            for (size_t k = 0; k < 9; ++k)
                gE_elec[9*n + k] = grd4[12*n + k];

            // M
            for (size_t k = 0; k < 3; ++k) {
                gE_elec[io  + k] += gamma1*grd4[12*n + 9 + k]; // O
                gE_elec[ih1 + k] += gamma2*grd4[12*n + 9 + k]; // H
                gE_elec[ih2 + k] += gamma2*grd4[12*n + 9 + k]; // H
            }

            for (size_t k = 0; k < 3; ++k) {
                gE_elec[ih1 + k] += GRDQ(0, 0, k)*phi[4*n + 1]  // phi(h1)
                                  + GRDQ(0, 1, k)*phi[4*n + 2]  // phi(h2)
                                  + GRDQ(0, 2, k)*phi[4*n + 3]; // phi(M)

                gE_elec[ih2 + k] += GRDQ(1, 0, k)*phi[4*n + 1]  // phi(h1)
                                  + GRDQ(1, 1, k)*phi[4*n + 2]  // phi(h2)
                                  + GRDQ(1, 2, k)*phi[4*n + 3]; // phi(M)

                gE_elec[io + k] += GRDQ(2, 0, k)*phi[4*n + 1]  // phi(h1)
                                 + GRDQ(2, 1, k)*phi[4*n + 2]  // phi(h2)
                                 + GRDQ(2, 2, k)*phi[4*n + 3]; // phi(M)
            }
        }
    } // if (gE_elec)

    if (gE_ind == 0)
        return;

    // charge-dipole interactions

    std::fill(phi, phi + natom, 0.0);
    std::fill(grd4, grd4 + natom3, 0.0);

    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;

        for (size_t j = 0; j < natom; ++j) {
            const size_t j3 = 3*j;
            const double qj = charge[j];

            const bool skip_ij = (i/4 == j/4);

            if (skip_ij)
                continue; // skip this (i, j) pair

            double Rij[3], Rsq(0), diR(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = xyz[i3 + k] - xyz[j3 + k];
                Rsq += Rij[k]*Rij[k];
                diR += dipole[i3 + k]*Rij[k];
            }

            double ts1, ts2;
            ttm4::smear2(std::sqrt(Rsq), AA[i%4]*AA[j%4], aCD, ts1, ts2);

            for (size_t k = 0; k < 3; ++k) {
                const double derij =
                    qj*(3*ts2*diR*Rij[k] - ts1*dipole[i3 + k]);

                grd4[i3 + k] += derij;
                grd4[j3 + k] -= derij;
            }

            phi[j] -= ts1*diR;
        }
    }

    // dipole-dipole interactions

    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;
        for (size_t j = i + 1; j < natom; ++j) {
            const size_t j3 = 3*j;

            double Rij[3], Rsq(0), diR(0), djR(0), didj(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = xyz[i3 + k] - xyz[j3 + k];
                Rsq += Rij[k]*Rij[k];
                diR += dipole[i3 + k]*Rij[k];
                djR += dipole[j3 + k]*Rij[k];
                didj += dipole[i3 + k]*dipole[j3 + k];
            }

            const bool ij_from_same_water = (i/4 == j/4);
            const double aDD = ij_from_same_water ?
                (i%4 == 0 ? aDD_intra_OH : aDD_intra_HH) : aDD_inter;

            double ts1, ts2, ts3;
            ttm4::smear3(std::sqrt(Rsq), AA[i%4]*AA[j%4], aDD, ts1, ts2, ts3);

            for (size_t k = 0; k < 3; ++k) {
                const double derij =
                   - 3*ts2*(didj*Rij[k] + djR*dipole[i3 + k]
                                        + diR*dipole[j3 + k])
                   + 15*ts3*diR*djR*Rij[k];

                grd4[i3 + k] += derij;
                grd4[j3 + k] -= derij;
            }
        }
    }

    // grd4 ==>> gE_ind

    for (size_t n = 0; n < nw; ++n) {

        const size_t io  = 9*n + 0;
        const size_t ih1 = 9*n + 3;
        const size_t ih2 = 9*n + 6;
        // O H H

        for (size_t k = 0; k < 9; ++k)
            gE_ind[9*n + k] = grd4[12*n + k];

        // M
        for (size_t k = 0; k < 3; ++k) {
            gE_ind[io  + k] += gamma1*grd4[12*n + 9 + k]; // O
            gE_ind[ih1 + k] += gamma2*grd4[12*n + 9 + k]; // H
            gE_ind[ih2 + k] += gamma2*grd4[12*n + 9 + k]; // H
        }

        for (size_t k = 0; k < 3; ++k) {
            gE_ind[ih1 + k] += GRDQ(0, 0, k)*phi[4*n + 1]  // phi(h1)
                             + GRDQ(0, 1, k)*phi[4*n + 2]  // phi(h2)
                             + GRDQ(0, 2, k)*phi[4*n + 3]; // phi(M)

            gE_ind[ih2 + k] += GRDQ(1, 0, k)*phi[4*n + 1]  // phi(h1)
                             + GRDQ(1, 1, k)*phi[4*n + 2]  // phi(h2)
                             + GRDQ(1, 2, k)*phi[4*n + 3]; // phi(M)

            gE_ind[io + k] += GRDQ(2, 0, k)*phi[4*n + 1]  // phi(h1)
                            + GRDQ(2, 1, k)*phi[4*n + 2]  // phi(h2)
                            + GRDQ(2, 2, k)*phi[4*n + 3]; // phi(M)
        }
    }

}

//----------------------------------------------------------------------------//

} // namespace ttm

////////////////////////////////////////////////////////////////////////////////
