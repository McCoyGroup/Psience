#include <cmath>
#include <cassert>

#include <iostream>

#include "gammq.h"
#include "macros.h"

////////////////////////////////////////////////////////////////////////////////

namespace ttm4 {

//----------------------------------------------------------------------------//

void smear01(const double& RESTRICT r12,
             const double& RESTRICT AA, // (polarA*polarB)^(1/6)
             const double& RESTRICT a,
             double& RESTRICT ts0,
             double& RESTRICT ts1)
{
    if (!(r12 > 0.0)) {
        std::cerr << "smear01: r12 = " << r12 << std::endl;
    }

    assert(r12 > 0.0);

    const double dri = 1.0/r12;
    const double drsqi = dri*dri;

    const double g34 = std::exp(ttm::gammln(3.0/4.0));

    const double rA = r12/AA;
    const double rA4 = std::pow(rA, 4);
    const double exp1 = std::exp(-a*rA4);
    const double a_mrt = std::pow(a, 1.0/4.0);

    ts0 = (1.0 - exp1 + a_mrt*rA*g34*ttm::gammq(3.0/4.0, a*rA4))*dri;
    ts1 = (1.0 - exp1)*dri*drsqi;
}

//----------------------------------------------------------------------------//

void smear2(const double& RESTRICT r12,
            const double& RESTRICT AA, // (polarA*polarB)^(1/6)
            const double& RESTRICT a,
            double& RESTRICT ts1,
            double& RESTRICT ts2)
{
    const double dri = 1.0/r12;
    const double drsqi = dri*dri;

    const double rA = r12/AA;
    const double rA4 = std::pow(rA, 4);
    const double exp1 = std::exp(-a*rA4);

    ts1 = (1.0 - exp1)*dri*drsqi;
    ts2 = (ts1 - (4.0/3.0)*a*exp1*rA4*dri*drsqi)*drsqi;
}

//----------------------------------------------------------------------------//

void smear3(const double& RESTRICT r12,
            const double& RESTRICT AA, // (polarA*polarB)^(1/6)
            const double& RESTRICT a,
            double& RESTRICT ts1,
            double& RESTRICT ts2,
            double& RESTRICT ts3)
{
    const double dri = 1.0/r12;
    const double drsqi = dri*dri;

    const double AA4 = std::pow(AA, 4);
    const double rA = r12/AA;
    const double rA4 = std::pow(rA, 4);
    const double exp1 = std::exp(-a*rA4);

    ts1 = (1.0 - exp1)*dri*drsqi;
    ts2 = (ts1 - (4.0/3.0)*a*exp1*rA4*dri*drsqi)*drsqi;
    ts3 = (ts2 - (4.0/15.0)*dri*a*(4*a*rA4 - 1.0)*exp1/AA4)*drsqi;
}

//----------------------------------------------------------------------------//

} // namespace ttm4

////////////////////////////////////////////////////////////////////////////////
