#include <cmath>

#include "tang-toennies.h"
#include "x2b-dispersion.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

inline double x6(const double& C6, const double& d6,
                 const double* p1, const double* p2)
{
    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    const double dz = p1[2] - p2[2];

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);

    const double tt6 = x2o::tang_toennies(6, d6*r);

    const double inv_rsq = 1.0/rsq;
    const double inv_r6 = inv_rsq*inv_rsq*inv_rsq;

    return - C6*tt6*inv_r6;
}

//----------------------------------------------------------------------------//

const double if6 = 1.0/x2o::factorial<6>();

//----------------------------------------------------------------------------//

inline double x6(const double& C6, const double& d6,
                 const double* p1, const double* p2,
                       double* g1,       double* g2)
{
    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    const double dz = p1[2] - p2[2];

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);

    const double d6r = d6*r;
    const double tt6 = x2o::tang_toennies(6, d6r);

    const double inv_rsq = 1.0/rsq;
    const double inv_r6 = inv_rsq*inv_rsq*inv_rsq;

    const double e6 = C6*tt6*inv_r6;

    const double grd = 6*e6*inv_rsq - C6*std::pow(d6, 7)*if6*std::exp(-d6r)/r;

    g1[0] += dx*grd;
    g2[0] -= dx*grd;

    g1[1] += dy*grd;
    g2[1] -= dy*grd;

    g1[2] += dz*grd;
    g2[2] -= dz*grd;

    return - e6;
}

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace x2o {

////////////////////////////////////////////////////////////////////////////////

const double x2b_dispersion::C6_HH = 2.009358600184719e+01; // kcal/mol * A^(-6)
const double x2b_dispersion::C6_OH = 8.349556669872743e+01; // kcal/mol * A^(-6)
const double x2b_dispersion::C6_OO = 2.373212214147944e+02; // kcal/mol * A^(-6)

const double x2b_dispersion::d6_HH =  9.406475169954112e+00; // A^(-1)
const double x2b_dispersion::d6_OH =  9.775202425217957e+00; // A^(-1)
const double x2b_dispersion::d6_OO =  9.295485815062264e+00; // A^(-1)

//----------------------------------------------------------------------------//

double x2b_dispersion::eval(const double* w1, const double* w2)
{
    const double* Oa  = w1;
    const double* Ha1 = w1 + 3;
    const double* Ha2 = w1 + 6;

    const double* Ob  = w2;
    const double* Hb1 = w2 + 3;
    const double* Hb2 = w2 + 6;

    const double HH6 =
        x6(C6_HH, d6_HH, Ha1, Hb1)
      + x6(C6_HH, d6_HH, Ha1, Hb2)
      + x6(C6_HH, d6_HH, Ha2, Hb1)
      + x6(C6_HH, d6_HH, Ha2, Hb2);

    const double OH6 =
        x6(C6_OH, d6_OH, Oa, Hb1)
      + x6(C6_OH, d6_OH, Oa, Hb2)
      + x6(C6_OH, d6_OH, Ob, Ha1)
      + x6(C6_OH, d6_OH, Ob, Ha2);

    const double OO6 =
        x6(C6_OO, d6_OO, Oa, Ob);

    return HH6 + OH6 + OO6;
}

//----------------------------------------------------------------------------//

double x2b_dispersion::eval
    (const double* w1, const double* w2, double* g1, double* g2)
{
    const double* Oa  = w1;
    const double* Ha1 = w1 + 3;
    const double* Ha2 = w1 + 6;

    const double* Ob  = w2;
    const double* Hb1 = w2 + 3;
    const double* Hb2 = w2 + 6;

    double* gOa  = g1;
    double* gHa1 = g1 + 3;
    double* gHa2 = g1 + 6;

    double* gOb  = g2;
    double* gHb1 = g2 + 3;
    double* gHb2 = g2 + 6;

    const double HH6 =
        x6(C6_HH, d6_HH, Ha1, Hb1, gHa1, gHb1)
      + x6(C6_HH, d6_HH, Ha1, Hb2, gHa1, gHb2)
      + x6(C6_HH, d6_HH, Ha2, Hb1, gHa2, gHb1)
      + x6(C6_HH, d6_HH, Ha2, Hb2, gHa2, gHb2);

    const double OH6 =
        x6(C6_OH, d6_OH, Oa, Hb1, gOa, gHb1)
      + x6(C6_OH, d6_OH, Oa, Hb2, gOa, gHb2)
      + x6(C6_OH, d6_OH, Ob, Ha1, gOb, gHa1)
      + x6(C6_OH, d6_OH, Ob, Ha2, gOb, gHa2);

    const double OO6 =
        x6(C6_OO, d6_OO, Oa, Ob, gOa, gOb);

    return HH6 + OH6 + OO6;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace x2o

////////////////////////////////////////////////////////////////////////////////
