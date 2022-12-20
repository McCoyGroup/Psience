#ifndef X2B_DISPERSION_H
#define X2B_DISPERSION_H

namespace x2o {

//
//  TT6*C6/R^6
//

struct x2b_dispersion {

    static const double C6_HH;
    static const double C6_OH;
    static const double C6_OO;

    static const double d6_HH;
    static const double d6_OH;
    static const double d6_OO;

    static double eval(const double*, const double*);

    // including gradients

    static double eval(const double*, const double*, double*, double*);
};

} // namespace x2o

#endif // X2B_DISPERSION_H
