#ifndef MBPOL_H
#define MBPOL_H

#include "ttm4-es.h"

namespace x2o {

struct mbpol {

    double operator()(size_t nw, const double*) const;
    double operator()(size_t nw, const double*, double*) const;

protected:
    mutable ttm::ttm4_es m_ttm4;
};

} // namespace x2o

#endif // MBPOL_H
