#ifndef TTM4_ES_H
#define TTM4_ES_H

#include <cstddef>

namespace ttm {

//
// class incapsulating TTM4 electrostatics
//

struct ttm4_es {

    // computed for vibrationally averaged (CC-pol) monomer

    static const double molecular_polarizability; // in A^3
    static const double molecular_dipole; // in D

    ttm4_es();
    ~ttm4_es();

    void operator()(size_t nw,
                    const double*, // O H H O H H
                    double& E_elec, double* gE_elec,
                    double& E_ind,  double* gE_ind);

private:
    size_t  m_nw;
    double* m_mem;

    void allocate(size_t);
};

} // namespace ttm

#endif // TTM4_ES_H
