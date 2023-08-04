#ifndef TTM4_SMEAR_H
#define TTM4_SMEAR_H

//
// smeared interactions for TTM4
//

namespace ttm4 {

void smear01(const double& r12, const double& AA,
             const double& a, double& ts0, double& ts1);

void smear2(const double& r12, const double& AA,
            const double& a, double& ts1, double& ts2);

void smear3(const double& r12, const double& AA,
            const double& a, double& ts1, double& ts2, double& ts3);

} // namespace ttm4

#endif // TTM4_SMEAR_H
