

//#include "Python.h"
#include "plzffi/FFIModule.hpp"
//#include "plzffi/FFIDynamicLibrary.hpp"
#include <numeric>
#include <string>
#include <vector>

using namespace plzffi;

// Put all module functions in here
// Declare our linkage
extern "C" {
    double calcpot_(int*, double*, const double*);// comes from libmbpol.so
    double calcpotg_(int *nw, double *Vpot, const double *x, double *g);// comes from libmbpol.so
}

namespace LegacyMBPol {

    // all functions should take a single argument of type "Arguments &"
    // data can be extracted by using the `.value` method with the name of the parameter
    double mbpol(Arguments &params) {

        auto nwaters = params.value<int>("nwaters");
        auto coords = params.value<double*>("coords");

        double pot_val;
        calcpot_(&nwaters, &pot_val, coords);
        return pot_val / 627.5094740631;

    }

    double test_pot(Arguments &params) {

        return 0.0;

    }

    double test_val(Arguments &params) {

        auto val = params.value<double>("val");
        return val;

    }

//    unsigned char test_uchar(Arguments &params) {
//
//        auto val = params.value<unsigned char>("val");
//        return val;
//
//    }

    std::vector<double> mbpol_vec(Arguments &params) {

        auto nwaters = params.value<int>("nwaters");
        auto coords = params.value<double*>("coords");
        auto coords_shape = params.shape("coords"); // we can infer this...

        std::vector<double> energies(coords_shape[0]);
        auto block_size = std::accumulate(coords_shape.begin()+1, coords_shape.end(), 1, std::multiplies<>());

        double pot_val;
        for (size_t w = 0; w < coords_shape[0]; w++) {
            calcpot_(&nwaters, &pot_val, coords + (block_size*w));
            energies[w] = pot_val / 627.5094740631;
        }

//        printf("    > sucessfully got energy %f %f %f!\n", energies[0], energies[1], energies[2]);

        return energies;
    }

    FFICompoundType energy_grad_type {
            {"energy", "grad"},
            {FFIType::Double, FFIType::Double},
            {{}, {0, 3}}
    };

//    FFICompoundType energy_grad_type {
//            {"energy", "grad"},
//            {FFIType::Double, FFIType::Double},
//            {{}, {}}
//    };

//    FFICompoundType energy_grad_type {
//            {"energy"},//, "grad"},
//            {FFIType::Double}//, FFIType::Double},
//    };

//    FFICompoundType energy_grad_type {
//    };

    FFICompoundReturn mbpol_grad(Arguments &params) {

        FFICompoundReturn res(energy_grad_type);


//        printf("...0?\n");

        auto nwaters = params.value<int>("nwaters");

//        printf(".10?\n");

        auto coords = params.value<double*>("coords");

//        printf("...1?\n");

        std::vector<double> grad(nwaters*9);
        double pot_val;

//        printf("...2?\n");

        calcpotg_(&nwaters, &pot_val, coords, grad.data());
        pot_val = pot_val / 627.5094740631;
        for (size_t i = 0; i < grad.size(); i++) {
            grad[i] = grad[i] / 627.5094740631; // Convert to Hartree
        }

//        printf("...3?\n");

        res.set<double>("energy", pot_val);
        res.set<double>("grad", std::move(grad), {(size_t)nwaters, 3, 3});

        return res;

    }

    FFICompoundReturn mbpol_grad_vec(Arguments &params) {

        FFICompoundReturn res(energy_grad_type);

        auto nwaters = params.value<int>("nwaters");
        auto coords = params.value<double*>("coords");
        auto coords_shape = params.shape("coords");

        std::vector<double> energies(coords_shape[0]);
        std::vector<double> grad(coords_shape[0]*nwaters*9);

        auto block_size = std::accumulate(coords_shape.begin()+1, coords_shape.end(), 1, std::multiplies<>());
        size_t grad_size = nwaters*9;
        size_t grad_offset = 0;
        for (size_t w = 0; w < coords_shape[0]; w++) {
            grad_offset = w * grad_size;
            calcpotg_(&nwaters, energies.data()+w, coords +(block_size*w), grad.data()+grad_offset);
            energies[w] /= 627.5094740631;
            for (size_t i = 0; i < grad_size; i++) {
                grad[grad_offset+i] /= 627.5094740631; // Convert to Hartree
            }
        }

        res.set<double>("energy", std::move(energies));
        res.set<double>("grad", std::move(grad), {coords_shape[0], (size_t)nwaters, 3, 3});

//        printf("    > sucessfully got energy %f %f %f!\n", energies[0], energies[1], energies[2]);

        return res;
    }

    bool mbpol_grad_vec_buffered(Arguments &params) {

        auto nwaters = params.value<int>("nwaters");
        auto coords = params.value<double*>("coords");
        auto coords_shape = params.shape("coords");

        auto energies = params.value<double*>("energies");
        auto grad = params.value<double*>("gradients");

//        std::vector<double> energies(coords_shape[0]);
//        std::vector<double> grad(coords_shape[0]*nwaters*9);

        auto block_size = std::accumulate(coords_shape.begin()+1, coords_shape.end(), 1, std::multiplies<>());
        size_t grad_size = nwaters*9;
        size_t grad_offset = 0;
        for (size_t w = 0; w < coords_shape[0]; w++) {
            grad_offset = w * grad_size;
            calcpotg_(&nwaters, energies+w, coords + (block_size*w), grad+grad_offset);
            energies[w] /= 627.5094740631;
            for (size_t i = 0; i < grad_size; i++) {
                grad[grad_offset+i] /= 627.5094740631; // Convert to Hartree
            }
        }

        return true;
    }

        // need a load function that can be called in PYMODINIT
    void load(FFIModule *mod) {
        // load modules and return python def

//        // add data for test pots
//        mod->add<unsigned char>(
//                "test_uchar",
//                {
//                },
//                test_uchar
//                );

        // add data for test pots
        mod->add<double>(
                "test_pot",
                {
                },
                test_pot
                );

        // add data for test pots
        mod->add<double>(
                "test_val",
                {
                    {"val", FFIType::Double, {}},
                },
                test_val
                );

        // add data for first obj
        mod->add<double>(
                "get_pot",
                {
                        {"nwaters", FFIType::Int, {}},
                        {"coords", FFIType::Double, {0, 3, 3}},
                },
                mbpol
                );

        // add data for version with gradient
        mod->add(
                "get_pot_grad",
                {
                        {"nwaters", FFIType::Int, {}},
                        {"coords", FFIType::Double, {0, 3, 3}},
                },
                energy_grad_type,
                mbpol_grad
        );

        // add data for vectorized version
        mod->add<double>(
                "get_pot_vec",
                {
                        {"nwaters", FFIType::Int, {}},
                        {"coords", FFIType::Double, {0, 0, 3, 3}},
                },
                mbpol_vec
        );

        // add data for vectorized version with gradient
        mod->add(
                "get_pot_grad_vec",
                {
                        {"nwaters", FFIType::Int, {}},
                        {"coords", FFIType::Double, {0, 0, 3, 3}},
                },
                energy_grad_type,
                mbpol_grad_vec
        );

        // add data for buffered/vectorized version with gradient
        mod->add<bool>(
                "get_pot_grad_vec_buffered",
                {
                        {"nwaters", FFIType::Int, {}},
                        {"coords", FFIType::Double, {0, 0, 3, 3}},
                        {"energies", FFIType::Double, {0}},
                        {"gradients", FFIType::Double, {0, 0, 3, 3}},
                },
                mbpol_grad_vec_buffered
        );

    }

    static FFIModule Data(
        "LegacyMBPol",
        "provides linkage for legacy version of MB-Pol",
        load
        );
}

PyMODINIT_FUNC PyInit_LegacyMBPol(void) {
    return LegacyMBPol::Data.attach();
}