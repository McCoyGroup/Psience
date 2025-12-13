import McUtils.Coordinerds as coordops
import McUtils.Numputils as nput
import McUtils.Devutils as dev
from McUtils.Data import AtomData, UnitsData
import numpy as np

from . import paths

def load_acetylene():
    acet_coords =  [[0., -0.00000001, -1.22523364],
                     [0., 0.00000002, 1.22523343],
                     [0., 0.00000003, -3.61388469],
                     [0., -0.00000005, 3.61388491]]
    # acet_zmat = coordops.functionalized_zmatrix(
    #         2,
    #         single_atoms=[0, 1]
    #     )
    specs_1 = [
        (0, 1), (0, 2), (1, 3), (1, 0, 2), (0, 1, 3), (2, 0, 1, 3)
    ]
    return acet_coords, specs_1

def load_vinylidene():
    minimum_2 = [[ 0.,          0.,  1.53899513],
                 [ 0.,          0., -0.89818466],
                 [-1.78011952,  0., -1.92243141],
                 [ 1.78011952,  0., -1.92243141]]
    # zm_2 = coordops.functionalized_zmatrix(
    #     2,
    #     ethyl_positions=[1]
    # )
    # specs_2 = coordops.extract_zmatrix_internals(zm_2)
    specs_2 = [
        (0, 1), (1, 2), (1, 3), (0, 1, 2), (0, 1, 3), (2, 0, 1, 3)
    ]

    return minimum_2, specs_2

def apply_smooth_acet_interp(minimum_1, specs_1, minimum_2, specs_2,
                             interp_type='sig'
                             ):

    ics_11 = nput.internal_coordinate_tensors(
        minimum_1,
        specs_1,
        order=0
    )[0]

    ics_21 = nput.internal_coordinate_tensors(
        minimum_2,
        specs_1,
        order=0
    )[0]

    ics_22 = nput.internal_coordinate_tensors(
        minimum_2,
        specs_2,
        order=0
    )[0]
    ics_12 = nput.internal_coordinate_tensors(
        minimum_1,
        specs_2,
        order=0
    )[0]

    specs_3 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 0, 1, 3)]
    ics_23 = nput.internal_coordinate_tensors(
        minimum_2,
        specs_3,
        order=0
    )[0]
    ics_13 = nput.internal_coordinate_tensors(
        minimum_1,
        specs_3,
        order=0
    )[0]

    # print(specs_1)
    # print(ics_11)
    # print(specs_2)
    # print(ics_21)
    # print(ics_22)

    # _, dist_conv2 = get_internal_distance_conversion(specs_2)
    # print(nput.distance_matrix(minimum_2, return_triu=True))
    # print(dist_conv2(ics_22))
    # return

    iterp_x = np.linspace(0, 1, 100)
    ic_interp_1 = ics_11[np.newaxis, :] * (1 - iterp_x[:, np.newaxis]) + ics_21[np.newaxis, :] * iterp_x[:, np.newaxis]
    ic_interp_2 = ics_12[np.newaxis, :] * (1 - iterp_x[:, np.newaxis]) + ics_22[np.newaxis, :] * iterp_x[:, np.newaxis]
    ic_interp_3 = ics_13[np.newaxis, :] * (1 - iterp_x[:, np.newaxis]) + ics_23[np.newaxis, :] * iterp_x[:, np.newaxis]

    d_ic_1 = ic_interp_1 - ics_11[np.newaxis, :]
    d_ic_2 = ic_interp_2 - ics_22[np.newaxis, :]

    interp_norms = np.array([
        np.linalg.norm(d_ic_1, axis=1),
        np.linalg.norm(d_ic_2, axis=1)
    ])
    # which_interp = np.argmin(interp_norms, axis=0)
    if interp_type == 'exp':
        percentages = np.exp(-interp_norms[0] / (interp_norms[1] + 1e-12))
    else:
        c = interp_norms[1][0] / (interp_norms[0][-1] + interp_norms[1][0])
        s = 10
        k = -s / c if c < .5 else -s / (1 - c)
        sigmoid_weight = lambda p: (
            1 / (1 + np.exp(-k * (c - p)))
        )
        percentages = 1 - sigmoid_weight(iterp_x)
    # print(percentages)
    # return
    dist_specs, conv_1 = coordops.get_internal_distance_conversion(specs_1)
    _, conv_2 = coordops.get_internal_distance_conversion(specs_2)
    _, conv_3 = coordops.get_internal_distance_conversion(specs_3)

    dists1 = conv_1(ic_interp_1)
    dists2 = conv_2(ic_interp_2)
    dists3 = conv_3(ic_interp_3)
    dists12 = percentages[:, np.newaxis] * dists1 + (1 - percentages[:, np.newaxis]) * dists2
    geoms1 = nput.points_from_distance_matrix(dists1, use_triu=True, target_dim=3)
    geoms2 = nput.points_from_distance_matrix(dists2, use_triu=True, target_dim=3)
    geoms3 = nput.points_from_distance_matrix(dists3, use_triu=True, target_dim=3)
    geoms12 = nput.points_from_distance_matrix(dists12, use_triu=True, target_dim=3)

    d_1 = nput.internal_coordinate_tensors(geoms1, dist_specs, order=1)[1:]
    d_2 = nput.internal_coordinate_tensors(geoms2, dist_specs, order=1)[1:]
    d_13 = nput.internal_coordinate_tensors(geoms1, specs_3, order=1)[1:]
    d_23 = nput.internal_coordinate_tensors(geoms2, specs_3, order=1)[1:]
    d_3 = nput.internal_coordinate_tensors(geoms3, specs_3, order=1)[1:]
    # d_12 = percentages[:, np.newaxis, np.newaxis] * d_1 + (1 - percentages[:, np.newaxis, np.newaxis]) * d_2
    d_12 = nput.internal_coordinate_tensors(geoms12, dist_specs, order=1)[1:]
    d_123 = nput.internal_coordinate_tensors(geoms12, specs_3, order=1)[1:]
    # d_121 = nput.internal_coordinate_tensors(geoms12, specs_1, order=1)[1:]
    # d_122 = nput.internal_coordinate_tensors(geoms12, specs_2, order=1)[1:]
    m_h = np.array([AtomData[a, "Mass"] for a in ["C", "C", "H", "H"]]) * UnitsData.convert("AtomicMassUnits",
                                                                                            "ElectronMass")
    m_d = np.array([AtomData[a, "Mass"] for a in ["C", "C", "D", "H"]]) * UnitsData.convert("AtomicMassUnits",
                                                                                            "ElectronMass")
    # print(m_h, m_d)
    # g121_h = nput.metric_tensor(d_121, m_h)
    # g122_h = nput.metric_tensor(d_122, m_h)
    g12_h = nput.metric_tensor(d_12, m_h)
    # g121_d = nput.metric_tensor(d_121, m_d)
    # g122_d = nput.metric_tensor(d_122, m_d)
    g12_d = nput.metric_tensor(d_12, m_d)
    g123_d = nput.metric_tensor(d_123, m_d)
    g123_h = nput.metric_tensor(d_123, m_h)
    g2_d = nput.metric_tensor(d_2, m_d)
    g2_3_d = nput.metric_tensor(d_23, m_d)
    g1_d = nput.metric_tensor(d_1, m_d)
    g1_3_d = nput.metric_tensor(d_13, m_d)
    g3_d = nput.metric_tensor(d_3, m_d)
    g2_h = nput.metric_tensor(d_2, m_h)
    g2_3_h = nput.metric_tensor(d_23, m_h)
    g1_3_h = nput.metric_tensor(d_13, m_h)
    g3_h = nput.metric_tensor(d_3, m_h)
    g1_h = nput.metric_tensor(d_1, m_h)

    # print()
    # print(ics_12)
    # print(ics_22)
    # print(
    #     nput.internal_coordinate_tensors(
    #         geoms1,
    #         specs_2,
    #         order=0
    #     )[0]
    # )
    # geoms12 = np.concatenate([
    #         geoms1[which_interp == 0],
    #         geoms2[which_interp == 1]
    #     ], axis=0)
    # geoms = conv_1(ic_interp_1)
    return {
        "smooth": geoms12.tolist(),
        "acet": geoms1.tolist(),
        "vinny": geoms2.tolist(),
        "merge": geoms3.tolist(),
        "perc": percentages.tolist(),
        "g12": g12_h.tolist(),
        "g12_CCHD": g12_d.tolist(),
        "g123_CCHD": g123_d.tolist(),
        "g123": g123_h.tolist(),
        # "g121":g121_h.tolist(),
        # "g121_CCHD":g121_d.tolist(),
        # "g122":g122_h.tolist(),
        # "g122_CCHD":g122_d.tolist(),
        "g1_CCHD": g1_d.tolist(),
        "g1": g1_h.tolist(),
        "g1_3_CCHD": g1_3_d.tolist(),
        "g1_3": g1_3_h.tolist(),
        # "g11_CCHD":g11_d.tolist(),
        "g2_CCHD": g2_d.tolist(),
        "g2_3_CCHD": g2_3_d.tolist(),
        "g3_CCHD": g3_d.tolist(),
        "g2": g2_h.tolist(),
        "g2_3": g2_3_h.tolist(),
        "g3": g3_h.tolist(),
        # "g22_CCHD":g22_d.tolist()
    }

def write_acet_interp_file(acet_struct=None, vin_struct=None,
                           filename='acet_interp_structs.json',
                           interp_type='sig'
                           ):
    min_1, specs_1 = load_acetylene()
    min_2, specs_2 = load_vinylidene()
    if acet_struct is None:
        acet_struct = min_1
    if vin_struct is None:
        vin_struct = min_2
    output = paths.path(filename)
    data = apply_smooth_acet_interp(
        acet_struct, specs_1,
        vin_struct, specs_2,
        interp_type=interp_type
    )
    dev.write_json(output, data)
    return output

def load_acet_interp_file(filename='acet_interp_structs.json'):
    return dev.read_json(paths.path(filename))

def get_emb_geom_energies(evaluator, emb_geom_data, key, cache=None):
    if cache is None:
        cache = {}
    if key not in cache:
        cache[key] = evaluator.calculate_energy(coords=emb_geom_data[key])
    return cache[key]
