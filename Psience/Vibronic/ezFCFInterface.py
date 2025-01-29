import os
import uuid
import numpy as np
import weakref
import tempfile as tf
import subprocess
import io
from McUtils.Jupyter import JHTML
from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput
import McUtils.Iterators as itut
from McUtils.ExternalPrograms import ExternalProgramRunner

__all__ = [
    "ezXML",
    "ezFCFInterface"
]

class ezXML:
    class ezTag(JHTML.XML.DeclarativeElement):
        all_remapped_attrs = weakref.WeakKeyDictionary()
        def __init__(self, *elems, flag=None, **attrs):
            preserved_attrs = {}
            if flag is not None:
                preserved_attrs["flag"] = self.sanitize_value(flag)
            remapped_attrs = {
                "attr-"+str(uuid.uuid4()): av
                for av in attrs.items()
            }
            self.all_remapped_attrs[self] = {k:v[0] for k,v in remapped_attrs.items()}
            super().__init__(
                *elems,
                **preserved_attrs,
                **{
                    k:self.sanitize_value(v[1]) for k,v in remapped_attrs.items()
                }
            )
        @classmethod
        def sanitize_value(cls, val):
            if isinstance(val, int):
                return str(val).lower()
            else:
                return str(val)
        def tostring(self, renormalize_padding=True):
            base = super().tostring(short_empty_elements=False, prettify=True)
            base = base.replace(
                "&#10;", "\n"
            )
            indent_marker = '---indent-marker---'
            pad = "\n" + indent_marker if renormalize_padding else ""
            for remapping in self.all_remapped_attrs.values():
                for k,v in remapping.items():
                    base = base.replace(k, pad+v)
            this_pad = pad if renormalize_padding else "\n"
            base = base.replace("></", ">" + this_pad + "</")
            # base = base.replace("\" />", "\"" + this_pad + "/>")
            if renormalize_padding:
                new_lines = []
                lines = base.splitlines()
                indent = ""
                for l in lines:
                    if l.strip().startswith("<") and not l.strip().startswith("</"):
                        indent = l[:l.index("<")] + "  "
                    elif l.startswith(indent_marker):
                        l = indent + l[len(indent_marker):]
                    new_lines.append(l)
                base = "\n".join(new_lines)
            return base

    class job_parameters(ezTag):  known_opts = {"temperature", "spectrum_intensity_threshold"}

    class do_not_excite_subspace(ezTag):  known_opts = {"size", "normal_modes"}

    class initial_state(ezTag):  known_opts = {}

    class target_state(ezTag):  known_opts = {}

    class the_only_initial_state(ezTag):  known_opts = {}

    class energy_thresholds(ezTag):  known_opts = {"units"}

    class print_franck_condon_matrices(ezTag):  known_opts = {"flag"}
    class print_normal_modes(ezTag):  known_opts = {"flag"}

    class parallel_approximation(ezTag):  known_opts = {"max_vibr_excitations_in_initial_el_state",
                                                        "max_vibr_excitations_in_target_el_state", "combination_bands",
                                                        "use_normal_coordinates_of_target_states"}

    class max_vibr_to_store(ezTag):  known_opts = {"target_el_state"}

    class single_excitation(ezTag):  known_opts = {"ini", "targ"}

    class dushinsky_rotations(ezTag):  known_opts = {"target_state", "max_vibr_excitations_in_initial_el_state",
                                                     "max_vibr_excitations_in_target_el_state"}

    class geometry(ezTag):
        known_opts = {"number_of_atoms", "linear", "units", "text"}
        def __init__(self,
                     *,
                     number_of_atoms,
                     text,
                     linear="false",
                     units="angstr",
                     **opts
        ):
            super().__init__(
                number_of_atoms=number_of_atoms,
                linear=linear,
                units=units,
                **opts,
                text=text
            )

    class manual_atoms_reordering(ezTag):  known_opts = {"new_order"}

    class manual_coordinates_transformation(ezTag):
        known_opts = {
            "shift_along_x", "shift_along_y", "shift_along_z",
            "rotate_around_x", "rotate_around_y", "rotate_around_z"
        }

    class normal_modes(ezTag):
        known_opts = {"if_mass_weighted", "text", "atoms"}
        def __init__(self,
                     *,
                     atoms,
                     if_mass_weighted,
                     text,
                     **opts
        ):
            super().__init__(
                atoms=atoms,
                if_mass_weighted=if_mass_weighted,
                **opts,
                text=text
            )

    class frequencies(ezTag):  known_opts = {"text"}

    class excitation_energy(ezTag):  known_opts = {"units"}

    class manual_normal_modes_reordering(ezTag):  known_opts = {"new_order"}

    class input(ezTag):  known_opts = {"job"}


class ezFCFInterface:
    def __init__(self, atoms, gs_nms, es_nms, excitations,
                 masses=None,
                 ground_states=None,
                 include_rotation=True, rotation_order='gs',
                 rotation_blocks=None,
                 rotation_center=None, logger=None, # ignored
                 mode_reordering=None,
                 rotation_method='duschinsky',
                 mass_weight=False,
                 dimensionless=False,
                 always_run_parallel=True,
                 print_all=True,
                 embed=False):
        self.atoms, self.atom_map = self.prep_masses(atoms, masses)
        self.excitation = excitations
        self.ground_states = ground_states
        if embed:
            raise NotImplementedError("need to add embedding to this branch too...")
        if rotation_method is not None and rotation_method != 'duschinsky':
            raise ValueError(f"ezFCF only supports Duschinsky rotations (got {rotation_method})")
        if dimensionless:
            raise ValueError("ezFCF doesn't support dimensionless normal modes")
        if mass_weight:
            gs_nms = gs_nms.make_mass_weighted()
            es_nms = es_nms.make_mass_weighted()
        self.gs_nms = gs_nms
        self.es_nms = es_nms
        self.include_rotation = include_rotation
        self.rotation_order = rotation_order
        self.rotation_blocks = rotation_blocks
        self.mode_reordering = mode_reordering
        self.always_run_parallel = always_run_parallel
        self.print_all = print_all

    @classmethod
    def _check_listable(cls, excitations):
        try:
            return all(
                isinstance(exc, np.ndarray) or all(
                    nput.is_numeric(e) or all(nput.is_numeric(ee) for ee in e)
                    for e in exc
                )
                for exc in excitations
            )
        except TypeError:
            return False


    def format(self, job_type="harmonic_pes", temperature=0, spectrum_intensity_threshold=1e-8):
        els = [
        ]
        if self.print_all:
            els.append(ezXML.print_normal_modes(flag="true"))
            els.append(ezXML.print_normal_modes(flag="true"))

        els.append(
            ezXML.job_parameters(
                temperature=temperature,
                spectrum_intensity_threshold=spectrum_intensity_threshold
            )
        )

        excitation_options, excitation_elems = self.prep_excitations(self.es_nms, self.ground_states, self.excitation)
        include_rotation = self.include_rotation
        if not include_rotation or self.always_run_parallel:
            els.append(
                self.prep_parallel(
                    *excitation_elems,
                    **excitation_options
                )
            )
        if include_rotation:
            els.append(
                self.prep_duschinsky(
                    *excitation_elems,
                    **excitation_options
                )
            )

        mode_reordering = self.mode_reordering
        return ezXML.input(
            *els,
            self.prep_initial_state(self.atoms, self.gs_nms),
            self.prep_target_state(self.atoms, self.es_nms, mode_reordering=mode_reordering),
            job=job_type
        )

    def prep_masses(self, atoms, masses, units=None):
        if masses is None:
            masses = [AtomData[a]["Mass"] for a in atoms]
        if units is None:
            if masses[0] > 500:
                units = "ElectronMass"
            else:
                units = "AtomicMassUnits"
        masses = np.array(masses) * UnitsData.convert(units, "AtomicMassUnits")

        if atoms is None:
            atoms = [""] * len(masses)
            rem = np.arange(len(masses))
            for key, val in AtomData.data.items():
                diffs = val["Mass"] - masses[rem,]
                matches = np.where(abs(diffs) < .05)
                if len(matches) > 0 and len(matches[0]) > 0:
                    match_pos = rem[matches]
                    for i in match_pos:
                        atoms[i] = val["Symbol"]
                    rem = np.delete(rem, match_pos[0])
                    if len(rem) == 0:
                        break
            else:
                for i in rem: atoms[i] = "X"

        atom_map = {}
        new_atoms = []
        ndig = int(1 + np.floor(np.log10(len(atoms))))
        for i,(atom, mass) in enumerate(zip(atoms, masses)):
            if abs(mass - atom_map.get(atom, mass)) > 1e-4:
                atom = atom + f"00{i:0>{ndig}d}"
                while atom in atom_map:
                    atom = atom + f"0{i}"
            atom_map[atom] = mass
            new_atoms.append(atom)
        return new_atoms, atom_map

    def format_masses_file(self, atom_map):
        return "<masses>\n" + "\n".join(f"  <{at}> {mass} </{at}>" for at, mass in atom_map.items()) + "\n</masses>"

    class _write_dir:
        def __init__(self, dir=None, dir_prefix=None, dir_suffix=None, delete=True):
            self.dir = dir
            self._temp_dirs = []
            self.delete = delete
            self.opts = {'prefix':dir_prefix, 'suffix':dir_suffix}
        def __enter__(self):
            if self.dir is None:
                td = tf.TemporaryDirectory(**self.opts)
                self._temp_dirs.append(td)
                return td.__enter__()
            else:
                return self.dir
        def __exit__(self, exc_type, exc_val, exc_tb):
            if self.dir is None:
                td = self._temp_dirs.pop()
                td.__exit__(exc_type, exc_val, exc_tb)

    class ezFCFRunner(ExternalProgramRunner):
        default_opts = dict(prefix='ezFCF-', suffix='.xml')
        def __init__(self, job, binary, **opts):
            self.job = job
            super().__init__(binary, **opts)
        def prep_dir(self, dir):
            with open(os.path.join(dir, 'atomicMasses.xml'), 'w+') as mass_file:
                mass_file.write(self.job.format_masses_file(self.job.atom_map) + "\n\n")

        text_file_extensions = ExternalProgramRunner.text_file_extensions + ['.spectrum_parallel', '.spectrum_duschinsky']
    def run(self, ezFCF_binary,
            dir=None, dir_prefix=None, dir_suffix=None, mode='w',
            prefix='ezFCF-', suffix='.xml', delete=True,
            raise_errors=True,
            **job_opts):
        return self.ezFCFRunner(
            self,
            ezFCF_binary,
            dir=dir, dir_prefix=dir_prefix, dir_suffix=dir_suffix, mode=mode,
            prefix=prefix, suffix=suffix, delete=delete,
            raise_errors=raise_errors
        ).run(
            self.format(**job_opts).tostring() + "\n\n",
            return_auxiliary_files={"parallel": "{name}.spectrum_parallel", "duschinsky": '{name}.spectrum_dushinsky'}
        )

        # results = {}
        # with self._write_dir(dir=dir, dir_prefix=dir_prefix, dir_suffix=dir_suffix) as dir:
        #     with tf.NamedTemporaryFile(dir=dir, mode=mode, prefix=prefix, suffix=suffix, delete=False) as inp:
        #         results['input_file'] = self.format(**job_opts).tostring() + "\n\n"
        #         inp.write(results['input_file'])
        #     try:
        #         if os.path.isfile(ezFCF_binary):
        #             ezFCF_binary = os.path.abspath(ezFCF_binary)
        #         res = subprocess.run([ezFCF_binary, inp.name], capture_output=True, cwd=dir)
        #         results['process_output'] = res
        #         for key,ext in {"parallel":".spectrum_parallel", "duschinsky":'.spectrum_dushinsky'}.items():
        #             test = inp.name + ext
        #             if os.path.isfile(test):
        #                 with open(test) as strm:
        #                     results[key] = strm.read()
        #         err = res.stderr.decode().strip()
        #         if raise_errors and len(err) > 0:
        #             raise IOError(err)
        #     finally:
        #         if delete:
        #             try:
        #                 os.remove(inp.name)
        #             except OSError:
        #                 pass
        #
        # if not delete:
        #     return inp, results
        # else:
        #     return results
        # # if dir is not None:

    _exc_aliases = {
        'max_quanta':'max_vibr_excitations_in_target_el_state'
    }
    _gsstate_aliases = {
        'max_quanta':'max_vibr_excitations_in_initial_el_state'
    }
    @classmethod
    def canonicalize_excitation_options(cls, nms, threshold=None, fixed_modes=None, ground_states=None, **opts):
        base_opts = {cls._exc_aliases.get(k, k):v for k,v in opts.items()}
        elems = []

        gs_threshold = None
        if ground_states is not None:
            if isinstance(ground_states, dict):
                if 'threshold' in ground_states:
                    ground_states = ground_states.copy()
                    gs_threshold = ground_states['threshold']
                    del ground_states['threshold']
                base_opts.update(
                    {
                        cls._gsstate_aliases.get(k, k): v for k, v in ground_states.items()
                        if k not in {'threshold'}
                    }
                )
            else:
                elems.append(
                    ezXML.the_only_initial_state(cls.prep_state_str(nms, ground_states))
                )

        if fixed_modes is not None:
            if nput.is_numeric(fixed_modes): fixed_modes = [fixed_modes]
            elems.append(
                ezXML.do_not_excite_subspace(
                    size=len(fixed_modes),
                    norm_modes=" ".join(str(n) for n in fixed_modes)
                )
            )
        if threshold is not None or gs_threshold is not None:
            subelems = []
            if gs_threshold is not None:
                subelems.append(
                    ezXML.initial_state(
                        gs_threshold * UnitsData.convert("Hartrees", "Wavenumbers"),
                        units='cm-1'
                    )
                )
            if threshold is not None:
                subelems.append(
                    ezXML.target_state(
                        threshold * UnitsData.convert("Hartrees", "Wavenumbers"),
                        units='cm-1'
                    )
                )
            elems.append(
                ezXML.energy_thresholds(*subelems)
            )
        return base_opts, elems


    @classmethod
    def prep_excitations(cls, nms, ground_states, excitations):
        """
        Dispatcher to get appropriate state spaces
        :param excitations:
        :param check:
        :return:
        """

        opts = {}
        elems = []
        states = []

        try:
            exc0 = excitations[0][0]
        except (IndexError, KeyError, TypeError):
            refilter = True
        else:
            refilter = (
                    not isinstance(exc0, int)
                    or any(len(e) < len(nms.freqs) for e in excitations)
            )

        if refilter:
            if isinstance(excitations, dict):
                new_opts, new_elems = cls.canonicalize_excitation_options(nms, **excitations)
                opts.update(new_opts)
                elems.extend(new_elems)
            elif isinstance(excitations, (int, np.integer)) or isinstance(excitations[0], (int, np.integer)):
                opts.update({'max_vibr_excitations_in_target_el_state': excitations})
            elif cls._check_listable(excitations):
                from ..BasisReps import StateMaker
                state = StateMaker(len(nms.freqs), mode='normal')
                states.extend(state(*exc) for exc in states)
            else:
                raise ValueError("can't build a state space from {}".format(excitations))
        else:
            states.extend(excitations)

        gs_states = []
        if ground_states is not None:
            # gs_states = []

            try:
                exc0 = excitations[0][0]
            except (IndexError, KeyError, TypeError):
                refilter = True
            else:
                refilter = (
                        not isinstance(exc0, int)
                        or any(len(e) < nms.freqs for e in excitations)
                )

            if refilter:
                if isinstance(excitations, dict):
                    new_opts, new_elems = cls.canonicalize_excitation_options(excitations)
                    opts.update(new_opts)
                    elems.extend(new_elems)
                elif isinstance(excitations, (int, np.integer)) or isinstance(excitations[0], (int, np.integer)):
                    opts.update({'max_vibr_excitations_in_target_el_state': excitations})
                elif cls._check_listable(excitations):
                    from ..BasisReps import StateMaker
                    state = StateMaker(len(nms.freqs), mode='normal')
                    gs_states.extend(state(*exc) for exc in states)
                else:
                    raise ValueError("can't build a state space from {}".format(excitations))
            else:
                gs_states.extend(excitations)

            if len(gs_states) == 1:
                elems.append(
                    ezXML.the_only_initial_state(cls.prep_state_str(nms, ground_states))
                )
            # elif len(gs_states) > 2:
            #     raise ValueError("ezFCF doesn't support state-by-state corrections, need to find a new approach")
        else:
            opts['max_vibr_excitations_in_initial_el_state'] = 0

        if len(states) > 0:
            if len(gs_states) == 0:
                gs_states = [[0]*len(nms.freqs)]
            starts = [
                cls.prep_state_str(nms, s) for s in gs_states
            ]
            ends = [
                cls.prep_state_str(nms, s) for s in states
            ]
            for i in starts:
                for j in ends:
                    elems.append(
                        ezXML.single_excitation(
                            ini=i,
                            targ=j
                        )
                    )
            # raise ValueError("ezFCF doesn't support state-by-state corrections, need to find a new approach")
        return opts, elems

    @classmethod
    def prep_state_str(cls, nms, state):
        if isinstance(state, str):
            return state
        else:
            if not len(state) == len(nms.freqs) or not nput.is_numeric(state[0]):
                from ..BasisReps import StateMaker
                state = StateMaker(len(nms.freqs), mode='normal')(*state)
            nzp = np.nonzero(state)
            if len(nzp) == 0 or len(nzp[0]) == 0: return "0"
            return ",".join(f"{i}v{state[i]}" for i in nzp[0])

    def prep_parallel(self, *elems, rotation_order='gs',
                      max_vibr_excitations_in_initial_el_state=0,
                      max_vibr_excitations_in_target_el_state=0,
                      **etc):
        elems = list(elems)
        opts = dict(
            dict(use_normal_coordinates_of_target_states=rotation_order == 'gs'),
            **etc
        )
        if self.print_all:
            elems.append(ezXML.print_franck_condon_matrices(flag=True))
            elems.append(ezXML.print_franck_condon_matrices(flag=True))
        return ezXML.parallel_approximation(
            *elems,
            max_vibr_excitations_in_initial_el_state=max_vibr_excitations_in_initial_el_state,
            max_vibr_excitations_in_target_el_state=max_vibr_excitations_in_target_el_state,
            **opts
        )

    def prep_duschinsky(self, *elems, rotation_order='gs',
                        max_vibr_excitations_in_initial_el_state=0,
                        max_vibr_excitations_in_target_el_state=0,
                        **etc):
        elems = list(elems)
        opts = dict(
            dict(target_state=0 if rotation_order == 'es' else 1),
            **etc
        )
        if self.print_all:
            elems.append(ezXML.print_franck_condon_matrices(flag=True))
        return ezXML.dushinsky_rotations(
            *elems,
            max_vibr_excitations_in_initial_el_state=max_vibr_excitations_in_initial_el_state,
            max_vibr_excitations_in_target_el_state=max_vibr_excitations_in_target_el_state,
            **opts
        )

    @classmethod
    def prep_initial_state(cls, atoms, nms, autoembed=False, excitation_energy_units="cm-1"):
        els = [
            # ezXML.excitation_energy(units=excitation_energy_units),
        ]
        if not autoembed:
            els.append(
                ezXML.manual_coordinates_transformation(
                    shift_along_x=0, shift_along_y=0, shift_along_z=0,
                    rotate_around_x=0, rotate_around_y=0, rotate_around_z=0
                )
            )
        return cls.prep_state(ezXML.initial_state,
                              atoms, nms,
                              *els#excitation_energy_units=excitation_energy_units
                              )

    @classmethod
    def prep_target_state(cls, atoms, nms, autoembed=False, mode_reordering=None, excitation_energy_units="cm-1"):
        """

        :param atoms:
        :param nms:
        :param excitation_energy_units:
        :return:
        """
        els = [
            ezXML.excitation_energy(units=excitation_energy_units),
        ]
        if not autoembed:
            els.append(
                ezXML.manual_coordinates_transformation(
                    shift_along_x=0, shift_along_y=0, shift_along_z=0,
                    rotate_around_x=0, rotate_around_y=0, rotate_around_z=0
                )
            )
        if mode_reordering is not None:
            els.append(
                ezXML.manual_normal_modes_reordering(
                    new_order=" ".join(str(int(x)) for x in mode_reordering)
                )
            )
        return cls.prep_state(
            ezXML.target_state,
            atoms, nms,
            *els,
        )

    @classmethod
    def prep_state(cls, tag, atoms, nms, *elems):
        return tag(
            *elems,
            cls.prep_geometry(atoms, nms),
            cls.prep_normal_modes(atoms, nms),
            cls.prep_frequencies(nms)
        )

    @classmethod
    def prep_geometry(cls, atoms, nms, linear=False, units="BohrRadius"):
        return ezXML.geometry(
            number_of_atoms=str(len(atoms)),
            linear=str(linear).lower(),
            units="angstr",
            text="\n" + "\n".join(
                f"{a:<8} {xyz[0]: >12.6f}  {xyz[1]: >12.6f}  {xyz[2]: >12.6f}"
                for a,xyz in zip(
                    atoms,
                    nms.remove_mass_weighting().origin * UnitsData.convert(units, "Angstroms")
                )
            ) + "\n"
        )

    @classmethod
    def format_modes_block(cls, modes):
        splits = itut.chunked(modes.T, 3)
        lines = []
        for block in splits:
            block = np.concatenate([b.reshape(-1, 3) for b in block], axis=1)
            block = [
                ['{:>10.6f}'.format(b) for b in l] + [""] * (9 - len(l))
                for l in block
            ]
            block_lines = [
                (
                    "   " + (" ".join(l[:3]))
                    + "   " + (" ".join(l[3:6]))
                    + "   " + (" ".join(l[6:9]))
                ).rstrip()
                for l in block
            ]
            lines.extend(block_lines)
            lines.append("")

        return "\n".join(lines)

    @classmethod
    def parse_modes_block(cls, modes_str, num_atoms=None):
        # blocks = [b.strip() for b in modes_str.splitlines()]
        # block_groups = itut.split(blocks, lambda l,_:len(l.strip()) == 0)
        modes = []
        if num_atoms is not None:
            base_array = np.loadtxt(io.StringIO(modes_str.replace("-", " -").replace("\n", " "))).flatten()
            nchunks = len(base_array) // (3 * 3 * num_atoms)
            blocks = np.array_split(base_array, nchunks)
            for block in blocks:
                modes.extend(b.T.flatten()[:, np.newaxis] for b in itut.chunked(block.reshape(num_atoms, -1).T, 3))
        else:
            blocks = modes_str.split("\n\n")
            for block in blocks:
                if len(block.strip()) > 0:
                    block = np.loadtxt(io.StringIO(block.replace("-", " -")), dtype=float)
                    modes.extend(b.T.flatten()[:, np.newaxis] for b in itut.chunked(block.T, 3))
        return np.concatenate(modes, axis=1)

    @classmethod
    def prep_normal_modes(cls, atoms, nms):
        modes = nms.modes_by_coords / UnitsData.convert("BohrRadius", "Angstroms")
        if nms.mass_weighted:
            modes = modes / np.sqrt(UnitsData.convert("ElectronMass", "AtomicMassUnits"))
        return ezXML.normal_modes(
            atoms=" ".join(atoms),
            text="\n"+cls.format_modes_block(modes) + "\n",
            if_mass_weighted=str(nms.mass_weighted).lower()
        )

    @classmethod
    def format_freqs_block(cls, freqs):
        splits = itut.chunked(freqs, 3)
        lines = []
        for l in splits:
            l = [" " + '{:>10.2f}'.format(b) for b in l] + [""] * (3 - len(l))
            lines.append("".join(l).strip())

        return "\n".join(lines)
    @classmethod
    def parse_freqs_block(cls, freqs):
        return np.loadtxt(io.StringIO(freqs.replace("\n", " ")))

    @classmethod
    def prep_frequencies(cls, nms):
        return ezXML.frequencies(
            text="\n"+cls.format_freqs_block(nms.freqs * UnitsData.convert("Hartrees", "Wavenumbers")) + "\n"
        )

    @classmethod
    def parse_state(self, state_xml):
        from ..Molecools import Molecule

        if isinstance(state_xml, str):
            state_xml = JHTML.XML.parse(state_xml)

        nms_block = state_xml.find("normal_modes")
        try:
            mw = nms_block["if_mass_weighted"]
        except AttributeError:
            mw = False
        else:
            mw = mw in {"true", "y"}
        ats = nms_block["atoms"].split()
        raw_modes = ezFCFInterface.parse_modes_block(
            nms_block["text"],
            num_atoms=len(ats)
        ) #/ UnitsData.convert("Angstroms", "BohrRadius")
        if mw:
            mode_norms = np.diag(raw_modes.T @ raw_modes)
            scaling = mode_norms[0]
            if scaling > 10:
                scaling = UnitsData.convert("Angstroms", "BohrRadius") * np.sqrt(
                    UnitsData.convert("AtomicMassUnits", "ElectronMass")
                )
            elif scaling > 2:
                scaling = UnitsData.convert("Angstroms", "BohrRadius")
            else:
                scaling = 1
            raw_modes = raw_modes * scaling

            masses = 1 / np.sqrt(
                    np.array([AtomData[s, "Mass"] for s in ats])
                    * UnitsData.convert("AtomicMassUnits", "ElectronMass")
            )
            masses = np.repeat(masses, 3)
            modes = np.diag(masses) @ raw_modes
            inverse = raw_modes.T @ np.diag(1/masses)
            # modes = raw_modes
        else:
            # TODO: shoud
            inverse = raw_modes.T
            modes = np.linalg.pinv(inverse)
        freqs = ezFCFInterface.parse_freqs_block(
            state_xml.find("frequencies")["text"]
        ) * UnitsData.convert("Wavenumbers", "Hartrees")

        return Molecule.from_string(
            state_xml.find("geometry")["text"],
            "xyz",
            units="Angstroms",
            normal_modes={
                "matrix": modes,
                "inverse": inverse,
                "freqs": freqs
            }
        )

    @classmethod
    def parse_fc_model(cls, input_xml, **opts):
        from .FCFs import FranckCondonModel

        if os.path.isfile(input_xml):
            with open(input_xml) as input_file:
                input_xml = input_file.read()

        if isinstance(input_xml, str):
            input_xml = JHTML.XML.parse(input_xml)

        gs = cls.parse_state(input_xml.find('initial_state'))
        es = cls.parse_state(input_xml.find('target_state'))

        include_rotation = input_xml.find('duschinsky_rotations') is not None
        mode_reordering = input_xml.find('target_state').find("manual_normal_modes_reordering")
        if mode_reordering is not None:
            mode_reordering = [int(s) for s in mode_reordering['new_order'].split()]

        return FranckCondonModel.from_mols(
            gs, es,
            include_rotation=include_rotation,
            mode_reordering=mode_reordering,
            remove_transrot=False,
            **opts
        )
