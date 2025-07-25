import numpy as np
import McUtils.Devutils as dev
import McUtils.Numputils as nput
from McUtils.ExternalPrograms import PILInterface
from McUtils.Plots import ColorPalette

__all__ = [
    "SpectrumExtractor"
]

class SpectrumExtractor:

    def __init__(self, image_data, color_space='rgb', operation_color_space='lab'):
        image_data = np.asanyarray(image_data)
        if color_space != operation_color_space:
            image_data = ColorPalette.color_convert(image_data, color_space, operation_color_space)
        self.img = image_data
        self.color_space = operation_color_space

    @classmethod
    def from_pil(cls, pil_image, color_space=None, **etc):
        if color_space is None:
            color_space = pil_image.mode.lower()
            if color_space == 'l':
                color_space = 'grayscale'
            elif color_space == 'rgba':
                color_space = 'rgb'

        image_data = np.moveaxis(np.asarray(pil_image), -1, 0)[:3]
        if color_space == 'rgb':
            image_data = np.clip(image_data, 0, 255)

        return cls(
            image_data,
            color_space=color_space,
            **etc
        )
    @classmethod
    def from_file(cls, file, color_space=None, **etc):
        pil = PILInterface.from_file(file)

        return cls.from_pil(pil, color_space=color_space, **etc)
    @classmethod
    def from_url(cls, file, color_space=None):
        pil = PILInterface.from_url(file)
        return cls.from_pil(pil, color_space=color_space)

    default_tolerances = [25, 25, 25]
    def find_pixels(
            self,
            color,
            distance_tol=None,
            tolerances=None,
            color_space='rgb',
            search_color_space=None,
            image_range=None
    ):
        if isinstance(color, str):
            color = ColorPalette.parse_color_string(color)

        if search_color_space is None:
            search_color_space = self.color_space
        color = ColorPalette.color_convert(color, color_space, search_color_space)

        image_data = self.img
        if search_color_space is not None and search_color_space != self.color_space:
            image_data = ColorPalette.color_convert(image_data, self.color_space, search_color_space)

        if image_range is not None:
            y,x = image_range
            if nput.is_numeric(x):
                y,x = None, [y,x]
            if x is not None:
                xm, xM = x
                if not nput.is_int(xm):
                    xm = int(image_data.shape[-2] * xm)
                if not nput.is_int(xM):
                    xM = int(image_data.shape[-2] * xM)
                image_data = image_data[:, xm:xM, :]
            if y is not None:
                ym, yM = y
                if not nput.is_int(ym):
                    ym = int(image_data.shape[-1] * ym)
                if not nput.is_int(yM):
                    yM = int(image_data.shape[-1] * yM)
                image_data = image_data[:, :, ym:yM]

        if distance_tol is not None:
            raise NotImplementedError("finding pixels from color distance not supported")
        else:
            if tolerances is None:
                tolerances = self.default_tolerances

            mask = np.full(image_data[0].shape, True)
            for channel, center, tol in zip(
                image_data,
                color,
                tolerances
            ):
                mask = np.logical_and(
                    mask,
                    np.logical_and(
                        channel >= center - tol,
                        channel <= center + tol
                    )
                )

        return np.where(mask)

    def find_spectrum_lines(self,
                            pixel_positions,
                            max_pixel_distance=.02,
                            min_line_cutoff=0.05,
                            smoothing=True,
                            line_split_cutoff=5,
                            allow_disjoint=False,
                            spectrum_direction='up'
                            ):

        if max_pixel_distance < 1:
            max_pixel_distance = max_pixel_distance * max(self.img.shape[-1:])

        lines = []
        for x,Y in zip(*nput.group_by(pixel_positions[1], pixel_positions[0])[0]):
            Y = np.sort(Y) # just in case
            # if smoothing:
            y_diffs = np.diff(Y)
            split_points = np.where(np.abs(y_diffs) > line_split_cutoff)
            if len(split_points) > 0 and len(split_points[0]) > 0:
                Y = [
                    g
                    for g in np.split(Y, split_points[0]+1)
                ]
            else:
                Y = [Y]

            if spectrum_direction.lower() == 'down':
                Y = reversed(Y)

            if len(lines) == 0:
                lines.extend([[x, y]] for y in Y)
            else:
                add_pos = set()
                for y in Y:
                    min_y = y[0]
                    max_y = y[-1]
                    for n,l in enumerate(lines):
                        if n in add_pos: continue

                        xx, yy = l[-1]
                        dx = x-xx
                        dy = min([abs(min_y-yy[0]), abs(max_y - yy[-1])])
                        if abs(dx)+abs(dy) < max_pixel_distance:
                            l.append([x, y])
                            if not allow_disjoint:
                                add_pos.add(n)
                            break
                    else:
                        lines.append([[x, y]])


        if smoothing:
            lines = [
                [[x, np.average(y)] for x, y in l]
                for l in lines
            ]
        else:
            lines = [
                sum(
                    ([[x, yy] for yy in y] for x, y in l),
                    []
                )
                for l in lines
            ]

        if min_line_cutoff < 1:
            min_line_cutoff = len(pixel_positions[0]) * min_line_cutoff
        return [
            np.array(l).T
            for l in sorted(
                (ll for ll in lines if len(ll) > min_line_cutoff),
                key=lambda x:-len(x)
            )
        ]

    default_merge_tolerances = [10, 10, 10]
    def get_dominant_colors(self, bins=255, color_space='lab', min_counts=2, merge_tolerances=None):
        image_data = self.img
        if color_space != self.color_space:
            image_data = ColorPalette.color_convert(image_data, self.color_space, color_space)
        dists, (x, y, z) = np.histogramdd(np.moveaxis(image_data, 0, -1).reshape(-1, 3), bins=bins)
        take = np.where(dists > min_counts)
        base_colors = {
            (
                (x[i] + x[i+1])/2,
                (y[j] + y[j+1])/2,
                (z[k] + z[k+1])/2,
            ):int(dists[i,j,k])
            for i,j,k in zip(*take)
        }

        base_colors = {
            k:base_colors[k]
            for k in sorted(base_colors, key=lambda x:-base_colors[x])
        }
        if merge_tolerances is None:
            merge_tolerances = self.default_merge_tolerances

        test_keys = list(base_colors.keys())
        merges = {}
        for n,c in enumerate(base_colors.keys()):
            for k in test_keys[:max(n-1, 0)]:
                if c == k or k in merges: continue
                if all(abs(k[i] - c[i]) <= merge_tolerances[i] for i in range(3)):
                    merges[c] = k
                    break

        for c,k in merges.items():
            base_colors[k] += base_colors[c]
            del base_colors[c]

        # print(
        #     sum(v for v in base_colors.values()) / (image_data.shape[-1] * image_data.shape[-2])
        # )

        return {
            k:base_colors[k]
            for k in sorted(base_colors, key=lambda x:-base_colors[x])
        }


        #TOOD: move this to OpenCV
        # all_pixels = np.moveaxis(self.img, 0, -1).reshape(-1, 3)
        # pix = {}
        # if tolerances is None:
        #     tolerances = self.default_tolerances
        # tolerances = np.asanyarray(tolerances)
        # m = 0
        # for p in all_pixels:
        #     old_pix = None
        #     new_pix = {}
        #     m += 1
        #     for k in pix:
        #         if all(k[i] - p[i] <= tolerances[i] for i in range(3)):
        #             n = pix[k]
        #             p = tuple((n*k[i] + p[i]) / (n+1) for i in range(3))
        #             new_pix[p] = n + 1
        #             old_pix = k
        #             break
        #     else:
        #         p = tuple(p)
        #         new_pix[p] = 1
        #     if old_pix is not None:
        #         del pix[old_pix]
        #     pix.update(new_pix)
        #
        # return pix


    def extract_spectra(self,
                        color=None,
                        use_exact_color=False,
                        image_range=None,
                        max_dominant_percentage=.2,
                        allow_grayscale=False,
                        color_space='rgb',
                        dominant_color_space='lab',
                        dominant_bins=255,
                        min_dominant_component=50,
                        extract_lines=True,
                        spectrum_direction='up',
                        x_range=(0, 1),
                        y_range=(0, 1),
                        preserve_x_range=True,
                        preserve_y_range=False,
                        **opts
                        ):
        if isinstance(color, str):
            color = ColorPalette.parse_color_string(color)[:3]
            color_space = 'rgb'
        if color is None or not use_exact_color:
            opts = dev.OptionsSet(opts)
            dom_opts, opts = opts.split(None,
                                        {
                                            'bins',
                                            'merge_tolerances'
                                        }
                                        )
            doms = self.get_dominant_colors(dominant_bins,
                                            min_counts=min_dominant_component,
                                            color_space=dominant_color_space,
                                            **dom_opts)
            npixels = self.img.shape[-2] * self.img.shape[-1]
            max_npix = npixels * max_dominant_percentage
            if color is None:
                color = next((
                    d for d in doms
                    if (
                        doms[d] < max_npix
                        and (
                                allow_grayscale
                                or ColorPalette.color_convert(d, dominant_color_space, 'lch')[1] > .1
                        )
                    )
                ))
            else:
                c = ColorPalette.color_convert(color, color_space, dominant_color_space)
                color = next((
                    d for d in
                    sorted(doms,
                           key=lambda k: (np.linalg.norm(c - np.array(k))**2) * (doms[k] / npixels)
                           )
                ))
            color_space = dominant_color_space
        opts = dev.OptionsSet(opts)
        line_opts, pixel_opts = opts.split(None,
                                           {
                                               'max_pixel_distance',
                                               'min_line_percentage',
                                               'smoothing',
                                               'line_split_cutoff',
                                               'allow_disjoint'
                                           }
                                           )
        y, x = self.find_pixels(color, image_range=image_range, color_space=color_space, **pixel_opts)
        if len(x) == 0 or len(y) == 0:
            return color, []

        y = self.img.shape[1] - y
        if extract_lines:
            lines = self.find_spectrum_lines((x, y), spectrum_direction=spectrum_direction, **line_opts)
        else:
            lines = [np.array([x, y])]


        if x_range is not None:
            if nput.is_int(x_range):
                x_range = [0, x_range]

            if preserve_x_range:
                x_mins = [
                    np.min(xx)
                    for xx, yy in lines
                    if len(xx) > 0
                ]
                if len(x_mins) == 0:
                    x_mins = [np.min(x)]

                x_maxes = [
                    np.max(xx)
                    for xx, yy in lines
                    if len(xx) > 0
                ]
                if len(x_maxes) == 0:
                    x_maxes = [np.max(x)]

                x_span = [
                    min(x_mins),
                    max(x_maxes)
                ]
                x_scale = (x_range[1] - x_range[0]) / (x_span[1] - x_span[0])
                lines = [
                    [
                        (x - x_span[0]) * x_scale + x_range[0],
                        y
                    ]
                    for x,y in lines
                ]
            else:
                lines = [
                    [
                        (x - np.min(x)) * (x_range[1] - x_range[0]) / (np.max(x) - np.min(x)) + x_range[0],
                        y
                    ]
                    for x, y in lines
                ]

        if y_range is not None:
            if nput.is_int(y_range):
                y_range = [0, y_range]

            if preserve_y_range:
                y_span = [
                    min(np.min(y) for x,y in lines),
                    max(np.max(y) for x,y in lines)
                ]
                y_scale = (y_range[1] - y_range[0]) / max([(y_span[1] - y_span[0]), 1])
                lines = [
                    [
                        x,
                        (y - y_span[0]) * y_scale + y_range[0],
                    ]
                    for x,y in lines
                ]
            else:
                lines = [
                    [
                        x,
                        (y - np.min(y)) * (y_range[1] - y_range[0]) / (np.max(y) - np.min(y)) + y_range[0]
                    ]
                    for x, y in lines
                ]

        return color, lines

