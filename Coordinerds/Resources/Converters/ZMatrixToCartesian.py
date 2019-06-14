from .CoordinateSystemConverter import CoordinateSystemConverter
from .CommonCoordinateSystems import CartesianCoordinates3D, ZMatrixCoordinates
from ..CoordinateTransformations.TransformationUtilities import vec_norms, vec_angles, pts_dihedrals, vec_normalize, vec_crosses
from ..CoordinateTransformations.TransformationUtilities import affine_matrix, merge_transformation_matrix, rotation_matrix, affine_multiply
import numpy as np
# this import gets bound at load time, so unfortunately PyCharm can't know just yet
# what properties its class will have and will try to claim that the files don't exist

class ZMatrixToCartesianConverter(CoordinateSystemConverter):
    """A converter class for going from ZMatrix coordinates to Cartesian coordinates

    """

    @property
    def types(self):
        return (ZMatrixCoordinates, CartesianCoordinates3D)

    @staticmethod
    def zmatrix_affine_transforms(centers, vecs1, vecs2, angles, dihedrals):
        """Builds a single set of affine transformation matrices to apply to the vecs1 to get the next set of points

        :param refs1:
        :type refs1:
        :param refs2:
        :type refs2:
        :param refs3:
        :type refs3:
        :param angles:
        :type angles:
        :param dihedrals:
        :type dihedrals:
        :return:
        :rtype:
        """
        crosses = vec_crosses(vecs2, vecs1)
        rot_mats_1 = rotation_matrix(crosses, angles) # not entirely
        if dihedrals is not None:
            # this is where things break down
            # looks like I don't get the correct rotation into the dihedral frame
            rot_mats_2 = rotation_matrix(vecs1, -dihedrals)
            # raise Exception(dihedrals, vecs2, rot_mats_2)
            rot_mat = np.matmul(rot_mats_2, rot_mats_1)
        else:
            rot_mat = rot_mats_1
        transfs = affine_matrix(rot_mat, centers)
        # raise Exception(transfs, rot_mats_1, crosses, angles, vecs1)
        return transfs

    def build_next_points(self, refs1, dists, refs2, angles, refs3, dihedrals):
        vecs1 = dists*vec_normalize(refs2 - refs1)
        if dihedrals is not None:
            vecs2 = refs3 - refs1
        else:
            # vecs2 = np.broadcast_to([0., 1., 0.], (len(refs1), 3))
            # if no dihedral for now we'll let our new axis be some random things in the x-y plane
            vecs2 = np.concatenate((np.random.uniform(.5, 1, (len(refs1), 2)), np.zeros((len(refs1), 1))), axis=-1)
        transfs = self.zmatrix_affine_transforms(refs1, vecs1, vecs2, angles, dihedrals)
        newstuff = affine_multiply(transfs, vecs1)
        return newstuff


    def convert_many(self, coordlist, origins=None, axes=None, use_rad=True, **kw):
        """Expects to get a list of configurations
        These will look like:
            [
                [anchor, dist, ref, angle, plane, dihedral ]
                ...
            ]
        **For efficiency it is assumed that all configurations have the same length**

        :param coordlist:
        :type coordlist:
        :param origins:
        :type origins:
        :param axes:
        :type axes:
        :param use_rad:
        :type use_rad:
        :param kw:
        :type kw:
        """
        #still need to handle the n < 4 case...
        sysnum = len(coordlist)
        coordnum = len(coordlist[0])
        total_points = np.empty((sysnum, coordnum+1, 3))

        # gotta build the first three points by hand but everything else can be constructed iteratively
        if origins is None:
            origins = [0, 0, 0]
        if not isinstance(origins, np.ndarray):
            origins = np.array(origins)
        if len(origins.shape) < 2:
            origins = np.broadcast_to(origins, (sysnum, 3))
        total_points[:, 0] = origins

        # set up the next points by just setting them along the x-axis by default
        if axes is None:
            axes = [1, 0, 0]
        if not isinstance(axes, np.ndarray):
            axes = np.array(axes)
        if len(axes.shape) < 2:
            axes = np.broadcast_to(axes, (sysnum, 3))
        axes = vec_normalize(axes)
        dists = coordlist[:, 0, 1]
        ref_points_1 = np.reshape(dists, (sysnum, 1)) * axes
        ref_points_1 += origins
        total_points[:, 1] = ref_points_1

        # iteratively build the rest of the coords with one special cases for n=2
        for i in range(1, coordnum):

            # Get the distances away
            # raise Exception(coordlist[:, i, [0, 2, 4]])
            ref_coords1 = coordlist[:, i, 0].astype(np.int8) - 1 # reference atom numbers for first coordinate
            refs1 = total_points[np.arange(sysnum), ref_coords1] # get the actual reference coordinates
            dists = np.reshape(coordlist[:, i, 1], (sysnum, 1)) # pull the requisite distances

            ref_coords2 = coordlist[:, i, 2].astype(np.int8) - 1 # reference atom numbers for second coordinate
            refs2 = total_points[np.arange(sysnum), ref_coords2] # get the actual reference coordinates for the angle
            angle = coordlist[:, i, 3] # pull the requisite angle values
            if not use_rad:
                angle = np.deg2rad(angle)

            if i == 1:
                refs3 = None
                dihed = None
            else:
                ref_coords3 = coordlist[:, i, 4].astype(np.int8) - 1 # reference atom numbers for dihedral ref coordinate
                refs3 = total_points[np.arange(sysnum), ref_coords3] # get the actual reference coordinates for the dihed
                dihed = coordlist[:, i, 5] # pull proper dihedral values
                if not use_rad:
                    dihed = np.deg2rad(dihed)

            # I FUCKING CAN"T GENERATE A FUCKING PIECE OF SHIT RANDOM Z MATRIX AND THAT IS ALL THAT IS FUCKING ME RIGHT NOW
            # raise Exception((coordlist[:, :, [0, 2]], coordlist[:, i, [0, 2]], ref_coords1-ref_coords2))
            ref_points = self.build_next_points(refs1, dists, refs2, angle, refs3, dihed) # iteratively build z-mat
            total_points[:, i+1] = ref_points

        return total_points




    def convert(self, coords, **kw):
        """dipatches to convert_many but only pulls the first"""
        return self.convert_many(np.reshape(coords, (1,)+coords.shape), **kw)[0]


__converters__ = [ ZMatrixToCartesianConverter() ]