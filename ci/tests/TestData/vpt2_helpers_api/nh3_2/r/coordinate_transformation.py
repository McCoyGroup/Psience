import numpy as np

def conversion(r, t, f, **kwargs):
    cp1 = np.cos(f[..., 3])  # skip three for embedding
    ct1 = np.cos(t[..., 2])  # skip two for embedding
    ct2 = np.cos(t[..., 3])
    st1 = np.sin(t[..., 2])
    st2 = np.sin(t[..., 3])
    f[..., 3] = np.arccos(st1 * st2 * cp1 + ct1 * ct2)
    return np.array([r, t, f])

def inverse(r, t, f, **kwargs):
    cp1 = np.cos(f[..., 3])
    ct1 = np.cos(t[..., 2])
    ct2 = np.cos(t[..., 3])
    st1 = np.sin(t[..., 2])
    st2 = np.sin(t[..., 3])
    f[..., 3] = np.arccos((cp1 - ct1 * ct2) / (st1 * st2))
    return np.array([r, t, f])