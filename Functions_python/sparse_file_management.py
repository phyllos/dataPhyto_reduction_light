import h5py
import numpy as np
from scipy import sparse


def read_vector(h5, name):
    return np.asarray(h5[name][()]).squeeze()


def write_vector(h5, name, data):
    arr = np.asarray(data).squeeze()
    h5.create_dataset(name, data=arr)


def write_text(h5, name, text):
    h5.create_dataset(name, data=np.bytes_(text))


def read_sparse(h5, name):
    """Read either MATLAB v7.3 sparse group or native scipy CSC group."""
    g = h5[name]
    if "MATLAB_sparse" in g.attrs:
        n_rows = int(g.attrs["MATLAB_sparse"])
        jc = np.asarray(g["jc"][()], dtype=np.int64)
        ir = np.asarray(g["ir"][()], dtype=np.int64)
        data = np.asarray(g["data"][()], dtype=np.float64)
        n_cols = len(jc) - 1
        return sparse.csc_matrix((data, ir, jc), shape=(n_rows, n_cols))
    data = np.asarray(g["data"][()], dtype=np.float64)
    indices = np.asarray(g["indices"][()], dtype=np.int64)
    indptr = np.asarray(g["indptr"][()], dtype=np.int64)
    shape = tuple(np.asarray(g["shape"][()], dtype=np.int64))
    return sparse.csc_matrix((data, indices, indptr), shape=shape)


def write_sparse(h5, name, mat):
    """Write scipy sparse matrix in simple native HDF5 CSC format."""
    mat = mat.tocsc()
    g = h5.create_group(name)
    g.attrs["format"] = "csc"
    g.create_dataset("data", data=mat.data.astype(np.float64, copy=False))
    g.create_dataset("indices", data=mat.indices.astype(np.int64, copy=False))
    g.create_dataset("indptr", data=mat.indptr.astype(np.int64, copy=False))
    g.create_dataset("shape", data=np.asarray(mat.shape, dtype=np.int64))


def copy_if_exists(fin, fout, names):
    for name in names:
        if name in fin:
            obj = fin[name]
            if isinstance(obj, h5py.Dataset):
                fout.create_dataset(name, data=obj[()])
