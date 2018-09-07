from builtins import range
from builtins import object
import os
import os.path as osp
import shutil as shu
try:
    import pyfits as pf
except Exception as e:
    pass

import re
import numpy as nm


def pack256(*li):
    rr = ""
    for l in li:
        rr += l + '\0' * (256 - len(l))
    return rr


def is_cldf(name):
    f = open(name)
    try:
        list(f.keys())
        return True
    except Exception as e:
        return False


def open(name, mode="r"):
    return File(name, mode)


_metadata = "_mdb"


class File(object):
    def __init__(self, name, mode="r"):
        self._mode = '+'
        if mode == "w":
            self._create(name)
            return
        if mode == "r" or mode == "r-":
            self._name = name
            self._mode = "-"
            return
        if mode == "r+":
            self._name = name

    def _parsemetadata(self, path=""):
        if not path:
            path = self._name
        f = file(osp.join(path, _metadata))
        dct = {}
        for l in f:
            if not l.strip():
                continue
            id0 = l.find(" ")
            key = l[:id0]
            id1 = l[id0 + 1:].find(" ") + id0 + 1
            typ = l[id0 + 1:id1]
            data = l[id1 + 1:-1]
            if typ == "int":
                dct[key] = int(data)
                continue
            if typ == "float":
                dct[key] = float(data)
                continue
            if typ == "str":
                dct[key] = data
                continue
            f.close()
            raise TypeError("unknown type '%s' for metadata '%s'" % (typ, key))
        f.close()
        return dct

    def _writemetadata(self, dct, path=""):
        if not path:
            path = self._name
        f = file(osp.join(path, _metadata), "w")
        for k, v in list(dct.items()):
            if type(v) == str:
                typ = "str"
                modi = "%s"
            elif type(v) in (bool, int, int, nm.int32, nm.int64):
                typ = "int"
                v = int(v)
                modi = "%d"
            elif type(v) in (float, nm.float32, nm.float64):
                typ = "float"
                modi = "%.10g"
            else:
                raise TypeError("bad type %s" % type(v))
            f.write(("%s %s " + modi + "\n") % (k, typ, v))
        f.close()

    def remove(self, name):
        if osp.exists(name):
            if osp.isdir(name):
                shu.rmtree(name)
            else:
                os.remove(name)
        else:
            dct = self._parsemetadata(osp.split(name)[0])
            if osp.split(name)[1] in list(dct.keys()):
                del dct[osp.split(name)[1]]
                self._writemetadata(dct, osp.split(name)[0])

    def _create(self, name):
        if osp.isdir(name):
            shu.rmtree(name)
        os.mkdir(name)
        f = file(osp.join(name, _metadata), "w")
        f.write("")
        f.close()
        self._name = name

    def __contains__(self, key):
        try:
            self[key]
        except Exception:
            return False
        return True

    def __getitem__(self, key):
        fkey = osp.join(self._name, key)
        if fkey[-1] == '/':
            fkey = fkey[:-1]
        if osp.exists(fkey):
            if osp.isdir(fkey):
                return File(fkey, "r" + self._mode)
            try:
                return pf.open(fkey)[0].data
            except Exception:
                value = file(fkey).read()
                if key + "__type__" in self and self[key + "__type__"] == "str_array":
                    rvalue = []
                    p0 = value.find("\n")
                    nv = int(value[:p0])
                    value = value[p0 + 1:]
                    for i in range(nv):
                        p1 = value.find("\n")
                        nc = int(value[:p1])
                        rvalue += [value[p1 + 1:p1 + 1 + nc]]
                        value = value[p1 + 1 + nc + 1:]
                    return rvalue
                return value

        dct = self._parsemetadata(osp.split(fkey)[0])
        return dct[osp.split(fkey)[1]]

    def __setitem__(self, key, value):
        assert self._mode == '+'
        fkey = osp.join(self._name, key)
        if fkey[-1] == '/':
            fkey = fkey[:-1]
        self.remove(fkey)
        if isinstance(value, File):

            shu.copytree(value._name, fkey)
            return
        if type(value) in (list, tuple, nm.ndarray):
            if isinstance(value[0], str):
                tvalue = "%d\n" % len(value)
                for v in value:
                    tvalue += "%d\n" % len(v) + v + "\n"
                f = file(fkey, "w")
                f.write(tvalue)
                f.close()
                self[key + "__type__"] = "str_array"
                return
            value = nm.array(value)
            if value.dtype == nm.int32:
                value = value.astype(nm.int64)
            # print key,fkey,value.dtype
            pf.PrimaryHDU(value).writeto(fkey)
            return
        if type(value) == str and ("\n" in value or "\0" in value or len(value) > 50):
            # print key,len(value)

            f = file(fkey, "w")
            f.write(value)
            f.close()
            return
        dct = self._parsemetadata(osp.split(fkey)[0])
        dct[osp.split(fkey)[1]] = value
        self._writemetadata(dct, osp.split(fkey)[0])

    def create_group(self, name):
        assert self._mode == '+'
        return File(osp.join(self._name, name), "w")

    def create_dataset(self, name, data=None):
        assert data != None
        self[name] = data

    def __delitem__(self, key):
        assert self._mode == '+'
        fkey = osp.join(self._name, key)
        if fkey[-1] == '/':
            fkey = fkey[:-1]

        if osp.exists(fkey):
            self.remove(fkey)
            return
        dct = self._parsemetadata(osp.split(fkey)[0])
        del dct[osp.split(fkey)[1]]
        self._writemetadata(dct, osp.split(fkey)[0])

    def copy(self, a, b, c=""):
        if not c:
            self[b] = self[a]
        else:
            b[c] = a

    @property
    def attrs(self):
        return self

    def keys(self):
        dct = self._parsemetadata(self._name)
        ls = [el for el in os.listdir(
            self._name) if el[0] != '.' and el != _metadata]
        return ls + list(dct.keys())

    def items(self):
        ks = list(self.keys())
        return [(k, self[k]) for k in ks]

    def close(self):
        pass  # nothing to do


try:
    import h5py

    def hdf2cldf_grp(hdf, fdf):
        # first the metadata
        for kk in list(hdf.attrs.keys()):
            vl = hdf.attrs[kk]

            # print kk,type(vl)
            if type(vl) == str:
                sz = h5py.h5a.get_info(hdf.id, kk).data_size
                rr = vl.ljust(sz, '\0')
                fdf[kk] = rr
            else:
                fdf[kk] = vl
        # then the group/data
        for kk in list(hdf.keys()):
            if kk == "external_data":
                dts = hdf[kk][:]
                install_path = osp.join(fdf._name, "_external")
                os.mkdir(install_path)
                f = file(osp.join(install_path, "data.tar"), "w")
                f.write(dts.tostring())
                f.close()
                assert os.system("cd %s;tar xvf data.tar" % install_path) == 0
                assert os.system("cd %s;rm -f data.tar" % install_path) == 0
                fdf["external_dir"] = "."
                continue
            god = hdf[kk]
            if isinstance(god, h5py.Group):
                if not hasattr(fdf, kk):
                    fdf.create_group(kk)
                hdf2cldf_grp(god, fdf[kk])
            else:
                r = god[:]
                # print r
                if len(r) == 1:
                    r = r[0]
                fdf[kk] = r

    def hdf2cldf(ffin, ffout):
        hdf = h5py.File(ffin, "r")
        fdf = File(ffout, "w")
        hdf2cldf_grp(hdf, fdf)
except ImportError as e:
    pass
