#!/usr/bin/env python3
"""Shared binary-database reader and Fortran backend for direct FSCK training."""
import array
import ctypes as ct
import hashlib
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from pressure_selection import selected_pressure_folders

ROOT = Path(__file__).resolve().parent
FORTRAN_DIR = ROOT.parent / "fortran"
DEFAULT_DB = Path('/Volumes/home/cnp10/CODES/radiationDatabase/FSKTableV4')
PATTERN = re.compile(r'P\.(.+)_CO2\.(.+)_H2O\.(.+)_CO\.(.+)_fv\.(.+)_k\.dat$')
SOURCES = ['commonData.f90', 'quadlib.f90', 'commonRoutines.f90',
           'FSKDBPath.f90', 'libJCai.f90', 'fskTable.f90', 'fskTableV4.f90']

def build(compiler='gfortran', debug=False):
    executable = shutil.which(compiler)
    if not executable:
        raise RuntimeError(f'Fortran compiler not found: {compiler}')
    folder = ROOT / 'build'
    folder.mkdir(exist_ok=True)
    original = (FORTRAN_DIR/'fskTableV4.f90').read_text()
    extension = (FORTRAN_DIR/'row_bridge.inc').read_text()
    if original.count('end module fskTableV4') != 1:
        raise RuntimeError('Unexpected V4 module layout')
    patched = original.replace('end module fskTableV4', extension+'\nend module fskTableV4')
    (folder/'fskTableV4.f90').write_text(patched)
    library = folder / ('libfsck.dylib' if sys.platform == 'darwin' else 'libfsck.so')
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if debug else ['-O2']
    command = [executable, '-shared', '-fPIC', '-ffree-line-length-none'] + flags
    command += [str(FORTRAN_DIR/f) for f in SOURCES[:-1]] + [str(folder/'fskTableV4.f90'), '-o', str(library)]
    subprocess.run(command, cwd=folder, check=True)
    return library

class Backend:
    def __init__(self, library):
        self.lib = ct.CDLL(str(library))
        self.call = self.lib.fsck_file
        ptr = ct.POINTER(ct.c_double)
        self.call.argtypes = [ct.c_char_p, ct.c_int, ptr, ptr, ptr, ptr, ct.POINTER(ct.c_int)]
        self.call.restype = None
        self.output = (ct.c_double*(784*32))()
        self.g = (ct.c_double*32)()
        self.w = (ct.c_double*32)()
    def evaluate(self, path, state):
        encoded = os.fsencode(path)
        code = ct.c_int()
        self.call(encoded, len(encoded), (ct.c_double*5)(*state), self.output, self.g, self.w, ct.byref(code))
        if code.value:
            raise RuntimeError(f'Fortran read failure {code.value}: {path}')
        if not all(math.isfinite(x) for x in self.output):
            raise RuntimeError(f'Nonfinite Fortran results: {path}')
        return self.output

def files(database, pressure=None):
    directories = selected_pressure_folders(database.iterdir(), pressure)
    for directory in directories:
        for path in sorted(directory.glob('*_k.dat')):
            if path.name.startswith('.'):
                continue
            match = PATTERN.fullmatch(path.name)
            if not match:
                raise ValueError(f'Cannot parse database filename: {path}')
            state = tuple(float(x) for x in match.groups())
            if state[0] != float(directory.name):
                raise ValueError(f'Pressure mismatch: {path}')
            yield path, state

def validate_raw(path):
    raw = array.array('f')
    raw.frombytes(path.read_bytes())
    if raw.itemsize != 4 or len(raw) != 784*33:
        raise ValueError(f'Expected 784 records of 33 float32 values: {path}')
    if any(not math.isfinite(x) or x < 0 for x in raw):
        raise ValueError(f'Invalid kappa data or wrong native byte order: {path}')
    return raw

def source_digest():
    h=hashlib.sha256()
    for p in sorted(FORTRAN_DIR.iterdir()):
        if not p.is_file() or p.name.startswith('.'):
            continue
        h.update(p.name.encode()); h.update(p.read_bytes())
    return h.hexdigest()
