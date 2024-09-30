from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup
import sysconfig

extra_compile_args = sysconfig.get_config_var('CFLAGS').split()
extra_compile_args += ["-O3", "-Wall"]

ext_modules = [
    Pybind11Extension(
        "_atmosphere_bgs",
        ["src/main.cpp"],
        include_dirs = ["extern/eigen", "extern/pybind11/include", "extern/fmt/include"],
        cxx_std = 20,
        extra_compile_args = extra_compile_args,
    ),
]

#ext_modules[0]._add_ldflags(["-lstdc++exp"])

setup(
    ext_modules=ext_modules,
    zip_safe=False
)
