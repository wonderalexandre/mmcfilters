import os
import platform
import subprocess
import sys
import sysconfig
from pathlib import Path

from typing import Any, Dict

from packaging.version import Version, InvalidVersion
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


def read_version() -> str:
    """Read the package version from python/version.py."""
    about: Dict[str, Any] = {}
    version_file = Path(__file__).resolve().parent / "python" / "version.py"
    exec(version_file.read_text(encoding="utf-8"), about)
    return about["__version__"]


class CMakeExtension(Extension):
    """CMake-backed extension."""

    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = Path(sourcedir).resolve()


class CMakeBuild(build_ext):
    """Custom build command that integrates CMake with setuptools."""

    def run(self) -> None:
        try:
            out = subprocess.check_output(["cmake", "--version"], text=True)
        except OSError as exc:
            raise RuntimeError(
                "CMake must be installed to build the extensions for this package."
            ) from exc

        version_token = out.split("version", maxsplit=1)[-1].strip().split()[0]
        try:
            cmake_version = Version(version_token)
        except InvalidVersion as exc:
            raise RuntimeError(
                "Unable to identify the installed CMake version."
            ) from exc

        if platform.system() == "Windows" and cmake_version < Version("3.14"):
            raise RuntimeError("CMake >= 3.14 is required on Windows.")

        super().run()

    def build_extension(self, ext: Extension) -> None:
        extdir = Path(self.get_ext_fullpath(ext.name)).parent.resolve()
        prefix = sysconfig.get_config_var("LIBDIR")

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
        ]

        if prefix:
            cmake_args.append(f"-DPYTHON_LIBRARY_DIR={prefix}")

        cmake_args += [
            "-DBUILD_PYBIND=ON",
            "-DBUILD_TESTS=OFF",
            "-DMMC_PYBINDS_DIR=pybinds",
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args.append(f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}")
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args.append(f"-DCMAKE_BUILD_TYPE={cfg}")
            build_args += ["--", f"-j{os.cpu_count() or 2}"]

        env = os.environ.copy()
        version = self.distribution.get_version()
        cxxflags = env.get("CXXFLAGS", "")
        env["CXXFLAGS"] = f"{cxxflags} -DVERSION_INFO=\\\"{version}\\\""

        if self.parallel:
            env.setdefault("CMAKE_BUILD_PARALLEL_LEVEL", str(self.parallel))

        build_temp = Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)

        subprocess.check_call(["cmake", str(ext.sourcedir)] + cmake_args, cwd=build_temp, env=env)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=build_temp, env=env)

        self._move_output(ext)

    def _move_output(self, ext: Extension) -> None:
        ext_path = Path(self.get_ext_fullpath(ext.name)).resolve()
        source_path = Path(self.build_lib).resolve() / self.get_ext_filename(ext.name)
        destination = ext_path.parent

        if not source_path.exists():
            raise RuntimeError(
                f"Output file {source_path} not found. Check the CMake build output."
            )

        destination.mkdir(parents=True, exist_ok=True)
        self.copy_file(str(source_path), str(ext_path))


def read_long_description() -> str:
    readme_path = Path(__file__).resolve().parent / "README.md"
    return readme_path.read_text(encoding="utf-8")


NATIVE_EXTENSIONS = {
    "Linux": "*.so",
    "Darwin": "*.so",
    "Windows": "*.dll",
}

system = platform.system()
if system not in NATIVE_EXTENSIONS:
    raise RuntimeError(f"Platform {system} is not supported!")

setup(
    name="mmcfilters",
    version=read_version(),
    description="Library for connected image filtering based on morphological trees",
    long_description=read_long_description(),
    long_description_content_type="text/markdown",
    author="Wonder Alexandre Luz Alves",
    author_email="worderalexandre@gmail.com",
    license_files=["LICENSE"],
    license="GPL-3.0",
    url="https://github.com/wonderalexandre/ComponentTreeLearn",
    project_urls={
        "Source": "https://github.com/wonderalexandre/ComponentTreeLearn",
        "Documentation": "https://github.com/wonderalexandre/ComponentTreeLearn",
    },
    keywords=[
        "morphological trees",
        "mathematical morphology",
        "image processing",
        "computer vision",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Image Processing",
        "Programming Language :: Python",
        "Programming Language :: C++",
    ],
    ext_modules=[CMakeExtension("mmcfilters")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
    packages=["maf"],
    package_dir={"maf": "python"},
    package_data={"maf": ["*.py", NATIVE_EXTENSIONS[system]]},
    include_package_data=True,
)
