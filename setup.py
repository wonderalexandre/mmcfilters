import os
import re
import sys
import sysconfig
import platform
import subprocess

from pathlib import Path 
from distutils.version import LooseVersion
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import setuptools

import shutil  

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extension: "
                               + ", ".join(e.name for e in self.extensions))

        cmake_version = LooseVersion(re.search(r"version\s*([\d.]+)", out.decode()).group(1))
        if platform.system() == "Windows" and cmake_version < "3.14":
            raise RuntimeError("CMake >= 3.14 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        prefix = sysconfig.get_config_var("LIBDIR")

        # Argumentos do CMake
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
            "-DPYTHON_LIBRARY_DIR={}".format(prefix)
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        # Configurações específicas para plataformas
        if platform.system() == "Windows":
            cmake_args += ["-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            # Configurações para Linux/Unix
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["--", "-j2"]

        # Adicionar variáveis de ambiente para a compilação
        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get("CXXFLAGS", ""), self.distribution.get_version()).replace('"', '\\"')

        # Criar o diretório de build se não existir
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        # Executar o CMake
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)

        # Mover o output gerado para o local correto
        self.move_output(ext)

    def move_output(self, ext):
        extdir = Path(self.build_lib).resolve()
        dest_path = Path(self.get_ext_fullpath(ext.name)).resolve()
        source_path = extdir / self.get_ext_filename(ext.name)
        dest_directory = dest_path.parent

        # Verificar se o arquivo foi gerado
        if not source_path.exists():
            raise RuntimeError(f"Arquivo de saída {source_path} não encontrado. Verifique se o CMake gerou corretamente o arquivo .so.")

        # Criar diretório de destino se necessário e mover o arquivo
        dest_directory.mkdir(parents=True, exist_ok=True)
        self.copy_file(source_path, dest_path)

        



        
# Verifica o sistema operacional atual
system = platform.system()

# Define a extensão do arquivo nativo com base no sistema operacional
if system == "Linux":
    native_ext = "*.so"
elif system == "Darwin":
    native_ext = "*.so"  # Para macOS, o Pybind gera .so, mas você pode usar *.dylib para bibliotecas dinâmicas
elif system == "Windows":
    native_ext = "*.dll"
else:
    raise RuntimeError(f"Plataforma {system} não suportada!")

setup(
    name="mmcfilters",
    version="0.0.9",
    description="A simple library for connected image filtering based on morphological trees",
    long_description="",
    author="Wonder Alexandre Luz Alves",
    author_email="worderalexandre@gmail.com",
    license="GPL-3.0",
    url="https://github.com/wonderalexandre/ComponentTreeLearn",
    keywords="machine learning, morphological trees, mathematical morphology",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Image Processing",
        "Programming Language :: Python",
        "Programming Language :: C++",
    ],
    # Certificar-se de que as dependências necessárias estão instaladas
    setup_requires=["setuptools", "wheel", "cmake>=3.14"],
    ext_modules=[CMakeExtension('mmcfilters')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    packages=["maf"],  # Definindo o pacote maf
    package_dir={"maf": "python"},  # Atribuindo o diretório "python" ao pacote "maf"
    package_data={
        "maf": ["*.py", native_ext],  # Incluindo todos os arquivos Python no pacote e os modelos c++/pybinds
    },
)


#send to pypi
#1. python setup.py sdist
#2. pipenv run twine upload dist/* -r pypi