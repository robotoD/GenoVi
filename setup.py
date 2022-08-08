from setuptools import setup, find_packages
import codecs
from os import path
import re

extra_test = [
  'pytest>=4',
  'pytest-cov>=2',
],

extra_dev = [
  *extra_test,
]

extra_ci = [
  *extra_test,
  'python-coveralls',
]

here = path.abspath(path.dirname(__file__))

# Single-sourcing the package version: Read from init
def read(*parts):
    with codecs.open(path.join(here, *parts), 'r') as fp:
        return fp.read()


def find_version(*file_paths):
    print(*file_paths)
    version_file = read(*file_paths)
    print(version_file)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

setup(
  name = 'genovi',         # How you named your package folder (MyLib)
  packages = find_packages(where='.'),   # Chose the same as "name"
  version = find_version('scripts', '__init__.py'),      # Start with a small number and increase it with every change you make
  license='BY-NC-SA Creative Commons License',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Generates full genome images from a gbff file',   # Give a short description about your library
  author = 'Cumsille A. et al.',                   # Type in your name
  author_email = 'vicente.saona@sansano.usm.cl',      # Type in your E-Mail
  url = 'https://github.com/robotoD/GenoVi',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/robotoD/GenoVi/archive/refs/tags/0.1.6.tar.gz',    # I explain this later on
  package_data = {'genovi': ['scripts/dataset/cog-20.def.tab']},
  include_package_data = True,
  keywords = ['bioinformatics', 'genomics', 'CIRCOS'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'cairosvg>=2.5.2',
          'matplotlib>=3.5.2',
          'numpy>=1.20.2',
          'pandas>=1.2.4',
          'BioPython>=1.79',
          'argparse',
          'deepnog>=1.2.3',
          'scikit-learn',
          'torch>=1.2.0',
          'tqdm>=4.35.0',
          'Pillow'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Visualization',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: MIT License',   # Again, pick a license
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
    'Natural Language :: English',

  ],
  entry_points = {
              'console_scripts': [
                  'genovi=scripts.GenoVi:main',
              ],
          },

  extras_require={
    'test': extra_test,
    'dev': extra_dev,
  },
)
