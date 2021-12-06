from setuptools import setup, find_packages
setup(
  name = 'GenoVi',         # How you named your package folder (MyLib)
  packages = find_packages(),   # Chose the same as "name"
  version = '0.1.6',      # Start with a small number and increase it with every change you make
  license='BY-NC-SA Creative Commons License',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Generates full genome images from a gbff file',   # Give a short description about your library
  author = 'Cumsille A. et al.',                   # Type in your name
  author_email = 'vicente.saona@sansano.usm.cl',      # Type in your E-Mail
  url = 'https://github.com/robotoD/GenoVi',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/robotoD/GenoVi/archive/refs/tags/0.1.6.tar.gz',    # I explain this later on
  keywords = ['bioinformatics', 'genomics', 'CIRCOS'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'cairosvg',
          'numpy',
          'pandas',
          'BioPython',
          'argparse'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: MIT License',   # Again, pick a license
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
  ],
)