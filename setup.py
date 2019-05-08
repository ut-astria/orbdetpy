from setuptools import setup
from distutils.util import convert_path
from os import path, environ


here = path.abspath(path.dirname(__file__))
version_dict = {}
version_path = convert_path('orbdetpy/_version.py')

# Get the version number from the version path.
with open(version_path, 'r') as ver_file:
    exec(ver_file.read(), version_dict)

# Get the long description from the README file.
with open(path.join(here, 'README.rst')) as a_file:
    long_description = a_file.read()

with open('requirements.txt', 'r') as a_file:
    requirements = a_file.read()

CLASSIFIERS = ['Development Status :: 5 - Production/Stable',
               'Environment :: Console',
               'Intended Audience :: Developers',
               'Intended Audience :: Science/Research',
               'License :: OSI Approved :: GNU General Public License (GPL)',
               'Natural Language :: English',
               'Operating System :: OS Independent',
               'Programming Language :: Python',
               'Programming Language :: Python :: 3',
               'Programming Language :: Python :: 3.3',
               'Programming Language :: Python :: 3.4',
               'Programming Language :: Python :: 3.5',
               'Programming Language :: Python :: 3.6',
               'Programming Language :: Java',
               'Topic :: Scientific/Engineering',
               'Topic :: Scientific/Engineering :: Astronomy',
               'Topic :: Utilities']
PROJECT_URLS = {'Documentation': 'https://github.com/ut-astria/orbdetpy',
                'Source': 'https://github.com/ut-astria/orbdetpy',
                'Tracker': 'https://github.com/ut-astria/orbdetpy/issues'}
PACKAGE_DATA = {'orbdetpy': ['data/*', 'data/**/*', 'java/*', 'java/**/*',
                             'lib/*', 'lib/**/*']}

setup(name='orbdetpy',
      packages=['orbdetpy'],
      version=version_dict['__version__'],
      description='Orbit determination routines for Python',
      long_description=long_description,
      long_description_content_type='text/x-rst',
      author='Shiva Iyer',
      author_email='shiva.iyer@utexas.edu',
      url='https://github.com/ut-astria/orbdetpy',
      project_urls=PROJECT_URLS,
      keywords=['orbit_determination utilities orbit space'],
      classifiers=CLASSIFIERS,
      package_data=PACKAGE_DATA,
      install_requires=requirements)
