from setuptools import setup, find_packages
from py_sc_fermi import __version__ as VERSION

readme = 'README.md'
long_description = open( readme ).read()

config = {
    'description': 'Self-consistent Fermi Analysis',
    'long_description': long_description,
    'long_description_content_type': 'text/markdown',
    'author': 'Benjamin J. Morgan',
    'author_email': 'b.j.morgan@bath.ac.uk',
    'url': 'https://github.com/bjmorgan/py-sc-fermi',
    'download_url': "https://github.com/bjmorgan/py-sc-fermi/archive/%s.tar.gz" % (VERSION),
    'author_email': 'b.j.morgan@bath.ac.uk',
    'version': VERSION,
    'install_requires': open('requirements.txt').read(),
    'python_requires': '>=3.5',
    'license': 'MIT',
    'packages': [ 'py_sc_fermi' ],
    'scripts': [],
    'name': 'py-sc-fermi'
}

setup(**config)
