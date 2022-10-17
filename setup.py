from setuptools import setup, find_packages
from py_sc_fermi import __version__ as VERSION

readme = "README.md"
with open(readme) as f:
    long_description = f.read()

scripts = ["sc_fermi_solve"]

setup(
    description="Self-consistent Fermi Analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Benjamin J. Morgan, Alex G. Squires",
    author_email="b.j.morgan@bath.ac.uk",
    url="https://github.com/bjmorgan/py-sc-fermi",
    download_url="https://github.com/bjmorgan/py-sc-fermi/archive/%s.tar.gz"
    % (VERSION),
    version=VERSION,
    install_requires=open("requirements.txt").read(),
    python_requires=">=3.8",
    license="MIT",
    packages=find_packages(exclude=["docs", "tests*"]),
    name="py-sc-fermi",
    entry_points={
        "console_scripts": [f"{s} = py_sc_fermi.cli.{s}:main" for s in scripts]
    },
)
