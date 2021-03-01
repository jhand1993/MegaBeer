from setuptools import setup, find_packages
import pathlib
from MegaBeer.__init__ import __version__

base_path = pathlib.Path(__file__).parent

readme = (base_path / 'README.md').read_text()

setup(
    name='MegaBeer',
    version=__version__,
    description='Open source brewing tools tools',
    url='https://github.com/jhand1993/MegaBeer',
    author='Jared Hand',
    author_email='jared.hand1993@gmail.com',
    license='GPL V3',
    packages=find_packages(exclude=('tests')),
    zip_safe=False,
    install_requires=['numpy', 'scipy'],
    long_description=readme,
    long_description_content_type="text/markdown",
    include_package_data=True
)