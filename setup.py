""" Setup for brawn """

import os
from setuptools import setup, find_packages

with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


test_requirements = {
    'testing': [
        'coverage',
        'mypy == 0.982',
        'pytest >= 7.2.0, < 8',
        'pylint >= 2.15.2',
    ],
}


def get_version(path: str) -> str:
    """ Grabs the version from the package in order to have a single source of truth """
    with open(path, encoding="utf-8") as handle:
        for line in handle.read().splitlines():
            if line.startswith("VERSION"):
                return line.split()[-1].strip('"')
    raise RuntimeError("Unable to find version string")


setup(
    name='brawn',
    python_requires='>=3.9',
    version=get_version("brawn/constants.py"),
    author='Simon Shaw',
    author_email='SJShaw@users.noreply.github.com',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
    ],
    description='A tool for handling repetitive insertions into sequence alignments',
    entry_points={
        'console_scripts': [
            'brawn=brawn.__main__:entrypoint',
        ],
    },
    extras_require=test_requirements,
    license='GPLv3+',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(exclude=["tests"]),
    tests_require=test_requirements["testing"],
    url='https://github.com/SJShaw/brawn',
)
