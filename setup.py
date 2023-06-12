#!/usr/bin/env python3
"""Category-wide association study (CWAS)

This is a data analytic tool to perform stringent association tests
to find non-coding loci associated with autism spectrum disorder (ASD).
"""
import setuptools

setuptools.setup(
    name="cwas",
    version="1.0.0",
    license="MIT",
    author="Minwoo Jeong",
    author_email="jeongmwj@gmail.com",
    description=__doc__,
    url="https://github.com/mwjjeong/cwas",
    packages=setuptools.find_packages(),
    package_data={'': ['*.r', '*.R']},
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
    ],
    python_requires='>=3.7',
    entry_points={
        'console_scripts': [
            'cwas = cwas.__main__:main'
        ]
    }
)
