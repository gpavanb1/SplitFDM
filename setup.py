from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="SplitFDM",
    version="0.1.3",
    description="1D Finite-Difference Split Newton Solver",
    url="https://github.com/gpavanb1/SplitFDM",
    author="gpavanb1",
    author_email="gpavanb@gmail.com",
    license="MIT",
    packages=["splitfdm", "splitfdm.equations"],
    install_requires=["numpy", "numdifftools", "matplotlib", "splitnewton"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Topic :: Scientific/Engineering :: Mathematics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords="amr newton python finite-difference armijo optimization pseudotransient splitting",
    project_urls={  # Optional
        "Bug Reports": "https://github.com/gpavanb1/SplitFDM/issues",
        "Source": "https://github.com/gpavanb1/SplitFDM/",
    },
    zip_safe=False,
)
