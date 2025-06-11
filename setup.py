from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="resolutiontree",
    version="0.1.1",
    author="Joe Hou",
    author_email="joseph.houjue@gmail..com",
    description="Systematic exploration of clustering resolutions in single-cell analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/joe-jhou2/resolutiontree",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
    python_requires=">=3.8",
    install_requires = [
        "numpy>=1.19.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "scanpy>=1.8.0",
        "anndata>=0.8.0",
        "igraph>=0.9.0",
        "networkx>=2.6.0",
        "leidenalg>=0.8.0",
    ], 
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
            "sphinx",
        ],
    },
    keywords="single-cell, clustering, resolution, scanpy, leiden, visualization",
)