from setuptools import setup, find_packages

setup(
    name="rocketchip",  
    version="2.1.0",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "rocketchip_pkg": [ 
            "rules/*",
            "tools/*",
            "tests/*",
            "tests/seq_data/*",
        ],
    },
    install_requires=[
        "pandas>=2.2.0",
        "pyyaml>=6.0.1",
        "snakemake>=7.32.4",
        "macs3>=3.0.0b3",
    ],
    entry_points={
        "console_scripts": [
            "rocketchip=rocketchip_pkg.rocketchip:main",
        ],
    },
    description="Rocketchip: chromatin immunoprecipitation data analysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Viktoria Haghani",
    author_email="vhaghani@health.ucdavis.edu",
    url="https://github.com/vhaghani26/rocketchip",
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
