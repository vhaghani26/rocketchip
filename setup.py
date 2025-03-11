from setuptools import setup, find_packages

setup(
    name="rocketchip",  
    version="2.0.1",
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
        "pyyaml>=6.0.1",
        "snakemake>=7.32.4",
        "macs3>=3.0.0b3",
    ],
    entry_points={
        "console_scripts": [
            "rocketchip=rocketchip_pkg.rocketchip:main",
        ],
    },
    description="A pip-installable package for Rocketchip",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Viktoria Haghani",
    author_email="vhaghani@ucdavis.edu",
    url="https://github.com/vhaghani26/rocketchip",
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)