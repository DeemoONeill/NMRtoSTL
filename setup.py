from setuptools import setup

setup(
    name="NMRtoSTL",
    version="0.0.1",
    py_modules=[
        "NMRtoSTL/__init__",
        "NMRtoSTL/__main__",
        "NMRtoSTL/NMRtoSTL",
        "NMRtoSTL/importNMR",
    ],
    authors=["Matthew O'Neill", "Andrew Hall"],
    author_email="34174223+DeemoONeill@users.noreply.github.com",
    description="""A library and CLI tool for stacked 1d and 2d NMR into STL files for 3D printing""",
    install_requires=[dep for dep in open("requirements.txt").readlines()],
)
