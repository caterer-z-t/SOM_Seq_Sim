from setuptools import setup, find_packages

setup(
    name="SOM Seq Sim",  # Replace with your package name
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "seaborn",
        "minisom",
    ],
    extras_require={"dev": ["pytest", "pytest-cov", "flake8"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.12",
    author="Zachary Caterer, Victoria Hurd, Maddy Pernat, Con Muangkod",
    author_email="ztcaterer@colorado.edu",
    description="Self-Organizing Maps for Genetic Sequencing Simulation",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/caterer-z-t/SOM_Sew_Sim",
    project_urls={
        "Bug Tracker": "https://github.com/caterer-z-t/SOM-Seq_Sim/issues",
    },
)
