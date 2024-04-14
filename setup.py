import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="owlspec",
    version="0.1",
    author="Julian Held",
    author_email="julian.held@umn.edu",
    license='MIT',
    platforms=['any'],
    description="Library for optical emission spectroscopy of low-temperature plasmas.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mimurrayy/owl",
    packages=setuptools.find_packages(),
    install_requires=[
          'numpy>=1.17.0', 'scipy>=1.10.0', 'mendeleev>=0.12.0', 'astroquery>=0.4.5', 'roman>=4.0'
      ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
)
