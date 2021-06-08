# Data Listing
This repository contains everything needed to list the data available for MCNP.

## Requirements
- Python 3.6+

### Third-Party Python Libraries
We utilize a few third-party libraries. These are all widely used and possibly already installed if you use a scientific Python distribution. 

- [Jupyter-lab](https://jupyterlab.readthedocs.io/en/stable/index.html)&mdash;[Project Jupyter's](https://jupyter.org) next-generation notebook interface.

   Even though jupyter-lab is "next-generation", it is completely stable (as of June 2019) and is the future interface for the Jupyter Project.

- [pandas](https://pandas.pydata.org) Data analysis and manipulation tool.
- [tqdm](https://tqdm.github.io/) A Fast, Extensible Progress Bar for Python and CLI
- [jupyterlab_widgets](https://pypi.org/project/jupyterlab-widgets/) Necessary to make nice looking progress bar in Jupyter Notebook. This one is not required.

If they are not already installed, you can easily install them using [`pip`](https://pip.pypa.io/en/stable/).

```shell

# pip install jupyterlab   # If you are seeing these instructions, you already have jupyterlab
pip install pandas
pip install tqdm
pip install jupyterlab_widgets
```

## Getting Started
See [Getting Started](GettingStarted.ipynb) to begin using this tool. After viewing that, see [Listing](Listing.md) to see what data is on your machine.

## Issues and Questions
If you have a question, please start a [discussion](https://github.com/NuclearData/DataListing/discussions/) with the developers.

If you discover a bug or deficiencies, please submit an [Issue](https://github.com/NuclearData/DataListing/issues).

## Contributing
We welcome [Pull Requests](https://github.com/NuclearData/DataListing/pulls) for improvements and bug fixes. Please start a [discussion](https://github.com/NuclearData/DataListing/discussions) with the developers before spending a lot of time on creating a new feature.

## Distribution
This package is distributed on GitHub at <https://github.com/NuclearData/DataListing>.
If you are concerned about what this does you are free to examine the code yourself.

## Copyright and License
This repository and associated code is governed by the [Copyright](Copyright) and [LICENSE](LICENSE) contained in the relevant files.
