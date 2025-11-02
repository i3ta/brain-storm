# CS 4641 Team 24

This is the repository for the Georgia Tech CS 4641 Fall 2025 project for team 24.

## Notes

This repository is currently a work in progress. Read below for more
information and installation instructions.

### Structure

Most of the data analysis pipelines are set up using Snakemake for its Slurm
integration capabilities, so a lot of the repository is structured around that.
All of the workflow-related scripts are contained within the `workflow/` folder,
and local development code are contained in the local package `daikon`, which is
in the `src/` folder.

In the `workflow/` folder, all the Snakemake workflows start from the
`Snakefile` file, which includes all the files inside the `rules/` folder.
Snakemake operates on rules that determine which scripts can produce which
files, and by chaining them together can produce the output files necessary. The
scripts that the rules call are all inside of the `scripts/` folder.

All of the results end up in the `results/` folder. The rules create a folder
inside of the `results/` folder corresponding to the name of each dataset that
is in the `<project root>/data/` folder. This folder is not commited to git
because it is too large. Locally, we can use a small subset of 50 images per
class for testing, and on PACE we can store our images in the `scratch` folder
(300GB capacity) and connect it to where the pipeline expects the folder to be
via a symlink.

> Note: The local development package is named `daikon` for no reason at all --
> I just decided that was a good name for it.

## Installation

### Setup environment

To create the conda environment, run the following command from the root of the
repository:

```sh
conda env create --file environments/environment.yml
```

Then, before running the scripts, run the following to activate the environment:

```sh
conda activate cs4641-project
```

This will install all dependencies, as well as set up the local package in edit
mode so you can edit code and import it in whatever scripts you want to run.

### Set up on PACE-ICE

Setting the repository to work on PACE is similar. Follow the below
instructions to get started. This setup takes a bit more time, but can be
beneficial if there is a lot of data to process.

1. Get access to PACE-ICE and use ssh to connect to PACE.
2. Activate Anaconda 3

   ```sh
   module load anaconda3
   ```

3. Download the dataset files into the scratch folder. This folder has a larger
   capacity than your home folder. Make sure the dataset folder is directly in
   the `scratch` folder (your dataset should be at
   `~/scratch/brain-tumor-dataset/`).
4. Clone this repository and `cd` into the root of the repository
5. Create a symlink to the scratch folder to serve as your data folder:

   ```sh
   ln -s ~/scratch/ data/
   ```

6. Create the PACE environment, which is different from the normal development
   environment:

   ```sh
   conda env create --file environment/pace.yml
   ```

7. Activate the environment

   ```sh
   conda activate pace
   ```

8. To run the snakemake pipeline, run snakemake as you normally would from the
   `workflow/` folder, but instead of calling `snakemake -c 1`, call
   `./setup_slurm.sh <output files>`

> For more information, check out [Intro to PACE ICE](https://github.com/guru-desh/Intro-To-PACE-ICE)

### Pages

To set up pages, you need to first install [bun](https://bun.sh/). Then, you
can run the following instructions to start website development:

1. Enter the `pages/` directory

   ```sh
   cd pages/
   ```

2. Install website dependencies

   ```sh
   bun ci
   ```

3. To run the website in development mode, run

   ```sh
   bun run dev
   ```

4. Once you have made your changes and want to push the website, run

   ```sh
   bun run deploy
   ```
