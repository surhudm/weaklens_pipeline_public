#!/bin/bash
# Create a conda environment for weaklens_pipeline, build the pipeline, and register
# the environment as a Jupyter kernel.
#
# Usage (from anywhere):
#   bash notebooks/setup_weaklens_env.sh            # environment named "weaklens_pipeline_tutorial"
#   bash notebooks/setup_weaklens_env.sh myenv      # custom environment name
#   bash notebooks/setup_weaklens_env.sh /path/env  # environment at a path (conda -p)
set -e

ENV="${1:-weaklens_pipeline_tutorial}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ "$ENV" == */* ]]; then
    ENV_OPT=(-p "$ENV")
else
    ENV_OPT=(-n "$ENV")
fi

# Make "conda activate" usable inside this script
eval "$(conda shell.bash hook)"

echo "=== 1/4 Creating conda environment: $ENV"
conda create -y "${ENV_OPT[@]}" -c conda-forge --override-channels \
    python=3.10 "numpy<2" scipy astropy pandas pyyaml matplotlib fitsio \
    "swig=3.0.12" gsl "libboost-devel=1.84" gnupg ipykernel pip
conda activate "$ENV"

echo "=== 2/4 Installing Python requirements"
pip install -r "$REPO/requirements.txt"

echo "=== 3/4 Building weaklens_pipeline"
cd "$REPO"
rm -rf build
CPPFLAGS="-I$CONDA_PREFIX/include -DBOOST_TIMER_ENABLE_DEPRECATED" \
LDFLAGS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib" \
python setup.py install
cd /
if ! python -c "import weaklens, tomography" 2>/dev/null; then
    # The SWIG-generated .py files are sometimes not installed on the first pass (see README)
    cd "$REPO"
    CPPFLAGS="-I$CONDA_PREFIX/include -DBOOST_TIMER_ENABLE_DEPRECATED" \
    LDFLAGS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib" \
    python setup.py install
    cd /
fi
python -c "import weaklens, tomography; print('weaklens OK:', weaklens.__file__)"

echo "=== 4/4 Registering Jupyter kernel \"Python (weaklens_pipeline_tutorial)\""
KERNEL_OPT=(--user)
if [[ -n "$KERNEL_PREFIX" ]]; then
    KERNEL_OPT=(--prefix "$KERNEL_PREFIX")
fi
"$CONDA_PREFIX/bin/python" -m ipykernel install "${KERNEL_OPT[@]}" --name weaklens_pipeline_tutorial --display-name "Python (weaklens_pipeline_tutorial)"

echo
echo "Done. Activate with:  conda activate $ENV"
echo "In Jupyter / VS Code, select the kernel \"Python (weaklens_pipeline_tutorial)\"."
