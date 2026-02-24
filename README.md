[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18752672.svg)](https://doi.org/10.5281/zenodo.18752672)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/psaxena2/GlazeBalogaPlumeCode_PythonSaxena/HEAD?labpath=VolcanicPlumeCodePython_GlazeBaloga.ipynb)

# Glaze-Baloga Volcanic Plume Model: Python Implementation

Python implementation of the Glaze-Baloga integral volcanic plume model for simulating buoyant rise of 4-component plumes in arbitrary planetary atmospheres.

## Description

Based on the original IDL procedure by Lori Glaze and Stephen Baloga, this model solves five coupled differential equations describing volcanic plume dynamics (density, radius, velocity, temperature) using the Runge-Kutta method.

**Key Features:**
- 4-component plumes (solids, dry air, water vapor, liquid water)
- Arbitrary atmospheric compositions (CO₂, N₂, Earth-like, etc.)
- Neutral buoyancy height (NBH) calculation using sign change method
- Water vapor condensation physics
- VEI 2-8 parameterization
- Reads GCM atmospheric profiles
- Interactive Jupyter notebook with documentation

## Installation

```bash
# Clone repository
git clone https://github.com/[YOUR-USERNAME]/glaze-baloga-plume-model.git
cd glaze-baloga-plume-model

# Install dependencies
pip install -r requirements.txt

# Launch Jupyter
jupyter notebook GlazeBalogaPlumeModel.ipynb
