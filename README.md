Currently includes a script for downloading XAS data from the materials project.

Usage:

    python download.py [API key] [prefix] < [material_ids]

This will create two csv files:

[prefix]\_spectra.csv: Simulated XANES spectra for each element

[prefix]\_labels.csv: Labels to be predicted from the spectra (currently: information on the coordination environment)
