import sys
import pandas as pd
from pandas import DataFrame
from pymatgen.ext.matproj import MPRester
from pymatgen.core.structure import Structure
from pymatgen.analysis.xas.spectrum import XAS
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy

api_key = sys.argv[1]
out_prefix = sys.argv[2]

downloader = MPRester(api_key)

def record_to_spectrum(record):
    '''Converts an XAS record downloaded from the materials project using
    get_xas_spectrum to an XAS object'''
    s = record["spectrum"]
    return XAS(s["x"], s["y"], Structure.from_dict(s["structure"]),
               s["absorption_specie"], s["edge"], "XANES")

def add_descriptors(structure):
    '''Add coordination environment data to a pymatgen structure object's site
    properties. Return True if successful, false otherwise.'''
    lgf = LocalGeometryFinder()
    lgf.setup_structure(structure=structure)
    # Structure environment
    se = lgf.compute_structure_environments(maximum_distance_factor=1.41, only_cations=False)
    strategy = SimplestChemenvStrategy(se, distance_cutoff=1.4, angle_cutoff=0.3)
    # Coordination environments
    ces = [strategy.get_site_coordination_environment(site) \
            for site in structure.sites]
    if any(ce is None for ce in ces):
        return False
    # Extract descriptors from the coordination environments
    symbols = [ce[0] for ce in ces]
    sizes = [ce[1]["scaling_factor"] for ce in ces]
    # Add them to the structure object
    structure.add_site_property("symbol", symbols)
    structure.add_site_property("size", sizes)
    return True

def structure_to_table(structure):
    '''Given a pymatgen crystal structure object, convert the site properties
    into a table, represented as a pandas dataframe'''
    site_properties = structure.site_properties
    site_properties["element"] = [element.symbol for element in structure.species]
    return DataFrame(site_properties)
    
def spectrum_to_table(spectrum):
    '''Given a pymatgen XAS spectrum object, convert to a table, represented as
    a pandas dataframe'''
    return DataFrame({"element": [spectrum.absorbing_element] * len(spectrum.x),
            "x": spectrum.x, "y": spectrum.y})

structure_list = []
spectrum_list = []
# Loop over materials ids, stripping whitespace (in particular line breaks)
for material_id in [i.strip() for i in sys.stdin]:
    # Download structure 
    try:
        structure = downloader.get_entry_by_material_id(\
                    material_id,
                    inc_structure = "final",
                    conventional_unit_cell = False).structure
        structure_download_success = True
    except:
        print("Could not retrieve a structure for {} from database".format(material_id))
        structure_download_success = False

    # Download spectrum for each site and add them to the table
    spectrum_download_success = True
    for element in [element.symbol for element in structure.species]:
        spectrum_buffer = []
        try:
            spectrum_record = downloader.get_xas_data(material_id, element)
        except:
            print("Could not retrieve an XAS spectrum for {} from database".format(material_id))
            spectrum_download_success = False
            break
        # Sometimes spectrum construction fails due to invalid values in the spectrum
        try:
            spectrum = record_to_spectrum(spectrum_record)
        except ValueError as err:
            print("Spectrum for {} element {} is invalid. Error:".\
                    format(material_id, element))
            print(err)
            spectrum_download_success = False
            break
        one_spectrum_table = spectrum_to_table(spectrum)
        one_spectrum_table["id"] = material_id
        spectrum_buffer.append(one_spectrum_table)
    if spectrum_download_success:
        spectrum_list += spectrum_buffer

    # Skip coordination environment if we couldn't get a spectrum, since it
    # takes a long time
    if structure_download_success and spectrum_download_success:
        # Annotate the structure with labels
        annotation_success = add_descriptors(structure)
        if annotation_success:
            # Add the structure to the table
            one_label_table = structure_to_table(structure)
            one_label_table["id"] = material_id
            structure_list.append(one_label_table)
        else:
            print("Could not calculate coordination environment for {}".format(material_id))


label_table = pd.concat(structure_list)
spectrum_table = pd.concat(spectrum_list)

label_table.to_csv(out_prefix + "_labels.csv", index = False)
spectrum_table.to_csv(out_prefix + "_spectra.csv", index = False)
