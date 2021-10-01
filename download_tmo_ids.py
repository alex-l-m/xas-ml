from pymatgen.ext.matproj import MPRester

api_key = "QfuInHqVEyO8BoWa"
downloader = MPRester(api_key)

transition_metal_oxides = downloader.get_entries({"elements":{"$in":["Co", "Cr", "Cu", "Fe", "Mn", "Ni", "Sc", "Ti", "V", "Zn"], "$all": ["O"]}, "nelements":2})

for i in transition_metal_oxides:
    print(i.entry_id)
