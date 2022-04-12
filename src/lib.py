import ssl

repo = "cachome/ensembl-paralogs"
module = "interactions.py"

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

# Organisms configured for WikiPathways caching
organisms = [
    "Unspecified",
    "Acetobacterium woodii",
    "Anopheles gambiae",
    "Arabidopsis thaliana",
    "Bacillus subtilis",
    "Beta vulgaris",
    "Brassica napus",
    "Bos taurus",
    "Caenorhabditis elegans",
    "Canis familiaris",
    "Clostridium thermocellum",
    "Danio rerio",
    "Daphnia magna",
    "Daphnia pulex",
    "Drosophila melanogaster",
    "Escherichia coli",
    # "Equus caballus",
    "Gallus gallus",
    # "Glycine max",
    # "Gibberella zeae",
    "Homo sapiens",
    # "Hordeum vulgare",
    "Mus musculus",
    # "Mycobacterium tuberculosis",
    # "Oryza sativa",
    # "Pan troglodytes",
    # "Populus trichocarpa",
    "Rattus norvegicus",
    "Saccharomyces cerevisiae",
    "Solanum lycopersicum",
    "Sus scrofa",
    # "Vitis vinifera",
    "Xenopus tropicalis",
    # "Zea mays",
    "Plasmodium falciparum",

    "Macaca mulatta",
    "Macaca fascicularis",
    "Canis lupus familiaris",
    "Felis catus"
]
