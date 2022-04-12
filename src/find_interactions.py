import argparse
import glob
import os
import re
from time import sleep
import json as ljson
import gzip
import csv
from typing import Set

import requests
from lxml import etree

from lib import repo, module, ctx, organisms

# # Enable importing local modules when directly calling as script
# if __name__ == "__main__":
#     cur_dir = os.path.join(os.path.dirname(__file__))
#     sys.path.append(cur_dir + "/..")

# from lib import download_gzip

def get_pathway_ids_and_names(organism):
    base_url = "https://webservice.wikipathways.org/listPathways"
    params = f"?organism={organism}&format=json"
    url = base_url + params
    response = requests.get(url)
    data = response.json()
    ids_and_names = [[pw['id'], pw['name']] for pw in data['pathways']]
    return ids_and_names

def fetch_pathway_genes(data_dir, organism):
    """List genes symbols that are also TextLabels in WikiPathways
    """
    pathway_genes = []

    genes = []
    # E.g. https://raw.githubusercontent.com/eweitz/ideogram/master/dist/data/cache/homo-sapiens-genes.tsv
    genes_url = (
        "https://raw.githubusercontent.com/eweitz/ideogram/"
        f"master/dist/data/cache/{slug(organism)}-genes.tsv"
    )
    tsv_string = requests.get(genes_url).content.decode('utf-8')
    reader = csv.reader(tsv_string.splitlines(), delimiter="\t")
    for row in reader:
        if len(row) < 4 or row[0][0] == '#': continue
        # row: [chr, start, length, slim_id, symbol]
        gene = row[4] # more formally, gene symbol
        genes.append(gene)

    labels = get_gpml_labels(data_dir)

    for gene in genes:
        if gene in labels:
            pathway_genes.append(gene)

    print(f"Found {len(pathway_genes)} {organism} genes in WikiPathways")
    return pathway_genes

def get_gpml_labels(data_dir):
    print("Get GPML labels")
    labels = set()
    for gpml_path in glob.glob(f'{data_dir}/**/*.xml.gz'):
        # print('gpml_path')
        # print(gpml_path)
        with gzip.open(gpml_path, 'rb') as f:
            xml = f.read()
        # xml = gzip.decompress(xml_gz)
        tree = etree.fromstring(xml)
        elements = tree.xpath('//*')
        # print('len(elements)')
        # print(len(elements))
        for el in elements:
            if "TextLabel" in el.attrib:
                labels.add(el.attrib["TextLabel"])

    print(f"Found {len(labels)} labels in compressed GPML")
    return labels

def slug(value):
    return value.lower().replace(" ", "-")

def lossy_optimize_paralogs(json_str, gene):
    json = ljson.loads(json_str)
    # print('json')
    # print(json)
    for result in json["result"]:
        # print('result')
        # print(result)
        del result["score"]
        del result["url"]
        del result["revision"]
        del result["fields"]["indexerId"]
        del result["fields"]["source"]
        for field in result["fields"]:
            del result["fields"][field]["name"]
    # print('json after')
    # print(json)
    # print('\n\n')
    return ljson.dumps(json)


class WikiPathwaysCache():

    def __init__(self, output_dir="data/ixns_by_gene/", reuse=False):
        self.output_dir = output_dir
        self.tmp_dir = f"tmp/"
        self.reuse = reuse

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)

    def fetch_paralogs(self, genes, gene_dir):

        prev_error_pwids = []
        error_pwids = []

        error_path = gene_dir + "error_pwids.csv"
        if os.path.exists(error_path):
            with open(error_path) as f:
                prev_error_pwids = f.read().split(",")
                error_pwids = prev_error_pwids

        for gene in genes:
            json_path = gene_dir + gene + ".json"

            if self.reuse:
                if os.path.exists(json_path):
                    print(f"Found cache; skip processing {gene}")
                    continue
                elif id in prev_error_pwids:
                    print(f"Found previous error; skip processing {gene}")
                    continue

            # url = f"https://www.wikipathways.org/index.php/Pathway:{id}?view=widget"
            # base_url = "https://www.wikipathways.org//wpi/wpi.php"
            # url = f"{base_url}?action=downloadFile&type=gpml&pwTitle=Pathway:{id}"

            url = (
                "https://webservice.wikipathways.org/findInteractions"
                f"?query={gene}&format=json"
            )

            try:
                sleep(0.5)
                interactions = requests.get(url).text
            except Exception as e:
                print(f"Encountered error when stringifying JSON for {gene}")
                error_pwids.append(gene)
                with open(error_path, "w") as f:
                    f.write(",".join(error_pwids))
                sleep(0.5)
                continue

            print("Preparing and writing " + json_path)

            with open(json_path, "w") as f:
                f.write(interactions)

    def optimize_paralogs(self, gene_dir):

        optimize_error_pwids = []

        for json_path in glob.glob(f'{gene_dir}*.json'):
        # for json_path in ["tmp/homo-sapiens/WP231.json"]: # debug

            # print('gpml')
            # print(gpml)
            original_name = json_path.split("/")[-1]
            gene = original_name.split(".json")[0]
            # pwid = re.search(r"WP\d+", name).group() # pathway ID
            optimized_json_path = self.output_dir + gene + ".json.gz"

            # repo_url = f"https://github.com/{repo}/tree/main/"
            # code_url = f"{repo_url}src/{module}"
            # data_url = f"{repo_url}{optimized_json_path}"
            # wp_url = f"https://www.wikipathways.org/index.php/Pathway:{pwid}"
            # provenance = "\n".join([
            #     "<!--",
            #     f"  WikiPathways page: {wp_url}",
            #     f"  URL for this compressed file: {data_url}",
            #     # f"  Uncompressed GPML file: {original_name}",
            #     # f"  From upstream ZIP archive: {url}",
            #     f"  Source code for compression: {code_url}",
            #     "-->"
            # ])

            with open(json_path, 'rb') as f:
                json = f.read()

            if json.decode("utf-8") == '{"result":[]}':
                # print(f"Gene found, but no interactions for {gene}")
                continue

            print(f"Optimizing to create: {optimized_json_path}")

            try:
                json = lossy_optimize_paralogs(json, gene).encode()
                # json = lossless_optimize_paralogs(json, gene)
                json = gzip.compress(json)

            except Exception as e:
                handled = "Encountered error converting XML for pathway"
                handled2 = "not well-formed"
                if handled in str(e) or handled2 in str(e):
                    # print('Handled an error')
                    print(e)
                    optimize_error_pwids.append(gene)
                    continue
                else:
                    print('Encountered fatal error')
                    print(e)
                    # raise Exception(e)
                    continue

            with open(optimized_json_path, "wb") as f:
                f.write(json)

            # with open(optimized_json_path, "w") as f:
            #     f.write(json)

        num_errors = len(optimize_error_pwids)
        if num_errors > 0:
            print(f"{num_errors} pathways had optimization errors:")
            print(",".join(optimize_error_pwids))

    def populate_by_org(self, organism):
        """Fill caches for a configured organism
        """
        gene_dir = self.tmp_dir + "gene/"
        if not os.path.exists(gene_dir):
            os.makedirs(gene_dir)

        genes = fetch_pathway_genes(self.output_dir, organism)
        # self.fetch_intera ctions(genes, gene_dir)
        self.optimize_paralogs(gene_dir)

    def populate(self):
        """Fill caches for all configured organisms

        Consider parallelizing this.
        """
        organisms = ["Homo sapiens"] # Comment out to use all
        for organism in organisms:
            self.populate_by_org(organism)

# Command-line handler
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--output-dir",
        help=(
            "Directory to put outcome data.  (default: %(default))"
        ),
        default="data/gene/"
    )
    parser.add_argument(
        "--reuse",
        help=(
            "Whether to use previously-downloaded raw GPML zip archives"
        ),
        action="store_true"
    )
    args = parser.parse_args()
    output_dir = args.output_dir
    reuse = args.reuse

    WikiPathwaysCache(output_dir, reuse).populate()
