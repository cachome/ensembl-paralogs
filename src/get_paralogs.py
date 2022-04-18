"""Fetch paralogs from Ensembl

Example usage:
python3 src/get_paralogs.py --reuse
"""

import argparse
import os
import re
from time import sleep
import json as ljson
import gzip
import csv

import requests

from lib import repo, module, ctx, organisms

# # Enable importing local modules when directly calling as script
# if __name__ == "__main__":
#     cur_dir = os.path.join(os.path.dirname(__file__))
#     sys.path.append(cur_dir + "/..")

# from lib import download_gzip

def fetch_genes(organism):
    """List genes symbols for an organism
    """

    ids_by_gene = {}
    genes = []
    prefix = ''
    # E.g. https://raw.githubusercontent.com/eweitz/ideogram/master/dist/data/cache/homo-sapiens-genes.tsv
    genes_url = (
        "https://raw.githubusercontent.com/eweitz/ideogram/"
        f"master/dist/data/cache/{slug(organism)}-genes.tsv.gz"
    )
    print('genes_url')
    print(genes_url)
    response = requests.get(genes_url)

    if response.status_code != 200:
        print(f"Status code {response.status_code} for {genes_url}")
        return [genes, prefix]

    tsv_string = gzip.decompress(response.content).decode('utf-8')
    # print('tsv_string')
    # print(tsv_string)
    reader = csv.reader(tsv_string.splitlines(), delimiter="\t")
    for row in reader:
        if row[0][0] == '#':
            if '## prefix' in row[0]:
                prefix = row[0].split('prefix: ')[1]
                continue
        if len(row) < 4: continue
        # row: [chr, start, length, slim_id, symbol, description]
        gene = row[4] # more formally, gene symbol
        id = row[3]
        ids_by_gene[gene] = id
        genes.append(gene)

    return [genes, prefix, ids_by_gene]

def slug(value):
    return value.lower().replace(" ", "-")

def lossy_optimize_paralogs(json_str):
    json = ljson.loads(json_str)

    trimmed_ids = []
    for h in json["data"][0]["homologies"]:
        trimmed_id = re.sub(r'[A-Za-z]+0+', '', h["id"])
        trimmed_ids.append(int(trimmed_id))
    return trimmed_ids


class EnsemblCache():

    def __init__(self, output_dir="data/", reuse=False):
        self.output_dir = output_dir
        self.tmp_dir = f"tmp/"
        self.reuse = reuse

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)

    def fetch_paralogs(self, organism, genes, gene_dir):
        """For given genes, retrieve paralogs from Ensembl; write each to file
        """

        prev_error_pwids = []
        error_pwids = []

        error_path = gene_dir + "errors.csv"
        if os.path.exists(error_path):
            with open(error_path) as f:
                prev_error_pwids = f.read().split(",")
                error_pwids = prev_error_pwids

        for gene in genes:
            json_path = gene_dir + gene + ".json"

            if "/" in gene:
                print(f"Skipping gene with slash in name: {gene}")
                error_pwids.append(gene)
                continue

            if self.reuse:
                if os.path.exists(json_path):
                    print(f"Found cache; skip processing {gene}")
                    continue
                elif id in prev_error_pwids:
                    print(f"Found previous error; skip processing {gene}")
                    continue

            url = (
                f"https://rest.ensembl.org"
                f"/homology/symbol/{slug(organism).replace('-', '_')}/{gene}"
                f"?format=condensed&type=paralogues"
                f"&content-type=application/json"
            )

            print('url')
            print(url)

            try:
                sleep(0.1)
                paralogs = requests.get(url).text
            except Exception as e:
                print(f"Encountered error when stringifying JSON for {gene}")
                error_pwids.append(gene)
                with open(error_path, "w") as f:
                    f.write(",".join(error_pwids))
                sleep(0.1)
                continue

            print("Preparing and writing " + json_path)

            with open(json_path, "w") as f:
                f.write(paralogs)

    def optimize_paralogs(self, genes, ids_by_gene, tmp_gene_dir, gene_dir):
        optimize_errors = []

        genes_by_paralogs = {}

        num_redundant_paralogs = 0

        rows = []
        for gene in genes:
        # for gene in genes[:50]: # Helpful for debugging

            # Disregard fusion genes and genes like "AC113554.1"
            if (
                "/" in gene or
                re.match(r'^(AC|AL|BX|AF)[0-9]+\.[0-9]$', gene) #or
                # re.match(r'^RNU[0-9]-', gene)
            ): continue

            # original_name = json_path.split("/")[-1]
            # gene = original_name.split(".json")[0]
            json_path = tmp_gene_dir + gene + '.json'
            # json_path = 'tmp/homo-sapiens/ACE2.json'
            if not os.path.exists(json_path):
                optimize_errors.append(gene)
                continue

            # The same genes are often capitalized differently in different
            # organisms.  We can leverage this to decrease cache size by
            # ~2x.  E.g. human "MTOR" and orthologous mouse "Mtor".
            gene = gene.upper()

            # pwid = re.search(r"WP\d+", name).group() # pathway ID
            optimized_tsv_path = gene_dir + gene + ".tsv"

            # repo_url = f"https://github.com/{repo}/tree/main/"
            # code_url = f"{repo_url}src/{module}"
            # data_url = f"{repo_url}{optimized_tsv_path}"
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

            if json.decode("utf-8") == '{"data": [homologies: []]}':
                print(f"Gene found, but no paralogs for {gene}")
                continue

            print(f"Optimizing to create: {optimized_tsv_path}")

            try:
                paralogs = lossy_optimize_paralogs(json)
                paralogs.sort()
                if len(paralogs) == 0:
                    continue

                gene_id = ids_by_gene[gene]
                paralogs.append(int(gene_id))
                tmp_ids = paralogs
                tmp_ids.sort()

                tsv = "\t".join([str(i) for i in paralogs])
                tmp_tsv = "\t".join([str(i) for i in tmp_ids])

                # Many genes share all their paralogs (minus themselves) with
                # other genes.  This accounts for that, and compresses the
                # combined paralog cache size by ~3x.
                if "\t" in tsv:
                    if tmp_tsv in genes_by_paralogs:
                        tsv = '_' + genes_by_paralogs[tmp_tsv]
                        num_redundant_paralogs += 1
                    else:
                        genes_by_paralogs[tmp_tsv] = gene

                # Omit current gene from its own paralog list
                split_tsv = tsv.split("\t")
                if (gene_id in split_tsv):
                    split_tsv.remove(gene_id)
                tsv = "\t".join(split_tsv)
                tsv = gene_id + "\t" + tsv

                tsv = tsv.encode()
            except Exception as e:
                handled = "Encountered error converting TSV for gene"
                handled2 = "not well-formed"
                if handled in str(e) or handled2 in str(e):
                    # print('Handled an error')
                    print(e)
                    optimize_errors.append(gene)
                    continue
                else:
                    print('Encountered fatal error')
                    print(e)
                    # raise Exception(e)
                    continue

            with open(optimized_tsv_path, "wb") as f:
                f.write(tsv)

            tsv = tsv.decode('utf-8')
            if len(tsv) > 0:
                rows.append(f"{gene}\t{tsv}")

        print('num_redundant_paralogs')
        print(num_redundant_paralogs)
        num_errors = len(optimize_errors)
        if num_errors > 0:
            print(f"{num_errors} pathways had optimization errors:")
            print(",".join(optimize_errors))

        return rows

    def populate_by_org(self, organism):
        """Fill caches for a configured organism
        """
        tmp_gene_dir = self.tmp_dir + slug(organism) + "/"
        if not os.path.exists(tmp_gene_dir):
            os.makedirs(tmp_gene_dir)

        gene_dir = self.output_dir  + slug(organism) + "/"
        if not os.path.exists(gene_dir):
            os.makedirs(gene_dir)

        # gpml_dir = self.output_dir +  + organism + /"
        # if not os.path.exists(gpml_dir):
        #     os.makedirs(gpml_dir)

        [genes, prefix, ids_by_gene] = fetch_genes(organism)
        # self.fetch_paralogs(organism, genes, tmp_gene_dir)
        print('len(genes)')
        print(len(genes))
        rows = self.optimize_paralogs(genes, ids_by_gene, tmp_gene_dir, gene_dir)

        combined_dir = self.output_dir + "combined/"
        if not os.path.exists(combined_dir):
            os.makedirs(combined_dir)
        combined_path = combined_dir + slug(organism) + "-paralogs.tsv.gz"
        rows.insert(0, f"## prefix: {prefix}")
        combined = "\n".join(rows).encode()
        compressed_combined = gzip.compress(combined)
        with open(combined_path, "wb") as f:
            f.write(compressed_combined)
        print(f"Wrote {combined_path}")

    def populate(self):
        """Fill caches for all configured organisms

        Consider parallelizing this.
        """
        # organisms = ["Homo sapiens", "Mus musculus"] # Comment out to use all
        organisms = ["Homo sapiens"] # Comment out to use all
        # organisms = ["Caenorhabditis elegans"] # Comment out to use all
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
        default="data/"
    )
    parser.add_argument(
        "--reuse",
        help=(
            "Whether to use previously-downloaded raw data"
        ),
        action="store_true"
    )
    args = parser.parse_args()
    output_dir = args.output_dir
    reuse = args.reuse

    EnsemblCache(output_dir, reuse).populate()
