#!/usr/bin/env python

import yaml 
import argparse

parser = argparse.ArgumentParser(description="Retrieve information from meta.yml of each module for methods description.", epilog="Example usage: python methods.py -m methods_output.yml -n 'params'")
parser.add_argument("-m","--methods_output", help="Path to the methods output file", required=True)
parser.add_argument("-n", "--params", help="list of params obtained in the main.nf of the pipeline", required=True)
args = parser.parse_args()


with open (args.methods_output, 'r') as file:
    methods_output = yaml.safe_load(file)

with open (args.params, 'r') as file:
    params = yaml.safe_load(file)


tools = methods_output['tools']
workflow = methods_output['workflow']
pipeline_name = list(workflow.keys())[1]

print(tools)


genome = params.get('genome', 'not specified').split('/')[-1] if params.get('genome') else 'not specified'
genome_name = params.get('genome_name', 'not specified')
mirnas_name = params.get('mirnas_name', 'not specified')
trimgalore_pars = params['progPars'].get('trimgalore', 'not specified')
filter_unique =   "<p> non-unique reads have been filtered out </p>" if params.get('filter_unique') == "YES" else ''
shorttracks_pars = params['progPars'].get('shorttracks', 'not specified')



methods_description=f"""
id: "{ pipeline_name }-methods-description"
description: "Suggested text and references to use when describing pipeline usage within the methods section of a publication."
section_name: "{ pipeline_name } Methods Description"
section_href: "https://github.com/biocorecrg/{ pipeline_name }"
plot_type: "html"
## TODO nf-core: Update the HTML below to your preferred methods description, e.g. add publication citation for this pipeline
## You inject any metadata in the Nextflow '${workflow}' object
data: |
  <h4>Methods</h4>
  <p>The quality of raw reads was assessed using tool <a href="{tools['fastqc']['homepage']}">FastQC</a> version {tools['fastqc']['version']}.</p>
  <p>Adapter sequences and low-quality bases were subsequently removed with <a href="{tools['trimgalore']['homepage']}">TrimGalore</a> version {tools['trimgalore']['version']}, employing the parameters {trimgalore_pars}.</p>
  <p>The processed reads were then aligned to the reference genome {genome} using ShortStack.</p>
  <p>Aligned reads were annotated against miRBase and {genome_name} databases.</p>
  {filter_unique}
  <p>Coverage profiles were generated using ShortTracks with the {shorttracks_pars} readgroup options.</p>
  <p>Finally, PCA analysis was performed using DESeq2 on read counts after applying a variance-stabilizing transformation.</p>
  <p>The pipeline was executed with <a href="https://doi.org/10.1038/nbt.3820">Nextflow</a> 2017 {workflow['Nextflow']} with the following command:</p>
  <pre><code>${{workflow.commandLine}}</code></pre>
  <h4>References</h4>
  <ul>
    <li>Babraham Bioinformatics - FastQC A Quality Control tool for High Throughput Sequence Data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/. Accedido 30 de septiembre de 2025.</li>
    <li>{tools['trimgalore']['citation']}</li>
    <li>{tools['shortstack']['citation']}</li>
    <li>Kozomara, Ana, et al. «miRBase: From microRNA Sequences to Function». Nucleic Acids Research, vol. 47, n.o D1, enero de 2019, pp. D155-62. PubMed, https://doi.org/10.1093/nar/gky1141.</li>
    <li>Michael Love, Simon Anders. DESeq2. Bioconductor, 2017. DOI.org (Datacite), https://doi.org/10.18129/B9.BIOC.DESEQ2.</li>
    <li>Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319. doi: <a href="https://doi.org/10.1038/nbt.3820">10.1038/nbt.3820</a></li>
    <li>Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M. U., Di Tommaso, P., & Nahnsen, S. (2020). The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology, 38(3), 276-278. doi: <a href="https://doi.org/10.1038/s41587-020-0439-x">10.1038/s41587-020-0439-x</a></li>
  </ul>
  <div class="alert alert-info">
    <h5>Notes:</h5>
  <p>The pipeline was executed with Nextflow {workflow['Nextflow']}  (<a href="https://doi.org/10.1038/nbt.3820">Di Tommaso <em>et al.</em>, 2017</a>) with the following command:</p>
  <pre><code>${{workflow.commandLine}}</code></pre>
    <ul>
      <li>The command above does not include parameters contained in any configs or profiles that may have been used. Ensure the config file is also uploaded with your publication!</li>
      <li>You should also cite all software used within this run. Check the "Software Versions" of this report to get version information.</li>
    </ul>
  </div>

  """


with open('methods_description_mqc.yml','w') as file:
    file.write(methods_description)
