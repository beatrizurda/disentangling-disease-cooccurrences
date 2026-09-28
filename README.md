<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

# Disentangling the genetic and non-genetic origin of disease co-occurrences
Beatriz Urda-García<a href="https://orcid.org/0000-0002-3845-5751">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" />
Davide Cirillo<a href="https://orcid.org/0000-0002-3845-5751">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" />
</a>, Alfonso Valencia<a href="https://orcid.org/0000-0002-8937-6789">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" />
</a>

medRxiv: <a href="https://https://doi.org/10.1101/2021.07.22.21260979">https://doi.org/10.1101/2025.10.16.25338165</a>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#manuscript">Manuscript</a>
    </li>
    <li><a href="#frequently-used-terms">Frequently used terms</a></li>
    <li><a href="#code">Code</a></li>
  </ol>
</details>




<!-- MANUSCRIPT INFORMATION -->
## Manuscript

Preprint available in medRxiv at <a href="https://https://doi.org/10.1101/2025.10.16.25338165">https://doi.org/10.1101/2025.10.16.25338165</a>

### Abstract
Numerous diseases co-occur more than expected by chance, likely due to a combination of genetic and environmental factors. However, the extent to which these influences shape disease relationships remain unclear. Here, we integrate large-scale RNA-seq data and heritability measures from human diseases with genomic data from the UK Biobank to disentangle the genetic and non-genetic origins of disease co-occurrences (DCs). Our findings show that gene expression not only recovers but also expands upon genomically explained DCs, capturing disease relationships beyond genetic variation. Approximately 60% of transcriptomically inferred DCs have a detectable genomic component, whereas the remaining 40% are not explained by known genomic layers, suggesting contributions from regulatory or environmental mechanisms. Consistent with this interpretation, the relative contributions of transcriptomics and genomics reconstruct disease etiology and correlate with comorbidity burden, revealing key aspects of disease mechanisms. Additionally, we find that diseases do not generally co-occur based on their heritability, except when sharing SNPs. However, highly heritable diseases tend to have genetically driven co-occurrences, even with lowly heritable diseases. In contrast, transcriptomics explains DCs regardless of heritability, at least partly due to non-heritable mechanisms, such as regulatory or environmental. Integrating transcriptomic and genomic data provides near-complete coverage of DCs among the analyzed diseases, with a considerable portion likely rooted in factors beyond DNA sequence and, therefore, potentially modifiable.</a>. 


## Frequently used terms
- **Disease co-occurrences (DCs):** pair of diseases that tend to co-occur in the same person more than expected by chance (comorbidities and multumorbidities). 
- **Disease Similarity Network (DSN):** disease-disease network. Two diseases are connected if their gene expression profiles correlate significantly (positively or negatively).
- **Meta-patients:** groups of patients from a given disease with a similar gene expression profile. 
- **Stratified Similarity Network (SSN):** extension of the DSN network that includes the meta-patients. Hence, the SSN contains three types of links: 
      (1) links between diseases
      (2) links between diseases and meta-patients
      (3) links between meta-patients.


<!-- CODE -->
## Code

### Data Preparation
First, uniformly processed gene counts were dowloaded from the <a href="http://www.ilincs.org/apps/grein/">GREIN platform</a> and Summarized Emperiment Objects (SE) from Bioconductor were constructed for each study. Finally, SE objects corresponding to the same disease were merged (<code>merge_same_disease_se_objects.R</code>).

### RNA-seq pipeline
Then, we applied an RNA-seq pipeline to each disease separately and in parallel (<code>run_rnaseq_pipeline_for_disease.R</code>). Then, we clustered diseases based on their significantly dysregulated pathways (<code>molecular_insight_heatmap.R</code>).

### DSN generation
1. First, we computed distances between diseases (<code>Network_building/build_disease_level_network.py</code>)
2. We used the generated distances to obtain the Disease Similarity Network (DSN) (<code>generating_networks.R</code>)

### DSN overlap
1. Obtain the SE object for each icd9 (<code>generating_SE_objects_icd9_level.R</code>)
2. Run the RNA-seq pipeline at the icd9 level (<code>run_rnaseq_pipeline_for_disease.R</code>)
3. Compute distances between diseases (<code>build_ICD_level_network.R</code>)
4. Use the obtained distances to generate the ICD9 level DSN network (<code>generating_networks.R</code>)
5. Compute the overlap with the epidemiological network from Hidalgo et al. and Dong et al. (<code>network_overlap_icd.R</code>)

### Meta-patient definition and characterization. 
1. We used PAM and WARD algorithms to define meta-patients for each disease (groups of patients with a similar expression profile) (<code>defining_meta_patients.R</code>)
2. Then, we applied the RNA-seq pipeline for each meta-patient (<code>DEanalysis_for_metapatients.R</code>)

### SSN generation, analysis and overlap computation
1. First, we computed distances between diseases (<code>Network_building/build_metapatient_disease_network.py</code>)
2. We used the generated distances to obtain the Disease Similarity Network (DSN) (<code>generating_networks.R</code>)
3. We computed the overlap of the DSN with the epidemiological network from Hidalgo et al. and Dong et al.(<code>network_overlap_SSN.R</code>)

