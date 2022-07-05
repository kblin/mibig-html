# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" MiBIG specific sideloading """

import json
from os import path
import os
from typing import Any, Dict, List, Optional
import logging

from xml.dom import minidom
import time

from antismash.common.errors import AntismashInputError
from antismash.common.module_results import DetectionResults
from antismash.common.secmet import CDSFeature, SubRegion, Record
from antismash.common.secmet.locations import (
    CompoundLocation,
    FeatureLocation,
    location_contains_other,
)

from urllib.error import HTTPError
from mibig.converters.read.top import Everything
from mibig_taxa import TaxonCache


class MibigAnnotations(DetectionResults):
    def __init__(self, record_id: str, area: SubRegion, data: Everything, cache_file: str, pubmed_cache_file: str) -> None:
        super().__init__(record_id)
        self.data = data  # holds the original annotation json data
        # save calculated loci (relative to record), not annotated ones
        self.record_id = record_id
        self.area = area

        cache = TaxonCache(cache_file)
        self.taxonomy = cache.get(int(data.cluster.ncbi_tax_id))

        self.pubmed_cache = PubmedCache(pubmed_cache_file)

    def get_predicted_subregions(self) -> List[SubRegion]:
        return [self.area]

    def to_json(self) -> Dict[str, Any]:
        # save only information critical for deciding reusability
        loci = self.data.cluster.loci
        annotations = []
        for annot in (self.data.cluster.genes.annotations if self.data.cluster.genes else []):
            annotations.append(annot.to_json())
        extra_genes = []
        for extra_gene in (self.data.cluster.genes.extra_genes if self.data.cluster.genes else []):
            extra_genes.append(extra_gene.to_json())
        return {
            "record_id": self.record_id,
            "genbank_accession": loci.accession,
            "coords": (loci.start or -1, loci.end or -1),
            "gene_annotations": annotations,
            "extra_genes": extra_genes,
        }

    @staticmethod
    def from_json(prev: Dict[str, Any], record: Record, annotations_file: str,
                  cache_file: str, pubmed_cache_file: str) -> Optional["MibigAnnotations"]:
        with open(annotations_file) as handle:
            raw = json.load(handle)
            data = Everything(raw)

        # compare old vs new annotation, decide if we can 'reuse'
        can_reuse = True
        loci = data.cluster.loci
        gene_annotations = data.cluster.genes.annotations if data.cluster.genes else []
        extra_genes = data.cluster.genes.extra_genes if data.cluster.genes else []
        if loci.accession != prev["genbank_accession"]:
            logging.debug("Previous result's genbank_accession is not the same as the new one")
            can_reuse = False
        elif record.id != prev["record_id"]:
            logging.debug("Previous result's record_id is not the same as the new one")
            can_reuse = False
        elif (loci.start or -1) != prev["coords"][0] or (loci.end or -1) != prev["coords"][1]:
            logging.debug("Previous result's start/end coordinate is not the same as the new one")
            can_reuse = False
        elif len(gene_annotations) != len(prev["gene_annotations"]):
            logging.debug("Gene annotations have changed")
            can_reuse = False
        elif len(extra_genes) != len(prev["extra_genes"]):
            logging.debug("Additional genes have changed")
            can_reuse = False

        # if we can't reuse, stop running antismash, because CDS annotations won't be correct
        if can_reuse:
            product = ", ".join(data.cluster.biosynthetic_class)
            loci_region = FeatureLocation(
                loci.start - 1 if loci.start else 0,
                loci.end or len(record.seq)
            )
            area = SubRegion(loci_region, tool="mibig", label=product)
            return MibigAnnotations(record.id, area, data, cache_file, pubmed_cache_file)
        else:
            logging.error("Can't reuse MIBiG annotation.")
            raise AntismashInputError(
                "Genbank record or gene annotations are updated, can't reuse result")


def mibig_loader(annotations_file: str, cache_file: str, pubmed_cache_file: str, record: Record) -> MibigAnnotations:
    """This method will be called only when not reusing data"""
    with open(annotations_file) as handle:
        raw = json.load(handle)
        data = Everything(raw)

    product = ", ".join(data.cluster.biosynthetic_class)
    loci = data.cluster.loci
    loci_region = FeatureLocation(
        loci.start - 1 if loci.start else 0,
        loci.end or len(record.seq)
    )
    area = SubRegion(loci_region, tool="mibig", label=product)

    # re-annotate CDSes in the Record
    if data.cluster.genes:
        # extra genes
        for gene in data.cluster.genes.extra_genes:
            if gene.id and gene.location:
                exons = [FeatureLocation(exon.start - 1, exon.end, strand=gene.location.strand)
                         for exon in gene.location.exons]
                location = CompoundLocation(exons) if len(exons) > 1 else exons[0]
                if not location_contains_other(area.location, location):
                    raise ValueError(f"additional gene {gene.id} lies outside cluster")
                translation = gene.translation or record.get_aa_translation_from_location(location)
                cds_feature = CDSFeature(
                    location=location, locus_tag=gene.id, translation=translation)
                record.add_cds_feature(cds_feature, auto_deduplicate=False)
                record.add_alteration(f"{cds_feature.get_name()} was added")
        # re-annotation
        for cds_feature in record.get_cds_features_within_location(area.location):
            locus_tag = cds_feature.locus_tag
            protein_id = cds_feature.protein_id
            name = cds_feature.gene
            for annot in (data.cluster.genes.annotations if data.cluster.genes else []):
                if locus_tag and annot.id == locus_tag:
                    pass
                elif protein_id and annot.id == protein_id:
                    pass
                elif name and annot.name == name:
                    pass
                else:
                    continue
                if annot.product:
                    cds_feature.product = annot.product

    existing = set()
    for cds in record.get_cds_features():
        for name in [cds.locus_tag, cds.protein_id, cds.gene]:
            if name is not None:
                existing.add(name)

    referenced = set()
    if data.cluster.genes and data.cluster.genes.annotations:
        for gene in data.cluster.genes.annotations:
            referenced.add(gene.id)
    if data.cluster.nrp:
        for nrps_gene in data.cluster.nrp.nrps_genes:
            referenced.add(nrps_gene.gene_id)
        for thio in data.cluster.nrp.thioesterases:
            referenced.add(thio.gene)
    if data.cluster.polyketide:
        for synthase in data.cluster.polyketide.synthases:
            referenced.update(set(synthase.genes))
            for module in synthase.modules:
                referenced.update(set(module.genes))
    if data.cluster.saccharide:
        for transferase in data.cluster.saccharide.glycosyltransferases:
            referenced.add(transferase.gene_id)
    for compound in data.cluster.compounds:
        for moiety in compound.chem_moieties:
            if moiety.subcluster:
                referenced.update(set(moiety.subcluster))

    missing = referenced.difference(existing)
    if missing:
        raise ValueError(
            f"{data.cluster.mibig_accession} refers to missing genes: {', '.join(sorted(missing))}")

    return MibigAnnotations(record.id, area, data, cache_file, pubmed_cache_file)


class PubmedEntry:
    def __init__(self, title: str, authors: List[str], year: str, journal: str, pmid: str) -> None:
        self.title = title
        self.authors = authors
        self.year = year
        self.journal = journal
        self.pmid = pmid

    @property
    def info(self) -> str:
        return f"{self.authors[0]} et al., {self.journal} ({self.year}) PMID:{self.pmid}"

    def to_json(self):
        return {
            "title": self.title,
            "authors": self.authors,
            "year": self.year,
            "journal": self.journal,
            "pmid": self.pmid,
        }


class PubmedCache:
    def __init__(self, cache_file):
        self.mappings: Dict[str, PubmedEntry] = {}
        if cache_file and os.path.exists(cache_file):
            with open(cache_file, 'r') as handle:
                entries = json.load(handle)
            for entry_id, entry_values in entries.items():
                self.mappings[entry_id] = PubmedEntry(
                    entry_values["title"], entry_values["authors"], entry_values["year"], entry_values["journal"], entry_values["pmid"])

    def add(self, title: str, authors: List[str], year: str, journal: str, pmid: str) -> None:
        self.mappings[pmid] = PubmedEntry(title, authors, year, journal, pmid)

    def get(self, pmid: str) -> PubmedEntry:
        return self.mappings[pmid]

    def get_missing(self, want: List[str]) -> List[str]:
        have = self.mappings.keys()
        return sorted(list(want - have))
