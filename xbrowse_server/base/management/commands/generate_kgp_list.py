from __future__ import print_function
from optparse import make_option
import sys
import os
from django.core.management.base import BaseCommand
import json
from xbrowse_server.base.models import Project
import logging
import hashlib
from django.shortcuts import get_object_or_404
from xbrowse_server.base.models import ProjectTag, VariantTag
from xbrowse_server.mall import get_datastore
from xbrowse_server.mall import get_reference
from xbrowse_server.api import utils as api_utils

logger = logging.getLogger()
    
    
    
    
class Command(BaseCommand):
    """
    Generate a report of variants in specific project/family/gene in seqr
    """
    __VERSION__ = 'v0.0.1'
    
    
    def add_arguments(self, parser):
        parser.add_argument('args', nargs='*')
        parser.add_argument('--family_list')
 
    def handle(self,*args,**options):
        """
        Starting point for script
        
        Args:
            none needed
            
        Returns:
            Outputs a report printed to stdout
        """
        
        list_of_families_to_process = [] 
        if options['family_list']:
            list_of_families_to_process=self.process_family_list_file(options['family_list'])
        #going with all projects since the list doesn't have project names properly written to match safely (cutnpaste errors)
        all_projects = Project.objects.all()
        for project in all_projects:
            for fam in project.get_families():
                if fam.family_id in list_of_families_to_process:
                    fam_details = self.process_family(fam,list_of_families_to_process[fam.family_id])
                    self.write_to_file(fam_details,'/Users/harindra/Desktop/kgp_populated.txt')


    def write_to_file(self,fam_details,out_file_name):
        '''
        Write the given family details to file
        
        Args:
            fam_details (list(dict)): list of family details in form of dictionaries
            out_file_name (str): a name for the output file
        '''
        
        with open(out_file_name,'w') as out:
            for fam in fam_details:
                line = "%s,%s,%s,%s,%s,%s,%s,%s\n" % (  fam['seqr_family_page_link'],
                                                        fam['chromosome'],
                                                        fam['start'],
                                                        fam['stop'],
                                                        fam['reference_allele'],
                                                        fam['alternate_allele'],
                                                        fam['hgvs_c'],
                                                        fam['hgvs_p'])
                out.write(line)
                print(line)
        out.close()

        
    def process_family(self,fam,input_dets_on_fam):
        '''
        Process a single family
        
        Args:
            A Family object
            
        Returns:
            ([dict]) representing a list of dicts of variants of that family
        '''
        fam_details=[]
        for indiv in fam.get_individuals():
            genotype_data_for_indiv = self.get_genotype_data_for_indiv(fam.project.project_id,
                                             fam.family_id,
                                             indiv.indiv_id)
            for i,entry in enumerate(genotype_data_for_indiv):
                if input_dets_on_fam['gene_name'] in entry['auxiliary']['gene_symbol']:
                    fam_details.append({
                                    'seqr_family_page_link':'https://seqr.broadinstitute.org/project/' + fam.project.project_id  +'/family/' + fam.family_id,
                                    'start' : entry['variant']['start'],
                                    'stop' : entry['variant']['stop'],
                                    'chromosome' : entry['variant']['chromosome'],
                                    'reference_allele' : entry['variant']['reference_allele'],
                                    'alternate_allele' : entry['variant']['alternate_allele'],
                                    'hgvs_c':entry['variant']['hgvs_c'],
                                    'hgvs_p':entry['variant']['hgvs_p'],
                                 })
                    
        return fam_details

        
        
        
    def get_genotype_data_for_indiv(self,project_id,family_id,indiv_id):
        '''
        Gets genotype data for this individual
        
        Args:
            family_id (str): id of family
            project_id (str): id of family
            indiv_id (str): id of individualk
            
        Returns:
            (JSON) A JSON object of data
        '''
        project = get_object_or_404(Project, project_id=project_id)
        genomic_features=[]
        variants=[]
        project_tags = ProjectTag.objects.filter(project__project_id=project_id)
        for project_tag in project_tags:
            variant_tags = VariantTag.objects.filter(project_tag=project_tag)
            for variant_tag in variant_tags:
                if variant_tag.family is not None and family_id == variant_tag.family.family_id:
                    variant = get_datastore(project).get_single_variant(
                            project.project_id,
                            variant_tag.family.family_id,
                            variant_tag.xpos,
                            variant_tag.ref,
                            variant_tag.alt,
                    )

                    if variant is None:
                        logging.info("Variant no longer called in this family (did the callset version change?)")
                        continue
                    api_utils.add_extra_info_to_variants_project(get_reference(), project, [variant], add_family_tags=False,add_populations=False)
                    variants.append({"variant": variant.toJSON(),
                                     "tag": project_tag.title,
                                     "family": variant_tag.family.toJSON(),
                                     "tag_name": variant_tag.project_tag.tag,
                                 })
        
        current_genome_assembly = self.find_genome_assembly(project)
        genomic_features=[]
        for variant in variants:     
            #now we have more than 1 gene associated to these VAR postions,
            #so we will associate that information to each gene symbol
            for i,gene_id in enumerate(variant['variant']['gene_ids']):
                annotation_set_to_use = variant['variant']['annotation']['worst_vep_annotation_index']
                genomic_feature = {}
                genomic_feature['gene'] ={"id": gene_id }
                genomic_feature['variant']={
                                            'assembly':current_genome_assembly,
                                            'reference_allele':variant['variant']['ref'],
                                            'alternate_allele':variant['variant']['alt'],
                                            'start':variant['variant']['pos'],
                                            'stop':int(variant['variant']['pos_end']),
                                            'chromosome':variant['variant']['chr'],
                                            'hgvs_c':variant['variant']['annotation']['vep_annotation'][annotation_set_to_use]['hgvsc'],
                                            'hgvs_p':variant['variant']['annotation']['vep_annotation'][annotation_set_to_use]['hgvsp']
                                            }
                genomic_feature['zygosity'] = variant['variant']['genotypes'][indiv_id]['num_alt']
                gene_symbol=""
                if gene_id != "":
                    gene = get_reference().get_gene(gene_id)
                    gene_symbol = gene['symbol']
    
                genomic_feature['auxiliary']={
                                              "tag_name":variant['tag_name'],
                                              "gene_symbol":gene_symbol
                                              }
                genomic_features.append(genomic_feature) 
        return genomic_features



    def process_family_list_file(self,file_of_projects):
        '''
        Process the input file of family id, internal_project_ids, and gene_names
        
        Args:
            the input file name (str)
            
        Returns:
            (dict) A dictionary (key being seqr family id) of objects representing the parsed file
                   {'family_id':'', 'internal_project_id':'', 'gene_name':''}
        '''
        to_process={}
        with open(file_of_projects,'r') as fi:
            for line in fi:
                fields=line.rstrip().split('\t')
                to_process[fields[0]] ={
                                'family_id':fields[0],
                                'internal_project_id':fields[1],
                                'gene_name':fields[2]
                                }
        fi.close()
        return to_process  


        
    def find_genome_assembly(self,project):
        """
        Find the genome assembly of this individual
        Args:
            project: This is a seqr.project object reprenting the project
        Returns:
        The genome assembly version
        """
        if project.genome_version:
            return 'GRCh' + project.genome_version
        return 'GRCh37'
      