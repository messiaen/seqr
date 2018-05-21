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
from xbrowse import Variant
from xbrowse_server.mall import get_datastore
from xbrowse_server.mall import get_reference
from xbrowse_server.api import utils as api_utils
from pprint import pprint
logger = logging.getLogger()
from xbrowse_server.base.lookups import get_variants_by_tag
    
    
    
class Command(BaseCommand):
    """
    Generate a report of variants in specific project/family/gene in seqr
    """
    __VERSION__ = 'v0.0.1'
    
    
    def add_arguments(self, parser):
        parser.add_argument('args', nargs='*')
        parser.add_argument('--project_list')
 
    def handle(self,*args,**options):
        """
        Starting point for script
        
        Args:
            none needed
            
        Returns:
            Outputs a report printed to stdout
        """
        if options['project_list']:
            projects_to_process=self.process_family_list_file(options['project_list'])
                
        skipped_out = open("/Users/harindra/Desktop/skipper.txt",'w')
        all_families=[]
        for project_id_from_input,input_details in projects_to_process.iteritems():
            for line in input_details:
                family_id_from_input=line['family_id']
                gene_id_from_input=line['gene_name']
                print("----processing: " + project_id_from_input+"\t"+family_id_from_input + "\t" +gene_id_from_input+"\n")
                try:
                    project = get_object_or_404(Project, project_id=project_id_from_input)
                except Exception as e:
                    skipped_out.write(project_id_from_input+"\t"+family_id_from_input + "\t" +gene_id_from_input+"\t:" + e + "\n")
                    continue            
                for i,e in enumerate(get_variants_by_tag(project, 'Known gene for phenotype', family_id=family_id_from_input)):
                    family_details = self.process_family(e.toJSON(),family_id_from_input,project_id_from_input,gene_id_from_input)
                    if family_details:
                        all_families.append(family_details)
        results_out =  open("/Users/harindra/Desktop/kgp_populated.txt",'w')
        self.write_to_file(all_families, results_out)
        skipped_out.close()
        results_out.close()


    def process_family(self,fam,family_id_from_input,project_id_from_input,target_gene_name):
        '''
        Process a single family
        
        Args:
            (obj) A variant object
            
        Returns:
            [dict] representing information of target variant of that family. None if not a target family
        '''
        annotation_set_to_use = fam['annotation']['worst_vep_annotation_index']
        vep_annotation = fam['annotation']['vep_annotation'][annotation_set_to_use]
        
        if vep_annotation['gene_symbol'] in target_gene_name:
            return {
                            'family_id':family_id_from_input,
                            'gene_name':target_gene_name,
                            'cmg_internal_project_id':project_id_from_input,
                            'seqr_family_page_link':'https://seqr.broadinstitute.org/project/' + project_id_from_input  +'/family/' + family_id_from_input,
                            'start' : fam['pos'],
                            'stop' : fam['pos_end'],
                            'chromosome' : fam['chr'],
                            'reference_allele' : fam['ref'],
                            'alternate_allele' : fam['alt'],
                            'hgvs_c':vep_annotation['hgvsc'],
                            'hgvs_p':vep_annotation['hgvsp'],
                         }

    def write_to_file(self,fam_details,out_file):
        '''
        Write the given family details to file
        
        Args:
            fam_details (list(dict)): list of family details in form of dictionaries
            out_file_name (str): a name for the output file
        '''
        headers="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                                                                    'Family ID (CollPrefix_ID)',
                                                                    'CMG Internal Project ID(s)',
                                                                    'Gene Name',
                                                                    'Seqr family page Link',
                                                                    'Chromosome',
                                                                    'Start',
                                                                    'Stop',
                                                                    'Reference allele',
                                                                    'Alternate allele',
                                                                    'HGVS.C',
                                                                    'HGVS.P')
        out_file.write(headers)
        for i,fam in enumerate(fam_details):
            line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                                                    fam['family_id'],
                                                    fam['cmg_internal_project_id'],
                                                    fam['gene_name'],
                                                    fam['seqr_family_page_link'],
                                                    fam['chromosome'],
                                                    fam['start'],
                                                    fam['stop'],
                                                    fam['reference_allele'],
                                                    fam['alternate_allele'],
                                                    self.format_seqr_hgvs_c(fam['hgvs_c']),
                                                    self.format_seqr_hgvs_p(fam['hgvs_p']))
            out_file.write(line)


        
    
    def get_gene_symbol(self,variant):
        '''
        A bug instability in data model has this mix of data
        Given a gene_id, find the symbol from the reference
        
        Args:
            a gene id (str)        
        Returns:
            (str) a gene symbol; returns empty string if symbol is not found
        ''' 
        #gene = get_reference().get_gene(gene_id)
        #gene_symbol = gene['symbol'] if gene else ""
        #return gene_symbol
        if 'symbol' in variant:
            return variant['symbol']
        if 'gene_symbol' in variant:
            return variant['gene_symbol']
        return ""
    
    
    
    def format_seqr_hgvs_c(self,hgvsc_from_seqr):
        '''
        Given a HGVS.c str from seqr, format it
        
        Args:
            (str) a hgvs.c string
            
        Returns:
            (str) formatted string or empty string ""
        '''
        if hgvsc_from_seqr and "c." in hgvsc_from_seqr:
            return 'c.' + hgvsc_from_seqr.split('c.')[1]
        else:
            return hgvsc_from_seqr if hgvsc_from_seqr else ""
        
        
    def format_seqr_hgvs_p(self,hgvsp_from_seqr):
        '''
        Given a HGVS.p str from seqr, format it
        
        Args:
            (str) a hgvs.p string
            
        Returns:
            (str) formatted string or empty string ""
        '''
        if hgvsp_from_seqr and "p." in hgvsp_from_seqr:
            return 'p.' + hgvsp_from_seqr.split('p.')[1]
        else:
            return hgvsp_from_seqr if hgvsp_from_seqr else ""
            



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
                if fields[1] in to_process:
                    to_process[fields[1]].append({
                                'family_id':fields[0],
                                'internal_project_id':fields[1],
                                'gene_name':fields[2]
                                })
                else:
                    to_process[fields[1]] =[{
                                    'family_id':fields[0],
                                    'internal_project_id':fields[1],
                                    'gene_name':fields[2]
                                    }]
        fi.close()
        return to_process  



      