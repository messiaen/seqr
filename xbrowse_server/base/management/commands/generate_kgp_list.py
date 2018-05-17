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
        
        families_to_process = [] 
        if options['family_list']:
            families_to_process=self.process_family_list_file(options['family_list'])

        #going with all projects since the list doesn't have project names properly written to match safely (cutnpaste errors)
        all_projects = Project.objects.all()
        all_families=[]
        so = open("/Users/harindra/Desktop/skipper.txt",'w')
        for project in all_projects:
            for fam in project.get_families():
                if fam.family_id in families_to_process:
                    all_families.extend(self.process_family(fam,families_to_process[fam.family_id],so))
        self.write_to_file(all_families,'/Users/harindra/Desktop/kgp_populated.txt')
        so.close()

    def write_to_file(self,fam_details,out_file_name):
        '''
        Write the given family details to file
        
        Args:
            fam_details (list(dict)): list of family details in form of dictionaries
            out_file_name (str): a name for the output file
        '''
        
        with open(out_file_name,'w') as out:
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
            out.write(headers)
            for fam in fam_details:
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
                                                        fam['hgvs_c'],
                                                        fam['hgvs_p'])
                out.write(line)
        out.close()

        
    def process_family(self,fam,input_dets_on_fam,so):
        '''
        Process a single family
        
        Args:
            (obj) A Family object
            (list of objs with each being having a gene name) the information input in input file
            
        Returns:
            ([dict]) representing a list of dicts of variants of that family
        '''
        target_gene_symbols=[]
        for input_fam in input_dets_on_fam:
            target_gene_symbols.append(input_fam['gene_name'])
        fam_details=[]
        genotype_data = self.get_genotype_data(fam.project.project_id,fam.family_id)
        if input_fam['gene_name'] in target_gene_symbols and genotype_data ==[]:
            so.write(str(fam))
            so.write('\t')
            so.write(str(target_gene_symbols))
            so.write('\n')
        for i,entry in enumerate(genotype_data):
            if entry['gene_symbol'] in target_gene_symbols:
                fam_details.append({
                                'family_id':fam.family_id,
                                'gene_name':entry['gene_symbol'],
                                'cmg_internal_project_id':input_dets_on_fam[0]['internal_project_id'],
                                'seqr_family_page_link':'https://seqr.broadinstitute.org/project/' + fam.project.project_id  +'/family/' + fam.family_id,
                                'start' : entry['start'],
                                'stop' : entry['stop'],
                                'chromosome' : entry['chromosome'],
                                'reference_allele' : entry['reference_allele'],
                                'alternate_allele' : entry['alternate_allele'],
                                'hgvs_c':entry['hgvs_c'],
                                'hgvs_p':entry['hgvs_p'],
                             })
        return fam_details

        
        
        
    def get_genotype_data(self,project_id,family_id):
        '''
        Gets genotype data for this family
        
        Args:
            family_id (str): id of family
            project_id (str): id of family
            
        Returns:
            (list of dicts) 
        '''
        project = get_object_or_404(Project, project_id=project_id)
        genomic_features=[]
        variants=[]
        project_tags = ProjectTag.objects.filter(project__project_id=project_id)
        already_added_variant={}
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
                        #continue
                        variant = Variant.fromJSON({
                            'xpos': variant_tag.xpos,
                            'ref': variant_tag.ref,
                            'alt': variant_tag.alt,
                            'genotypes': {},
                            'extras': {'project_id': project.project_id, 'family_id': family_id}
                            })
                    api_utils.add_extra_info_to_variants_project(get_reference(), project, [variant], add_family_tags=False,add_populations=False)
                    
                    #import pprint
                    #pprint.pprint(variant.toJSON())
                    
                    if (variant.toJSON()['pos'],variant.toJSON()['pos_end'],variant.toJSON()['chr']) not in already_added_variant:
                        variants.append({"variant": variant.toJSON(),
                                         "family": variant_tag.family.toJSON()
                                     })
                        already_added_variant[(variant.toJSON()['pos'],variant.toJSON()['pos_end'],variant.toJSON()['chr'])]=1
        
        current_genome_assembly = self.find_genome_assembly(project)
        genomic_features=[]

        for variant in variants:
            if variant['variant']['annotation']:
                annotation_set_to_use = variant['variant']['annotation']['worst_vep_annotation_index']
                genomic_features.append({
                                    'gene_symbol':self.get_gene_symbol(variant['variant']['annotation']['vep_annotation'][annotation_set_to_use]),
                                    'assembly':current_genome_assembly,
                                    'reference_allele':variant['variant']['ref'],
                                    'alternate_allele':variant['variant']['alt'],
                                    'start':variant['variant']['pos'],
                                    'stop':int(variant['variant']['pos_end']),
                                    'chromosome':variant['variant']['chr'],
                                    'hgvs_c': self.format_seqr_hgvs_c(variant['variant']['annotation']['vep_annotation'][annotation_set_to_use]['hgvsc']),
                                    'hgvs_p': self.format_seqr_hgvs_p(variant['variant']['annotation']['vep_annotation'][annotation_set_to_use]['hgvsp'])
                                })
        return genomic_features
    
    
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
                if fields[0] in to_process:
                    to_process[fields[0]].append({
                                'family_id':fields[0],
                                'internal_project_id':fields[1],
                                'gene_name':fields[2]
                                })
                else:
                    to_process[fields[0]] =[{
                                    'family_id':fields[0],
                                    'internal_project_id':fields[1],
                                    'gene_name':fields[2]
                                    }]
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
      