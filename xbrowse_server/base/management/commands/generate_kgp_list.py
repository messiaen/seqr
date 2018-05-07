from optparse import make_option
import sys
import os
from django.core.management.base import BaseCommand
from xbrowse_server.reports.utilities import fetch_project_individuals_data
import json
import time
import datetime



class Command(BaseCommand):
    __VERSION__= '0.0.1'

    def add_arguments(self, parser):
        parser.add_argument('args', nargs='*')
        parser.add_argument('--family_id',
                    dest='family_id',
                    help='Generate report for this family only.'
                    )
    def handle(self, *args, **options):
        '''
        '''
        pass
   