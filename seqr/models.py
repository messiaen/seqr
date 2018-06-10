from abc import abstractmethod
import os
import uuid


from django.contrib.auth.models import User, Group
from django.db import models
from django.utils import timezone
from django.utils.text import slugify as __slugify

from guardian.shortcuts import assign_perm

from seqr.utils.xpos_utils import get_chrom_pos, get_xpos
from reference_data.models import GENOME_VERSION_GRCh37, GENOME_VERSION_GRCh38, GENOME_VERSION_CHOICES
from django.conf import settings

CAN_VIEW = 'can_view'
CAN_EDIT = 'can_edit'
IS_OWNER = 'is_owner'

_SEQR_OBJECT_PERMISSIONS = (
    (CAN_VIEW, CAN_VIEW),
    (CAN_EDIT, CAN_EDIT),
    (IS_OWNER, IS_OWNER),
)


def _slugify(text):
    # using _ instead of - makes ids easier to select, and use without quotes in a wider set of contexts
    return __slugify(text).replace('-', '_')


class ModelWithGUID(models.Model):
    MAX_GUID_SIZE = 30

    guid = models.CharField(max_length=MAX_GUID_SIZE, db_index=True, unique=True)

    created_date = models.DateTimeField(default=timezone.now, db_index=True)
    created_by = models.ForeignKey(User, null=True, blank=True, related_name='+', on_delete=models.SET_NULL)

    # used for optimistic concurrent write protection (to detect concurrent changes)
    last_modified_date = models.DateTimeField(null=True, blank=True,  db_index=True)

    class Meta:
        abstract = True

    @abstractmethod
    def _compute_guid(self):
        """Returns a human-readable label (aka. slug) for this object with only alphanumeric
        chars, '-' and '_'. This label doesn't need to be globally unique by itself, but should not
        be null or blank, and should be globally unique when paired with this object's created-time
        in seconds.
        """

    def __unicode__(self):
        return self.guid

    def json(self):
        """Utility method that returns a json {field-name: value-as-string} mapping for all fields."""
        return {k: v for k, v in self.__dict__.items() if not k.startswith('_')}

    def save(self, *args, **kwargs):
        """Create a GUID at object creation time."""

        being_created = not self.pk
        if being_created and not self.created_date:
            self.created_date = timezone.now()
        else:
            self.last_modified_date = timezone.now()

        super(ModelWithGUID, self).save(*args, **kwargs)

        if being_created:
            self.guid = self._compute_guid()[:ModelWithGUID.MAX_GUID_SIZE]
            super(ModelWithGUID, self).save()


class Project(ModelWithGUID):
    DISEASE_AREA = [(da.lower().replace(" ", "_"), da) for da in (
        "Blood", "Cardio", "Kidney", "Muscle", "Neurodev", "Orphan Disease", "Retinal")
    ]

    name = models.TextField()  # human-readable project name
    description = models.TextField(null=True, blank=True)

    # user groups that allow Project permissions to be extended to other objects as long as
    # the user remains is in one of these groups.
    owners_group = models.ForeignKey(Group, related_name='+', on_delete=models.PROTECT)
    can_edit_group = models.ForeignKey(Group, related_name='+', on_delete=models.PROTECT)
    can_view_group = models.ForeignKey(Group, related_name='+', on_delete=models.PROTECT)

    #primary_investigator = models.ForeignKey(User, null=True, blank=True, related_name='+')

    is_phenotips_enabled = models.BooleanField(default=False)
    phenotips_user_id = models.CharField(max_length=100, null=True, blank=True, db_index=True)

    is_mme_enabled = models.BooleanField(default=True)
    mme_primary_data_owner = models.TextField(null=True, blank=True, default=settings.MME_DEFAULT_CONTACT_NAME)
    mme_contact_url = models.TextField(null=True, blank=True, default=settings.MME_DEFAULT_CONTACT_HREF)
    mme_contact_institution = models.TextField(null=True, blank=True, default=settings.MME_DEFAULT_CONTACT_INSTITUTION)

    is_functional_data_enabled = models.BooleanField(default=False)
    disease_area = models.CharField(max_length=20, null=True, blank=True, choices=DISEASE_AREA)

    # legacy
    custom_reference_populations = models.ManyToManyField('base.ReferencePopulation', blank=True, related_name='+')
    deprecated_last_accessed_date = models.DateTimeField(null=True, blank=True, db_index=True)
    deprecated_project_id = models.TextField(default="", blank=True, db_index=True)  # replace with model's 'id' field


    def __unicode__(self):
        return self.name.strip()

    def _compute_guid(self):
        label = (self.name or self.deprecated_project_id).strip()
        return 'R%04d_%s' % (self.id, _slugify(str(label)))

    def save(self, *args, **kwargs):
        """Override the save method and create user permissions groups + add the created_by user.

        This could be done with signals, but seems cleaner to do it this way.
        """
        being_created = not self.pk

        if being_created:
            # create user groups
            self.owners_group = Group.objects.create(name="%s_%s_%s" % (_slugify(self.name.strip())[:30], 'owners', uuid.uuid4()))
            self.can_edit_group = Group.objects.create(name="%s_%s_%s" % (_slugify(self.name.strip())[:30], 'can_edit', uuid.uuid4()))
            self.can_view_group = Group.objects.create(name="%s_%s_%s" % (_slugify(self.name.strip())[:30], 'can_view', uuid.uuid4()))

        super(Project, self).save(*args, **kwargs)

        if being_created:
            assign_perm(user_or_group=self.owners_group, perm=IS_OWNER, obj=self)
            assign_perm(user_or_group=self.owners_group, perm=CAN_EDIT, obj=self)
            assign_perm(user_or_group=self.owners_group, perm=CAN_VIEW, obj=self)

            assign_perm(user_or_group=self.can_edit_group, perm=CAN_EDIT, obj=self)
            assign_perm(user_or_group=self.can_edit_group, perm=CAN_VIEW, obj=self)

            assign_perm(user_or_group=self.can_view_group, perm=CAN_VIEW, obj=self)

            # add the user that created this Project to all permissions groups
            user = self.created_by
            if user and not user.is_staff:  # staff have access too all resources anyway
                user.groups.add(self.owners_group, self.can_edit_group, self.can_view_group)

    def delete(self, *args, **kwargs):
        """Override the delete method to also delete the project-specific user groups"""

        super(Project, self).delete(*args, **kwargs)

        self.owners_group.delete()
        self.can_edit_group.delete()
        self.can_view_group.delete()

    class Meta:
        permissions = _SEQR_OBJECT_PERMISSIONS


class ProjectCategory(ModelWithGUID):
    projects = models.ManyToManyField('Project')
    name = models.TextField(db_index=True)  # human-readable category name
    # color = models.CharField(max_length=20, default="#1f78b4")

    def __unicode__(self):
        return self.name.strip()

    def _compute_guid(self):
        return 'PC%06d_%s' % (self.id, _slugify(str(self)))


class Family(ModelWithGUID):
    ANALYSIS_STATUS_ANALYSIS_IN_PROGRESS='I'
    ANALYSIS_STATUS_WAITING_FOR_DATA='Q'
    ANALYSIS_STATUS_CHOICES = (
        ('S', 'Solved'),
        ('S_kgfp', 'Solved - known gene for phenotype'),
        ('S_kgdp', 'Solved - gene linked to different phenotype'),
        ('S_ng', 'Solved - novel gene'),
        ('Sc_kgfp', 'Strong candidate - known gene for phenotype'),
        ('Sc_kgdp', 'Strong candidate - gene linked to different phenotype'),
        ('Sc_ng', 'Strong candidate - novel gene'),
        ('Rcpc', 'Reviewed, currently pursuing candidates'),
        ('Rncc', 'Reviewed, no clear candidate'),
        ('I', 'Analysis in Progress'),
        ('Q', 'Waiting for data'),
    )

    CAUSAL_INHERITANCE_MODE_CHOICES = (
        ('r', 'recessive'),    # the actual inheritance model (the one in phenotips is the external inheritance model)
        ('u', 'unknown'),
        ('d', 'dominant'),
        ('x', 'x-linked recessive'),
        ('n', 'de novo'),


    )

    project = models.ForeignKey('Project', on_delete=models.PROTECT)

    # WARNING: family_id is unique within a project, but not necessarily unique globally.
    family_id = models.CharField(db_index=True, max_length=100)
    display_name = models.CharField(db_index=True, max_length=100, null=True, blank=True)  # human-readable name

    description = models.TextField(null=True, blank=True)

    pedigree_image = models.ImageField(null=True, blank=True, upload_to='pedigree_images')

    analysis_notes = models.TextField(null=True, blank=True)
    analysis_summary = models.TextField(null=True, blank=True)

    causal_inheritance_mode = models.CharField(max_length=20, default='u', choices=CAUSAL_INHERITANCE_MODE_CHOICES)

    coded_phenotype = models.TextField(null=True, blank=True)
    post_discovery_omim_number = models.TextField(null=True, blank=True)

    analysis_status = models.CharField(
        max_length=10,
        choices=[(s[0], s[1][0]) for s in ANALYSIS_STATUS_CHOICES],
        default="Q"
    )

    internal_analysis_status = models.CharField(
        max_length=10,
        choices=[(s[0], s[1][0]) for s in ANALYSIS_STATUS_CHOICES],
        null=True,
        blank=True
    )

    internal_case_review_notes = models.TextField(null=True, blank=True)
    internal_case_review_summary = models.TextField(null=True, blank=True)

    def __unicode__(self):
        return self.family_id.strip()

    def _compute_guid(self):
        return 'F%06d_%s' % (self.id, _slugify(str(self)))

    class Meta:
        unique_together = ('project', 'family_id')


class Individual(ModelWithGUID):
    SEX_CHOICES = (
        ('M', 'Male'),
        ('F', 'Female'),
        ('U', 'Unknown'),
    )

    AFFECTED_STATUS_CHOICES = (
        ('A', 'Affected'),
        ('N', 'Unaffected'),
        ('U', 'Unknown'),
    )

    CASE_REVIEW_STATUS_IN_REVIEW = "I"
    CASE_REVIEW_STATUS_CHOICES = (
        ('N', 'Not In Review'),
        ('I', 'In Review'),
        ('U', 'Uncertain'),
        ('A', 'Accepted'),
        ('R', 'Not Accepted'),
        ('Q', 'More Info Needed'),
        ('P', 'Pending Results and Records'),
        ('W', 'Waitlist'),
        ('WD', 'Withdrew'),
        ('IE', 'Ineligible'),
        ('DP', 'Declined to Participate'),
    )

    CASE_REVIEW_STATUS_ACCEPTED_FOR_OPTIONS = (
        ('S', 'Store DNA'),
        ('A', 'Array'),   # allow multiple-select. No selection = Platform Uncertain
        ('E', 'Exome'),
        ('G', 'Genome'),
        ('R', 'RNA-seq'),
        ('P', 'Reprocess'),

    )

    SEX_LOOKUP = dict(SEX_CHOICES)
    AFFECTED_STATUS_LOOKUP = dict(AFFECTED_STATUS_CHOICES)
    CASE_REVIEW_STATUS_LOOKUP = dict(CASE_REVIEW_STATUS_CHOICES)
    CASE_REVIEW_STATUS_REVERSE_LOOKUP = {name.lower(): key for key, name in CASE_REVIEW_STATUS_CHOICES}

    family = models.ForeignKey(Family, on_delete=models.PROTECT)

    # WARNING: individual_id is unique within a family, but not necessarily unique globally
    individual_id = models.TextField(db_index=True)
    maternal_id = models.TextField(null=True, blank=True, db_index=True)  # individual_id of mother
    paternal_id = models.TextField(null=True, blank=True, db_index=True)  # individual_id of father
    # add ForeignKeys for mother Individual & father Individual?

    sex = models.CharField(max_length=1, choices=SEX_CHOICES, default='U')
    affected = models.CharField(max_length=1, choices=AFFECTED_STATUS_CHOICES, default='U')

    display_name = models.TextField(default="", blank=True)

    notes = models.TextField(blank=True, null=True)

    case_review_status = models.CharField(max_length=2, choices=CASE_REVIEW_STATUS_CHOICES, null=True, blank=True)
    case_review_status_accepted_for = models.CharField(max_length=10, null=True, blank=True)
    case_review_status_last_modified_date = models.DateTimeField(null=True, blank=True, db_index=True)
    case_review_status_last_modified_by = models.ForeignKey(User, null=True, blank=True, related_name='+', on_delete=models.SET_NULL)
    case_review_discussion = models.TextField(null=True, blank=True)

    phenotips_patient_id = models.CharField(max_length=30, null=True, blank=True, db_index=True)    # PhenoTips internal id
    phenotips_eid = models.CharField(max_length=165, null=True, blank=True)  # PhenoTips external id
    phenotips_data = models.TextField(null=True, blank=True)

    mme_id = models.CharField(max_length=50, null=True, blank=True)
    mme_submitted_data = models.TextField(null=True, blank=True)

    def __unicode__(self):
        return self.individual_id.strip()

    def _compute_guid(self):
        return 'I%06d_%s' % (self.id, _slugify(str(self)))

    class Meta:
        unique_together = ('family', 'individual_id')


class UploadedFileForFamily(models.Model):
    family = models.ForeignKey(Family, on_delete=models.PROTECT)
    name = models.TextField()
    uploaded_file = models.FileField(upload_to="uploaded_family_files", max_length=200)
    uploaded_by = models.ForeignKey(User, null=True, on_delete=models.SET_NULL)
    uploaded_date = models.DateTimeField(null=True, blank=True)


class UploadedFileForIndividual(models.Model):
    individual = models.ForeignKey(Individual, on_delete=models.PROTECT)
    name = models.TextField()
    uploaded_file = models.FileField(upload_to="uploaded_individual_files", max_length=200)
    uploaded_by = models.ForeignKey(User, null=True, on_delete=models.SET_NULL)
    uploaded_date = models.DateTimeField(null=True, blank=True)


class ProjectLastAccessedDate(models.Model):
    """Used to provide a user-specific 'last_accessed' column in the project table"""
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    last_accessed_date = models.DateTimeField(auto_now=True, db_index=True)


class Sample(ModelWithGUID):
    """Represents a single biological sample.

    A sample can have be used to generate multiple types of analysis results, depending on the
    sample type. For example, an exome, genome or rna sample can be used to generate an aligned bam,
    a variant callset, CNV callset, etc., and an rna sample can also yield ASE, and splice junction
    data.

    For now, all meta-data on these analysis results is stored in the sample record, but
    if versioning is needed for analysis results, it'll be necessary to create a separate table
    for each analysis type, where records have a many(analysis-versions)-to-one(sample) relationship with this table.
    """

    SAMPLE_STATUS_CHOICES = (
        ('S', 'In Sequencing'),
        ('I', 'Interim'),    # needs additional sequencing to reach expected (95x) coverage
        ('C', 'Complete'),   # sample sequencing complete and achieved expected coverage
        ('A', 'Abandoned'),  # sample failed sequencing
    )

    SAMPLE_TYPE_WES = 'WES'
    SAMPLE_TYPE_WGS = 'WGS'
    SAMPLE_TYPE_RNA = 'RNA'
    SAMPLE_TYPE_ARRAY = 'ARRAY'
    SAMPLE_TYPE_CHOICES = (
        (SAMPLE_TYPE_WES, 'Exome'),
        (SAMPLE_TYPE_WGS, 'Whole Genome'),
        (SAMPLE_TYPE_RNA, 'RNA'),
        (SAMPLE_TYPE_ARRAY, 'ARRAY'),
        # ('ILLUMINA_INFINIUM_250K', ),
    )

    sample_type = models.CharField(max_length=3, choices=SAMPLE_TYPE_CHOICES)

    individual = models.ForeignKey('Individual', on_delete=models.PROTECT, null=True)

    # This sample_id should be used for looking up this sample in the underlying dataset (for
    # example, for variant callsets, it should be the VCF sample id). It is not a ForeignKey
    # into another table.
    sample_id = models.TextField(db_index=True)

    # biological sample status
    sample_status = models.CharField(max_length=1, choices=SAMPLE_STATUS_CHOICES, default='S')

    #funding_source = models.CharField(max_length=20, null=True)

    is_external_data = models.BooleanField(default=False)


    #sample_batch = models.ForeignKey('SampleBatch', on_delete=models.PROTECT, null=True)

    # reference back to xbrowse base_project is a temporary work-around to support merging of
    # different projects into one - including those that contain different types of callsets
    # such as exome and genome
    deprecated_base_project = models.ForeignKey('base.Project', null=True, on_delete=models.SET_NULL)

    def __unicode__(self):
        return self.sample_id.strip()

    def _compute_guid(self):
        return 'S%06d_%s' % (self.id, _slugify(str(self)))

    #class Meta:
    #    unique_together = ('sample_batch', 'sample_id')


class Dataset(ModelWithGUID):
    """Dataset represents a single analysis type and result file"""

    ANALYSIS_TYPE_ALIGNMENT = 'ALIGN'
    ANALYSIS_TYPE_VARIANT_CALLS = 'VARIANTS'
    ANALYSIS_TYPE_SV = 'SV'
    ANALYSIS_TYPE_BREAKPOINTS = 'BREAK'
    ANALYSIS_TYPE_SPLICE_JUNCTIONS = 'SPLICE'
    ANALYSIS_TYPE_ASE = 'ASE'
    ANALYSIS_TYPE_CHOICES = (
        (ANALYSIS_TYPE_ALIGNMENT, 'Alignment'),
        (ANALYSIS_TYPE_VARIANT_CALLS, 'Variant Calls'),
        (ANALYSIS_TYPE_SV, 'SV Calls'),
        (ANALYSIS_TYPE_BREAKPOINTS, 'Breakpoints'),
        (ANALYSIS_TYPE_SPLICE_JUNCTIONS, 'Splice Junction Calls'),
        (ANALYSIS_TYPE_ASE, 'Allele Specific Expression'),
    )

    DATASET_STATUS_QUEUED = 'Q'
    DATASET_STATUS_LOADING = 'G'
    DATASET_STATUS_LOADED = 'L'
    DATASET_STATUS_FAILED = 'F'
    DATASET_STATUS_CHOICES = (
        (DATASET_STATUS_QUEUED, 'Queued'),
        (DATASET_STATUS_LOADING, 'Loading'),
        (DATASET_STATUS_LOADED, 'Loaded'),
        (DATASET_STATUS_FAILED, 'Failed'),
    )

    name = models.TextField(null=True, blank=True)
    description = models.TextField(null=True, blank=True)
    status = models.CharField(null=True, blank=True, db_index=True, max_length=1, choices=DATASET_STATUS_CHOICES)

    genome_version = models.CharField(max_length=10, choices=GENOME_VERSION_CHOICES, default=GENOME_VERSION_GRCh37)

    # When a dataset is copied from source_file_path to an internal seqr database or directory,
    # the dataset_id should be the pointer used to query this data. Although global uniqueness
    # is not enforced, the dataset_id value should avoid collisions, and should be derived only from
    # properties of the dataset itself (eg. creation date, size, or md5) so that if a dataset
    # is added a second time, it would be assigned the same dataset id as before.
    # This will allow datasets to be processed and loaded only once, but shared between projects if
    # needed by using the same dataset_id in the Dataset records of different projects.
    dataset_id = models.TextField(null=True, blank=True, db_index=True)   # elasticsearch index

    dataset_location = models.TextField(null=True, blank=True, db_index=True)

    analysis_type = models.CharField(max_length=10, choices=ANALYSIS_TYPE_CHOICES)

    source_file_path = models.TextField(db_index=True)

    is_loaded = models.BooleanField(default=False)
    loaded_date = models.DateTimeField(null=True, blank=True)

    samples = models.ManyToManyField('Sample')

    # for convenience, add a pointer to the project that this dataset and samples belong to
    project = models.ForeignKey('Project', null=True, on_delete=models.CASCADE)

    #tool = models.TextField(null=True, blank=True)
    #tool_version = models.TextField(null=True, blank=True)

    def __unicode__(self):
        return self.guid

    def _compute_guid(self):
        filename = os.path.basename(self.source_file_path).split(".")[0]
        return 'D%06d_%s_%s' % (self.id, self.analysis_type[0:3], filename)


# TODO AliasFields work for lookups, but save/update doesn't work?
class AliasField(models.Field):
    def contribute_to_class(self, cls, name, private_only=False):
        super(AliasField, self).contribute_to_class(cls, name, private_only=True)
        setattr(cls, name, self)

    def __get__(self, instance, instance_type=None):
        return getattr(instance, self.db_column)


#class SampleBatch(ModelWithGUID):
#    """Represents a set of biological samples that were processed together."""
#
#    notes = models.TextField(null=True, blank=True)
#
#    def __unicode__(self):
#        return self.name.strip()
#
#    def _compute_guid(self):
#        return 'D%05d_%s' % (self.id, _slugify(str(self)))


class VariantTagType(ModelWithGUID):
    """
    Previous color choices:
        '#1f78b4',
        '#a6cee3',
        '#b2df8a',
        '#33a02c',
        '#fdbf6f',
        '#ff7f00',
        '#ff0000',
        '#cab2d6',
        '#6a3d9a',
        '#8F754F',
        '#383838',
    """
    project = models.ForeignKey('Project', on_delete=models.CASCADE)

    name = models.TextField()
    category = models.TextField(null=True, blank=True)
    description = models.TextField(null=True, blank=True)
    color = models.CharField(max_length=20, default="#1f78b4")
    order = models.FloatField(null=True)
    is_built_in = models.BooleanField(default=False)  # built-in tags (eg. "Pathogenic") can't be modified by users through the UI

    def __unicode__(self):
        return self.name.strip()

    def _compute_guid(self):
        return 'VTT%05d_%s' % (self.id, _slugify(str(self)))

    class Meta:
        unique_together = ('project', 'name', 'color')


class VariantTag(ModelWithGUID):
    variant_tag_type = models.ForeignKey('VariantTagType', on_delete=models.CASCADE)

    genome_version = models.CharField(max_length=10, choices=GENOME_VERSION_CHOICES, default=GENOME_VERSION_GRCh37)
    xpos_start = models.BigIntegerField()
    xpos_end = models.BigIntegerField(null=True)
    xpos = AliasField(db_column="xpos_start")
    ref = models.TextField()
    alt = models.TextField()

    lifted_over_genome_version = models.CharField(max_length=10, null=True, blank=True, choices=GENOME_VERSION_CHOICES)
    lifted_over_xpos_start = models.BigIntegerField(null=True)
    lifted_over_xpos = AliasField(db_column="lifted_over_xpos_start")

    # Cache genotypes and annotations for the variant as gene id and consequence - in case the dataset gets deleted, etc.
    saved_variant_json = models.TextField(null=True, blank=True)

    # context in which a variant tag was saved
    family = models.ForeignKey('Family', null=True, blank=True, on_delete=models.SET_NULL)
    search_parameters = models.TextField(null=True, blank=True)  # aka. search url

    def __unicode__(self):
        chrom, pos = get_chrom_pos(self.xpos_start)
        return "%s:%s: %s" % (chrom, pos, self.variant_tag_type.name)

    def _compute_guid(self):
        return 'VT%07d_%s' % (self.id, _slugify(str(self)))

    class Meta:
        index_together = ('xpos_start', 'ref', 'alt', 'genome_version')

        unique_together = ('variant_tag_type', 'genome_version', 'xpos_start', 'xpos_end', 'ref', 'alt', 'family')


class VariantNote(ModelWithGUID):
    project = models.ForeignKey('Project', null=True, on_delete=models.SET_NULL)

    note = models.TextField(null=True, blank=True)
    submit_to_clinvar = models.BooleanField(default=False)

    genome_version = models.CharField(max_length=10, choices=GENOME_VERSION_CHOICES, default=GENOME_VERSION_GRCh37)
    xpos_start = models.BigIntegerField()
    xpos_end = models.BigIntegerField(null=True)
    xpos = AliasField(db_column="xpos_start")
    ref = models.TextField()
    alt = models.TextField()

    lifted_over_genome_version = models.CharField(max_length=10, null=True, blank=True, choices=GENOME_VERSION_CHOICES)
    lifted_over_xpos_start = models.BigIntegerField(null=True)
    lifted_over_xpos = AliasField(db_column="lifted_over_xpos_start")

    # Cache genotypes and annotations for the variant as gene id and consequence - in case the dataset gets deleted, etc.
    saved_variant_json = models.TextField(null=True, blank=True)

    # these are for context - if note was saved for a family or an individual
    family = models.ForeignKey('Family', null=True, blank=True, on_delete=models.SET_NULL)
    search_parameters = models.TextField(null=True, blank=True)  # aka. search url

    def __unicode__(self):
        chrom, pos = get_chrom_pos(self.xpos_start)
        return "%s:%s: %s" % (chrom, pos, (self.note or "")[:20])

    def _compute_guid(self):
        return 'VT%07d_%s' % (self.id, _slugify(str(self)))


class LocusList(ModelWithGUID):
    """List of gene ids or regions"""

    name = models.TextField(db_index=True)
    description = models.TextField(null=True, blank=True)

    is_public = models.BooleanField(default=False)

    def __unicode__(self):
        return self.name.strip()

    def _compute_guid(self):
        return 'LL%05d_%s' % (self.id, _slugify(str(self)))

    class Meta:
        permissions = _SEQR_OBJECT_PERMISSIONS


class LocusListGene(ModelWithGUID):
    locus_list = models.ForeignKey('LocusList', on_delete=models.CASCADE)

    gene_id = models.TextField(db_index=True)

    description = models.TextField(null=True, blank=True)

    def __unicode__(self):
        return "%s:%s" % (self.locus_list, self.gene_id)

    def _compute_guid(self):
        return 'LLG%07d_%s' % (self.id, _slugify(str(self)))

    class Meta:
        unique_together = ('locus_list', 'gene_id')


class LocusListInterval(ModelWithGUID):
    locus_list = models.ForeignKey('LocusList', on_delete=models.CASCADE)

    genome_version = models.CharField(max_length=10, choices=GENOME_VERSION_CHOICES, default=GENOME_VERSION_GRCh37)
    chrom = models.CharField(max_length=2)
    start = models.IntegerField()
    end = models.IntegerField()

    description = models.TextField(null=True, blank=True)

    def __unicode__(self):
        return "%s:%s:%s-%s" % (self.locus_list, self.chrom, self.start, self.end)

    def _compute_guid(self):
        return 'LLI%07d_%s' % (self.id, _slugify(str(self)))

    class Meta:
        unique_together = ('locus_list', 'genome_version', 'chrom', 'start', 'end')


"""
class FamilyGroup(ModelWithGUID):
    project = models.ForeignKey(Project, on_delete=models.CASCADE)

    name = models.TextField()
    description = models.TextField(null=True, blank=True)

    families = models.ManyToManyField(Family)

    def __unicode__(self):
        return self.name
"""

