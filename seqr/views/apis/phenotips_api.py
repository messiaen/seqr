"""
seqr integrates with the `PhenoTips <http://phenotips.org>`_ UI so that users can enter
detailed phenotype information for individuals.
The PhenoTips server is installed locally so that it's not visible to the outside network, and
seqr then acts as a proxy for all HTTP requests between the PhenoTips web-based UI (which is
running in users' browser), and the PhenoTips server (running on the internal network).

This proxy setup allows seqr to check authentication and authorization before allowing users to
access patients in PhenoTips, and is similar to how seqr manages access to the SQL database and
other internal systems.

This module implements the proxy functionality + methods for making requests to PhenoTips HTTP APIs.

PhenoTips API docs are at:

https://phenotips.org/DevGuide/API
https://phenotips.org/DevGuide/RESTfulAPI
https://phenotips.org/DevGuide/PermissionsRESTfulAPI
"""

import json
import logging
#import pickle
import re
import requests
from requests.auth import HTTPBasicAuth

import settings

from django.contrib.auth.decorators import login_required
from django.views.decorators.csrf import csrf_exempt
from django.http import HttpResponse
from django.core.exceptions import PermissionDenied, ObjectDoesNotExist

from reference_data.models import HumanPhenotypeOntology
from seqr.models import Project, CAN_EDIT, CAN_VIEW, Individual
from seqr.views.apis.auth_api import API_LOGIN_REQUIRED_URL
from seqr.views.utils.permissions_utils import check_permissions

logger = logging.getLogger(__name__)


DEBUG = False

PHENOTIPS_QUICK_SAVE_URL_REGEX = "/bin/preview/data/(P[0-9]{1,20})"

DO_NOT_PROXY_URL_KEYWORDS = [
    '/delete',
    '/logout',
    '/login',
    '/admin',
    '/CreatePatientRecord',
    '/bin/PatientAccessRightsManagement',
    '/ForgotUsername',
    '/ForgotPassword',
]


def create_patient(project, patient_eid, patient_json=None):
    """Create a new PhenoTips patient record with the given patient id.

    Args:
        project (Model): PhenoTips permissions will be set to only allow access from this seqr project.
        patient_eid (string): external id (eg. "NA12878") to use for the patient. Must be globally unique across all PhenoTips patients.
        patient_json (dict): optional - phenotips patient record like the object returned by get_patient_data(..).
    Raises:
        PhenotipsException: if unable to create patient record
    """
    url = '/bin/PhenoTips/OpenPatientRecord?create=true&eid=%(patient_eid)s' % locals()
    auth_tuple = _get_phenotips_uname_and_pwd_for_project(project.phenotips_user_id)
    _make_api_call('GET', url, auth_tuple=auth_tuple, verbose=False, parse_json_resonse=False)

    #url = '/rest/patients/eid/%(patient_eid)s' % locals()

    #data = {}
    #if patient_json is not None:
    #    data = json.dumps(patient_json)

    #auth_tuple = _get_phenotips_uname_and_pwd_for_project(project.phenotips_user_id)
    #_make_api_call('POST', url, data=data, auth_tuple=auth_tuple, verbose=False, parse_json_resonse=False)

    patient_data = get_patient_data(project, patient_eid, is_external_id=True)
    patient_id = patient_data['id']

    username, _ = _get_phenotips_uname_and_pwd_for_project(project.phenotips_user_id)
    response = add_user_to_patient(username, patient_id, allow_edit=True)
    logger.info("Added PhenoTips user %(username)s to %(patient_id)s: %(response)s" % locals())
    username_read_only, _ = _get_phenotips_uname_and_pwd_for_project(project.phenotips_user_id, read_only=True)
    add_user_to_patient(username_read_only, patient_id, allow_edit=False)
    logger.info("Added PhenoTips user %(username)s to %(patient_id)s: %(response)s" % locals())

    return patient_data


def get_patient_data(project, patient_id, is_external_id=False):
    """Retrieves patient data from PhenoTips and returns a json obj.

    Args:
        project (Model): used to retrieve PhenoTips credentials
        patient_id (string): PhenoTips patient id (either internal eg. "P000001" or external eg. "NA12878")
        is_external_id (bool): whether the provided id is an external id
    Returns:
        dict: json dictionary containing all PhenoTips information for this patient
    Raises:
        PhenotipsException: if unable to retrieve data from PhenoTips
    """

    if is_external_id:  url = '/rest/patients/eid/%(patient_id)s' % locals()
    else:               url = '/rest/patients/%(patient_id)s' % locals()

    auth_tuple = _get_phenotips_uname_and_pwd_for_project(project.phenotips_user_id)
    return _make_api_call('GET', url, auth_tuple=auth_tuple, verbose=False)


def update_patient_data(project, patient_id, patient_json, is_external_id=False):
    """Updates patient data in PhenoTips to the values in patient_json.

    Args:
        project (Model): used to retrieve PhenoTips credentials
        patient_id (string): PhenoTips patient id (either internal eg. "P000001" or external eg. "NA12878")
        patient_json (dict): phenotips patient record like the object returned by get_patient_data(..).
        is_external_id (bool): whether the provided id is an external id
    Raises:
        PhenotipsException: if api call fails
    """
    if not patient_json:
        raise ValueError("patient_json arg is empty")

    if is_external_id:  url = '/rest/patients/eid/%(patient_id)s' % locals()
    else:               url = '/rest/patients/%(patient_id)s' % locals()

    auth_tuple = _get_phenotips_uname_and_pwd_for_project(project.phenotips_user_id, read_only=False)
    return _make_api_call('PUT', url, data=json.dumps(patient_json), auth_tuple=auth_tuple, expected_status_code=204)


def delete_patient(project, patient_id, is_external_id=False):
    """Deletes patient from PhenoTips for the given patient_id.

    Args:
        project (Model): used to retrieve PhenoTips credentials
        patient_id (string): PhenoTips patient id (either internal eg. "P000001" or external eg. "NA12878")
        is_external_id (bool): whether the provided id is an external id
    Raises:
        PhenotipsException: if api call fails
    """

    if is_external_id:  url = '/rest/patients/eid/%(patient_id)s' % locals()
    else:               url = '/rest/patients/%(patient_id)s' % locals()

    auth_tuple = _get_phenotips_uname_and_pwd_for_project(project.phenotips_user_id, read_only=False)
    return _make_api_call('DELETE', url, auth_tuple=auth_tuple, expected_status_code=204)


def update_patient_field_value(project, patient_id, field_name, field_value, is_external_id=False):
    """ Utility method for updating one field in the patient record, while keeping other fields
    the same. For field descriptions, see https://phenotips.org/DevGuide/JSONExport1.3

    Args:
        project (Model): used to retrieve PhenoTips credentials
        patient_id (string): PhenoTips patient id (either internal eg. "P000001" or external eg. "NA12878")
        field_name (string): PhenoTips patient field name (eg. "family_history").
        field_value (string or dict): PhenoTips HPO terms.
        is_external_id (bool): whether the provided patient id is an external id
    Raises:
        PhenotipsException: if api call fails
    """
    if field_name not in set([
        "allergies",
        "apgar",
        "clinicalStatus",
        "date_of_birth",
        "date_of_death",
        "disorders",
        "ethnicity",
        "family_history",
        "features",
        "genes",
        "global_age_of_onset",
        "global_mode_of_inheritance",
        "life_status",
        "nonstandard_features",
        "notes",
        "prenatal_perinatal_history",
        "sex",
        "solved",
        "specificity",
        "variants",
    ]):
        raise ValueError("Unexpected field_name: %s" % (field_name, ))

    patient_json = get_patient_data(project, patient_id, is_external_id=is_external_id)

    patient_json[field_name] = field_value

    update_patient_data(project, patient_id, patient_json, is_external_id=is_external_id)


def set_patient_hpo_terms(project, patient_id, hpo_terms, is_external_id=False):
    """Utility method for specifying a list of HPO IDs for a patient.

    Args:
        project (Model): used to retrieve PhenoTips credentials
        patient_id (string): PhenoTips patient id (either internal eg. "P000001" or external eg. "NA12878")
        hpo_terms (list): list of HPO ID strings (eg. ["HP:00012345", "HP:0012346", ...])
        is_external_id (bool): whether the provided id is an external id
    Raises:
        PhenotipsException: if api call fails
    """
    field_name = "features"
    field_value = [{"id": hpo_term} for hpo_term in hpo_terms]

    update_patient_field_value(project, patient_id, field_name, field_value, is_external_id=is_external_id)


def add_user_to_patient(username, patient_id, allow_edit=True):
    """Grant a PhenoTips user access to the given patient.

    Args:
        username (string): PhenoTips username to grant access to.
        patient_id (string): PhenoTips internal patient id.
    """
    headers = { "Content-Type": "application/x-www-form-urlencoded" }
    data = {
        'collaborator': 'XWiki.' + str(username),
        'patient': patient_id,
        'accessLevel': 'edit' if allow_edit else 'view',
        'xaction': 'update',
        'submit': 'Update'
    }

    url = '/bin/get/PhenoTips/PatientAccessRightsManagement?outputSyntax=plain'
    return _make_api_call(
        'POST',
        url,
        http_headers=headers,
        data=data,
        auth_tuple=(settings.PHENOTIPS_ADMIN_UNAME, settings.PHENOTIPS_ADMIN_PWD),
        expected_status_code=204
    )


def create_phenotips_user(username, password):
    """Creates a new user in PhenoTips"""

    headers = { "Content-Type": "application/x-www-form-urlencoded" }
    data = { 'parent': 'XWiki.XWikiUsers' }

    url = '/rest/wikis/xwiki/spaces/XWiki/pages/%(username)s' % locals()
    _make_api_call(
        'PUT',
        url,
        http_headers=headers,
        data=data,
        auth_tuple=(settings.PHENOTIPS_ADMIN_UNAME, settings.PHENOTIPS_ADMIN_PWD),
        parse_json_resonse=False,
        expected_status_code=[201, 202],
    )

    data = {
        'className': 'XWiki.XWikiUsers',
        'property#password': password,
        #'property#first_name': first_name,
        #'property#last_name': last_name,
        #'property#email': email_address,
    }

    url = '/rest/wikis/xwiki/spaces/XWiki/pages/%(username)s/objects' % locals()
    return _make_api_call(
        'POST',
        url,
        data=data,
        auth_tuple=(settings.PHENOTIPS_ADMIN_UNAME, settings.PHENOTIPS_ADMIN_PWD),
        parse_json_resonse=False,
        expected_status_code=201,
    )


@login_required
@csrf_exempt
def phenotips_pdf_handler(request, project_guid, patient_id):
    """Requests the PhenoTips PDF for the given patient_id, and forwards PhenoTips' response to the client.

    Args:
        request: Django HTTP request object
        project_guid (string): project GUID for the seqr project containing this individual
        patient_id (string): PhenoTips internal patient id
    """

    url = "/bin/export/data/%(patient_id)s?format=pdf&pdfcover=0&pdftoc=0&pdftemplate=PhenoTips.PatientSheetCode" % locals()
    project = Project.objects.get(guid=project_guid)

    check_permissions(project, request.user, CAN_VIEW)

    auth_tuple = _get_phenotips_username_and_password(request.user, project, permissions_level=CAN_VIEW)

    return _send_request_to_phenotips('GET', url, scheme=request.scheme, auth_tuple=auth_tuple)



@login_required
@csrf_exempt
def phenotips_edit_handler(request, project_guid, patient_id):
    """Request the PhenoTips Edit page for the given patient_id, and forwards PhenoTips' response to the client.

    Args:
        request: Django HTTP request object
        project_guid (string): project GUID for the seqr project containing this individual
        patient_id (string): PhenoTips internal patient id
    """

    # query string forwarding needed for PedigreeEditor button
    query_string = request.META["QUERY_STRING"]
    url = "/bin/edit/data/%(patient_id)s?%(query_string)s" % locals()

    project = Project.objects.get(guid=project_guid)

    check_permissions(project, request.user, CAN_EDIT)

    auth_tuple = _get_phenotips_username_and_password(request.user, project, permissions_level=CAN_EDIT)

    #if 'current_phenotips_session' not in request.session:
    #    phenotips_session = requests.Session()
    #    request.session['current_phenotips_session'] = pickle.dumps(phenotips_session)
    #else:
    #    phenotips_session = pickle.loads(request.session['current_phenotips_session'])

    return _send_request_to_phenotips('GET', url, scheme=request.scheme, auth_tuple=auth_tuple) #, session=phenotips_session)


@login_required(login_url=API_LOGIN_REQUIRED_URL)
@csrf_exempt
def proxy_to_phenotips_handler(request):
    """This django view accepts GET and POST requests and forwards them to PhenoTips"""

    url = request.get_full_path()
    if any([k for k in DO_NOT_PROXY_URL_KEYWORDS if k.lower() in url.lower()]):
        logger.warn("Blocked proxy url: " + str(url))
        return HttpResponse(status=204)
    logger.info("Proxying url: " + str(url))

    # forward the request to PhenoTips, and then the PhenoTips response back to seqr
    http_headers = _convert_django_META_to_http_headers(request.META)
    http_headers = {key: value for key, value in http_headers.items() if key.lower() not in HTTP_REQUEST_HEADERS_TO_NOT_PROXY}

    #if 'current_phenotips_session' not in request.session:
    #    phenotips_session = requests.Session()
    #    request.session['current_phenotips_session'] = pickle.dumps(phenotips_session)
    #else:
    #    phenotips_session = pickle.loads(request.session['current_phenotips_session'])

    # Some PhenoTips endpoints that use HTTP redirects lose the phenotips JSESSION auth cookie
    # along the way, and don't proxy correctly. Using a Session object as below to store the cookies
    # provides a work-around.
    phenotips_session = requests.Session()
    for key, value in request.COOKIES.items():
        phenotips_session.cookies.set(key, value)

    http_response = _send_request_to_phenotips(
        request.method,
        url,
        scheme=request.scheme,
        http_headers=http_headers,
        data=request.body,
        session=phenotips_session,
    )

    # if this is the 'Quick Save' request, also save a copy of phenotips data in the seqr SQL db.
    match = re.match(PHENOTIPS_QUICK_SAVE_URL_REGEX, url)
    if match:
        _handle_phenotips_save_request(patient_id=match.group(1), http_headers=http_headers)

    return http_response


def _make_api_call(
        method,
        url,
        http_headers=None,
        data=None,
        auth_tuple=None,
        expected_status_code=200,
        parse_json_resonse=True,
        verbose=False):
    """Utility method for making an API call and then parsing & returning the json response.

    Args:
        method (string): 'GET' or 'POST'
        url (string): url path, starting with '/' (eg. '/bin/edit/data/P0000001')
        data (string): request body - used for POST, PUT, and other such requests.
        auth_tuple (tuple): ("username", "password") pair
        expected_status_code (int or list): expected server response code
        parse_json_resonse (bool): whether to parse and return the json response
        verbose (bool): whether to print details about the request & response
    Returns:
        json object or None if response content is empty
    """

    try:
        response = _send_request_to_phenotips(method, url, http_headers=http_headers, data=data, auth_tuple=auth_tuple, verbose=verbose)
    except requests.exceptions.RequestException as e:
        raise PhenotipsException(e.message)
    if (isinstance(expected_status_code, int) and response.status_code != expected_status_code) or (
        isinstance(expected_status_code, list) and response.status_code not in expected_status_code):
        raise PhenotipsException("Unable to retrieve %s. response code = %s: %s" % (
            url, response.status_code, response.reason_phrase))

    if parse_json_resonse:
        if not response.content:
            return {}

        try:
            return json.loads(response.content)
        except ValueError as e:
            logger.error("Unable to parse PhenoTips response for %s request to %s" % (method, url))
            raise PhenotipsException("Unable to parse response for %s:\n%s" % (url, e))


def _send_request_to_phenotips(method, url, scheme="http", http_headers=None, data=None, auth_tuple=None, session=None, verbose=False):
    """Send an HTTP request to a PhenoTips server.
    (see PhenoTips API docs: https://phenotips.org/DevGuide/RESTfulAPI)

    Args:
        method (string): 'GET' or 'POST'
        url (string): url path, starting with '/' (eg. '/bin/edit/data/P0000001')
        scheme: request scheme (typically "http" or "https")
        http_headers: (dict): HTTP headers to send
        data (bytes): body of a POST request
        auth_tuple: ("username", "password") pair

    Returns:
        HttpResponse from the PhenoTips server.
    """

    if http_headers:
        http_headers['Host'] = settings.PHENOTIPS_SERVER

    r = requests
    if session is not None:
        r = session
    else:
        r = requests.Session()

    if method == "GET":
        method_impl = r.get
    elif method == "POST":
        method_impl = r.post
    elif method == "PUT":
        method_impl = r.put
    elif method == "HEAD":
        method_impl = r.head
    elif method == "DELETE":
        method_impl = r.delete
    else:
        raise ValueError("Unexpected HTTP method: %s. %s" % (method, url))

    auth = HTTPBasicAuth(*auth_tuple) if auth_tuple is not None else None

    if not url.startswith("http"):
        if not url.startswith("/"):
            raise ValueError("%s url doesn't start with /" % url)

        url = "%s://%s%s" % (scheme, settings.PHENOTIPS_SERVER, url)

    if verbose or DEBUG:
        logger.info("Sending %(method)s request to %(url)s" % locals())
        if http_headers:
            logger.info("  headers:")
            for key, value in sorted(http_headers.items(), key=lambda i: i[0]):
                logger.info("---> %(key)s: %(value)s" % locals())
        if data:
            logger.info("  data: %(data)s" % locals())
        if auth:
            logger.info("  auth: %(auth_tuple)s" % locals())

    response = method_impl(url, headers=http_headers, data=data, auth=auth)

    http_response = HttpResponse(
        status=response.status_code,
        content=response.content,
        reason=response.reason,
        charset=response.encoding
    )
    if verbose or DEBUG:
        from requests_toolbelt.utils import dump
        data = dump.dump_all(response)
        logger.info("===> dump - phenotips_api:\n" + str(data))


        logger.info("  response: <Response: %s> %s" % (response.status_code, response.reason))
        logger.info("  response-headers:")
        for key, value in sorted(response.headers.items(), key=lambda i: i[0]):
            logger.info("<--- %(key)s: %(value)s" % locals())

    for header_key, header_value in response.headers.items():
        if header_key.lower() not in HTTP_RESPONSE_HEADERS_TO_NOT_PROXY:
            http_response[header_key] = header_value

    return http_response


def _handle_phenotips_save_request(patient_id, http_headers):
    """Update the seqr SQL database record for this patient with the just-saved phenotype data."""

    url = '/rest/patients/%s' % patient_id

    http_headers = {key: value for key, value in http_headers.items() if key == 'Cookie'}
    response = _send_request_to_phenotips('GET', url, http_headers=http_headers)
    if response.status_code != 200:
        logger.error("ERROR: unable to retrieve patient json. %s %s %s" % (
            url, response.status_code, response.reason_phrase))
        return

    patient_json = json.loads(response.content)

    try:
        if patient_json.get('external_id'):
            # prefer to use the external id for legacy reasons: some projects shared phenotips
            # records by sharing the phenotips internal id, so in rare cases, the
            # Individual.objects.get(phenotips_patient_id=...) may match multiple Individual records
            individual = Individual.objects.get(phenotips_eid=patient_json['external_id'])
        else:
            individual = Individual.objects.get(phenotips_patient_id=patient_json['id'])

    except ObjectDoesNotExist as e:
        logger.error("ERROR: PhenoTips patient id %s not found in seqr Individuals." % patient_json['id'])
        return

    _update_individual_phenotips_data(individual, patient_json)


def _update_individual_phenotips_data(individual, patient_json):
    """Process and store the given patient_json in the given Individual model.

    Args:
        individual (Individual): Django Individual model
        patient_json (json): json dict representing the patient record in PhenoTips
    """

    # for each HPO term, get the top level HPO category (eg. Musculoskeletal)
    for feature in patient_json.get('features', []):
        hpo_id = feature['id']
        try:
            feature['category'] = HumanPhenotypeOntology.objects.get(hpo_id=hpo_id).category_id
        except ObjectDoesNotExist as e:
            logger.error("ERROR: PhenoTips HPO id %s not found in seqr HumanPhenotypeOntology table." % hpo_id)

    individual.phenotips_data = json.dumps(patient_json)
    individual.phenotips_patient_id = patient_json['id']  # phenotips internal id
    individual.phenotips_eid = patient_json.get('external_id')  # phenotips external id
    individual.save()


def _get_phenotips_uname_and_pwd_for_project(phenotips_user_id, read_only=False):
    """Return the PhenoTips username and password for this seqr project"""
    if not phenotips_user_id:
        raise ValueError("Invalid phenotips_user_id: " + str(phenotips_user_id))

    uname = phenotips_user_id + ('_view' if read_only else '')
    pwd = phenotips_user_id + phenotips_user_id

    return uname, pwd


def _get_phenotips_username_and_password(user, project, permissions_level):
    """Checks if user has permission to access the given project, and raises an exception if not.

    Args:
        user (User): the django user object
        project(Model): Project model
        permissions_level (string): 'edit' or 'view'
    Raises:
        PermissionDenied: if user doesn't have permission to access this project.
    Returns:
        2-tuple: PhenoTips username, password that can be used to access patients in this project.
    """
    if permissions_level == CAN_EDIT:
        uname, pwd = _get_phenotips_uname_and_pwd_for_project(project.phenotips_user_id, read_only=False)
    elif permissions_level == CAN_VIEW:
        uname, pwd = _get_phenotips_uname_and_pwd_for_project(project.phenotips_user_id, read_only=True)
    else:
        raise ValueError("Unexpected auth_permissions value: %s" % permissions_level)

    auth_tuple = (uname, pwd)

    return auth_tuple


def _convert_django_META_to_http_headers(meta_dict):
    """Converts django request.META dictionary into a dictionary of HTTP headers"""
    def convert_key(key):
        key = key.replace("HTTP_", "")
        tokens = key.split("_")
        capitalized_tokens = map(lambda x: x.capitalize(), tokens)
        return "-".join(capitalized_tokens)

    http_headers = {
        convert_key(key): str(value)
        for key, value in meta_dict.items()
        if key.startswith("HTTP_") or (key in ('CONTENT_LENGTH', 'CONTENT_TYPE') and value)
    }

    return http_headers


class PhenotipsException(Exception):
    pass

HTTP_REQUEST_HEADERS_TO_NOT_PROXY = { k.lower() for k in [
    'Connection',
    'X-Real-Ip',
    'X-Forwarded-Host',
]}

HTTP_RESPONSE_HEADERS_TO_NOT_PROXY = { k.lower() for k in [
    'Connection',
    'Keep-Alive',
    'Proxy-Authenticate',
    'Proxy-Authorization',
    'TE',
    'Trailers',
    'Transfer-Encoding',
    'Upgrade',
]}
