import elasticsearch
import json
import settings
from seqr.utils.redis_utils import get_redis_client

VARIANT_DOC_TYPE = 'variant'


# make encoded values as human-readable as possible
ES_FIELD_NAME_ESCAPE_CHAR = '$'
ES_FIELD_NAME_BAD_LEADING_CHARS = set(['_', '-', '+', ES_FIELD_NAME_ESCAPE_CHAR])
ES_FIELD_NAME_SPECIAL_CHAR_MAP = {
    '.': '_$dot$_',
    ',': '_$comma$_',
    '#': '_$hash$_',
    '*': '_$star$_',
    '(': '_$lp$_',
    ')': '_$rp$_',
    '[': '_$lsb$_',
    ']': '_$rsb$_',
    '{': '_$lcb$_',
    '}': '_$rcb$_',
}


def get_es_client():
    return elasticsearch.Elasticsearch(host=settings.ELASTICSEARCH_SERVICE_HOSTNAME)


def _encode_name(s):
    """Applies a reversable encoding to the special chars in the given name or id string, and returns the result.
    Among other things, the encoded string is a valid elasticsearch or mongodb field name.

    See:
    https://discuss.elastic.co/t/special-characters-in-field-names/10658/2
    https://discuss.elastic.co/t/illegal-characters-in-elasticsearch-field-names/17196/2
    """
    s = s.replace(ES_FIELD_NAME_ESCAPE_CHAR, 2 * ES_FIELD_NAME_ESCAPE_CHAR)
    for original_value, encoded in ES_FIELD_NAME_SPECIAL_CHAR_MAP.items():
        s = s.replace(original_value, encoded)
    if s[0] in ES_FIELD_NAME_BAD_LEADING_CHARS:
        s = ES_FIELD_NAME_ESCAPE_CHAR + s
    return s


def _decode_name(s):
    """Decodes a name or id string that was encoded by #_encode_name(..) and returns the original string"""

    if s.startswith(ES_FIELD_NAME_ESCAPE_CHAR):
        s = s[1:]
    for original_value, encoded in ES_FIELD_NAME_SPECIAL_CHAR_MAP.items():
        s = s.replace(encoded, original_value)
    s = s.replace(2*ES_FIELD_NAME_ESCAPE_CHAR, ES_FIELD_NAME_ESCAPE_CHAR)
    return s


def get_mapping(index_pattern):
    es_client = get_es_client()
    redis_client = get_redis_client()

    mappings = None
    redis_key = "es_utils__get_mapping__{}".format(index_pattern)
    if redis_client:
        mappings_string = redis_client.get(redis_key)
        if mappings_string is not None:
            mappings = json.loads(mappings_string)

    if mappings is None:
        mappings = es_client.indices.get_mapping(index_pattern)

    if redis_client:
        mappings_string = json.dumps(mappings)
        redis_client.set(redis_key, mappings_string)

    return mappings


def find_index_that_contains_sample(index_pattern, sample_id):
    """Given an elasticsearch index pattern (eg. 'index_name_prefix*') and a sample id, return the specific index name
    that matches index_pattern and contains genotypes for this sample (eg. 'index_name_prefix_3').

    Args:
         index_pattern (string):
         sample_id (string):

    Return:
        list: index names that match the given pattern and contain the given sample, or an empty list if no matches
    """

    mappings = get_mapping(index_pattern)

    encoded_sample_id = _encode_name(sample_id)
    genotype_field = "{}_num_alt".format(encoded_sample_id)

    matching_indices = []
    for index_name, index_mapping in mappings.items():
        if genotype_field in index_mapping["mappings"]["variant"]["properties"]:
            matching_indices.append(index_name)

    return matching_indices


def get_index_meta(index_name, doc_type_name="*"):
    _meta = {}
    es_client = get_es_client()
    mappings = es_client.indices.get_mapping(index=index_name, doc_type=doc_type_name)
    for index_type_name, mapping in mappings.get(index_name, {}).get('mappings', {}).items():
        _meta.update(mapping.get("_meta", {}))

    return _meta


def set_index_meta(self, index_name, index_type_name, _meta):
    index_mapping = {
        "_meta": _meta,
    }

    es_client = get_es_client()
    es_client.indices.put_mapping(index=index_name, doc_type=index_type_name, body=index_mapping)
