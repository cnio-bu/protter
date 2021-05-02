from collections import OrderedDict
import re


def parse_sdrf_key_value_field(field):
    if field in ("not applicable","not available"):
        res = None
    else:
        subfields = re.split("\s*;\s*",field.strip('; '))
        res = OrderedDict([re.split("\s*=\s*",x,maxsplit=1)
                           for x in subfields])
    return res
