import json


def readJSON(json_file):
    d  = []
    with open(json_file) as fd:
        for l in fd:
            d.append(l.strip())
    return json.loads(''.join(d))


def dumpJSON(json_file, content_dictionary):
    fd = open(json_file, 'w')
    fd.write(json.dumps(content_dictionary,
                        separators = (',', ':'),
                        indent     = 2).encode('utf8'))
    fd.close()
