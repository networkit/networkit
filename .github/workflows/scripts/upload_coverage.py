import requests
import json
import os

URL = os.getenv('COVERALLS_ENDPOINT', 'https://coveralls.io') + "/api/v1/jobs"

# Append additional data to the coveralls-API request
with open(os.getenv('COVERALLS_DUMP_FILE')) as dump:
    dumpdata = json.load(dump)
    dumpdata['service_name'] = 'Github Actions'
    dumpdata['service_number'] = os.getenv('COVERALLS_SVC_NUM')
    response = requests.post(URL, files={'json_file': json.dumps(dumpdata)},
                             verify=True)
    try:
        result = response.json()
    except ValueError:
        result = {'error': 'Failure to submit data. '
                  'Response [%(status)s]: %(text)s' % {
                      'status': response.status_code,
                      'text': response.text}}
    print(result)