# coding: utf-8
#
# Copyright (C) 2012 André Panisson
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# -------------------------------------------------------------------
#
# This file has been modified by porting it to Python3.
#
# -------------------------------------------------------------------
# Changed type of autoflush attribute to integer, so autoflush every x
# send requests is possible.
#
"""
Allow a Python script to communicate with Gephi using the Gephi Graph Streaming protocol and plugin.
"""

__author__ = 'panisson@gmail.com'

import urllib.request, urllib.error, urllib.parse
import json
import time

class JSONClient(object):

    def __init__(self, autoflush=0, enable_timestamps=False, process_event_hook=None):
        self.data = ""
        self.autoflush = autoflush
        self.unflushedDumps = 0
        self.enable_timestamps = enable_timestamps

        if enable_timestamps:
            def default_peh(event):
                event['t'] = int(time.time())
                return event
        else:
            default_peh = lambda e: e
        if process_event_hook is None:
            self.peh = default_peh
        else:
            self.peh = lambda e: default_peh(process_event_hook(e))

    def flush(self):
        if len(self.data) > 0:
            self._send(self.data)
            self.data = ""
        
    def incrementUnflushedDumps(self):
        self.unflushedDumps = self.unflushedDumps + 1
        if self.unflushedDumps > self.autoflush:
            self.flush()
            self.unflushedDumps = 0
    def _send(self, data):
        print('passing')
        pass

    def add_node(self, id, flush=True, **attributes):
        self.data += json.dumps(self.peh({"an":{id:attributes}})) + '\r\n'
        self.incrementUnflushedDumps()

    def change_node(self, id, flush=True, **attributes):
        self.data += json.dumps(self.peh({"cn":{id:attributes}})) + '\r\n'
        self.incrementUnflushedDumps()

    def delete_node(self, id):
        self._send(json.dumps(self.peh({"dn":{id:{}}})) + '\r\n')

    def add_edge(self, id, source, target, directed=True, **attributes):
        attributes['source'] = source
        attributes['target'] = target
        attributes['directed'] = directed
        self.data += json.dumps(self.peh({"ae":{id:attributes}})) + '\r\n'
        self.incrementUnflushedDumps()

    def change_edge(self, id, source, target, directed=True, **attributes):
        attributes['source'] = source
        attributes['target'] = target
        attributes['directed'] = directed
        self.data += json.dumps(self.peh({"ce":{id:attributes}})) + '\r\n'
        self.incrementUnflushedDumps()

    def delete_edge(self, id):
        self._send(json.dumps(self.peh({"de":{id:{}}})) + '\r\n')

    def clean(self):
        self._send(json.dumps(self.peh({"dn":{"filter":"ALL"}})) + '\r\n')

class GephiClient(JSONClient):

    def __init__(self, url='http://127.0.0.1:8080/workspace0', autoflush=False):
        JSONClient.__init__(self, autoflush)
        self.url = url

    def _send(self, data):
        conn = urllib.request.urlopen(self.url+ '?operation=updateGraph', data.encode('utf-8'))
        return conn.read()

class GephiFileHandler(JSONClient):

    def __init__(self, out, **params):
        params['autoflush'] = True
        JSONClient.__init__(self, **params)
        self.out = out

    def _send(self, data):
        self.out.write(data)
