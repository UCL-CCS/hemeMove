#!/usr/bin/env python
# encoding: utf-8
"""
report.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import os

class Report(object):
    def __init__(self,config,graphs):
        self.name=config['name']
        self.graphs={label:graphs[label].specialise(config['defaults'],specialisation) for label,specialisation in config['graphs'].iteritems()}

    def prepare(self,results):
        for graph in self.graphs.values():
            graph.prepare(results)
    
    def write(self,report_path):
        for label,graph in self.graphs.iteritems():
            graph.write_data(open(os.path.join(report_path,label),'w'))