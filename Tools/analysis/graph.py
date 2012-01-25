#!/usr/bin/env python
# encoding: utf-8
"""
graph.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import copy
import itertools
import csv

class Graph(object):
    def __init__(self,config):
        self.select={}
        for prop in ['name','select','curves','dependent','independent']:
            if prop in config:
                setattr(self,prop,config[prop])
            
    def specialise(self,*specialisations):
        """Return a copy of this graph, with the configuration modified by some additional dicts of dicts.
        """
        result=copy.deepcopy(self)
        for specialisation in specialisations:
            if not specialisation:
                continue
            if 'select' in specialisation:
                result.select.update(specialisation['select'])
            for prop in ['curves','dependent','independent']:
                if prop in specialisation:
                    getattr(result,prop).append(specialisation[prop])
        return result
    def prepare(self,results):
        def filtration(result):
            return all([getattr(result,prop)==value for prop,value in self.select.iteritems()])
        self.filtered_results=filter(filtration,results)
        def curve_key(result):
            return tuple([getattr(result,prop) for prop in self.curves])
        self.curve_results={key:list(group) for key,group in itertools.groupby(sorted(self.filtered_results,key=curve_key),curve_key)}
        # for now, support only a single dependent and a single independent variable on each curve
        self.curve_data={curve:[ (getattr(result,self.dependent[0]),getattr(result,self.independent[0])) for result in results ] for curve,results in self.curve_results.iteritems()}
    def write_data(self,csv_file):
        writer=csv.writer(csv_file)
        fieldnames=[self.dependent[0],self.independent[0]]
        csv_file.write("#CSV File produced by hemelb reporter\n")
        csv_file.write("#Groups are separated by "+str(tuple(self.curves))+"\n#")
        writer.writerow(fieldnames)
        for curve,results in sorted(self.curve_data.iteritems()):
            writer.writerow([])
            writer.writerow([])
            csv_file.write("#"+str(curve)+"\n")
            writer.writerows(results)
            