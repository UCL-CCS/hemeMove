#!/usr/bin/env python
# encoding: utf-8
"""
test_report.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import os
import unittest
import tempfile

from ..report import Report
from ..graph import Graph
from ..results_collection import ResultsCollection
import fixtures

class TestReport(unittest.TestCase):
    def setUp(self):
        self.graphs={'performance_versus_cores': Graph(fixtures.GraphConfig('performance_versus_cores'))}
        self.r=Report(fixtures.ReportConfig('planck_performance'),self.graphs)
        self.results=ResultsCollection(fixtures.Results('cylinders').path,fixtures.ResultsConfig('example'))
    def test_construct(self):
        r=self.r
        self.assertEqual('Performance on planck',r.name)
        self.assertEqual({'type':'hemelb','machine':'planck'},r.graphs['performance_versus_cores'].select)
        self.assertEqual(['cores'],r.graphs['performance_versus_cores'].curves)
        self.assertEqual(['total'],r.graphs['performance_versus_cores'].dependent)
        self.assertEqual(['steps'],r.graphs['performance_versus_cores'].independent)
    def test_prepare(self):
        self.r.prepare(self.results.results)
    def test_write(self):
        self.r.prepare(self.results.results)
        self.r.write(fixtures.ReportOutput('planck_performance').path)
        