#!/usr/bin/env python
# encoding: utf-8
"""
test_report.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import os
import unittest

from ..report import Report
from ..graph import Graph
import fixtures

class TestReport(unittest.TestCase):
    def setUp(self):
        self.graphs={'performance_versus_cores': Graph(fixtures.GraphConfig('performance_versus_cores'))}
    def test_construct(self):
        r=Report(fixtures.ReportConfig('planck_performance'),self.graphs)
        self.assertEqual('Performance on planck',r.name)
        self.assertEqual({'type':'hemelb','machine':'planck'},r.graphs['performance_versus_cores'].select)
        self.assertEqual(['cores'],r.graphs['performance_versus_cores'].curves)
        self.assertEqual(['total'],r.graphs['performance_versus_cores'].dependent)
        self.assertEqual(['steps'],r.graphs['performance_versus_cores'].independent)
