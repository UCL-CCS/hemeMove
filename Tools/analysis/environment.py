#!/usr/bin/env python
# encoding: utf-8
"""
environment.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 University College London. All rights reserved.
"""
import yaml
import os

localroot=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
#Load and invoke the default non-machine specific config JSON dictionaries.
config=yaml.load(open(os.path.join(localroot,'Tools','analysis','config_defaults.yml')))
config.update(yaml.load(open(os.path.join(localroot,'Tools','analysis','config.yml'))))

