#!/usr/bin/env python
# encoding: utf-8
"""
curate.py

Created by James Hetherington on 2012-01-27.
Copyright (c) 2012 UCL. All rights reserved.
"""

import sys
import os
import argparse
import shutil
import socket
import time
from config import config, config_user
from remote_hemelb import RemoteHemeLB

class Driver(object):
    """Gather a group of results, and manipulate them in some way"""
    def __init__(self,clargs):
        # By default, an unsupplied argument does not create a result property
        self.parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
        self.parser.add_argument("--retry",action='store_true',default=False)
        self.define_args()
        options,extra=self.parser.parse_known_args(clargs)

        remote=extra[1]
        self.options=vars(options)

        if remote in config:
            config.update(config.get(remote,{}))
            config.update(config_user)
            config.update(config.get(remote,{}))
        else:
            config['address']=remote

        self.port=config['port']
        self.steering_id=config['steering_id']
        self.address=config['address']
        
        # This Remote will have now connected, but not attempted to send/receive
        if  self.options['retry']:
            while True:
                try:
                    self.hemelb=RemoteHemeLB(port=self.port,address=self.address,steering_id=self.steering_id)
                    break
                except socket.error as (errno,message):
                    if errno!=61:
                        raise
                    else:
                        print("# No steering server, will retry")
                        time.sleep(10)
        else:
            self.hemelb=RemoteHemeLB(port=self.port,address=self.address,steering_id=self.steering_id)

    def define_args(self):
       pass
       
    def act(self):
        self.hemelb.step()

def main():
    Intervention(sys.argv).act()

if __name__ == '__main__':
    main()
