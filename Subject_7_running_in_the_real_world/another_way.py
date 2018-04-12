#!/usr/bin/env python

import subprocess as sp

## Another way to run!
output = sp.check_output('ls', shell=True)
output_list = output.split('\n')
for cur_line in output_list:
    if len (cur_line) > 0:
        print "wc",cur_line
