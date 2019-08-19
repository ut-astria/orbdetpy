# testodet.py - Example program for orbit determination.
# Copyright (C) 2018-2019 University of Texas
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import time
import argparse
from orbdetpy import determineOrbit

def main(args):
    print("OD start : %s" % time.strftime("%Y-%m-%d %H:%M:%S"),
          flush=True)
    determineOrbit(args.config, args.input, output_file = args.output)
    print("OD end   : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Test script for orbit determination.')
    parser.add_argument('config', help='Path to config file.')
    parser.add_argument('input', help='Path to input file.')
    parser.add_argument('output', help='Path to output file.')
    args = parser.parse_args()
    main(args)
