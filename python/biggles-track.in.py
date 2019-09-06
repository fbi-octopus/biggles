#!@PYTHON_EXECUTABLE@

import os
import sys

# Dodgy hack to automagically set the path for this script to allow for re-locatable tarball packages. If this script
# is at /foo/bar/bin/biggles, set the path to include /foo/bar/@PYTHON_SITE_PACKAGES@.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '@PYTHON_SITE_PACKAGES@')))

# Wrapper script for biggles command
if __name__ == '__main__':
    from biggles.tracking import main
    main()
