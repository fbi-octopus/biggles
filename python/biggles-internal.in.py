# Equivalent to invoking the biggles command but uses the biggles module in biggles _build_ path.
import sys
sys.path = '@BIGGLES_INTERNAL_PYTHONPATH@'.split(':') + sys.path

# Wrapper script for biggles command
if __name__ == '__main__':
    from biggles.tools import main
    main()
