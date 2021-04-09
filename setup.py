from setuptools import setup
import re
import os
VSRE = r"^__VERSION__ = ['\']([^'\']*)['\']"
with open(os.path.join('BFEE2', 'version.py')) as version_file:
    version = re.search(VSRE, version_file.read(), re.M).group(1)
setup(version=version)
