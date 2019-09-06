'''
the file descriptor returns a dictionary containing username, hostname, a uuid4 string and the utc time
'''

__all__ = ["filedescriptor"]
__version__ = "0.0.1"

import os
import datetime
import uuid

def filedescriptor():
    result = {}
    result["user"] = os.getenv("USER")
    result["host"] = os.uname()[1]
    result["utc"]  = str(datetime.datetime.utcnow())
    result["uuid4"] = str(uuid.uuid4())
    return result
