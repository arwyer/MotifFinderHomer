import os
import uuid
newpath = os.path.join("config", str(uuid.uuid4()))
print(newpath)
