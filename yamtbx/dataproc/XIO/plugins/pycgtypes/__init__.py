from __future__ import absolute_import
from __future__ import unicode_literals
# Initialization

# Import types into the cgtypes namespace
from .vec3 import vec3
from .vec4 import vec4
from .mat3 import mat3
from .mat4 import mat4
from .quat import quat

# CGKit version
# Applies to the whole kit and is a cheap way for applications to check
# which version is installed
_CGKit_version = "1.2.0"
