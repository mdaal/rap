### external imports
#import numpy as np
print("rap level init")


# Bring sweep() up to package top level so that  one creates a swp object as:
#     import rap
#     swp = rap.sweep()
# instead of:
#     import rap
#     swp = rap.sweeps.sweep.sweep()
from .sweeps.sweep import sweep
from .pulses.pulse import pulse
