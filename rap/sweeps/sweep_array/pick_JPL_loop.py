import numpy as np

def pick_JPL_loop(metadata, loop, Sweep_Array):
    freq, z = metadata.Loop_Data_Column_Names
    loop.index = 0
    # self.loop.normalization = None
    loop.z = Sweep_Array[0][z]
    loop.freq = Sweep_Array[0][freq]