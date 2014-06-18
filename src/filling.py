import numpy as np
powers = np.linspace(0,7,8)

max_len = 640*2**(-powers)
min_len = max_len/2.0
lambda_values = 0.01 * (640.0/max_len)**(3.0/2.0)

for i in range(8):
    run="./snap2fill.x /home/springel/universe/AnalysisIllustris/AddGrads/L75n455FP/snapdir_136/snap_136 ../data/gas_web/gasweb_snap_136 ../data/filling/ILL3_136 %f %f %f"%(min_len[i], max_len[i], lambda_values[i])
    print run
