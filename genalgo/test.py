import numpy as np
import bladecoordinates
import STLgen
import math
import matplotlib.pyplot as plt
lead1=[-0.012054612,0.196694492,0]
lead2=[0.050003115,0.200079,0]
trail1=[0.036167406,0.209908331,0]
trail2=[0.1,0.200079,0]
target1=bladecoordinates.blade1
target2=bladecoordinates.blade2
target1=target1.transpose()
target2=target2.transpose()
print(target1[0])
target1[0] -= trail1[0]
target2[0] -= lead2[0]
target1[1] -= trail1[1]
target2[1] -= lead2[1]
slope1 = math.degrees(math.atan(abs((lead1[1] - trail1[1]) / (lead1[0] - trail1[0]))))
slope2 = math.degrees(math.atan(abs((lead2[1] - trail2[1]) / (lead2[0] - trail2[0]))))
target1=target1.transpose()
target2=target2.transpose()
target1 = STLgen.rotate(target1,slope1, [0, 0, 1])
target2 = STLgen.rotate(target2,slope2, [0, 0, 1])
target1 = STLgen.rotate(target1,-19, [0, 0, 1])
target2 = STLgen.rotate(target2,19, [0, 0, 1])
target1=target1.transpose()
target2=target2.transpose()
target1[0] -= 0.01
target1=target1.transpose()
target2=target2.transpose()
plt.plot(target1.transpose()[0],target1.transpose()[1])
plt.plot(target2.transpose()[0],target2.transpose()[1])
plt.axes().set_aspect('equal')
plt.show()
