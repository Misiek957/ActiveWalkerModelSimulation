import numpy as np
import math

grad = np.pi * 3 / 5
curDir1 = 1.94
alpha = np.pi / 4
theta_max = 0.52


def forbidden_angle(curDir):
    global grad, alpha, theta_max
    # if abs((((grad-curDir)+np.pi/2) % np.pi)-np.pi/2) < alpha:
    #     if ((((grad - curDir) + np.pi / 2) % np.pi) - np.pi / 2) >= 0:  # Current direction more than 0, less than alpha, within the forbidden zone
    #         # NewcurDir = curDir - (((grad - curDir) + (np.pi / 2)) % np.pi/2) - (np.pi / 2) + alpha

    grad_max = math.tan(theta_max)
    # print(self.mod)
    # if self.mod < grad_max:
    #     # Ignore gradient processing if maximum gradient is not exceeded anywhere
    #     return curDir
    # else: # Calc forbidden angle alpha

    # convert to 2pi
    if curDir < 0:
        curDir = curDir % (2*np.pi)

    thetaForbMax = (grad + alpha) % (2*np.pi)
    thetaForbMin = (grad - alpha) % (2*np.pi)

    if (curDir % np.pi) <= thetaForbMin or (curDir % np.pi) >= thetaForbMax:
        return curDir # returns curDir unchanged as it is not in the forbidden zone

    check1 = (thetaForbMax - curDir) % (2*np.pi) # push left (AC)
    check2 = ((curDir%np.pi) - thetaForbMin) % (2*np.pi) # push right (C)

    if (thetaForbMax - curDir) % (2*np.pi) <= alpha:
        # Case 1: curDir lies closer to thetaForbMax, pushes left of gradient, (add), takes priority if curDir = theta
        dif = (thetaForbMax - curDir) % (2*np.pi)  # abs to be removed?
        curDir = (curDir + dif) % np.pi  # Modulus in case of overlapping the 2pi line
    elif ((curDir%np.pi) - thetaForbMin) % (2*np.pi) < alpha:
        # Case 2: curDir lies closer to thetaForbMin, pushes right of gradient, (subtract)
        dif = (curDir - thetaForbMin) % (2*np.pi)
        curDir = curDir - dif
    # curDir is [0,2pi]

    # convert to pi
    if 0 <= curDir <= np.pi:
        return curDir
    elif np.pi < curDir < 2 * np.pi:
        curDir = curDir - 2 * np.pi
        return curDir

curDir2 = forbidden_angle(curDir1)
print(curDir2)
