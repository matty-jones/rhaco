import numpy as np
import matplotlib.pyplot as plt

def loadAtoms(fileName):
    with open(fileName, 'r') as fileHandle:
        data = fileHandle.readlines()
    atomPosns = []
    for line in data:
        if 'HETATM' in line:
            atomPosn = [_ for _ in line.split(' ') if len(_) > 0][5:8]
            atomPosns.append(list(map(float, atomPosn)))
    return atomPosns


def calculateAngles(atomPosns):
    # Calculate angles for the following atom index vectors
    # There's a 'special' atom that's across a periodic cell image
    specialAtom = np.array(atomPosns[12]) + np.array([0.0, 26.738, 0.0])
    specialVectorPairs = [[[atomPosns[6], specialAtom], [atomPosns[6], atomPosns[114]]],
                          [[specialAtom, atomPosns[6]], [specialAtom, atomPosns[114]]]]
    angles = []
    calcAngles = [[[5, 6], [5, 153]],
                  [[6, 5], [6, 153]],
                  [[18, 15], [18, 62]],
                  [[15, 18], [15, 62]],
                  [[15, 43], [15, 59]],
                  [[43, 15], [43, 59]],
                  [[43, 5], [43, 129]],
                  [[5, 43], [5, 129]]]
    for vectorPair in calcAngles:
        firstVector = np.array(atomPosns[vectorPair[0][1]]) - np.array(atomPosns[vectorPair[0][0]])
        secondVector = np.array(atomPosns[vectorPair[1][1]]) - np.array(atomPosns[vectorPair[1][0]])
        angle = np.arccos(np.dot(firstVector, secondVector) / (np.linalg.norm(firstVector) * np.linalg.norm(secondVector)))
        angles.append(angle)
    for vectorPair in specialVectorPairs:
        firstVector = np.array(vectorPair[0][1]) - np.array(vectorPair[0][0])
        secondVector = np.array(vectorPair[1][1]) - np.array(vectorPair[1][0])
        angles.append(np.arccos(np.dot(firstVector, secondVector) / (np.linalg.norm(firstVector) * np.linalg.norm(secondVector))))
    #plt.figure()
    #plt.hist(angles)
    #plt.show()
    return np.mean(angles), [specialAtom, np.array(atomPosns[18]), np.array(atomPosns[38])]


def getNewAtomPosn(planeAtoms, angle):
    # Oxygen atom position is defined as:
    # 1) Half way between the two bonded metal atoms
    # 2) In the same plane as the two bonded metal atoms and the central ionically bonded species
    # 3) Offset from the line between the two bonded metal atoms towards the ionically bonded species,
    #    such that the angle between the vector describing the bond between the two bonded metal atoms
    #    and the vector describing the bond between one of the metals and the oxygen atom is an angle 
    #    of 0.274 rad.
    # planeAtoms = [metalAtom1Posn, metalAtom2Posn, ionicCentralAtomPosn]
    newAtomCoords = 0.5 * (np.array(planeAtoms[0]) + np.array(planeAtoms[1]))
    # Coordinate transformation to set the origin as planeAtoms[0]
    [x, y, z] = newAtomCoords
    # Use MorphCT's master rotation matrix infrastructure
    # Rotation matrix calculations from: http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
    # The array that describes the 3D rotation of (x, y, z) around the point (a, b, c) through
    # the unit axis <u, v, w> by the angle theta is given by:
    # [ (a(v^2 + w^2) - u(bv + cw - ux - vy - wz))(1 - cos(theta)) + x*cos(theta) + (-cv + bw - wy + vz)sin(theta),
    #   (b(u^2 + w^2) - v(au + cw - ux - vy - wz))(1 - cos(theta)) + y*cos(theta) + (cu - aw + wx - uz)sin(theta),
    #   (c(u^2 + v^2) - w(au + bv - ux - vy - wz))(1 - cos(theta)) + z*cos(theta) + (-bu + av - vx + uy)sin(theta) ]
    [a, b, c] = np.array(planeAtoms[0])
    vec1 = (np.array(planeAtoms[1]) - np.array(planeAtoms[0])) / np.linalg.norm(np.array(planeAtoms[1]) - np.array(planeAtoms[0]))
    vec2 = (np.array(planeAtoms[2]) - np.array(planeAtoms[0])) / np.linalg.norm(np.array(planeAtoms[2]) - np.array(planeAtoms[0]))
    [u, v, w] = np.cross(vec1, vec2) / np.linalg.norm(np.cross(vec1, vec2))
    theta = angle
    newPosition = np.array([(a * (v**2 + w**2) - u * ((b * v) + (c * w) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (x * np.cos(theta)) + ((-(c * v) + (b * w) - (w * y) + (v * z)) * np.sin(theta)),
                   (b * (u**2 + w**2) - v * ((a * u) + (c * w) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (y * np.cos(theta)) + (((c * u)  - (a * w) + (w * x) - (u * z)) * np.sin(theta)),
                   (c * (u**2 + v**2) - w * ((a * u) + (b * v) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (z * np.cos(theta)) + ((-(b * u) + (a * v) - (v * x) + (u * y)) * np.sin(theta))])
    print("New Atom Coords =", newPosition)



if __name__ == "__main__":
    atomPosns = loadAtoms('templateM1.pdb')
    averageAngle, planeAtoms = calculateAngles(atomPosns)
    print("Mean angle =", averageAngle)
    getNewAtomPosn(planeAtoms, averageAngle)

