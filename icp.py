import numpy as np
from scipy.spatial.distance import cdist

def best_fit_transform(A, B):
    '''
    Input:
      A: Nx3 numpy array of corresponding 3D points
      B: Nx3 numpy array of corresponding 3D points
    Returns:
      T: 4x4 homogeneous transformation matrix
      R: 3x3 rotation matrix
      t: 3x1 column vector
    '''

    assert len(A) == len(B)

    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    AA = A - centroid_A
    BB = B - centroid_B

    H = np.dot(AA.T, BB)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)

    if np.linalg.det(R) < 0:
       Vt[2,:] *= -1
       R = np.dot(Vt.T, U.T)

    t = centroid_B.T - np.dot(R,centroid_A.T)

    T = np.identity(4)
    T[0:3, 0:3] = R
    T[0:3, 3] = t

    return T, R, t

def nearest_neighbor(src, dst):
    '''
    Input:
        src: Nx3 array of points
        dst: Nx3 array of points
    Output:
        distances: Euclidean distances of the nearest neighbor
        indices: dst indices of the nearest neighbor
    '''

    all_dists = cdist(src, dst, 'euclidean')
    indices = all_dists.argmin(axis=1)
    distances = all_dists[np.arange(all_dists.shape[0]), indices]
    return distances, indices

def icp(A, B, init_pose=None, max_iterations=20, tolerance=0.001):
    '''
    Input:
        A: Nx3 numpy array of source 3D points
        B: Nx3 numpy array of destination 3D point
        init_pose: 4x4 homogeneous transformation
        max_iterations: exit algorithm after max_iterations
        tolerance: convergence criteria
    Output:
        T: final homogeneous transformation
        distances: Euclidean distances (errors) of the nearest neighbor
    '''

    src = np.ones((4,A.shape[0]))
    dst = np.ones((4,B.shape[0]))
    src[0:3,:] = np.copy(A.T)
    dst[0:3,:] = np.copy(B.T)

    if init_pose is not None:
        src = np.dot(init_pose, src)

    prev_error = 0

    for i in range(max_iterations):
        # find the nearest neighbours between the current source and destination points
        distances, indices = nearest_neighbor(src[0:3,:].T, dst[0:3,:].T)

        # compute the transformation between the current source and nearest destination points
        T,_,_ = best_fit_transform(src[0:3,:].T, dst[0:3,indices].T)

        # update the current source
        src = np.dot(T, src)

        # check error
        mean_error = np.sum(distances) / distances.size
        if abs(prev_error-mean_error) < tolerance:
            break
        prev_error = mean_error

    # calculate final transformation
    T,_,_ = best_fit_transform(A, src[0:3,:].T)

    return T, distances

#testing for point set matching
A = np.zeros((8,3))
B = np.zeros((11,3))

A[0][0] = 43.89
A[0][1] = -5.88
A[0][2] = 106.99
A[1][0] = 42.02
A[1][1] = 20.52
A[1][2] = 112.52
A[2][0] = 42.01
A[2][1] = 25.39
A[2][2] = 113.25
A[3][0] = 44.95
A[3][1] = 4.69
A[3][2] = 112.60
A[4][0] = 44.12
A[4][1] = 17.96
A[4][2] = 115.15
A[5][0] = 48.26
A[5][1] = -1.37
A[5][2] = 113.59
A[6][0] = 46.28
A[6][1] = 7.03
A[6][2] = 114.58
A[7][0] = 47.0
A[7][1] = 18.52
A[7][2] = 117.65

B[0] = np.array([72.78,7.12,146.10])
B[1] = np.array([70.19,24.80,148.67])
B[2] = np.array([76.21,18.28,147.20])
B[3] = np.array([72.71,17.69,148.09])
B[4] = np.array([70.67,4.62,145.95])
B[5] = np.array([64.38,-5.40,143.59])
B[6] = np.array([72.47,-1.16,143.85])
B[7] = np.array([69.82,19.81,148.32])
B[8] = np.array([77.0,25.0,150.0])
B[9] = np.array([80.0,-10.0,140.0])
B[10] = np.array([83.0,30.0,145.0])

T,d = icp(A,B)
print T

#also get the final positions of the points after the transformation
R_f = T[0:3, 0:3]
t_f = T[0:3, 3]

A_f = np.dot(R_f, A.T).T + t_f

print "A: "
print A_f
